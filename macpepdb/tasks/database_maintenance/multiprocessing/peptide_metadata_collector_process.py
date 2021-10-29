# std imports
from multiprocessing import Queue, Event, Array
from multiprocessing.connection import Connection as ProcessConnection
from psycopg2.extras import execute_batch
from queue import Empty as EmptyQueueError
from typing import ClassVar

# external imports
import psycopg2

# internal imports
from macpepdb.models.peptide import Peptide
from macpepdb.models.peptide_metadata import PeptideMetadata
from macpepdb.utilities.generic_process import GenericProcess


class PeptideMetadataCollectorProcess(GenericProcess):
    """
    Collects and update metadata for peptides.

    Parameters
    ----------
    termination_event : Event
        Event for terminating the process
    id : int
        ID for identification in logs
    database_url : str
        Database URL, e.g. postgres://username:password@host:port/database
    peptides_queue : Queue
        Queue for peptides
    update_counter : Array
        Counter for updates
    log_connection : multiprocessing.connection.Connection
        Connection to the log process
    empty_queue_and_stop_flag : Event
        Event which indicates that the process can stop as soon as the queue is empty.
    """


    DECLARED_PEPTIDE_METADATA_UPSERT_NAME: ClassVar[str] = "update_peptide_metadata"
    """Name of the prepared statement for peptide metadata insert/update (upsert)
    """

    DECLARED_PEPTIDE_METADATA_UPSERT: ClassVar[str] = (
        f"PREPARE {DECLARED_PEPTIDE_METADATA_UPSERT_NAME} AS "
        f"INSERT INTO {PeptideMetadata.TABLE_NAME} (partition, mass, sequence, is_swiss_prot, is_trembl, taxonomy_ids, unique_taxonomy_ids, proteome_ids) VALUES ($1, $2, $3, $4, $5, $6, $7, $8) "
        "ON CONFLICT (partition, mass, sequence) DO "
        f"UPDATE SET is_swiss_prot = $4, is_trembl = $5, taxonomy_ids = $6, unique_taxonomy_ids = $7, proteome_ids = $8 WHERE {PeptideMetadata.TABLE_NAME}.partition = $1 AND {PeptideMetadata.TABLE_NAME}.mass = $2 AND {PeptideMetadata.TABLE_NAME}.sequence = $3;"
    )
    """Query to prepare statement for prepared statement for peptide metadata insert/update (upsert).  Order of placeholder: partition, mass, sequence, is_swiss_prot, is_trembl, taxonomy_ids, unique_taxonomy_ids, proteome_ids
    """

    DECLARED_PEPTIDE_UPDATE_METADATA_STATUS_NAME: ClassVar[str] = "update_peptide_metadata_update_flag"
    """Name of the prepared statement for flagging the peptide as updated.
    """

    DECLARED_PEPTIDE_UPDATE_METADATA_STATUS: ClassVar[str] = (
        f"PREPARE {DECLARED_PEPTIDE_UPDATE_METADATA_STATUS_NAME} AS "
        f"UPDATE {Peptide.TABLE_NAME} SET is_metadata_up_to_date = true WHERE partition = $1 AND mass = $2 AND sequence = $3;"
    )
    """Query to prepare statement for flagging the peptide as updated. Order of placeholder: partition, mass, sequence
    """

    def __init__(self, termination_event: Event, id: int, database_url: str, peptides_queue: Queue, update_counter: Array, log_connection: ProcessConnection, empty_queue_and_stop_flag: Event):
        super().__init__(termination_event)
        self.__id = id
        self.__database_url = database_url
        self.__peptides_queue = peptides_queue
        self.__update_counter = update_counter
        self.__log_connection = log_connection
        self.__empty_queue_and_stop_flag = empty_queue_and_stop_flag

    def run(self):
        """
        Collects and updates the information from the referenced peptides und set update flag to false.
        """
        self.activate_signal_handling()
        self.__log_connection.send(f"peptide update worker {self.__id} is online")
        database_connection = None

        # Let the process run until empty_queue_and_stop_flag is true and peptides_queue is empty or database_maintenance_stop_event is true.
        while (not self.__empty_queue_and_stop_flag.is_set() or not self.__peptides_queue.empty()) and not self.termination_event.is_set():
            try:
                # Open/reopen database connection
                if not database_connection or (database_connection and database_connection.closed != 0):
                    database_connection = psycopg2.connect(self.__database_url)
                    with database_connection:
                        with database_connection.cursor() as database_cursor:
                            database_cursor.execute(self.__class__.DECLARED_PEPTIDE_METADATA_UPSERT)
                            database_cursor.execute(self.__class__.DECLARED_PEPTIDE_UPDATE_METADATA_STATUS)
                # Wait 5 seconds for a new set of peptides
                peptides = self.__peptides_queue.get(True, 5)
                while True:
                    try:
                        with database_connection:
                            with database_connection.cursor() as database_cursor:
                                # Make sure each peptides contains its metadata
                                for peptide in peptides:
                                    peptide.fetch_metadata_from_proteins(database_cursor)
                                PeptideMetadata.bulk_insert(
                                    database_cursor,
                                    peptides,
                                    # Override statement with the prepared statement
                                    statement=f"EXECUTE {self.__class__.DECLARED_PEPTIDE_METADATA_UPSERT_NAME} (%s, %s, %s ,%s ,%s, %s, %s, %s);",
                                    is_prepared_statement=True
                                )
                                execute_batch(
                                    database_cursor,
                                    f"EXECUTE {self.__class__.DECLARED_PEPTIDE_UPDATE_METADATA_STATUS_NAME} (%s, %s, %s)",
                                    [(peptide.partition, peptide.mass, peptide.sequence) for peptide in peptides],
                                    page_size=len(peptides)
                                )
                        with self.__update_counter:
                            self.__update_counter[0] += len(peptides)
                        break
                    except BaseException as error:
                        self.__log_connection.send(f"error occured, see:\n{error}\ntry again")
            # Catch errors which occure during database connect
            except psycopg2.Error as error:
                self.__log_connection.send("Error when opening the database connection, see:\n{}".format(error))
            # Catch queue.Empty which is thrown when protein_queue.get() timed out
            except EmptyQueueError:
                pass

        self.__log_connection.send(f"peptide update worker {self.__id} is stopping")
        if database_connection and database_connection.closed != 0:
            with database_connection:
                with database_connection.cursor() as database_cursor:
                    database_cursor.execute(f"DEALLOCATE {self.__class__.DECLARED_PEPTIDE_METADATA_UPSERT_NAME};")
                    database_cursor.execute(f"DEALLOCATE {self.__class__.DECLARED_PEPTIDE_UPDATE_METADATA_STATUS_NAME};")
            database_connection.close()
        self.__log_connection.close()
