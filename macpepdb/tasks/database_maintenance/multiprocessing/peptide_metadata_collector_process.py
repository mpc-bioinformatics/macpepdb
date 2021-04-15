import psycopg2
import traceback

from multiprocessing import Queue, Event, Array
from multiprocessing.connection import Connection as ProcessConnection
from psycopg2.extras import execute_batch
from queue import Empty as EmptyQueueError


from ....models.peptide import Peptide
from ....utilities.generic_process import GenericProcess


class PeptideMetadataCollectorProcess(GenericProcess):
    def __init__(self, id: int, database_url: str, peptide_sequence_queue: Queue, update_counter: Array, log_connection: ProcessConnection, empty_queue_and_stop_flag: Event, immediate_stop_event: Event):
        super().__init__()
        self.__id = id
        self.__database_url = database_url
        self.__peptide_sequence_queue = peptide_sequence_queue
        self.__update_counter = update_counter
        self.__log_connection = log_connection
        self.__empty_queue_and_stop_flag = empty_queue_and_stop_flag
        self.__immediate_stop_event = immediate_stop_event

    def run(self):
        """
        Collects the information from the referenced peptides und set update flag to false.
        """
        self.__log_connection.send(f"peptide update worker {id} is online")
        database_connection = None
        PREPARED_STATEMENT_NAME = "updatepeptide_metadata"
        PREPARE_STATEMENT_QUERY = f"PREPARE {PREPARED_STATEMENT_NAME} AS UPDATE {Peptide.TABLE_NAME} SET is_metadata_up_to_date = true, is_swiss_prot = $1, is_trembl = $2, taxonomy_ids = $3, unique_taxonomy_ids = $4, proteome_ids = $5 WHERE sequence = $6;"

        # Let the process run until empty_queue_and_stop_flag is true and peptide_sequence_queue is empty or immediate_stop_event is true.
        while (not self.__empty_queue_and_stop_flag.is_set() or not self.__peptide_sequence_queue.empty()) and not self.__immediate_stop_event.is_set():
            try:
                # Open/reopen database connection
                if not database_connection or (database_connection and database_connection.closed != 0):
                    database_connection = psycopg2.connect(self.__database_url)
                    with database_connection:
                        with database_connection.cursor() as database_cursor:
                            database_cursor.execute(PREPARE_STATEMENT_QUERY)
                # Wait 5 seconds for a new set of IDs
                peptide_sequences = self.__peptide_sequence_queue.get(True, 5)
                while True:
                    try:
                        with database_connection:
                            with database_connection.cursor() as database_cursor:
                                execute_batch(
                                    database_cursor,
                                    f"EXECUTE {PREPARED_STATEMENT_NAME} (%s, %s, %s ,%s ,%s, %s);",
                                    [Peptide.generate_metadata_update_values(database_cursor, peptide_sequence) for peptide_sequence in peptide_sequences],
                                    page_size=len(peptide_sequences)
                                )
                        with self.__update_counter:
                            self.__update_counter[0] += len(peptide_sequences)
                        break
                    except BaseException as error:
                        self.__log_connection.send(f"error occured, see:\n{error}\ntry again")
            # Catch errors which occure during database connect
            except psycopg2.Error as error:
                self.__log_connection.send("Error when opening the database connection, see:\n{}".format(error))
            # Catch queue.Empty which is thrown when protein_queue.get() timed out
            except EmptyQueueError:
                pass

        self.__log_connection.send(f"peptide update worker {id} is stopping")
        if database_connection and database_connection.closed != 0:
            with database_connection:
                with database_connection.cursor() as database_cursor:
                    database_cursor.execute(f"DEALLOCATE {PREPARED_STATEMENT_NAME};")
            database_connection.close()
        self.__log_connection.close()
