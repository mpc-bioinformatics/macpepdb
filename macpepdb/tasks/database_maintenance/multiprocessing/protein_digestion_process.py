# std imports
import random
import time
from multiprocessing import Event, Queue, Array
from multiprocessing.connection import Connection as ProcessConnection
from queue import Empty as EmptyQueueError

# external imports
import psycopg2

# internal imports
from macpepdb.models.protein import Protein
from macpepdb.proteomics.enzymes.digest_enzyme import DigestEnzyme
from macpepdb.utilities.generic_process import GenericProcess


class ProteinDigestionProcess(GenericProcess):
    UNSOLVEABLE_ERROR_FACTOR_LIMIT = 2
    """
    Sequentially digests proteins from the given queue and inserts them and their proteins into the given database.
    """

    def __init__(self, termination_event: Event, id: int, database_url: str, protein_queue: Queue, enzyme: DigestEnzyme, general_log: ProcessConnection, unprocessible_protein_log: ProcessConnection, statistics: Array, finish_event: Event):
        """
        termination_event : Event
            Event for terminating the process
        id : int
            ID for identification in logs
        database_url : str
            Database URL, e.g. postgres://username:password@host:port/database
        protein_queue : Queue
            Protein queue
        enzyme : DigestEnzyme
            Digest enzyme
        general_log : multiprocessing.connection.Connection
            Connection to the log process
        unprocessible_protein_log : multiprocessing.connection.Connection
            connection to log process which logs unprocessible proteins
        statistics : Array
            Shared array to collect statistics [inserted proteins, inserted peptides, errors]
        finish_event : Event
            Event which indicates that the process can stop as soon as the queue is empty.
        """
        super().__init__(termination_event)
        self.__id = id
        self.__database_url = database_url
        self.__protein_queue = protein_queue
        self.__general_log = general_log
        self.__unprocessible_protein_log = unprocessible_protein_log
        self.__statistics = statistics
        self.__enzyme = enzyme
        self.__finish_event = finish_event

    def run(self):
        """
        Starts the process which digests proteins and inserts the peptides into the database.
        """
        self.activate_signal_handling()
        self.__general_log.send("digest worker {} is online".format(self.__id))
        database_connection = None

        # Let the process run until finish_event is true and protein_queue is empty or termination_event is true.
        while (not self.__finish_event.is_set() or not self.__protein_queue.empty()) and not self.termination_event.is_set():
            try:
                # Open/reopen database connection
                if not database_connection or (database_connection and database_connection.closed != 0):
                    database_connection = psycopg2.connect(self.__database_url)

                # Try to get a protein from the queue, timeout is 2 seconds
                protein = self.__protein_queue.get(True, 5)
                
                # Variables for loop control
                unsolvable_error_factor = 0
                try_transaction_again = True
                while try_transaction_again:
                    number_of_new_peptides = 0
                    error = None
                    try:
                        count_protein = False
                        number_of_new_peptides = 0
                        with database_connection:
                            with database_connection.cursor() as database_cursor:
                                existing_protein = Protein.select(database_cursor, ("accession = %s", (protein.accession,)))
                                if existing_protein is not None:
                                    number_of_new_peptides = existing_protein.update(database_cursor, protein, self.__enzyme)
                                else:
                                    number_of_new_peptides = Protein.create(database_cursor, protein, self.__enzyme)
                                    count_protein = True

                        # Commit was successfully stop while-loop and add statistics
                        try_transaction_again = False
                        self.__statistics.acquire()
                        if count_protein:
                            self.__statistics[0] += 1
                        self.__statistics[1] += number_of_new_peptides
                        self.__statistics.release()
                    # Rollback is done implcit by `with database_connection`
                    # Each error increases the unsolveable error factor differently. If the factor reaches UNSOLVEABLE_ERROR_FACTOR_LIMIT the protein is logged as unprocessible
                    ## Catch violation of unique constraints. Usually a peptide which is already inserted by another transaction.
                    except psycopg2.errors.UniqueViolation as unique_violation_error:
                        error = unique_violation_error
                        if unsolvable_error_factor < self.__class__.UNSOLVEABLE_ERROR_FACTOR_LIMIT:
                            unsolvable_error_factor += 0.2
                    ## Catch deadlocks between transactions. This occures usually when 2 transactions try to insert the same peptides 
                    except psycopg2.errors.DeadlockDetected as deadlock_detected_error:
                        error = deadlock_detected_error
                        # Try again after 5 (first try) and 10 (second try) + a random number between 0 and 5 (both including) seconds maybe some blocking transactions can pass so this transaction will successfully finish on the next try.
                        if unsolvable_error_factor < self.__class__.UNSOLVEABLE_ERROR_FACTOR_LIMIT:
                            unsolvable_error_factor += 1
                            time.sleep(5 * unsolvable_error_factor + random.randint(0, 5))
                    ## Catch other errors.
                    except psycopg2.Error as base_error:
                        unsolvable_error_factor += self.__class__.UNSOLVEABLE_ERROR_FACTOR_LIMIT
                        error = base_error
                    finally:
                        # Log the last error if the unsolvable_error_factor exceeds the limit
                        if unsolvable_error_factor >= self.__class__.UNSOLVEABLE_ERROR_FACTOR_LIMIT:
                            self.__general_log.send("Exception on protein {}, see:\n{}".format(protein.accession, error))
                            self.__unprocessible_protein_log.send(protein.to_embl_entry())
                            self.__statistics.acquire()
                            self.__statistics[2] += 1
                            self.__statistics.release()
                            try_transaction_again = False
            # Catch errors which occure during database connect
            except psycopg2.Error as error:
                self.__general_log.send("Error when opening the database connection, see:\n{}".format(error))
            # Catch queue.Empty which is thrown when protein_queue.get() timed out
            except EmptyQueueError:
                pass
        # Close database connection
        if database_connection and database_connection.closed == 0:
            database_connection.close()
        self.__general_log.send("digest worker {} is stopping".format(self.__id))
        self.__general_log.close()
        self.__unprocessible_protein_log.close()