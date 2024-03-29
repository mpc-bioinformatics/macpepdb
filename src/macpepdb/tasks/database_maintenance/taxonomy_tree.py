# std imports
import pathlib
import time
import random
import signal
from multiprocessing import Process, Event, Queue, Pipe, get_context as get_processing_context
from multiprocessing.connection import Connection, wait
from queue import Full as QueueFullError, Empty as QueueEmptyError
from typing import Tuple

# external imports
import psycopg2

# internal imports
from macpepdb.tasks.database_maintenance.multiprocessing.logger_process import LoggerProcess
from macpepdb.models.taxonomy import Taxonomy, TaxonomyRank
from macpepdb.models.taxonomy_merge import TaxonomyMerge

class TaxonomyTree:
    """
    Controls the creation and updates of the taxonomy tree.

    Parameters
    ----------
    termination_event: Event
        Event for process termination
    nodes_dmp_path: pathlib.Path
        Path the nodes.dmp within the work directory
    names_dmp_path: pathlib.Path
        Path the names.dmp within the work directory
    merge_dmp_path: pathlib.Path
        Path the merge.dmp within the work directory
    delete_dmp_path: pathlib.Path
        Path the delnodes.dmp within the work directory
    database_url: str
        Database URL, e.g. postgres://username:password@host:port/database
    log_file: pathlib.Path
        Path to the log directory within the work directory
    number_of_threads: int
        Number of threads / processes
    """

    def __init__(self, termination_event: Event, nodes_dmp_path: pathlib.Path, names_dmp_path: pathlib.Path, merge_dmp_path: pathlib.Path, delete_dmp_path: pathlib.Path, database_url: str, log_file: pathlib.Path, number_of_threads: int):
        self.__termination_event = termination_event
        self.__nodes_dmp_path = nodes_dmp_path
        self.__names_dmp_path = names_dmp_path
        self.__merge_dmp_path = merge_dmp_path
        self.__delete_dmp_path = delete_dmp_path
        self.__database_url = database_url
        self.__log_file = log_file
        self.__number_of_threads = number_of_threads

    def maintain(self):
        """
        Updates the taxonomy tables.
        """
        # Initialize signal handler for TERM and INT
        signal.signal(signal.SIGTERM, self.__termination_event_handler)
        signal.signal(signal.SIGINT, self.__termination_event_handler)

        self.__build_taxonomies()
        self.__merge_taxonomies()
        self.__delete_taxonomies()

    def __termination_event_handler(self, signal_number, frame):
        """
        Sets the termination event to true. Callback for signal handling.

        Paremters
        ---------
        signal_number
            Signal ID
        frame
            Signal frame
        """
        self.__termination_event.set()

    def __start_work_and_logger_processes(self, work_function, logger_write_mode) -> Tuple[list, Process, Queue, Event]:
        """
        This will create working processes and one logger processes. The worker processes are connected with the logger process.
        Make sure the worker function has the following argument list: id (int), database_url (str), stop_flag (Event) and log_connection (Connection).
        
        Parameters
        ----------
        work_function : Callable[int, str, Event, Connection]
            Function which is executed by the worker processes
        logger_write_mode : str
            Write mode for the logger (see write mode for Pythons `open()`-function)

        Returns
        -------
        Tuple (array of worker processes, the logger process, the work queue, and a stop flag for the worker processes)
        """
        # Build processing context
        processing_context = get_processing_context("spawn")

        # Initialize queue with 3 times the thread number as limit to prevent RAM overflow.
        queue = processing_context.Queue(self.__number_of_threads * 3)
        # Initialize stop flag
        stop_flag = processing_context.Event()

        # Array to collect procs
        procs = []
        # Array to collect 
        log_connections = []

        for id in range(self.__number_of_threads):
            log_connection_read, log_connection_write = processing_context.Pipe(duplex=False)
            log_connections.append(log_connection_read)
            proc = processing_context.Process(target=work_function, args=(id, self.__database_url, stop_flag, queue, log_connection_write,))
            proc.start()
            # Close this handle of the writeable end to ensure only the process has one
            log_connection_write.close()
            procs.append(proc)

        logger = LoggerProcess(self.__termination_event, self.__log_file, logger_write_mode, log_connections)
        logger.start()

        return procs, logger, queue, stop_flag


    def __build_taxonomies(self):
        """
        Inserts nodes into the taxonomy tree.
        """
        if self.__nodes_dmp_path.exists() and self.__names_dmp_path.exists():
            if self.__nodes_dmp_path.is_dir():
                raise AssertionError(f"'{self.__nodes_dmp_path}' is a directory, not a file.")
            if self.__names_dmp_path.is_dir():
                raise AssertionError(f"'{self.__names_dmp_path}' is a directory, not a file.")

            procs, logger, queue, stop_flag = self.__start_work_and_logger_processes(self.process_new_taxonomies, "w")

            print("Build taxonomy tree ...")
            # Map taxonomy_id -> taxonomy
            taxonomies = {}
            # Open node file
            with self.__nodes_dmp_path.open("r") as nodes_file:
                for line in nodes_file:
                    # Parse line and create taxonomy without name
                    id, parent_id, rank = self.parse_node_line(line)
                    taxonomies[id] = Taxonomy(id, parent_id, "", rank)
            
            print("Read taxonomy names and insert into database ...")
            # Open name file
            with self.__names_dmp_path.open("r") as names_file:
                taxonomies_chunk = []
                for line in names_file:
                    # Parse line
                    id, name, name_class = self.parse_name_line(line)
                    # Check if name is scientific
                    if name_class == 'scientific name':
                        # Take remove taxonomy from dict
                        taxonomy = taxonomies.pop(id, None)
                        if taxonomy:
                            # Assigne name
                            taxonomy.name = name
                            taxonomies_chunk.append(taxonomy)
                    if len(taxonomies_chunk) == 1000:
                        # Put into queue
                        while True:
                            try:
                                queue.put(
                                    taxonomies_chunk,
                                    True,
                                    2
                                )
                                taxonomies_chunk = []
                                break
                            except QueueFullError:
                                continue
                if len(taxonomies_chunk) > 0:
                    # Put into queue
                    while True:
                        try:
                            queue.put(
                                taxonomies_chunk,
                                True,
                                2
                            )
                            taxonomies_chunk = []
                            break
                        except QueueFullError:
                            continue

            # Set stop flag
            stop_flag.set()

            # Wait for processing worker to stop
            for proc in procs:
                proc.join()

            # Wait for logger to stop
            logger.join()

            if len(taxonomies):
                with self.__log_file().open('a+') as log_file:
                    log_file.write(f"\n## No name was found for this taxonomies\n")
                    log_file.write(f"## id\t|\tparent_id\t|\trank\n")
                    for taxonomy in taxonomies.values():
                        log_file.write(f"{taxonomy.id}\t|\t{taxonomy.parent_id}\t|\t{str(taxonomy.rank)}\n")
                    log_file.write(f"####\n")
    
    def __merge_taxonomies(self):
        """
        Inserts merges into the taxonomy tree
        """
        if self.__merge_dmp_path.exists():
            if self.__merge_dmp_path.is_dir():
                raise AssertionError(f"'{self.__merge_dmp_path}' is a directory, not a file.")

            procs, logger, queue, stop_flag = self.__start_work_and_logger_processes(self.process_taxonomy_merge, "a+")

            print("Process taxonomy merges ...")
            # Open file with node merges.
            with self.__merge_dmp_path.open("r") as merge_file:
                taxonomy_merges_chunk = []
                # Read line line by line. Each line contains the old id and the new id
                for merge_line in merge_file:
                    source_id, target_id = self.parse_merge_line(merge_line)
                    taxonomy_merges_chunk.append(TaxonomyMerge(source_id, target_id))
                    if len(taxonomy_merges_chunk) == 1000:
                        # Try to put Taxonomy into processing queue until there is a slot for it
                        while True:
                            try:
                                queue.put(
                                    taxonomy_merges_chunk,
                                    True,
                                    2
                                )
                                taxonomy_merges_chunk = []
                                break
                            except QueueFullError:
                                continue
                if len(taxonomy_merges_chunk) > 0:
                    # Try to put Taxonomy into processing queue until there is a slot for it
                    while True:
                        try:
                            queue.put(
                                taxonomy_merges_chunk,
                                True,
                                2
                            )
                            taxonomy_merges_chunk = []
                            break
                        except QueueFullError:
                            continue

            # Set stop flag
            stop_flag.set()

            # Wait for processing worker to stop
            for proc in procs:
                proc.join()

            # Wait for logger to stop
            logger.join()

    def __delete_taxonomies(self):
        """
        Deletes taxonomy from the tree
        """
        if self.__delete_dmp_path.exists():
            if self.__delete_dmp_path.is_dir():
                raise AssertionError(f"'{self.__delete_dmp_path}' is a directory, not a file.")
            procs, logger, queue, stop_flag = self.__start_work_and_logger_processes(self.process_taxonomy_deletion, "a+")

            print("Process taxonomy deletions ...")
            # Open file
            with self.__delete_dmp_path.open("r") as delete_file:
                # Read line by line. Each line contains the id for removal
                for delete_line in delete_file:
                    taxonomy_id = self.parse_delete_line(delete_line)
                    # Try to put Taxonomy into processing queue until there is a slot for it
                    while True:
                        try:
                            queue.put(
                                taxonomy_id,
                                True,
                                2
                            )
                            break
                        except QueueFullError:
                            continue

            # Set stop flag
            stop_flag.set()

            # Wait for processing worker to stop
            for proc in procs:
                proc.join()

            # Wait for logger to stop
            logger.join()

    @staticmethod
    def process_new_taxonomies(id: int, database_url: str, stop_flag: Event, queue: Queue, log_connection: Connection):
        """
        Inserts new taxonomies into the tree

        Parameters
        ----------
        id : int
            Just a id for identifying the process in the logs
        database_url : str
            Database URL, e.g. postgres://username:password@host:port/database
        stop_flag : Event
            Multiprocessing event to stop the processes
        queue : Queue
            Multiprocessing queue for queuing the new nodes
        log_connection : Connection
            Multiprocessing connection to send logs to the logger
        """
        log_connection.send(f"insert worker {id} is online")
        database_connection = None

        while not stop_flag.is_set() or not queue.empty():
            if not database_connection or (database_connection and database_connection.closed != 0):
                database_connection = psycopg2.connect(database_url)
            try:
                # Take taxonomies from queue
                taxonomies = queue.get(True, 2)
                commit_errors = 0
                while True:
                    try:
                        with database_connection:
                            with database_connection.cursor() as database_cursor:
                                Taxonomy.bulk_insert(database_cursor, taxonomies)
                        break
                    except BaseException as error:
                        commit_errors += 1
                        # If there are tries left, sleep between 2 and 5 second and try again
                        if commit_errors < 3:
                            time.sleep(random.randint(0, 5))
                        else:
                            # Otherwise log the error and proceed with next taxonomy
                            commaseperated_taxonomy_ids = ", ".join([str(taxonomy.id) for taxonomy in taxonomies])
                            log_connection.send(f"taxonomy {commaseperated_taxonomy_ids} raises error:\n{error}\n")
                            break
            except QueueEmptyError:
                # Catch queue empty error (thrown on timeout)
                continue
        if database_connection and database_connection.closed == 0:
            database_connection.close()
        log_connection.send(f"insert worker {id} is stopping")
        # Close connection to logger
        log_connection.close()


    @staticmethod
    def process_taxonomy_merge(id: int, database_url: str, stop_flag: Event, queue: Queue, log_connection: Connection):
        """
        Processes taxonomy merges from a queue.

        Parameters
        ----------
        id : int
            Just a id for identifying the process in the logs
        database_url : str
            Database URL, e.g. postgres://username:password@host:port/database
        stop_flag : Event
            Multiprocessing event to stop the processes
        queue : Queue
            Multiprocessing queue for queuing the merges
        log_connection : Connection
            Multiprocessing connection to send logs to the logger
        """
        log_connection.send(f"merge worker {id} is online")
        database_connection = None
        while not stop_flag.is_set() or not queue.empty():
            if not database_connection or (database_connection and database_connection.closed != 0):
                database_connection = psycopg2.connect(database_url)
            try:
                # Take a taxonomy from queue
                taxonomy_merges = queue.get(True, 2)
                commit_errors = 0
                # Try to insert it into the database
                while True:
                    try:
                        with database_connection:
                            with database_connection.cursor() as database_cursor:
                                TaxonomyMerge.bulk_insert(database_cursor, taxonomy_merges)
                        break
                    except BaseException as error:
                        commit_errors += 1
                        # If there are tries left, sleep between 2 and 5 second and try again
                        if commit_errors < 3:
                            time.sleep(random.randint(2, 5))
                        else:
                            # Otherwise log the error and proceed with next taxonomy
                            commaseperated_taxonmy_merges = ", ".join([f"({taxonomy_merge.source_id},{taxonomy_merge.target_id})" for taxonomy_merge in taxonomy_merges])
                            log_connection.send(f"taxonomy merge {commaseperated_taxonmy_merges} raises error:\n{error}\n")
                            break
            except QueueEmptyError:
                # Catch queue empty error (thrown on timeout)
                continue
        if database_connection and database_connection.closed == 0:
            database_connection.close()
        log_connection.send(f"merge worker {id} is stopping")
        # Close connection to logger
        log_connection.close()

    @staticmethod
    def process_taxonomy_deletion(id: int, database_url: str, stop_flag: Event, queue: Queue, log_connection: Connection):
        """
        Processes taxonomy merges from a queue.

        Parameters
        ----------
        id : int
            Just a id for identifying the process in the logs
        database_url : str
            Database URL, e.g. postgres://username:password@host:port/database
        stop_flag : Event
            Multiprocessing event to stop the processes
        queue : Queue
            Multiprocessing queue for queuing the merges
        log_connection : Connection
            Multiprocessing connection to send logs to the logger
        """
        log_connection.send(f"deletion worker {id} is online")
        database_connection = None
        while not stop_flag.is_set() or not queue.empty():
            if not database_connection or (database_connection and database_connection.closed != 0):
                database_connection = psycopg2.connect(database_url)
            try:
                # Take a taxonomy id from queue
                taxonomy_id = queue.get(True, 2)
                with database_connection:
                    with database_connection.cursor() as database_cursor:
                        commit_errors = 0
                        # Try to insert it into the database
                        while True:
                            try:
                                database_cursor.execute("DELETE FROM taxonomies WHERE id = %s;", (taxonomy_id,))
                                database_cursor.execute("DELETE FROM taxonomy_merges WHERE source_id = %s OR target_id = %s;", (taxonomy_id, taxonomy_id))
                                database_connection.commit()
                                break
                            except BaseException as error:
                                database_connection.rollback()
                                commit_errors += 1
                                if commit_errors < 3:
                                    time.sleep(random.randint(2, 5))
                                else:
                                    log_connection.send(f" taxonomy deletion {taxonomy_id} raises error:\n{error}\n")
                                    break
            except QueueEmptyError:
                # Catch queue empty error (thrown on timeout)
                continue
        if database_connection and database_connection.closed == 0:
            database_connection.close()
        log_connection.send(f"deletion worker {id} is stopping")
        # Close connection to logger
        log_connection.close()


    @classmethod
    def parse_node_line(cls, line: str) -> tuple:
        """
        Parses node lines and returns id, parent id and rank
        
        Parameters
        ----------
        line : str
            Line from nodes.dmp

        Returns
        -------
        Tuple (id, parent id, rank)
        """
        node_attributes = cls.split_dmp_file_line(line)
        return int(node_attributes[0]), int(node_attributes[1]), TaxonomyRank.from_string(node_attributes[2])

    @classmethod
    def parse_name_line(cls, line: str) -> tuple:
        """
        Parses names line and returns id, name and name class

        Parameters
        ----------
        line : str
            Line from names.dmp

        Returns
        -------
        Tuple (id, name, name class)
        """
        node_attributes = cls.split_dmp_file_line(line)
        return  int(node_attributes[0]), node_attributes[1], node_attributes[3]

    @classmethod
    def parse_merge_line(cls, line: str) -> tuple:
        """
        Parses merge line and returns source id and target id
        
        Parameters
        ----------
        line : str
            Line from merged.dmp

        Returns
        -------
        Tuple (source_id, target_id )
        """
        node_attributes = cls.split_dmp_file_line(line)
        return int(node_attributes[0]), int(node_attributes[1])


    @classmethod
    def parse_delete_line(cls, line: str) -> int:
        """
        Parses delnodes line and returns taxonomy id

        Parameters
        ----------
        line : str
            Line from delnodes.dmp

        Returns
        -------
        Taxonomy id
        """
        return int(cls.split_dmp_file_line(line)[0])

    @classmethod
    def split_dmp_file_line(cls, line: str) -> list:
        """
        Splits line from dmp line. Separator is '|'
        Each element of the line will be stripped from prepending and appending whitespaces.
        
        Parameters
        ----------
        line : str
            Line from dmp-file

        Returns
        -------
        List
        """
        return [elem.strip() for elem in line.split("|")[:-1]]
