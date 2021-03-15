# Python modules
import pathlib
import time
import random
from multiprocessing import Process, Event, Queue, Pipe, get_context as get_processing_context
from multiprocessing.connection import Connection, wait
from queue import Full as QueueFullError, Empty as QueueEmptyError

# 3rd party modules
import psycopg2

# Internal imports
from .multiprocessing.logger_process import LoggerProcess
from ...models.taxonomy import Taxonomy, TaxonomyRank
from ...models.taxonomy_merge import TaxonomyMerge
from ...models.protein import Protein

class TaxonomyTree:

    def __init__(self, nodes_dmp_path: pathlib.Path, names_dmp_path: pathlib.Path, merge_dmp_path: pathlib.Path, delete_dmp_path: pathlib.Path, database_url: str, log_file: pathlib.Path, number_of_threads: int):
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
        self.__build_taxonomies()
        self.__merge_taxonomies()
        self.__delete_taxonomies()


    # This will create working processes and one logger processes. The worker processes are connected with the logger process.
    # Make sure the worker function has the following argument list: id (int), database_url (str), stop_flag (Event) and log_connection (Connection).
    # Returns (array of worker processes, the logger process, the work queue, and a stop flag for the worker processes)
    def __start_work_and_logger_processes(self, work_function, logger_write_mode) -> (list, Process, Queue, Event):
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

        logger = LoggerProcess(self.__log_file, logger_write_mode, processing_context)
        logger.start(log_connections)

        return procs, logger, queue, stop_flag


    def __build_taxonomies(self):
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
                            # Put into queue
                            while True:
                                try:
                                    queue.put(
                                        taxonomy,
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

            if len(taxonomies):
                with self.__log_file().open('a+') as log_file:
                    log_file.write(f"\n## No name was found for this taxonomies\n")
                    log_file.write(f"## id\t|\tparent_id\t|\trank\n")
                    for taxonomy in taxonomies.values():
                        log_file.write(f"{taxonomy.id}\t|\t{taxonomy.parent_id}\t|\t{str(taxonomy.rank)}\n")
                    log_file.write(f"####\n")
    
    def __merge_taxonomies(self):
        if self.__merge_dmp_path.exists():
            if self.__merge_dmp_path.is_dir():
                raise AssertionError(f"'{self.__merge_dmp_path}' is a directory, not a file.")

            procs, logger, queue, stop_flag = self.__start_work_and_logger_processes(self.process_taxonomy_merge, "a+")

            print("Process taxonomy merges ...")
            # Open file with node merges.
            with self.__merge_dmp_path.open("r") as merge_file:
                # Read line line by line. Each line contains the old id and the new id
                for merge_line in merge_file:
                    source_id, target_id = self.parse_merge_line(merge_line)
                    # Try to put Taxonomy into processing queue until there is a slot for it
                    while True:
                        try:
                            queue.put(
                                TaxonomyMerge(source_id, target_id),
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

    def __delete_taxonomies(self):
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

    # Inserts new taxonomies into the database
    @staticmethod
    def process_new_taxonomies(id: int, database_url: str, stop_flag: Event, queue: Queue, log_connection: Connection):
        log_connection.send(f"insert worker {id} is online")
        database_connection = None

        while not stop_flag.is_set() or not queue.empty():
            if not database_connection or (database_connection and database_connection.closed != 0):
                database_connection = psycopg2.connect(database_url)
            try:
                # Take a taxonomy from queue
                taxonomy = queue.get(True, 2)
                with database_connection:
                    with database_connection.cursor() as database_cursor:
                        database_cursor.execute("SELECT id FROM taxonomies WHERE id = %s;", (taxonomy.id,))
                        taxonomy_row = database_cursor.fetchone()
                        if not taxonomy_row:
                            commit_errors = 0
                            # Try to insert it into the database
                            while True:
                                try:
                                    database_cursor.execute("INSERT INTO taxonomies (id, parent_id, name, rank) VALUES (%s, %s, %s, %s);", (taxonomy.id, taxonomy.parent_id, taxonomy.name, taxonomy.rank.value))
                                    database_connection.commit()
                                    break
                                except BaseException as error:
                                    # On error, do a rollback
                                    database_connection.rollback()
                                    commit_errors += 1
                                    # If there are tries left, sleep between 2 and 5 second and try again
                                    if commit_errors < 3:
                                        time.sleep(random.randint(0, 5))
                                    else:
                                        # Otherwise log the error and proceed with next taxonomy
                                        log_connection.send(f" taxonomy {taxonomy.id}|{taxonomy.parent_id}|{taxonomy.name} raises error:\n{error}\n")
                                        break
            except QueueEmptyError:
                # Catch queue empty error (thrown on timeout)
                continue
        database_connection.close()
        log_connection.send(f"insert worker {id} is stopping")
        # Close connection to logger
        log_connection.close()


    @staticmethod
    def process_taxonomy_merge(id: int, database_url: str, stop_flag: Event, queue: Queue, log_connection: Connection):
        log_connection.send(f"merge worker {id} is online")
        database_connection = None
        while not stop_flag.is_set() or not queue.empty():
            if not database_connection or (database_connection and database_connection.closed != 0):
                database_connection = psycopg2.connect(database_url)
            try:
                # Take a taxonomy from queue
                taxonomy_merge = queue.get(True, 2)
                with database_connection:
                    with database_connection.cursor() as database_cursor:
                        database_cursor.execute("SELECT source_id, target_id FROM taxonomy_merges WHERE source_id = %s AND target_id = %s;", (taxonomy_merge.source_id, taxonomy_merge.target_id))
                        taxonomy_merge_row = database_cursor.fetchone()
                        if not taxonomy_merge_row:
                            commit_errors = 0
                            # Try to insert it into the database
                            while True:
                                try:
                                    database_cursor.execute("INSERT INTO taxonomy_merges (source_id, target_id) VALUES (%s, %s);", (taxonomy_merge.source_id, taxonomy_merge.target_id))
                                    database_cursor.execute("UPDATE proteins SET taxonomy_id = %s WHERE taxonomy_id = %s;", (taxonomy_merge.source_id, taxonomy_merge.target_id))
                                    database_cursor.execute("UPDATE peptides SET taxonomy_ids = (SELECT array_agg(taxonomy_id) FROM proteins WHERE id = ANY(protein_ids)) WHERE taxonomy_ids && %s;", ([taxonomy_merge.source_id],))
                                    database_connection.commit()
                                    break
                                except BaseException as error:
                                    # Do a rollback
                                    database_connection.rollback()
                                    commit_errors += 1
                                    # If there are tries left, sleep between 2 and 5 second and try again
                                    if commit_errors < 3:
                                        time.sleep(random.randint(2, 5))
                                    else:
                                        # Otherwise log the error and proceed with next taxonomy
                                        log_connection.send(f" taxonomy merge {taxonomy_merge.source_id}|{taxonomy_merge.target_id} raises error:\n{error}\n")
                                        break
            except QueueEmptyError:
                # Catch queue empty error (thrown on timeout)
                continue
        database_connection.close()
        log_connection.send(f"merge worker {id} is stopping")
        # Close connection to logger
        log_connection.close()

    @staticmethod
    def process_taxonomy_deletion(id: int, database_url: str, stop_flag: Event, queue: Queue, log_connection: Connection):
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
        database_connection.close()
        log_connection.send(f"deletion worker {id} is stopping")
        # Close connection to logger
        log_connection.close()


    @classmethod
    def parse_node_line(cls, line: str) -> tuple:
        """
        Parse node lines and returns id, parent id and rank
        @param line Line from nodes.dmp
        @return tupel (id, parent id, rank)
        """
        node_attributes = cls.split_dmp_file_line(line)
        return int(node_attributes[0]), int(node_attributes[1]), TaxonomyRank.from_string(node_attributes[2])

    @classmethod
    def parse_name_line(cls, line: str) -> tuple:
        """
        Parse names line and returns id, name and name class
        @param line Line from names.dmp
        @return tupel (id, name, name class)
        """
        node_attributes = cls.split_dmp_file_line(line)
        return  int(node_attributes[0]), node_attributes[1], node_attributes[3]

    @classmethod
    def parse_merge_line(cls, line: str) -> tuple:
        """
        Parse merge line and returns source id and target id
        @param line Line from merged.dmp
        @return tupel (source_id, target_id )
        """
        node_attributes = cls.split_dmp_file_line(line)
        return int(node_attributes[0]), int(node_attributes[1])


    @classmethod
    def parse_delete_line(cls, line: str) -> int:
        """
        Parse delnodes line and returns taxonomy id
        @param line Line from delnodes.dmp
        @return int Taxonomy id
        """
        return int(cls.split_dmp_file_line(line)[0])

    @classmethod
    def split_dmp_file_line(cls, line: str) -> list:
        """
        Splits line from dmp line. Separator is '|'
        Each element of the line will be stripped from prepending and appending whitespaces.
        @param line Line from dmp-file
        @return list
        """
        return [elem.strip() for elem in line.split("|")[:-1]]
