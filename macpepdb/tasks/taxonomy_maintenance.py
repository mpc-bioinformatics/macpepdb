import pathlib
import time
import random


from multiprocessing import Process, Event, Queue, Pipe, get_context as get_processing_context
from multiprocessing.connection import Connection, wait
from queue import Full as QueueFullError, Empty as QueueEmptyError

from sqlalchemy import create_engine, or_
from sqlalchemy.orm import sessionmaker

from ..models.taxonomy import Taxonomy, TaxonomyRank
from ..models.taxonomy_merge import TaxonomyMerge
from ..models.protein import Protein

class TaxonomyMaintenance():

    def __init__(self, nodes_dmp_path: pathlib.Path, names_dmp_path: pathlib.Path, merge_dmp_path: pathlib.Path, delete_dmp_path: pathlib.Path, database_url: str, log_file: pathlib.Path, number_of_threads: int):
        self.__nodes_dmp_path = nodes_dmp_path
        self.__names_dmp_path = names_dmp_path
        self.__merge_dmp_path = merge_dmp_path
        self.__delete_dmp_path = delete_dmp_path
        self.__database_url = database_url
        self.__log_file = log_file
        self.__number_of_threads = number_of_threads

    def start(self):
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

        logger = Process(target=self.logging_worker, args=(log_connections, self.__log_file, logger_write_mode))
        logger.start()

        return procs, logger, queue, stop_flag

        


    def __build_taxonomies(self):
        procs, logger, queue, stop_flag = self.__start_work_and_logger_processes(self.process_new_taxonomies, "w")

        # Open node file
        with self.__nodes_dmp_path.open("r") as nodes_file:
            # Open names file
            with self.__names_dmp_path.open("r") as names_file:
                name_attributes =  {
                    "id": -1,
                    "name": None,
                    "name_class": ""
                }

                # Read node file line by line, each line contains one node.
                # Nodes are sorted by id
                for node_line in nodes_file:
                    node_attributes = self.parse_node_line(node_line)
                    taxonomy_name = None
                    # Search names file for name of the current node.
                    # Name file is sorted by node id too.
                    while True:
                        # If ids matches, and name_class is "scientific name" name of node is found
                        if name_attributes["id"] == node_attributes["id"] and name_attributes["name_class"] == "scientific name":
                            taxonomy_name = name_attributes["name"]
                            break
                        # If id of name file is greater then node id, no name is found (insert into SQL will throw an error, due to name is None)
                        elif name_attributes["id"] > node_attributes["id"]:
                            break
                        try:
                            name_attributes = self.parse_name_line(next(names_file))
                        except StopIteration:
                            break

                    # Try to put Taxonomy into processing queue until there is a slot for it
                    while True:
                        try:
                            queue.put(
                                Taxonomy(node_attributes["id"], node_attributes["parent_id"], taxonomy_name, node_attributes["rank"]),
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
    
    def __merge_taxonomies(self):
        procs, logger, queue, stop_flag = self.__start_work_and_logger_processes(self.process_taxonomy_merge, "a+")

        # Open file with node merges.
        with self.__merge_dmp_path.open("r") as merge_file:
            # Read line line by line. Each line contains the old id and the new id
            for merge_line in merge_file:
                merge_attributes = self.parse_merge_line(merge_line)
                # Try to put Taxonomy into processing queue until there is a slot for it
                while True:
                    try:
                        queue.put(
                            TaxonomyMerge(merge_attributes["source_id"], merge_attributes["target_id"]),
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
        procs, logger, queue, stop_flag = self.__start_work_and_logger_processes(self.process_taxonomy_deletion, "a+")

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
        # Create database connection
        engine = create_engine(database_url, pool_size = 1, max_overflow = 0, pool_timeout = 3600)
        SessionClass = sessionmaker(bind = engine, autoflush=False)
        while not stop_flag.is_set() or not queue.empty():
            # Open new session
            session = SessionClass()
            try:
                # Take a taxonomy from queue
                taxonomy = queue.get(True, 2)
                # Check if the taxonomy already exists
                existing_taxonomy = session.query(Taxonomy).filter(Taxonomy.id == taxonomy.id).one_or_none()
                if not existing_taxonomy:
                    commit_errors = 0
                    # Try to insert it into the database
                    while True:
                        try:
                            session.add(taxonomy)
                            session.commit()
                            break
                        except BaseException as error:
                            # On error, do a rollback
                            session.rollback()
                            commit_errors += 1
                            # If there are tries left, sleep between 2 and 5 second and try again
                            if commit_errors < 3:
                                time.sleep(random.randint(2, 5))
                            else:
                                # Otherwise log the error and proceed with next taxonomy
                                log_connection.send(f" taxonomy {taxonomy.id}|{taxonomy.parent_id}|{taxonomy.name} raises error:\n{error}\n")
                                break
            except QueueEmptyError:
                # Catch queue empty error (thrown on timeout)
                continue
            # Close session
            session.close()
        # Dispose all connection within the engine
        engine.dispose()
        log_connection.send(f"insert worker {id} is stopping")
        # Close connection to logger
        log_connection.close()


    @staticmethod
    def process_taxonomy_merge(id: int, database_url: str, stop_flag: Event, queue: Queue, log_connection: Connection):
        log_connection.send(f"merge worker {id} is online")
        # Create database connection
        engine = create_engine(database_url, pool_size = 1, max_overflow = 0, pool_timeout = 3600)
        SessionClass = sessionmaker(bind = engine, autoflush=False)
        while not stop_flag.is_set() or not queue.empty():
            # Open new session
            session = SessionClass()
            try:
                # Take a taxonomy from queue
                taxonomy_merge = queue.get(True, 2)
                # Check if the merge already exists
                existing_taxonomy_merge = session.query(TaxonomyMerge).filter(TaxonomyMerge.source_id == taxonomy_merge.source_id, TaxonomyMerge.target_id == taxonomy_merge.target_id).one_or_none()
                if not existing_taxonomy_merge:
                    commit_errors = 0
                    # Try to insert it into the database
                    while True:
                        try:
                            # Insert merge
                            session.add(taxonomy_merge)
                            # Update all taxonomy id of proteins which have the old id.
                            session.query(Protein).filter(Protein.taxonomy_id == taxonomy_merge.source_id).update({Protein.taxonomy_id: taxonomy_merge.target_id})
                            session.commit()
                            break
                        except BaseException as error:
                            # Do a rollback
                            session.rollback()
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
            # Close session
            session.close()
        # Dispose all connection within the engine
        engine.dispose()
        log_connection.send(f"merge worker {id} is stopping")
        # Close connection to logger
        log_connection.close()

    @staticmethod
    def process_taxonomy_deletion(id: int, database_url: str, stop_flag: Event, queue: Queue, log_connection: Connection):
        log_connection.send(f"deletion worker {id} is online")
        # Create a database connection
        engine = create_engine(database_url, pool_size = 1, max_overflow = 0, pool_timeout = 3600)
        SessionClass = sessionmaker(bind = engine, autoflush=False)
        while not stop_flag.is_set() or not queue.empty():
            # Create a session
            session = SessionClass()
            try:
                # Take a taxonomy id from queue
                taxonomy_id = queue.get(True, 2)
                commit_errors = 0
                # Try to insert it into the database
                while True:
                    try:
                        # Delete the Taxonomy
                        session.query(Taxonomy).filter(Taxonomy.id == taxonomy_id).delete()
                        # Delete the merges which have the deleted taxonomy as source or target
                        session.query(TaxonomyMerge).filter(or_(TaxonomyMerge.source_id == taxonomy_id, TaxonomyMerge.target_id == taxonomy_id)).delete()
                        session.commit()
                        break
                    except BaseException as error:
                        session.rollback()
                        commit_errors += 1
                        if commit_errors < 3:
                            time.sleep(random.randint(2, 5))
                        else:
                            log_connection.send(f" taxonomy deletion {taxonomy_id} raises error:\n{error}\n")
                            break
            except QueueEmptyError:
                # Catch queue empty error (thrown on timeout)
                continue
            # Close session
            session.close()
        # Dispose all connection within the engine
        engine.dispose()
        log_connection.send(f"deletion worker {id} is stopping")
        # Close connection to logger
        log_connection.close()


    @classmethod
    def parse_node_line(cls, line: str):
        node_attributes = cls.split_dmp_file_line(line)
        return {
            "id": int(node_attributes[0]),
            "parent_id": int(node_attributes[1]),
            "rank": TaxonomyRank.from_string(node_attributes[2])
        }

    @classmethod
    def parse_name_line(cls, line: str):
        node_attributes = cls.split_dmp_file_line(line)
        return {
            "id": int(node_attributes[0]),
            "name": node_attributes[1],
            "name_class": node_attributes[3]
        }

    @classmethod
    def parse_merge_line(cls, line: str):
        node_attributes = cls.split_dmp_file_line(line)
        return {
            "source_id": int(node_attributes[0]),
            "target_id": int(node_attributes[1])
        }

    @classmethod
    def parse_delete_line(cls, line: str):
        return int(cls.split_dmp_file_line(line)[0])

    @classmethod
    def split_dmp_file_line(cls, line: str):
        return [elem.strip() for elem in line.split("|")[:-1]]

    @staticmethod
    def logging_worker(process_connections: list, log_file_path: pathlib.Path, write_mode: str):
        with log_file_path.open(write_mode) as log_file:
            log_file.write("error logger is online\n")
            log_file.flush()
            while process_connections:
                for conn in wait(process_connections):
                    try:
                        message = conn.recv()
                    except EOFError:
                        process_connections.remove(conn)
                    else:
                        log_file.write(f"{message}\n")
                        log_file.flush()
            log_file.write("error logger is stopping\n")
            # will be flushed on file close

    @classmethod
    def maintain_from_command_line_interface(cls, args):
        maintenance = cls(
            pathlib.Path(args.nodes_dmp_path),
            pathlib.Path(args.names_dmp_path),
            pathlib.Path(args.merge_dmp_path),
            pathlib.Path(args.delete_dmp_path),
            args.database_url,
            pathlib.Path(args.log_file),
            args.threads
        )
        maintenance.start()


    @classmethod
    def comand_line_arguments(cls, subparsers):
        parser = subparsers.add_parser('taxonomy-maintenance', help="Builds and updates the taxonomy tree.")
        parser.add_argument("--nodes-dmp-path", "-t", type=str, required=True, help="Dump file with taxonomies")
        parser.add_argument("--names-dmp-path", "-n", type=str, required=True, help="Dump file with taxonomy names")
        parser.add_argument("--merge-dmp-path", "-m", type=str, required=True, help="Dump file with merged taxonomies")
        parser.add_argument("--delete-dmp-path", "-r", type=str, required=True, help="Dump file with removed taxonomies")
        parser.add_argument("--database-url", "-d", type=str, required=True, help="Database url")
        parser.add_argument("--log-file", "-l", type=str, required=True, help="Log file")
        parser.add_argument("--threads", type=int, required=True, help="Number of threads")
        parser.set_defaults(func=cls.maintain_from_command_line_interface)