import pathlib
import time
import sys
import re
import sys
import csv
import random
import signal
import os
import datetime

from multiprocessing import Process, Value, Queue, Array, Event, Pipe, get_context as get_process_context
from multiprocessing.connection import Connection, wait
from ctypes import c_bool, c_ulonglong
from queue import Empty, Full
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.orm.attributes import flag_modified
from sqlalchemy.sql import exists
from sqlalchemy.exc import DBAPIError
from sqlalchemy.orm.exc import NoResultFound

from ..proteomics.enzymes.digest_enzyme import DigestEnzyme
from ..proteomics.file_reader.uniprot_text_reader import UniprotTextReader
from ..models.protein import Protein
from ..models.peptide import Peptide
from ..models.protein_merge import ProteinMerge
from ..models.maintenance_information import MaintenanceInformation

class Digestion:
    STATSTIC_FILE_HEADER = ["seconds", "inserted_proteins", "inserted_peptides", "unsolvable errors", "protein_insert_rate", "peptide_insert_rate", "error_rate"]
    STATSTIC_PRINT_MIN_COLUMN_WIDTH = 12
    

    def __init__(self, input_file_paths: list, log_file_path: pathlib.Path, unprocessible_proteins_file_path: pathlib.Path, statistics_csv_file_path: pathlib.Path, statistics_write_period: int, thread_count: int, enzyme_name: str, maximum_number_of_missed_cleavages: int, minimum_peptide_length: int, maximum_peptide_length: int):
        self.__log_file_path = log_file_path
        self.__unprocessible_proteins_file_path = unprocessible_proteins_file_path
        self.__thread_count = thread_count
        self.__maximum_number_of_missed_cleavages = maximum_number_of_missed_cleavages
        self.__minimum_peptide_length = minimum_peptide_length
        self.__maximum_peptide_length = maximum_peptide_length
        self.__enzyme_name = enzyme_name
        EnzymeClass = DigestEnzyme.get_enzyme_by_name(enzyme_name)
        self.__enzyme = EnzymeClass(self.__maximum_number_of_missed_cleavages, self.__minimum_peptide_length, self.__maximum_peptide_length)
        self.__input_file_paths = input_file_paths
        self.__statistics_write_period = statistics_write_period
        self.__statistics_csv_file_path = statistics_csv_file_path
        self.__stop_signal = False

    def digest_to_database(self, database_url: str):
        self.__set_databse_status(database_url, True)
        self.__load_or_set_digestion_informations(database_url)

        # Very important 
        process_context = get_process_context("spawn")

        # Events to control processes
        empty_queue_and_stop_flag = process_context.Event()
        logger_stop_flag = process_context.Event()
        stop_immediately_flag = process_context.Event()
        # Providing a maximum size prevents overflowing RAM and makes sure every process has enough work to do.
        protein_queue = process_context.Queue(3 * self.__thread_count)
        # Array for statistics [created_proteins, failed_proteins, created_peptides, processed_peptides]
        statistics = process_context.Array(c_ulonglong, 3)

        # Initialize signal handler for TERM and INT
        signal.signal(signal.SIGTERM, self.__stop_signal_handler)
        signal.signal(signal.SIGINT, self.__stop_signal_handler)

        unprocessible_channels = []
        log_connections = []

        # Log connection for the main process
        log_connection_read, log_connection = process_context.Pipe(duplex=False)
        log_connections.append(log_connection_read)

        print(f"to stop the digestion gracefully, send TERM or INT signal to process {os.getpid()}")
        log_connection.send(f"process id = {os.getpid()}")

        # Start digest worker
        digest_processes = []
        for worker_id in range(0, self.__thread_count):
            unprocessible_connection_read, unprocessible_connection_write = process_context.Pipe(duplex=False)
            unprocessible_channels.append(unprocessible_connection_read)
            log_connection_read, log_connection_write = process_context.Pipe(duplex=False)
            log_connections.append(log_connection_read)
            digest_worker = process_context.Process(target=self.digest_worker, args=(worker_id, empty_queue_and_stop_flag, stop_immediately_flag, protein_queue, log_connection_write, unprocessible_connection_write, statistics, self.__enzyme, database_url,))
            digest_worker.start()
            digest_processes.append(digest_worker)
            # Close these copies of the writeable site of the channels, so only the processes have working copies.
            unprocessible_connection_write.close()
            log_connection_write.close()

        # Start logger for unprocessible proteins
        log_connection_read, log_connection_write = process_context.Pipe(duplex=False)
        unprocessible_proteins_process = process_context.Process(target=self.unprocessible_proteins_worker, args=(unprocessible_channels, self.__unprocessible_proteins_file_path, log_connection_write,))
        unprocessible_proteins_process.start()
        log_connections.append(log_connection_read)
        # Close this copy
        log_connection_write.close()

        # Statistic writer
        log_connection_read, log_connection_write = process_context.Pipe(duplex=False)
        statistics_logger_process = process_context.Process(target=self.statistics_worker, args=(logger_stop_flag, statistics, self.__statistics_write_period, self.__statistics_csv_file_path, log_connection_write,))
        statistics_logger_process.start()
        log_connections.append(log_connection_read)
        # Close this copy
        log_connection_write.close()

        logger_process = process_context.Process(target=self.logging_worker, args=(log_connections, self.__log_file_path,))
        logger_process.start()

        for input_file_path in self.__input_file_paths:
            # Break file loop
            if self.__stop_signal:
                break
            with input_file_path.open("r") as input_file:
                input_file_reader = UniprotTextReader(input_file)

                # Enqueue proteins
                for protein, protein_merges in input_file_reader:
                    # Break protein read loop
                    if self.__stop_signal:
                        break
                    # Retry lopp
                    while True:
                        # Break retry loop
                        if self.__stop_signal:
                            break
                        try:
                            # Try to queue a protein for digestion, wait 5 seconds for a free slot
                            protein_queue.put((protein, protein_merges), True, 5)
                            # Continue with the next protein from the file
                            break
                        # Catch timouts of protein_queue.put()
                        except Full:
                            # Try again
                            pass

        # Signalling all processes to 
        if self.__stop_signal:
            stop_immediately_flag.set()
        empty_queue_and_stop_flag.set()


        # Wait for all digest workers to stop
        for digest_worker in digest_processes:
            digest_worker.join()

        # It is possible that the last statistic is not entirely transfered when the logger stops. So wait a second to make sure all information are transfered.
        time.sleep(1)

        # Stop logger after each digest worker is finished
        logger_stop_flag.set()
        
        # Wait for the logger processes to stop
        statistics_logger_process.join()
        unprocessible_proteins_process.join()

        log_connection.close()
        logger_process.join()

        now = datetime.datetime.utcnow()
        epoch = datetime.datetime(1970,1,1)
        last_update_timestamp = (now - epoch).total_seconds()
        self.__set_databse_status(database_url, False, last_update_timestamp)


    def __load_or_set_digestion_informations(self, database_url: str):
        """
        Loads the digestion information from previous digestions from the databsae and overrides the given ones.
        If no digestion information were found, it will save the current one.
        @param database_url Database URL
        """
        engine = create_engine(database_url, pool_size = 1, max_overflow = 0, pool_timeout = 3600)
        session_factory = sessionmaker(bind = engine, autoflush=False)
        session = session_factory()
        digestion_information = session.query(MaintenanceInformation).filter(MaintenanceInformation.key == MaintenanceInformation.DIGESTION_PARAMTERS_KEY).one_or_none()
        if digestion_information:
            self.__maximum_number_of_missed_cleavages = digestion_information.values['maximum_number_of_missed_cleavages']
            self.__minimum_peptide_length = digestion_information.values['minimum_peptide_length']
            self.__maximum_peptide_length = digestion_information.values['maximum_peptide_length']
            self.__enzyme_name = digestion_information.values['enzyme_name']
            EnzymeClass = DigestEnzyme.get_enzyme_by_name(self.__enzyme_name)
            self.__enzyme = EnzymeClass(self.__maximum_number_of_missed_cleavages, self.__minimum_peptide_length, self.__maximum_peptide_length)
        else:
            digestion_information_values = {
                'enzyme_name': self.__enzyme_name,
                'maximum_number_of_missed_cleavages': self.__maximum_number_of_missed_cleavages,
                'minimum_peptide_length': self.__minimum_peptide_length,
                'maximum_peptide_length': self.__maximum_peptide_length
            }
            digestion_information = MaintenanceInformation(MaintenanceInformation.DIGESTION_PARAMTERS_KEY, digestion_information_values)
            session.add(digestion_information)
            session.commit()
        session.close()

    def __set_databse_status(self, database_url: str, is_maintenance_mode: bool, last_update_timestamp: int = None):
        """
        Sets the database status
        @param is_maintenance_mode Indicate if set to maintenance mode
        @param last_update_timestamp UTC timestamp in seconds
        """
        engine = create_engine(database_url, pool_size = 1, max_overflow = 0, pool_timeout = 3600)
        session_factory = sessionmaker(bind = engine, autoflush=False)
        session = session_factory()
        database_status = session.query(MaintenanceInformation).filter(MaintenanceInformation.key == MaintenanceInformation.DATABASE_STATUS_KEY).one_or_none()
        if database_status:
            database_status.values['maintenance_mode'] = is_maintenance_mode
            if last_update_timestamp:
                database_status.values['last_update'] = last_update_timestamp
            # SQL-Achemy does not register changes to a JSON-column, unless we flag it as modified
            flag_modified(database_status, 'values')
        else:
            digestion_information_values = {
                'maintenance_mode': is_maintenance_mode
            }
            if last_update_timestamp:
                digestion_information_values['last_update'] = last_update_timestamp
            else:
                digestion_information_values['last_update'] = 0
            digestion_information = MaintenanceInformation(MaintenanceInformation.DATABASE_STATUS_KEY, digestion_information_values)
            session.add(digestion_information)
        session.commit()
        session.close()


    @staticmethod
    def digest_worker(id: int, empty_queue_and_stop_flag: Event, stop_immediately_flag: Event, protein_queue: Queue, log_connection: Connection, unprocessible_log_connection: Connection, statistics: Array, enzyme: DigestEnzyme, database_url: str):
        log_connection.send("digest worker {} is online".format(id))
        engine = create_engine(database_url, pool_size = 1, max_overflow = 0, pool_timeout = 3600)
        SessionClass = sessionmaker(bind = engine, autoflush=False)
        # Let the process run until empty_queue_and_stop_flag is true and protein_queue is empty or stop_immediately_flag is true.
        while (not empty_queue_and_stop_flag.is_set() or not protein_queue.empty()) and not stop_immediately_flag.is_set():
            try:
                # Try to get a protein from the queue, timeout is 2 seconds
                protein, protein_merges = protein_queue.get(True, 5)
                
                # Variables for loop control
                unsolvable_errors = 0
                try_transaction_again = True
                while try_transaction_again:
                    session = SessionClass()
                    number_of_new_peptides = 0
                    try:
                        count_protein = False
                        # Check if the Protein exists by its accession or secondary accessions
                        accessions = [protein.accession] + [protein_merge.source_accession for protein_merge in protein_merges]
                        existing_protein = session.query(Protein).filter(Protein.accession.in_(accessions)).first()
                        if existing_protein:
                            need_commit, number_of_new_peptides = Digestion.update_protein(session, protein, existing_protein, protein_merges, enzyme)
                            if need_commit:
                                session.commit()
                            else:
                                session.rollback()
                        else:
                            number_of_new_peptides = Digestion.create_protein(session, protein, protein_merges, enzyme)
                            session.commit()
                            count_protein = True

                        # Commit was successfully stop while-loop and add statistics
                        try_transaction_again = False
                        if count_protein:
                            statistics[0] += 1
                        statistics[1] += number_of_new_peptides
                    except:
                        # Rollback session
                        session.rollback()
                        # Remove all peptides from protein
                        protein.peptides = []
                        # Try again after 5 (first try) and 10 (second try) + a random number between 0 and 5 (both including) seconds maybe some blocking transactions can pass so this transaction will successfully finish on the next try.
                        # If this is the third time an unsolvable error occures give up and log the error.
                        if unsolvable_errors < 2:
                            unsolvable_errors += 1
                            time.sleep(5 * unsolvable_errors + random.randint(0, 5))
                        # Log the error on the third try and put the protein in unprocessible queue
                        else:
                            exception = sys.exc_info()[1]
                            log_connection.send("Exception on protein {}, see:\n{}".format(protein.accession, exception))
                            unprocessible_log_connection.send(protein.to_embl_entry(protein_merges))
                            statistics[2] += 1
                            try_transaction_again = False
                    finally:
                        session.close()
            # Catch queue.Empty which is thrown when protein_queue.get() timed out
            except Empty:
                pass
        log_connection.send("digest worker {} is stopping".format(id))
        log_connection.close()
        unprocessible_log_connection.close()

    @staticmethod
    def create_protein(session, protein: Protein, protein_merges: list, enzyme: DigestEnzyme) -> int:
        """
        Creates a new protein. Does not check if the protein already exists.
        @param session Open database session, with open transaction.
        @param protein Protein to digest
        @param protein_merges List of protein merges for the given protein
        @param enzyme Digest enzym
        @return int Number of newly inserted peptides
        """
        session.add(protein)
        peptides = enzyme.digest(protein)

        for merge in protein_merges:
            session.add(merge)

        existing_peptide_query = session.query(Peptide).filter(Peptide.sequence.in_([peptide.sequence for peptide in peptides]))
        existing_peptides_dict = {peptide.sequence: peptide for peptide in existing_peptide_query.all()}
        for peptide in peptides:
            if peptide.sequence in existing_peptides_dict:
                protein.peptides.append(existing_peptides_dict[peptide.sequence])
            else:
                protein.peptides.append(peptide)

        return len(peptides) - len(existing_peptides_dict)

    @staticmethod
    def update_protein(session, updated_protein: Protein, existing_protein: Protein, protein_merges: list, enzyme: DigestEnzyme) -> (bool, int):
        """
        Creates a new protein. Does not check if the protein already exists.
        @param session Open database session
        @param updated_protein Protein with updates from file
        @param existing_protein Protein from database
        @param protein_merges List of protein merges for the updated_protein
        @param enzyme Digest enzym
        @return tuple (bool, int) First element indicates if the session needs to commit, seconds is the number of newly inserted peptides
        """
        needs_commit = False
        new_peptides_count = 0
        if updated_protein.accession != existing_protein.accession:
            # Replace the target accession with the new protein accession
            session.execute(ProteinMerge.__table__.update().where(ProteinMerge.target_accession == existing_protein.accession).values(target_accession = updated_protein.accession))
            # Create all missing protein merges
            for protein_merge in protein_merges:
                protein_merge_query = session.query(ProteinMerge).filter(ProteinMerge.source_accession == protein_merge.source_accession, ProteinMerge.target_accession == protein_merge.target_accession)
                if not protein_merge_query.one_or_none():
                    session.add(protein_merge)
            existing_protein.accession = updated_protein.accession
            needs_commit = True
        if updated_protein.taxonomy_id != existing_protein.taxonomy_id:
            existing_protein.taxonomy_id = updated_protein.taxonomy_id
            needs_commit = True
        if updated_protein.proteome_id != existing_protein.proteome_id:
            existing_protein.proteome_id = updated_protein.proteome_id
            needs_commit = True
        if updated_protein.sequence != existing_protein.sequence:
            existing_protein.sequence = updated_protein.sequence
            currently_referenced_peptides = existing_protein.peptides.all()
            new_peptide_map = {peptide.sequence: peptide for peptide in enzyme.digest(existing_protein)}

            # Check if the referenced peptides are in the new peptide map
            for peptide in currently_referenced_peptides:
                if not peptide.sequence in new_peptide_map:
                    # If not, delete reference, because the peptide is not longer in the sequence
                    existing_protein.peptides.remove(peptide)
                else:
                    # If contained remove them from the new peptides
                    new_peptide_map.pop(peptide.sequence)


            existing_peptide_query = session.query(Peptide).filter(Peptide.sequence.in_([peptide.sequence for peptide in new_peptide_map.values()]))
            existing_peptides_dict = {peptide.sequence: peptide for peptide in existing_peptide_query.all()}
            # Iterate over the new peptides 
            for peptide in new_peptide_map.values():
                if peptide.sequence in existing_peptides_dict:
                    # If the peptide exists, create a reference
                    existing_protein.peptides.append(existing_peptides_dict[peptide.sequence])
                else:
                    # If not create it and add reference
                    existing_protein.peptides.append(peptide)
            new_peptides_count = len(new_peptide_map) - len(existing_peptides_dict)
            needs_commit = True

        return needs_commit, new_peptides_count


    @staticmethod
    def logging_worker(process_connections: list, log_file_path: pathlib.Path):
        with log_file_path.open("w") as log_file:
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
            log_file.write("error logger is stopping")
            # will be flushed on file close
    
    @staticmethod
    def unprocessible_proteins_worker(process_connections: list, unprocessible_proteins_file_path: pathlib.Path, log_connection: Connection):
        log_connection.send("unprocessible proteins logger is online")
        with unprocessible_proteins_file_path.open("w") as unprocessible_proteins_file:
            while process_connections:
                for conn in wait(process_connections):
                    try:
                        embl_entry = conn.recv()
                    except EOFError:
                        process_connections.remove(conn)
                    else:
                        unprocessible_proteins_file.write(embl_entry + "\n")
                        unprocessible_proteins_file.flush()
        log_connection.send("unprocessible proteins logger is stopping")
        log_connection.close()

    @staticmethod
    def statistics_worker(stop_flag: Event, statistics: Array, write_period: int, statistics_file_path: pathlib.Path, log_connection: Connection):
        log_connection.send("statistics logger is online")
        # Snapshot of last written statistic to calculate difference
        last_statistics = [0] * len(statistics)
        start_at = time.time()
        Digestion.print_statistic_row(Digestion.STATSTIC_FILE_HEADER)
        with statistics_file_path.open("w") as statistics_file:
            statistics_writer = csv.writer(statistics_file)
            statistics_writer.writerow(Digestion.STATSTIC_FILE_HEADER)
            statistics_file.flush()
            while not stop_flag.is_set():
                # Prepare array for next write
                current_statistics = []
                # Wait for next write
                stop_flag.wait(write_period)
                # Calculate seconds after start
                current_time = int(time.time() - start_at)
                # Get current statistics
                for value in statistics:
                    current_statistics.append(value)
                # Initialize csv row
                csv_row = [current_time]
                # Assign current statistics to csv row
                for value in current_statistics:
                    csv_row.append(value)
                # Calulate differences to last statistics (= rates)
                for idx in range(0, len(current_statistics)):
                    csv_row.append(current_statistics[idx] - last_statistics[idx])
                # Write csv row to csv file
                statistics_writer.writerow(csv_row)
                statistics_file.flush()
                # Output to console
                Digestion.print_statistic_row(csv_row)
                # Assign new 'snapshot'
                last_statistics = current_statistics
        log_connection.send("statistics logger is offline")
        log_connection.close()
    
    @staticmethod
    def print_statistic_row(row: list):
        output = []
        for idx in range(0, len(Digestion.STATSTIC_FILE_HEADER)):
            column_width = max(
                Digestion.STATSTIC_PRINT_MIN_COLUMN_WIDTH,
                len(Digestion.STATSTIC_FILE_HEADER[idx])
            )
            output.append("{value:>{column_width}}".format(column_width=column_width, value=row[idx]))
        print("\t".join(output))

    def __stop_signal_handler(self, signal_number, frame):
        self.__stop_signal = True


    @classmethod
    def digest_from_command_line(cls, args):
        digestion = cls(
            [pathlib.Path(input_file_path) for input_file_path in args.input_file_paths],
            pathlib.Path(args.log_file_path),
            pathlib.Path(args.unprocessible_proteins_file_path),
            pathlib.Path(args.statistics_csv_file_path),
            args.statistics_write_period,
            args.thread_count,
            args.enzyme_name,
            args.maximum_number_of_missed_cleavages,
            args.minimum_peptide_length,
            args.maximum_peptide_length
        )
        digestion.digest_to_database(args.database_url)

    @classmethod
    def comand_line_arguments(cls, subparsers):
        parser = subparsers.add_parser('digestion', help="Digest proteins from a file and store the resulting peptides in a PostgreSQL-database.")
        parser.add_argument("--input-file-paths", "-i", required=True, help="Protein file (UniProt text format. Parameter can be used multiple times)", action='append')
        parser.add_argument("--log-file-path", "-l", type=str, required=True, help="Log file")
        parser.add_argument("--unprocessible-proteins-file-path", "-u", type=str, required=True, help="File for unprocessible proteins")
        parser.add_argument("--statistics-csv-file-path", "-s", type=str, required=True, help="File for digest statistics")
        parser.add_argument("--statistics-write-period", type=int, help="Seconds between writes to the statistics file (default: 900)", default=900)
        parser.add_argument("--thread-count", "-t", type=int, required=True, help="Number of concurrent digestions and inserts to the database")
        parser.add_argument("--enzyme-name", "-e", type=str, required=True, help="The name of the enzyme", choices=DigestEnzyme.get_known_enzymes())
        parser.add_argument("--maximum-number-of-missed-cleavages", "-c", type=int, help="Maximum number of missed cleavages (default: 2)", default="2")
        parser.add_argument("--minimum-peptide-length", type=int, help="Minimum peptide length (default: 5)", default="5")
        parser.add_argument("--maximum-peptide-length", type=int, help="Maximum peptide length (default: 60)", default="60")
        parser.add_argument("--database-url", "-d", type=str, required=True, help="Database url for postgres")
        parser.set_defaults(func=cls.digest_from_command_line)


