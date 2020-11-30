import pathlib
import time
import sys
import re
import sys
import csv
import random
import signal
import os

from multiprocessing import Process, Value, Queue, Array, Event, Pipe, get_context as get_process_context
from multiprocessing.connection import Connection, wait
from ctypes import c_bool, c_ulonglong
from queue import Empty, Full
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.sql import exists
from sqlalchemy.exc import DBAPIError
from sqlalchemy.orm.exc import NoResultFound

from ..proteomics.enzymes.digest_enzyme import DigestEnzyme
from ..proteomics.file_reader.file_reader import FileReader
from ..models.protein import Protein
from ..models.peptide import Peptide
from ..models.protein_merge import ProteinMerge

class Digestion:
    STATSTIC_FILE_HEADER = ["seconds", "inserted_proteins", "inserted_peptides", "unsolvable errors", "protein_insert_rate", "peptide_insert_rate", "error_rate"]
    STATSTIC_PRINT_MIN_COLUMN_WIDTH = 12
    

    def __init__(self, input_file_paths: list, input_format: str, log_file_path: pathlib.Path, unprocessible_proteins_fasta_file_path: pathlib.Path, statistics_csv_file_path: pathlib.Path, statistics_write_period: int, thread_count: int, enzyme_name: str, maximum_number_of_missed_cleavages: int, minimum_peptide_length: int, maximum_peptide_length: int, start_with_protein: int):
        self.__log_file_path = log_file_path
        self.__unprocessible_proteins_fasta_file_path = unprocessible_proteins_fasta_file_path
        self.__thread_count = thread_count
        EnzymeClass = DigestEnzyme.get_enzyme_by_name(enzyme_name)
        self.__enzyme = EnzymeClass(maximum_number_of_missed_cleavages, minimum_peptide_length, maximum_peptide_length)
        self.__file_reader_class = FileReader.get_reader_by_file_format(input_format)
        self.__input_file_paths = input_file_paths
        self.__statistics_write_period = statistics_write_period
        self.__statistics_csv_file_path = statistics_csv_file_path
        self.__stop_signal = False
        self.__start_with_protein = start_with_protein

    def digest_to_database(self, database_url: str):
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
        unprocessible_proteins_process = process_context.Process(target=self.unprocessible_proteins_worker, args=(unprocessible_channels, self.__unprocessible_proteins_fasta_file_path, log_connection_write,))
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

        current_input_file = ""
        # Resets which each input file
        read_proteins_of_current_file = 0
        queued_proteins_of_current_file = 0

        for input_file_path in self.__input_file_paths:
            read_proteins_of_current_file = 0
            current_input_file = str(input_file_path.resolve())
            # Break file loop
            if self.__stop_signal:
                break
            with input_file_path.open("r") as input_file:
                input_file_reader = self.__file_reader_class(input_file)

                # Enqueue proteins
                for protein, protein_merges in input_file_reader:
                    read_proteins_of_current_file += 1
                    # Break protein read loop
                    if self.__stop_signal:
                        break
                    if read_proteins_of_current_file < self.__start_with_protein:
                        continue
                    # Retry lopp
                    while True:
                        # Break retry loop
                        if self.__stop_signal:
                            break
                        try:
                            # Try to queue a protein for digestion, wait 5 seconds for a free slot
                            protein_queue.put((protein, protein_merges), True, 5)
                            # Protein is successfully read, if it is queued
                            queued_proteins_of_current_file += 1
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

        # Report the protein number
        if self.__stop_signal:
            log_connection.send(f"receive stop signal at protein {queued_proteins_of_current_file - protein_queue.qsize()} in file {current_input_file}")

        # It is possible that the last statistic is not entirely transfered when the logger stops. So wait a second to make sure all information are transfered.
        time.sleep(1)

        # Stop logger after each digest worker is finished
        logger_stop_flag.set()
        
        # Wait for the logger processes to stop
        statistics_logger_process.join()
        unprocessible_proteins_process.join()

        log_connection.close()
        logger_process.join()



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
                # Create peptides
                peptides = enzyme.digest(protein)

                # Variables for loop control
                unsolvable_errors = 0
                try_transaction_again = True
                while try_transaction_again:
                    session = SessionClass()
                    try:
                        count_protein = False
                        # Check if the Protein exists
                        try:
                            protein = session.query(Protein).filter(Protein.accession == protein.accession).one()
                        except NoResultFound: 
                            # Add the new protein to the session
                            session.add(protein)
                            count_protein = True

                        existing_protein_merge_source_accessions = [row[0] for row in session.query(ProteinMerge.source_accession).filter(ProteinMerge.target_accession == protein.accession).all()]
                        for merge in protein_merges:
                            if not merge.source_accession in existing_protein_merge_source_accessions:
                                session.add(merge)


                        existing_peptides = session.query(Peptide).filter(Peptide.sequence.in_([peptide.sequence for peptide in peptides])).all()
                        existing_peptides_dict = {}
                        for peptide in existing_peptides:
                            existing_peptides_dict[peptide.sequence] = peptide
                        for peptide in peptides:
                            if peptide.sequence in existing_peptides_dict:
                                protein.peptides.append(existing_peptides_dict[peptide.sequence])
                            else:
                                protein.peptides.append(peptide)

                        session.commit()

                        # Commit was successfully stop while-loop and add statistics
                        try_transaction_again = False
                        if count_protein:
                            statistics[0] += 1
                        statistics[1] += len(peptides) - len(existing_peptides)
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
                            unprocessible_log_connection.send(protein.to_fasta_entry(protein_merges))
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
    def unprocessible_proteins_worker(process_connections: list, unprocessible_proteins_fasta_path: pathlib.Path, log_connection: Connection):
        log_connection.send("unprocessible proteins logger is online")
        with unprocessible_proteins_fasta_path.open("w") as fasta_file:
            while process_connections:
                for conn in wait(process_connections):
                    try:
                        fasta_entry = conn.recv()
                    except EOFError:
                        process_connections.remove(conn)
                    else:
                        fasta_file.write(fasta_entry + "\n")
                        fasta_file.flush()
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
            args.input_format,
            pathlib.Path(args.log_file_path),
            pathlib.Path(args.unprocessible_proteins_fasta_file_path),
            pathlib.Path(args.statistics_csv_file_path),
            args.statistics_write_period,
            args.thread_count,
            args.enzyme_name,
            args.maximum_number_of_missed_cleavages,
            args.minimum_peptide_length,
            args.maximum_peptide_length,
            args.start_with_protein
        )
        digestion.digest_to_database(args.database_url)

    @classmethod
    def comand_line_arguments(cls, subparsers):
        parser = subparsers.add_parser('digestion', help="Digest proteins from a file and store the resulting peptides in a PostgreSQL-database.")
        parser.add_argument("--input-file-paths", "-i", required=True, help="Protein file (can be used multiple times)", action='append')
        parser.add_argument("--input-format",  "-f", type=str, required=True, help="Format of the input files", choices=["fasta", "text"])
        parser.add_argument("--log-file-path", "-l", type=str, required=True, help="Log file")
        parser.add_argument("--unprocessible-proteins-fasta-file-path", "-u", type=str, required=True, help="File for unprocessible proteins")
        parser.add_argument("--statistics-csv-file-path", "-s", type=str, required=True, help="File for digest statistics")
        parser.add_argument("--statistics-write-period", type=int, help="Seconds between writes to the statistics file (default: 900)", default=900)
        parser.add_argument("--thread-count", "-t", type=int, required=True, help="Number of concurrent digestions and inserts to the database")
        parser.add_argument("--enzyme-name", "-e", type=str, required=True, help="The name of the enzyme", choices=DigestEnzyme.get_known_enzymes())
        parser.add_argument("--maximum-number-of-missed-cleavages", "-c", type=int, help="Maximum number of missed cleavages (default: 2)", default="2")
        parser.add_argument("--minimum-peptide-length", type=int, help="Minimum peptide length (default: 5)", default="5")
        parser.add_argument("--maximum-peptide-length", type=int, help="Maximum peptide length (default: 60)", default="60")
        parser.add_argument("--database-url", "-d", type=str, required=True, help="Database url for postgres")
        parser.add_argument("--start-with-protein", type=int, default=0, help="In case a previous digestion was gracefully stopped e.g. by shutdown or yourself. You can start where the previous digestion stopped. Look into the log file of the previous digestion to get a clue where it stopped.")
        parser.set_defaults(func=cls.digest_from_command_line)


