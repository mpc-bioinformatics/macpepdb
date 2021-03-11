import pathlib
import time
import sys
import re
import sys
import csv
import random
import signal
import os
import json

import psycopg2
from psycopg2.extras import execute_values

from multiprocessing import Process, Value, Queue, Array, Event, Pipe, get_context as get_process_context
from multiprocessing.connection import Connection, wait
from ctypes import c_bool, c_ulonglong
from queue import Empty, Full

from ...proteomics.enzymes.digest_enzyme import DigestEnzyme
from ...proteomics.file_reader.uniprot_text_reader import UniprotTextReader
from ...models.protein import Protein
from ...models.protein_peptide_association import ProteinPeptideAssociation
from ...models.peptide import Peptide
from ...models.maintenance_information import MaintenanceInformation
from .logger_process import LoggerProcess

class ProteinDigestion:
    STATSTIC_FILE_HEADER = ["seconds", "inserted_proteins", "inserted_peptides", "unsolvable errors", "protein_insert_rate", "peptide_insert_rate", "error_rate"]
    STATSTIC_PRINT_MIN_COLUMN_WIDTH = 12
    

    def __init__(self, protein_data_dir: pathlib.Path, log_dir_path: pathlib.Path, statistics_write_period: int, number_of_threads: int, enzyme_name: str, maximum_number_of_missed_cleavages: int, minimum_peptide_length: int, maximum_peptide_length: int, run_count: int):
        self.__log_file_path = log_dir_path.joinpath(f"digest_{run_count}.log")
        self.__unprocessible_proteins_fasta_file_path = log_dir_path.joinpath(f"unprocessible_proteins_{run_count}.txt")
        self.__number_of_threads = number_of_threads
        EnzymeClass = DigestEnzyme.get_enzyme_by_name(enzyme_name)
        self.__enzyme = EnzymeClass(maximum_number_of_missed_cleavages, minimum_peptide_length, maximum_peptide_length)
        self.__input_file_paths = [pathlib.Path(path) for path in protein_data_dir.glob('*.txt')]
        self.__input_file_paths += [pathlib.Path(path) for path in protein_data_dir.glob('*.dat')]
        self.__statistics_write_period = statistics_write_period
        self.__statistics_csv_file_path = log_dir_path.joinpath(f"statistics_{run_count}.csv")
        self.__stop_signal = False

    def run(self, database_url: str) -> int:
        """
        Reads the protein files in `work_dir/protein_data` and updates the database.
        @param database_url Datebase
        @return int Number of errors
        """
        self.__load_or_set_digestion_informations(database_url)

        # Very important 
        process_context = get_process_context("spawn")

        # Events to control processes
        empty_queue_and_stop_flag = process_context.Event()
        logger_stop_flag = process_context.Event()
        stop_immediately_flag = process_context.Event()
        # Providing a maximum size prevents overflowing RAM and makes sure every process has enough work to do.
        protein_queue = process_context.Queue(3 * self.__number_of_threads)
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
        for worker_id in range(0, self.__number_of_threads):
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

        logger_process = LoggerProcess(self.__log_file_path, "w", process_context)
        logger_process.start(log_connections)

        for input_file_path in self.__input_file_paths:
            current_input_file = str(input_file_path.resolve())
            # Break file loop
            if self.__stop_signal:
                break
            with input_file_path.open("r") as input_file:
                input_file_reader = UniprotTextReader(input_file)
                # Enqueue proteins
                for protein in input_file_reader:
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
                            protein_queue.put(protein, True, 5)
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

        return statistics[2]


    def __load_or_set_digestion_informations(self, database_url: str):
        """
        Loads the digestion information from previous digestions from the databsae and overrides the given ones.
        If no digestion information were found, it will save the current one.
        @param database_url Database URL
        """
        database_connection = psycopg2.connect(database_url)
        with database_connection.cursor() as database_cursor:
            database_cursor.execute("SELECT values FROM maintenance_information WHERE key = %s;", (MaintenanceInformation.DIGESTION_PARAMTERS_KEY,))
            digestion_information_row = database_cursor.fetchone()
            if digestion_information_row:
                digestion_information_values = digestion_information_row[0]
                DigestEnzymeClass = DigestEnzyme.get_enzyme_by_name(digestion_information_values['enzyme_name'])
                self.__enzyme = DigestEnzymeClass(
                    digestion_information_values['maximum_number_of_missed_cleavages'],
                    digestion_information_values['minimum_peptide_length'], 
                    digestion_information_values['maximum_peptide_length']
                )
            else:
                digestion_information_values = {
                    'enzyme_name': self.__enzyme.NAME,
                    'maximum_number_of_missed_cleavages': self.__enzyme.max_number_of_missed_cleavages,
                    'minimum_peptide_length': self.__enzyme.minimum_peptide_length,
                    'maximum_peptide_length': self.__enzyme.maximum_peptide_length
                }
                database_cursor.execute("INSERT INTO maintenance_information (key, values) VALUES (%s, %s);", (MaintenanceInformation.DIGESTION_PARAMTERS_KEY, json.dumps(digestion_information_values)))
                database_connection.commit()
        database_connection.close()

    @staticmethod
    def digest_worker(id: int, empty_queue_and_stop_flag: Event, stop_immediately_flag: Event, protein_queue: Queue, log_connection: Connection, unprocessible_log_connection: Connection, statistics: Array, enzyme: DigestEnzyme, database_url: str):
        log_connection.send("digest worker {} is online".format(id))
        database_connection = None

        # Let the process run until empty_queue_and_stop_flag is true and protein_queue is empty or stop_immediately_flag is true.
        while (not empty_queue_and_stop_flag.is_set() or not protein_queue.empty()) and not stop_immediately_flag.is_set():
            try:
                # Open/reopen database connection
                if not database_connection or (database_connection and database_connection.closed != 0):
                    database_connection = psycopg2.connect(database_url)

                # Try to get a protein from the queue, timeout is 2 seconds
                protein = protein_queue.get(True, 5)
                
                # Variables for loop control
                unsolvable_errors = 0
                try_transaction_again = True
                while try_transaction_again:
                    number_of_new_peptides = 0
                    try:
                        count_protein = False
                        number_of_new_peptides = 0
                        with database_connection:
                            with database_connection.cursor() as database_cursor:
                                # Check if the Protein exists by its accession or secondary accessions
                                accessions = [protein.accession] + protein.secondary_accessions
                                stored_protein = Protein.select(database_cursor, ("accession = ANY(%s)", [accessions]))
                                if stored_protein:
                                    number_of_new_peptides = stored_protein.update(database_cursor, protein, enzyme)
                                else:
                                    number_of_new_peptides = Protein.create(database_cursor, protein, enzyme)
                                    count_protein = True

                        # Commit was successfully stop while-loop and add statistics
                        try_transaction_again = False
                        statistics.acquire()
                        if count_protein:
                            statistics[0] += 1
                        statistics[1] += number_of_new_peptides
                        statistics.release()
                    # Catch all transaction errors
                    except psycopg2.Error as error:
                        # Rollback is done implcit by `with database_connection`
                        # Remove all peptides from protein
                        protein.peptides = []
                        # Try again after 5 (first try) and 10 (second try) + a random number between 0 and 5 (both including) seconds maybe some blocking transactions can pass so this transaction will successfully finish on the next try.
                        # If this is the third time an unsolvable error occures give up and log the error.
                        if unsolvable_errors < 2:
                            unsolvable_errors += 1
                            time.sleep(5 * unsolvable_errors + random.randint(0, 5))
                        # Log the error on the third try and put the protein in unprocessible queue
                        else:
                            log_connection.send("Exception on protein {}, see:\n{}".format(protein.accession, error))
                            unprocessible_log_connection.send(protein.to_embl_entry())
                            statistics.acquire()
                            statistics[2] += 1
                            statistics.release()
                            try_transaction_again = False
            # Catch errors which occure during databse connect
            except psycopg2.Error as error:
                log_connection.send("Error when opening the database connection, see:\n{}".format(error))
            # Catch queue.Empty which is thrown when protein_queue.get() timed out
            except Empty:
                pass
        # Close database connection
        if database_connection and database_connection.closed == 0:
            database_connection.close()
        log_connection.send("digest worker {} is stopping".format(id))
        log_connection.close()
        unprocessible_log_connection.close()

    @staticmethod
    def unprocessible_proteins_worker(process_connections: list, unprocessible_proteins_fasta_path: pathlib.Path, log_connection: Connection):
        log_connection.send("unprocessible proteins logger is online")
        with unprocessible_proteins_fasta_path.open("w") as unprocessible_proteins_file:
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
        ProteinDigestion.print_statistic_row(ProteinDigestion.STATSTIC_FILE_HEADER)
        with statistics_file_path.open("w") as statistics_file:
            statistics_writer = csv.writer(statistics_file)
            statistics_writer.writerow(ProteinDigestion.STATSTIC_FILE_HEADER)
            statistics_file.flush()
            while not stop_flag.is_set():
                # Prepare array for next write
                current_statistics = []
                # Wait for next write
                stop_flag.wait(write_period)
                # Calculate seconds after start
                current_time = int(time.time() - start_at)
                # Get current statistics
                current_statistics = [value for value in statistics]
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
                ProteinDigestion.print_statistic_row(csv_row)
                # Assign new 'snapshot'
                last_statistics = current_statistics
        log_connection.send("statistics logger is offline")
        log_connection.close()
    
    @staticmethod
    def print_statistic_row(row: list):
        output = []
        for idx in range(0, len(ProteinDigestion.STATSTIC_FILE_HEADER)):
            column_width = max(
                ProteinDigestion.STATSTIC_PRINT_MIN_COLUMN_WIDTH,
                len(ProteinDigestion.STATSTIC_FILE_HEADER[idx])
            )
            output.append("{value:>{column_width}}".format(column_width=column_width, value=row[idx]))
        print("\t".join(output))

    def __stop_signal_handler(self, signal_number, frame):
        self.__stop_signal = True
