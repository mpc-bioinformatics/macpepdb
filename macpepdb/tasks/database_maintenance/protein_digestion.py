import pathlib
import time
import signal
import json

import psycopg2

from ctypes import c_ulonglong
from queue import Full as FullQueueError

from .multiprocessing.logger_process import LoggerProcess
from .multiprocessing.protein_digestion_process import ProteinDigestionProcess
from .multiprocessing.statistics_logger_process import StatisticsLoggerProcess
from .multiprocessing.unprocessable_protein_logger_process import UnprocessableProteinLoggerProcess
from ...proteomics.enzymes.digest_enzyme import DigestEnzyme
from ...proteomics.file_reader.uniprot_text_reader import UniprotTextReader
from ...models.maintenance_information import MaintenanceInformation
from ... import process_context

class ProteinDigestion:
    STATISTIC_FILE_HEADER = ["seconds", "inserted_proteins", "inserted_peptides", "unsolvable errors", "protein_insert_rate", "peptide_insert_rate", "error_rate"]
    STATISTIC_PRINT_MIN_COLUMN_WIDTH = 12
    

    def __init__(self, protein_data_dir: pathlib.Path, log_dir_path: pathlib.Path, statistics_write_period: int, number_of_threads: int, enzyme_name: str, maximum_number_of_missed_cleavages: int, minimum_peptide_length: int, maximum_peptide_length: int, run_count: int):
        self.__log_file_path = log_dir_path.joinpath(f"digest_{run_count}.log")
        self.__unprocessible_proteins_embl_file_path = log_dir_path.joinpath(f"unprocessible_proteins_{run_count}.txt")
        self.__number_of_threads = number_of_threads
        EnzymeClass = DigestEnzyme.get_enzyme_by_name(enzyme_name)
        self.__enzyme = EnzymeClass(maximum_number_of_missed_cleavages, minimum_peptide_length, maximum_peptide_length)
        self.__input_file_paths = [pathlib.Path(path) for path in protein_data_dir.glob('*.txt')]
        self.__input_file_paths += [pathlib.Path(path) for path in protein_data_dir.glob('*.dat')]
        self.__statistics_write_period = statistics_write_period
        self.__statistics_csv_file_path = log_dir_path.joinpath(f"statistics_{run_count}.csv")
        self.__stop_signal = False

    def run(self, database_url: str) -> (int, bool):
        """
        Reads the protein files in `work_dir/protein_data` and updates the database.
        @param database_url Datebase
        @return tupel First element isthe number of errors and the second element is the status of the stop flag
        """
        self.__load_or_set_digestion_informations(database_url)

        # Events to control processes
        clear_queue_and_stop_event = process_context.Event()
        logger_stop_flag = process_context.Event()
        immediate_stop_event = process_context.Event()
        # Providing a maximum size prevents overflowing RAM and makes sure every process has enough work to do.
        protein_queue = process_context.Queue(3 * self.__number_of_threads)
        # Array for statistics [created_proteins, failed_proteins, created_peptides, processed_peptides]
        statistics = process_context.Array(c_ulonglong, 3)

        # Initialize signal handler for TERM and INT
        signal.signal(signal.SIGTERM, self.__stop_signal_handler)
        signal.signal(signal.SIGINT, self.__stop_signal_handler)

        unprocessable_log_connection = []
        general_log_connections = []

        # Start digest worker
        digest_processes = []
        for worker_id in range(0, self.__number_of_threads):
            unprocessible_connection_read, unprocessible_connection_write = process_context.Pipe(duplex=False)
            unprocessable_log_connection.append(unprocessible_connection_read)
            general_log_connection_read, general_log_connection_write = process_context.Pipe(duplex=False)
            general_log_connections.append(general_log_connection_read)
            digest_worker = ProteinDigestionProcess(worker_id, database_url, protein_queue, self.__enzyme, general_log_connection_write, unprocessible_connection_write, statistics, clear_queue_and_stop_event, immediate_stop_event)
            digest_worker.start()
            digest_processes.append(digest_worker)
            # Close these copies of the writeable site of the channels, so only the processes have working copies.
            unprocessible_connection_write.close()
            general_log_connection_write.close()

        # Start logger for unprocessible proteins
        general_log_connection_read, general_log_connection_write = process_context.Pipe(duplex=False)
        unprocessible_proteins_process = UnprocessableProteinLoggerProcess(self.__unprocessible_proteins_embl_file_path, unprocessable_log_connection, general_log_connection_write)
        unprocessible_proteins_process.start()
        general_log_connections.append(general_log_connection_read)
        # Close this copy
        general_log_connection_write.close()

        # Statistic writer
        general_log_connection_read, general_log_connection_write = process_context.Pipe(duplex=False)
        statistics_logger_process = StatisticsLoggerProcess(statistics, self.__statistics_csv_file_path, self.__statistics_write_period, self.__class__.STATISTIC_FILE_HEADER, self.__class__.STATISTIC_PRINT_MIN_COLUMN_WIDTH, general_log_connection_write, logger_stop_flag)
        statistics_logger_process.start()
        general_log_connections.append(general_log_connection_read)
        # Close this copy
        general_log_connection_write.close()

        logger_process = LoggerProcess(self.__log_file_path, "w", general_log_connections)
        logger_process.start()

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
                        except FullQueueError:
                            # Try again
                            pass

        # Signalling all processes to 
        if self.__stop_signal:
            immediate_stop_event.set()
        clear_queue_and_stop_event.set()

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

        logger_process.join()

        return statistics[2], self.__stop_signal

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

    def __stop_signal_handler(self, signal_number, frame):
        self.__stop_signal = True
