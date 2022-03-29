# std imports
import pathlib
import time
import signal
import shutil
import json
from ctypes import c_ulonglong
from datetime import datetime
from multiprocessing import Event, Array, Queue

# external imports
import psycopg2

# internal imports
from macpepdb import process_context
from macpepdb.models.maintenance_information import MaintenanceInformation
from macpepdb.proteomics.enzymes.digest_enzyme import DigestEnzyme
from macpepdb.tasks.database_maintenance.multiprocessing.embl_file_reader_process import EmblFileReaderProcess
from macpepdb.tasks.database_maintenance.multiprocessing.logger_process import LoggerProcess
from macpepdb.tasks.database_maintenance.multiprocessing.protein_digestion_process import ProteinDigestionProcess
from macpepdb.tasks.database_maintenance.multiprocessing.statistics_logger_process import StatisticsLoggerProcess
from macpepdb.tasks.database_maintenance.multiprocessing.unprocessable_protein_logger_process import UnprocessableProteinLoggerProcess

class ProteinDigestion:
    """
    Controls the protein digestion.

    Parameters
    ----------
    termination_event: Event
        Termination event
    protein_data_dir: pathlib.Path
        Path to the protein data directory of the workdir
    log_dir_path: pathlib.Path
        Path to the logs directory of the workdir
    statistics_write_period: int
        Second between statistic logs
    enzyme_name: str
        Name of the digestion enzym
    maximum_number_of_missed_cleavages: int
        Maximum number of missed cleavages
    minimum_peptide_length: int
        Minimum peptide length
    maximum_peptide_length: int
        Maximum peptide length
    run_count: int
        Indicates the current digestion run. Multiple runs may be necessary to digest all data.
    """

    STATISTIC_FILE_HEADER = ["seconds", "inserted_proteins", "inserted_peptides", "unsolvable_errors", "protein_insert_rate", "peptide_insert_rate", "error_rate"]
    """CSV header for the log file
    """

    def __init__(self, termination_event: Event, protein_data_dir: pathlib.Path, log_dir_path: pathlib.Path, statistics_write_period: int, number_of_threads: int, enzyme_name: str, maximum_number_of_missed_cleavages: int, minimum_peptide_length: int, maximum_peptide_length: int, run_count: int):
        self.__log_file_path = log_dir_path.joinpath(f"digest_{run_count}.log")
        self.__unprocessible_proteins_embl_file_path = log_dir_path.joinpath(f"unprocessible_proteins_{run_count}.txt")
        self.__number_of_threads = number_of_threads
        self.__max_protein_queue_size = 3 * self.__number_of_threads
        EnzymeClass = DigestEnzyme.get_enzyme_by_name(enzyme_name)
        self.__enzyme = EnzymeClass(maximum_number_of_missed_cleavages, minimum_peptide_length, maximum_peptide_length)
        self.__input_file_paths = [pathlib.Path(path) for path in protein_data_dir.glob('*.txt')]
        self.__input_file_paths += [pathlib.Path(path) for path in protein_data_dir.glob('*.dat')]
        self.__statistics_write_period = statistics_write_period
        self.__statistics_csv_file_path = log_dir_path.joinpath(f"statistics_{run_count}.csv")
        self.__termination_event = termination_event

    def run(self, database_url: str) -> int:
        """
        Reads the protein files in `work_dir/protein_data` and updates the database.

        Parameters
        ----------
        database_url: str
            Database URL, e.g. postgres://username:password@host:port/database

        Returns
        -------
        Number of digestion errors.
        """
        # Bevore initalizing the whole multiprocessing stuff, just check if there is an actual datafile to process.
        if len(self.__input_file_paths) == 0:
            return 0

        # Initialize signal handler for TERM and INT
        signal.signal(signal.SIGTERM, self.__termination_event_handler)
        signal.signal(signal.SIGINT, self.__termination_event_handler)

        self.__load_or_set_digestion_informations(database_url)

        # Events to control processes
        finish_digestion_event = process_context.Event()
        stop_logging_event = process_context.Event()
        # Providing a maximum size prevents overflowing RAM and makes sure every process has enough work to do.
        protein_queue = process_context.Queue(self.__max_protein_queue_size)
        # Array for statistics [created_proteins, created_peptides, number_of_errors]
        statistics = process_context.Array(c_ulonglong, 3)

        unprocessable_log_connection = []
        general_log_connections = []

        # Start digest worker
        digest_processes = []
        for worker_id in range(0, self.__number_of_threads):
            unprocessible_connection_read, unprocessible_connection_write = process_context.Pipe(duplex=False)
            unprocessable_log_connection.append(unprocessible_connection_read)
            general_log_connection_read, general_log_connection_write = process_context.Pipe(duplex=False)
            general_log_connections.append(general_log_connection_read)
            digest_worker = ProteinDigestionProcess(self.__termination_event, worker_id, database_url, protein_queue, self.__enzyme, general_log_connection_write, unprocessible_connection_write, statistics, finish_digestion_event)
            digest_worker.start()
            digest_processes.append(digest_worker)
            # Close these copies of the writeable site of the channels, so only the processes have working copies.
            unprocessible_connection_write.close()
            general_log_connection_write.close()

        # Start logger for unprocessible proteins
        general_log_connection_read, general_log_connection_write = process_context.Pipe(duplex=False)
        unprocessible_proteins_process = UnprocessableProteinLoggerProcess(self.__termination_event, self.__unprocessible_proteins_embl_file_path, unprocessable_log_connection, general_log_connection_write)
        unprocessible_proteins_process.start()
        general_log_connections.append(general_log_connection_read)
        # Close this copy
        general_log_connection_write.close()

        # Start process for EMBL reading
        general_log_connection_read, general_log_connection_write = process_context.Pipe(duplex=False)
        embl_file_reader_process = EmblFileReaderProcess(self.__termination_event, self.__input_file_paths, protein_queue, general_log_connection_write)
        embl_file_reader_process.start()
        general_log_connections.append(general_log_connection_read)
        # Close this copy
        general_log_connection_write.close()

        # Statistic writer
        general_log_connection_read, general_log_connection_write = process_context.Pipe(duplex=False)
        statistics_logger_process = StatisticsLoggerProcess(self.__termination_event, statistics, self.__statistics_csv_file_path, self.__statistics_write_period, self.__class__.STATISTIC_FILE_HEADER, general_log_connection_write, stop_logging_event)
        statistics_logger_process.start()
        general_log_connections.append(general_log_connection_read)
        # Close this copy
        general_log_connection_write.close()

        logger_process = LoggerProcess(self.__termination_event, self.__log_file_path, "w", general_log_connections)
        logger_process.start()

        while embl_file_reader_process.is_alive():
            self.print_status(statistics, protein_queue)
            time.sleep(1)

        # Signal processes to finish work
        if not self.__termination_event.is_set():
            finish_digestion_event.set()

        while any([proc.is_alive() for proc in digest_processes]):
            self.print_status(statistics, protein_queue, "wait for digest processes to stop...")
            time.sleep(1)

        # Wait for all digest workers to stop
        for digest_worker in digest_processes:
            digest_worker.join()

        # It is possible that the last statistic is not entirely transfered when the logger stops. So wait a second to make sure all information are transfered.
        time.sleep(1)

        # Stop logger after each digest worker is finished
        stop_logging_event.set()

        self.print_status(statistics, protein_queue, status = "stopping statistics logger ...")
        
        # Wait for the logger processes to stop
        statistics_logger_process.join()
        unprocessible_proteins_process.join()

        logger_process.join()

        self.print_status(statistics, protein_queue, status = "stopping logger ...")

        self.print_status(statistics, protein_queue, is_updatable=False)

        return statistics[2]

    def __load_or_set_digestion_informations(self, database_url: str):
        """
        Loads the digestion information from previous digestions from the databsae and overrides the given ones.
        If no digestion information were found, it will save the current one.

        Parameters
        ----------
        database_url: str
            Database URL, e.g. postgres://username:password@host:port/database
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

    def print_status(self, statistics: Array, protein_queue: Queue, status: str = "", is_updatable: bool = True):
        """
        Prints a status line.

        Parameters
        ----------
        statistics: Array
            Multiprocessing array which contains the inserted protein, inserted peptides and errors.
        protein_queue: Queue
            Multiprocessing queue which contains the current queued proteins for digestion.
        status: str
            Additional status message (optional)
        is_updatable: bool = True
            Indicates if this line is replaceable.
        """
        console_width, _ = shutil.get_terminal_size()
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        message = f"\r{timestamp}> {statistics[0]:,} proteins; {statistics[1]:,} peptides; {statistics[2]:,} errors; {protein_queue.qsize()} / {self.__max_protein_queue_size} queue"
        if len(status):
            message += f"; {status}"
        message += ' ' * (console_width - len(message) - 1)
        print(message, end="")
        if not is_updatable:
            print("")