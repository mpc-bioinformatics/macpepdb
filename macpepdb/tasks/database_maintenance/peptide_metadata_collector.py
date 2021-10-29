# std imports
import pathlib
import signal
import shutil
import time
from ctypes import c_ulonglong
from datetime import datetime
from multiprocessing import Event, Value, Queue

# internal imports
from macpepdb import process_context
from macpepdb.tasks.database_maintenance.multiprocessing.logger_process import LoggerProcess
from macpepdb.tasks.database_maintenance.multiprocessing.peptide_metadata_collector_process import PeptideMetadataCollectorProcess
from macpepdb.tasks.database_maintenance.multiprocessing.statistics_logger_process import StatisticsLoggerProcess
from macpepdb.tasks.database_maintenance.multiprocessing.updatable_peptide_collector_process import UpdatablePeptideCollectorProcess


class PeptideMetadataCollector:
    """
    Controls the metadata collection.

    Parameter
    ---------
    termination_event: Event
        Event which controlls the termination of subprocesses.
    log_dir_path: pathlib.Path
        Path to the logs directory
    statistics_write_period: int
        Seconds between statistics logs
    number_of_threads: int
        Number of threads/processes
    """

    STATISTIC_FILE_HEADER = ["seconds", "updated_peptides", "peptide_update_rate"]
    """CSV header for log file
    """

    def __init__(self, termination_event: Event, log_dir_path: pathlib.Path, statistics_write_period: int, number_of_threads: int):
        self.__log_file_path = log_dir_path.joinpath(f"metadata_collector.log")
        self.__number_of_threads = number_of_threads
        self.__max_peptide_queue_size = 3 * self.__number_of_threads
        self.__statistics_write_period = statistics_write_period
        self.__statistics_csv_file_path = log_dir_path.joinpath(f"metadata_collector_statistics.csv")
        self.__termination_event = termination_event

    def run(self, database_url: str):
        """
        Iterates through all peptides which are flagged for updates, collectes the meta information from the referenced proteins (review status, taxonomy IDs and proteome IDs) and stores them in the peptide.

        Parameters
        ----------
        database_url : str
            E.g. postgres://username:password@host:port/database
        """
        # Initialize signal handler for TERM and INT
        signal.signal(signal.SIGTERM, self.__termination_event_handler)
        signal.signal(signal.SIGINT, self.__termination_event_handler)
        
        print("Update protein information on peptides. This may take a while.")

        # Events to control processes
        finish_update_event = process_context.Event()
        logger_stop_event = process_context.Event()

        # Providing a maximum size prevents overflowing RAM and makes sure every process has enough work to do.
        peptide_queue = process_context.Queue(self.__max_peptide_queue_size)
        # Counter for updates
        update_counter = process_context.Array(c_ulonglong, 1)

        log_connections = []

        # Start digest worker
        peptide_update_processes = []
        for worker_id in range(0, self.__number_of_threads):
            general_log_connection_read, general_log_connection_write = process_context.Pipe(duplex=False)
            log_connections.append(general_log_connection_read)
            update_peptide_worker = PeptideMetadataCollectorProcess(self.__termination_event, worker_id, database_url, peptide_queue, update_counter, general_log_connection_write, finish_update_event)
            update_peptide_worker.start()
            peptide_update_processes.append(update_peptide_worker)
            # Close these copies of the writeable site of the channel, so only the processes have working copies.
            general_log_connection_write.close()

        # Start logger for unprocessible proteins
        general_log_connection_read, general_log_connection_write = process_context.Pipe(duplex=False)
        log_connections.append(general_log_connection_read)
        # Close this copy
        general_log_connection_write.close()

        # Statistic writer
        general_log_connection_read, general_log_connection_write = process_context.Pipe(duplex=False)
        statistics_logger_process = StatisticsLoggerProcess(self.__termination_event, update_counter, self.__statistics_csv_file_path, self.__statistics_write_period, self.__class__.STATISTIC_FILE_HEADER, general_log_connection_write, logger_stop_event)
        statistics_logger_process.start()
        log_connections.append(general_log_connection_read)
        # Close this copy
        general_log_connection_write.close()

        # Updatable peptide collector
        general_log_connection_read, general_log_connection_write = process_context.Pipe(duplex=False)
        updatable_peptide_collector_process = UpdatablePeptideCollectorProcess(self.__termination_event, database_url, peptide_queue, general_log_connection_write)
        updatable_peptide_collector_process.start()
        log_connections.append(general_log_connection_read)
        # Close this copy
        general_log_connection_write.close()

        logger_process = LoggerProcess(self.__termination_event, self.__log_file_path, "w", log_connections)
        logger_process.start()

        while updatable_peptide_collector_process.is_alive():
            self.print_status(update_counter, peptide_queue)
            time.sleep(1)

        # Signalling all processes to finish
        if not self.__termination_event.is_set():
            finish_update_event.set()

        while any([proc.is_alive() for proc in peptide_update_processes]):
            self.print_status(update_counter, peptide_queue, status="waiting for update processes to stop ...")
            time.sleep(1)

        # Wait for all digest workers to stop
        for update_peptide_worker in peptide_update_processes:
            update_peptide_worker.join()

        logger_stop_event.set()
        self.print_status(update_counter, peptide_queue, status="stopping logger ...")

        statistics_logger_process.join()

        logger_process.join()

        self.print_status(update_counter, peptide_queue, is_updatable=False)

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

    def print_status(self, update_counter: Value, peptide_queue: Queue, status: str = "", is_updatable: bool = True):
        """
        Prints a status line.

        Parameters
        ----------
        update_counter: Value
            Multiprocessing value which contains the update counter
        peptide_queue: Queue
            Multiprocessing queue which contains the current queued peptide which get an update.
        status: str
            Additional status message (optional)
        is_updatable: bool = True
            Indicates if this line is replaceable.
        """
        console_width, _ = shutil.get_terminal_size()
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        message = f"\r{timestamp}> {update_counter[0]:,} updated peptides; {peptide_queue.qsize()} / {self.__max_peptide_queue_size} queue"
        if len(status):
            message += f"; {status}"
        message += ' ' * (console_width - len(message) - 1)
        print(message, end="")
        if not is_updatable:
            print("")