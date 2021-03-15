import os
import pathlib
import psycopg2
import signal

from queue import Full as FullQueueError
from ctypes import c_ulonglong

from .multiprocessing.logger_process import LoggerProcess
from .multiprocessing.peptide_metadata_collector_process import PeptideMetadataCollectorProcess
from .multiprocessing.statistics_logger_process import StatisticsLoggerProcess
from ...models.peptide import Peptide
from ... import process_context


class PeptideMetadataCollector:
    STATISTIC_FILE_HEADER = ["seconds", "updated_pepide", "peptide_update_rate"]
    STATISTIC_PRINT_MIN_COLUMN_WIDTH = 12

    def __init__(self, log_dir_path: pathlib.Path, statistics_write_period: int, number_of_threads: int):
        self.__log_file_path = log_dir_path.joinpath(f"metadata_collector.log")
        self.__number_of_threads = number_of_threads
        self.__statistics_write_period = statistics_write_period
        self.__statistics_csv_file_path = log_dir_path.joinpath(f"metadata_collector_statistics.csv")
        self.__stop_signal = False


    def run(self, database_url: str):
        """
        Iterates through all peptides which are flagged for updates, collectes the meta information from the referenced proteins (review status, taxonomy IDs and proteome IDs) and stores them in the peptide.
        @param database_url
        """
        print("Update protein information on peptides. This may take a while.")

        # Events to control processes
        empty_queue_and_stop_flag = process_context.Event()
        logger_stop_flag = process_context.Event()
        immediate_stop_event = process_context.Event()

        # Providing a maximum size prevents overflowing RAM and makes sure every process has enough work to do.
        peptide_queue = process_context.Queue(3 * self.__number_of_threads)
        # Counter for updates
        update_counter = process_context.Array(c_ulonglong, 1)

        # Initialize signal handler for TERM and INT
        signal.signal(signal.SIGTERM, self.__stop_signal_handler)
        signal.signal(signal.SIGINT, self.__stop_signal_handler)

        log_connections = []

        print(f"to stop the peptide update gracefully, send TERM or INT signal to process {os.getpid()}")

        # Start digest worker
        peptide_update_processes = []
        for worker_id in range(0, self.__number_of_threads):
            general_log_connection_read, general_log_connection_write = process_context.Pipe(duplex=False)
            log_connections.append(general_log_connection_read)
            update_peptide_worker = PeptideMetadataCollectorProcess(worker_id, database_url, peptide_queue, update_counter, general_log_connection_write, empty_queue_and_stop_flag, immediate_stop_event)
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
        statistics_logger_process = StatisticsLoggerProcess(update_counter, self.__statistics_csv_file_path, self.__statistics_write_period, self.__class__.STATISTIC_FILE_HEADER, self.__class__.STATISTIC_PRINT_MIN_COLUMN_WIDTH, general_log_connection_write, logger_stop_flag)
        statistics_logger_process.start()
        log_connections.append(general_log_connection_read)
        # Close this copy
        general_log_connection_write.close()


        logger_process = LoggerProcess(self.__log_file_path, "w", log_connections)
        logger_process.start()

        database_connection = psycopg2.connect(database_url)
        enqueued_all_updatedable_peptides = False
        while not enqueued_all_updatedable_peptides:
            try:
                with database_connection:
                    with database_connection.cursor(name='peptide_update') as database_cursor:
                        database_cursor.execute(f"SELECT id FROM {Peptide.TABLE_NAME} WHERE is_metadata_up_to_date = false;")
                        while True:
                            # Fetch 1000 peptide IDs
                            peptides_chunk = [row[0] for row in database_cursor.fetchmany(1000)]
                            # Break while loop if no more peptides were found
                            if not len(peptides_chunk):
                                enqueued_all_updatedable_peptides = True
                                break
                            # Each peptide update process should update 100 peptides at a time, so enqueue them in list of 100 IDs
                            slice_start = 0
                            while slice_start < len(peptides_chunk):
                                while True:
                                    # Break retry loop
                                    if self.__stop_signal:
                                        break
                                    try:
                                        # Try to enqueue the ID slice, wait 5 seconds for a free slot
                                        peptide_queue.put(peptides_chunk[slice_start:slice_start+100], True, 5)
                                        # Continue with the next slice from the file
                                        break
                                    # Catch timouts of protein_queue.put()
                                    except FullQueueError:
                                        # Try again
                                        pass
                                slice_start += 100
            except psycopg2.Error as error:
                # Catch database errors and disconnects
                log_connection.send(f"error occured, see:\n{error}\ntry again")

        if database_connection and database_connection.closed != 0:
            database_connection.close()

        # Signalling all processes to 
        if self.__stop_signal:
            immediate_stop_event.set()
        empty_queue_and_stop_flag.set()

        # Wait for all digest workers to stop
        for update_peptide_worker in peptide_update_processes:
            update_peptide_worker.join()

        logger_stop_flag.set()
        statistics_logger_process.join()

        logger_process.join()

    def __stop_signal_handler(self, signal_number, frame):
        self.__stop_signal = True
