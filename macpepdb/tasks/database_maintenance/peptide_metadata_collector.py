import csv
import os
import pathlib
import psycopg2
import signal
import time

from multiprocessing import Process, Value, Queue, Array, Event, Pipe, get_context as get_process_context
from multiprocessing.connection import Connection
from queue import Empty, Full
from ctypes import c_ulonglong

from .logger_process import LoggerProcess
from ...models.peptide import Peptide


class PeptideMetadataCollector:
    STATISTIC_FILE_HEADER = ["seconds", "updated_pepide", "peptide_update_rate"]
    STATSTIC_PRINT_MIN_COLUMN_WIDTH = 12

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
        # Very important 
        process_context = get_process_context("spawn")

        # Events to control processes
        empty_queue_and_stop_flag = process_context.Event()
        logger_stop_flag = process_context.Event()
        stop_immediately_flag = process_context.Event()

        # Providing a maximum size prevents overflowing RAM and makes sure every process has enough work to do.
        peptide_queue = process_context.Queue(3 * self.__number_of_threads)
        # Counter for updates
        update_counter = process_context.Value(c_ulonglong)

        # Initialize signal handler for TERM and INT
        signal.signal(signal.SIGTERM, self.__stop_signal_handler)
        signal.signal(signal.SIGINT, self.__stop_signal_handler)

        log_connections = []

        # Log connection for the main process
        log_connection_read, log_connection = process_context.Pipe(duplex=False)
        log_connections.append(log_connection_read)

        print(f"to stop the peptide update gracefully, send TERM or INT signal to process {os.getpid()}")
        log_connection.send(f"process id = {os.getpid()}")

        # Start digest worker
        peptide_update_processes = []
        for worker_id in range(0, self.__number_of_threads):
            log_connection_read, log_connection_write = process_context.Pipe(duplex=False)
            log_connections.append(log_connection_read)
            update_peptide_worker = process_context.Process(target=PeptideMetadataCollector.update_peptide_worker, args=(worker_id, empty_queue_and_stop_flag, stop_immediately_flag, peptide_queue, log_connection_write, update_counter, database_url,))
            update_peptide_worker.start()
            peptide_update_processes.append(update_peptide_worker)
            # Close these copies of the writeable site of the channel, so only the processes have working copies.
            log_connection_write.close()

        # Start logger for unprocessible proteins
        log_connection_read, log_connection_write = process_context.Pipe(duplex=False)
        log_connections.append(log_connection_read)
        # Close this copy
        log_connection_write.close()

        # Statistic writer
        log_connection_read, log_connection_write = process_context.Pipe(duplex=False)
        statistics_logger_process = process_context.Process(target=PeptideMetadataCollector.statistic_worker, args=(logger_stop_flag, update_counter, self.__statistics_write_period, self.__statistics_csv_file_path, log_connection_write,))
        statistics_logger_process.start()
        log_connections.append(log_connection_read)
        # Close this copy
        log_connection_write.close()


        logger_process = LoggerProcess(self.__log_file_path, "w", process_context)
        logger_process.start(log_connections)

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
                                    except Full:
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
            stop_immediately_flag.set()
        empty_queue_and_stop_flag.set()


        # Wait for all digest workers to stop
        for update_peptide_worker in peptide_update_processes:
            update_peptide_worker.join()

        logger_stop_flag.set()
        statistics_logger_process.join()

        log_connection.close()
        logger_process.join()


    def update_peptide_worker(id: int, empty_queue_and_stop_flag: Event, stop_immediately_flag: Event, peptide_id_queue: Queue, log_connection: Connection, update_counter: Value, database_url: str):
        """
        Collects the information from the referenced peptides und set update flag to false.
        @param database_cursor psycopg2 cursor with open transaction
        @param peptide_id Peptide ID
        """
        log_connection.send(f"peptide update worker {id} is online")
        database_connection = None

        # Let the process run until empty_queue_and_stop_flag is true and peptide_id_queue is empty or stop_immediately_flag is true.
        while (not empty_queue_and_stop_flag.is_set() or not peptide_id_queue.empty()) and not stop_immediately_flag.is_set():
            try:
                # Open/reopen database connection
                if not database_connection or (database_connection and database_connection.closed != 0):
                    database_connection = psycopg2.connect(database_url)
                # Wait 5 seconds for a new set of IDs
                peptide_ids = peptide_id_queue.get(True, 5)
                while True:
                    try:
                        with database_connection:
                            with database_connection.cursor() as database_cursor:
                                for peptide_id in peptide_ids:
                                    Peptide.update_metadata(database_cursor, peptide_id)
                        with update_counter.get_lock():
                            update_counter.value += len(peptide_ids)
                        break
                    except BaseException as error:
                        log_connection.send(f"error occured, see:\n{error}\ntry again")
            # Catch errors which occure during databse connect
            except psycopg2.Error as error:
                log_connection.send("Error when opening the database connection, see:\n{}".format(error))
            # Catch queue.Empty which is thrown when protein_queue.get() timed out
            except Empty:
                pass

        log_connection.send(f"peptide update worker {id} is stopping")
        if database_connection and database_connection.closed != 0:
            database_connection.close()

    def __stop_signal_handler(self, signal_number, frame):
        self.__stop_signal = True


    @staticmethod
    def statistic_worker(stop_flag: Event, update_counter: Value, write_period: int, statistics_file_path: pathlib.Path, log_connection: Connection):
        log_connection.send("peptide update statistics logger is online")
        # Snapshot of last written statistic to calculate difference
        last_update_counter = update_counter.value
        start_at = time.time()
        PeptideMetadataCollector.print_statistic_row(PeptideMetadataCollector.STATISTIC_FILE_HEADER)
        with statistics_file_path.open("w") as statistics_file:
            statistics_writer = csv.writer(statistics_file)
            statistics_writer.writerow(PeptideMetadataCollector.STATISTIC_FILE_HEADER)
            statistics_file.flush()
            while not stop_flag.is_set():
                # Wait for next write
                stop_flag.wait(write_period)
                # Calculate seconds after start
                current_time = int(time.time() - start_at)
                current_counter = update_counter.value
                # Initialize csv row
                csv_row = [current_time, current_counter]
                # Calulate differences to last statistics (= rates)
                csv_row.append(current_counter - last_update_counter)
                # Write csv row to csv file
                statistics_writer.writerow(csv_row)
                statistics_file.flush()
                # Output to console
                PeptideMetadataCollector.print_statistic_row(csv_row)
                # Assign new 'snapshot'
                last_update_counter = current_counter
        log_connection.send("peptide update statistics logger is offline")
        log_connection.close()

    @staticmethod
    def print_statistic_row(row: list):
        output = []
        for idx in range(0, len(PeptideMetadataCollector.STATISTIC_FILE_HEADER)):
            column_width = max(
                PeptideMetadataCollector.STATSTIC_PRINT_MIN_COLUMN_WIDTH,
                len(PeptideMetadataCollector.STATISTIC_FILE_HEADER[idx])
            )
            output.append("{value:>{column_width}}".format(column_width=column_width, value=row[idx]))
        print("\t".join(output))

    def __stop_signal_handler(self, signal_number, frame):
        self.__stop_signal = True
