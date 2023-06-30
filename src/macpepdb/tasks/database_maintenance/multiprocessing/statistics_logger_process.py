# std imports
import csv
import pathlib
import time
from multiprocessing import Array, Event
from multiprocessing.connection import Connection as ProcessConnection
from typing import List

# internal imports
from macpepdb.utilities.generic_process import GenericProcess

class StatisticsLoggerProcess(GenericProcess):
    """
    Logs statistics.
    """

    def __init__(self, termination_event: Event, statistics: Array, statistics_file_path: pathlib.Path, write_period: int, file_header: List[str], log_connection: ProcessConnection, stop_logging_event: Event):
        """
        Parameters
        ----------
        termination_event : Event
            Event for terminating the process
        statistics : Array
            Shared array for statistics
        statistics_file_path : pathlib.Path
            Path the log file
        write_period : int
            Seconds between logs
        file_header : List[str]
            CSV header should have the same length as the statistics array
        log_connection : multiprocessing.connection.Connection
            Connection to the logging process.
        stop_logging_event : Event
            Indicates the process to stop.
        """
        super().__init__(termination_event)
        self.__statistics = statistics
        self.__statistics_file_path = statistics_file_path
        self.__write_period = write_period
        self.__file_header = file_header
        self.__log_connection = log_connection
        self.__stop_logging_event = stop_logging_event

    def run(self):
        """
        Starts the process to log.
        """
        self.activate_signal_handling()
        self.__log_connection.send("statistics logger is online")
        # Snapshot of last written statistic to calculate difference
        last_statistics = [0] * len(self.__statistics)
        start_at = time.time()
        with self.__statistics_file_path.open("w") as statistics_file:
            statistics_writer = csv.writer(statistics_file)
            statistics_writer.writerow(self.__file_header)
            statistics_file.flush()
            while not self.__stop_logging_event.is_set():
                # Prepare array for next write
                current_statistics = []
                # Wait for next write
                self.__stop_logging_event.wait(self.__write_period)
                # Calculate seconds after start
                current_time = int(time.time() - start_at)
                # Get current statistics
                current_statistics = [value for value in self.__statistics]
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
                # Assign new 'snapshot'
                last_statistics = current_statistics
        self.__log_connection.send("statistics logger is offline")
        self.__log_connection.close()