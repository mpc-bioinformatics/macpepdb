import csv
import pathlib
import time
from multiprocessing import Array, Event
from multiprocessing.connection import Connection as ProcessConnection

from ....utilities.generic_process import GenericProcess

class StatisticsLoggerProcess(GenericProcess):
    def __init__(self, statistics: Array, statistics_file_path: pathlib.Path, write_period: int, file_header: str, output_min_column_width: int, log_connection: ProcessConnection, stop_flag: Event):
        super().__init__()
        self.__statistics = statistics
        self.__statistics_file_path = statistics_file_path
        self.__write_period = write_period
        self.__file_header = file_header
        self.__output_min_column_width = output_min_column_width
        self.__log_connection = log_connection
        self.__stop_flag = stop_flag

    def run(self):
        self.__log_connection.send("statistics logger is online")
        # Snapshot of last written statistic to calculate difference
        last_statistics = [0] * len(self.__statistics)
        start_at = time.time()
        self.__print_statistic_row(self.__file_header, False)
        with self.__statistics_file_path.open("w") as statistics_file:
            statistics_writer = csv.writer(statistics_file)
            statistics_writer.writerow(self.__file_header)
            statistics_file.flush()
            while not self.__stop_flag.is_set():
                # Prepare array for next write
                current_statistics = []
                # Wait for next write
                self.__stop_flag.wait(self.__write_period)
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
                # Output to console, when the stop flag is not set, make line replaceable
                self.__print_statistic_row(csv_row, not self.__stop_flag.is_set())
                # Assign new 'snapshot'
                last_statistics = current_statistics
        print("\n")
        self.__log_connection.send("statistics logger is offline")
        self.__log_connection.close()

    def __print_statistic_row(self, row: list, is_replaceable: bool):
        line_end = "\r" if is_replaceable else "\n"
        output = []
        for idx in range(0, len(self.__file_header)):
            column_width = max(
                self.__output_min_column_width,
                len(self.__file_header[idx])
            )
            output.append(f"{row[idx]:>{column_width}}")
        print("\t".join(output), end=line_end)