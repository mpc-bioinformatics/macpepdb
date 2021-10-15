# std imports
import pathlib 
from multiprocessing import Event
from multiprocessing.connection import wait

# internal imports
from macpepdb.utilities.generic_process import GenericProcess

class LoggerProcess(GenericProcess):
    def __init__(self, termination_event: Event, log_file_path: pathlib.Path, write_mode: str, process_connections: list):
        """
        Writes messages from process_connections to log file. Process is running until all process connections closed by the other end (`EOFError`).
        @param termination_event 
        @param log_file_path Path to logfile
        @param write_mode Write mode for the log file
        @param process_connections List of `multiprocessing.connection.Connection`
        """
        super().__init__(termination_event)
        self.__log_file_path = log_file_path
        self.__write_mode = write_mode
        self.__process_connections = process_connections

    def run(self):
        self.activate_signal_handling()
        with self.__log_file_path.open(self.__write_mode) as log_file:
            log_file.write("error logger is online\n")
            log_file.flush()
            while self.__process_connections:
                for conn in wait(self.__process_connections):
                    try:
                        message = conn.recv()
                    except EOFError:
                        self.__process_connections.remove(conn)
                    else:
                        log_file.write(f"{message}\n")
                        log_file.flush()
            log_file.write("error logger is stopping")
            # will be flushed on file close
