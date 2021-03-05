import pathlib
from multiprocessing import Process
from multiprocessing.connection import wait

class LoggerProcess:
    def __init__(self, log_file_path: pathlib.Path, write_mode: str, processing_context):
        self.__log_file_path = log_file_path
        self.__write_mode = write_mode
        self.__processing_context = processing_context
        self.__process = None

    def start(self, process_connections: list):
        """
        Starts the logging process
        @param process_connections List of `multiprocessing.connection.Connection`
        """
        self.__process = self.__processing_context.Process(target=self.logging_worker, args=(process_connections, self.__log_file_path, self.__write_mode))
        self.__process.start()

    def join(self):
        if self.__process:
            self.__process.join()

    @staticmethod
    def logging_worker(process_connections: list, log_file_path: pathlib.Path, write_mode: str):
        """
        Writes messages from process_connections to log file. Method is running until all process connections closed by the other end (`EOFError`).
        @param process_connections List of `multiprocessing.connection.Connection`
        @param log_file_path Path to logfile
        @param write_mode Write mode for the log file
        """
        with log_file_path.open(write_mode) as log_file:
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
