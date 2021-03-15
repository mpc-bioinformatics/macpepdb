import pathlib 
from multiprocessing.connection import wait, Connection as ProcessConnection

from ....utilities.generic_process import GenericProcess

class UnprocessableProteinLoggerProcess(GenericProcess):
    def __init__(self, unprocessible_proteins_fasta_path: pathlib.Path, process_connections: list, log_connection: ProcessConnection):
        """
        Writes embl entries from process_connections to log file. Process is running until all process connections closed by the other end (`EOFError`).
        @param unprocessible_proteins_fasta_path Path to logfile
        @param process_connections List of `multiprocessing.connection.Connection`
        @param log_connection Connection to LoggerProcess
        """
        super().__init__()
        self.__unprocessible_proteins_fasta_path = unprocessible_proteins_fasta_path
        self.__process_connections = process_connections
        self.__log_connection = log_connection

    def run(self):
        self.__log_connection.send("unprocessible proteins logger is online")
        with self.__unprocessible_proteins_fasta_path.open("w") as unprocessible_proteins_file:
            while self.__process_connections:
                for conn in wait(self.__process_connections):
                    try:
                        embl_entry = conn.recv()
                    except EOFError:
                        self.__process_connections.remove(conn)
                    else:
                        unprocessible_proteins_file.write(embl_entry + "\n")
                        unprocessible_proteins_file.flush()
        self.__log_connection.send("unprocessible proteins logger is stopping")
        self.__log_connection.close()
