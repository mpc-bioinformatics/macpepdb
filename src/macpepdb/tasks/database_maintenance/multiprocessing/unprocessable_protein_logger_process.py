# std imports
import pathlib 
from multiprocessing import Event
from multiprocessing.connection import wait, Connection as ProcessConnection

# internal imports
from macpepdb.utilities.generic_process import GenericProcess

class UnprocessableProteinLoggerProcess(GenericProcess):
    """
    Writes embl entries from process_connections to log file. Process is running until all process connections closed by the other end (`EOFError`).
    """
    def __init__(self, termination_event: Event, unprocessible_proteins_fasta_path: pathlib.Path, process_connections: list, log_connection: ProcessConnection):
        """
        Parameters
        ----------
        termination_event : Event
            Event for terminating the process
        unprocessible_proteins_fasta_path : pathlib.Path
            Path to logfile
        process_connections : List[]
            List of read only from other processes
        log_connection : multiprocessing.connection.Connection
            Connection to log process

        Returns
        -------
        UnprocessableProteinLoggerProcess
        """
        super().__init__(termination_event)
        self.__unprocessible_proteins_fasta_path = unprocessible_proteins_fasta_path
        self.__process_connections = process_connections
        self.__log_connection = log_connection

    def run(self):
        """
        Starts process and logging.
        """
        self.activate_signal_handling()
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
