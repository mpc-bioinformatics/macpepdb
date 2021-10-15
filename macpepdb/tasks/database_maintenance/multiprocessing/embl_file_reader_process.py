# std imports
import pathlib
from queue import Full as FullQueueError
from multiprocessing import Event, Queue
from multiprocessing.connection import Connection as ProcessConnection
from typing import List

# internal imports
from macpepdb.utilities.generic_process import GenericProcess
from macpepdb.proteomics.file_reader.uniprot_text_reader import UniprotTextReader


class EmblFileReaderProcess(GenericProcess):
    """
    Reads proteins from EMBL-Files and put them into a queue for further processing.
    """
    def __init__(self, termination_event: Event, input_file_paths: List[pathlib.Path], protein_queue: Queue, general_log: ProcessConnection):
        super().__init__(termination_event)
        self.__input_file_paths = input_file_paths
        self.__protein_queue = protein_queue
        self.__general_log = general_log

    def run(self):
        self.activate_signal_handling()
        self.__general_log.send("EMBL-file reader is online")

        for input_file_path in self.__input_file_paths:
            with input_file_path.open("r") as input_file:
                input_file_reader = UniprotTextReader(input_file)
                # Enqueue proteins
                for protein in input_file_reader:
                    # Retry lopp
                    while not self.termination_event.is_set():
                        try:
                            # Try to queue a protein for digestion, wait 5 seconds for a free slot
                            self.__protein_queue.put(protein, True, 5)
                            # Continue with the next protein from the file
                            break
                        # Catch timouts of protein_queue.put()
                        except FullQueueError:
                            # Try again
                            pass
                    # Break protein read loop
                    if self.termination_event.is_set():
                        break
            # Break embl files loop
            if self.termination_event.is_set():
                break
        self.__general_log.send("EMBL-file reader is offline")
        self.__general_log.close()