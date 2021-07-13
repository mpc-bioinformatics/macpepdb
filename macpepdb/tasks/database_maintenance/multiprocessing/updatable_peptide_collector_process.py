from os import terminal_size
import psycopg2

from multiprocessing import Event, Queue
from multiprocessing.connection import Connection as ProcessConnection
from queue import Full as FullQueueError
from typing import List
from macpepdb.models import peptide


from macpepdb.models.peptide import Peptide
from macpepdb.utilities.generic_process import GenericProcess

class UpdatablePeptideCollectorProcess(GenericProcess):
    """
    Reads updatable peptides from the database and put them into a queue for further processing in batches of 100 peptides
    """
    PEPTIDE_BATCH = 100

    def __init__(self, termination_event: Event, database_url: str, peptide_queue: Queue, general_log: ProcessConnection):
        super().__init__(termination_event)
        self.__database_url = database_url
        self.__peptide_queue = peptide_queue
        self.__general_log = general_log

    def __enque_peptides(self, peptides: List[Peptide]):
        while not self.termination_event.is_set():
            try:
                # Try to enqueue the ID slice, wait 5 seconds for a free slot
                self.__peptide_queue.put(peptides, True, 5)
                break
            # Catch timouts of protein_queue.put()
            except FullQueueError:
                # Try again
                pass


    def run(self):
        self.activate_signal_handling()

        self.__general_log.send("Start enqueuing updatable peptides.")
        database_connection = None
        # retry loop
        while not self.termination_event.is_set():
            try:
                if not database_connection or (database_connection and database_connection.closed != 0):
                    database_connection = psycopg2.connect(self.__database_url)
                with database_connection.cursor(name='updatable_peptide_collector') as database_cursor:
                    database_cursor.itersize = 1000
                    database_cursor.execute(f"SELECT sequence, number_of_missed_cleavages FROM {Peptide.TABLE_NAME} WHERE is_metadata_up_to_date = false;")
                    peptides = []
                    for peptide_row in database_cursor:
                        peptides.append(Peptide(peptide_row[0], peptide_row[1]))
                        if len(peptides) == self.__class__.PEPTIDE_BATCH:
                            self.__enque_peptides(peptides)
                            peptides = []
                        if self.termination_event.is_set():
                            break # breat cursor loop
                    # Enqueue last batch of peptides
                    if len(peptides):
                        self.__enque_peptides(peptides)
                        peptides = []
                    break # break retry loop
            except psycopg2.Error as error:
                # Catch database errors and disconnects
                self.__general_log.send(f"error occured, see:\n{error}\ntry again")
        self.__general_log.send("All updatable peptides enqueued.")
        self.__general_log.close()