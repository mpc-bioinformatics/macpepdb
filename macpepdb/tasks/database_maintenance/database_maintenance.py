import datetime
import json
import pathlib
import glob
import shutil
import psycopg2

from enum import unique, Enum

from ...models.maintenance_information import MaintenanceInformation
from ...proteomics.enzymes.digest_enzyme import DigestEnzyme
from .taxonomy_tree import TaxonomyTree
from .protein_digestion import ProteinDigestion
from .peptide_metadata_collector import PeptideMetadataCollector

@unique
class DatabaseStatus(Enum):
    DIGESTION = "digest proteins"
    METADATA_COLLECTION = "collect peptide metadata"
    READY = "ready for normal use"

class DatabaseMaintenance():
    TAXONOMY_FILES = ['nodes.dmp', 'names.dmp', 'merged.dmp', 'delete.dmp']

    def __init__(self, database_url: str, work_dir: pathlib.Path, number_of_threads: int, statistics_write_period: int, enzyme_name: str, maximum_number_of_missed_cleavages: int, minimum_peptide_length: int, maximum_peptide_length: int):
        self.__database_url = database_url
        self.__work_dir = work_dir
        self.__number_of_threads = number_of_threads
        self.__statistics_write_period = statistics_write_period
        self.__enzyme_name = enzyme_name
        self.__maximum_number_of_missed_cleavages = maximum_number_of_missed_cleavages
        self.__minimum_peptide_length = minimum_peptide_length
        self.__maximum_peptide_length = maximum_peptide_length

        self.__protein_data_path = work_dir.joinpath('protein_data/')
        self.__taxonomy_data_path = work_dir.joinpath('taxonomy_data/')
        self.__logs_path = work_dir.joinpath('logs/')
        self.__temporary_protein_data_path = self.__work_dir.joinpath('.temp_protein_data/')

    def start(self):
        """
        Maintain (create/update) the database.
        """
        self.__prepare_logs_directory()
        self.__validate_folders()
        self.__maintain_taxonomy_tree()
        self.__maintain_proteins_digestion()
    
    def __validate_folders(self):
        """
        Checks if the necessary folders are present to maintain the database.
        """
        if not self.__protein_data_path.exists():
            raise AssertionError(f"'{self.__protein_data_path}' does not exists.")
        if not self.__protein_data_path.is_dir():
            raise AssertionError(f"'{self.__protein_data_path}' is not a directory.")

        if not self.__taxonomy_data_path.exists():
            raise AssertionError(f"'{self.__taxonomy_data_path}' does not exists.")
        if not self.__taxonomy_data_path.is_dir():
            raise AssertionError(f"'{self.__taxonomy_data_path}' is not a directory.")

    def __prepare_logs_directory(self) -> pathlib.Path:
        """
        Creates the logs folder in the work directory if necessary.
        If the log directory already exists and contains some files/folders, it will be renamed to `log_<index>` (increasing index). And a new empty folder is created.
        """
        if self.__logs_path.is_dir():
            log_path_content = [path for path in self.__logs_path.iterdir()]
            if len(log_path_content):
                log_rotate_idx = 0
                log_rotate_path = self.__logs_path.parent.joinpath(f"logs_0")
                while log_rotate_path.is_dir():
                    log_rotate_idx += 1
                    log_rotate_path = self.__logs_path.parent.joinpath(f"logs_{log_rotate_idx}")
                self.__logs_path.rename(log_rotate_path)
        self.__logs_path.mkdir(exist_ok=True)



    def __maintain_taxonomy_tree(self):
        """
        Builds/Updates the taxonmy tree
        """
        tree = TaxonomyTree(
            self.__taxonomy_data_path.joinpath('nodes.dmp'),
            self.__taxonomy_data_path.joinpath('names.dmp'),
            self.__taxonomy_data_path.joinpath('merged.dmp'),
            self.__taxonomy_data_path.joinpath('delnodes.dmp'),
            self.__database_url,
            self.__logs_path.joinpath('taxonomy.log'),
            self.__number_of_threads
        )
        tree.maintain()

    def __maintain_proteins_digestion(self):
        """
        Build/update the protein digestion
        """
        self.__set_databse_status(self.__database_url, DatabaseStatus.DIGESTION, True)
        # Run 0, first run with user supplied data in protein_data directory
        run_count = 0
        protein_data_path = self.__protein_data_path
        number_of_threads = self.__number_of_threads
        while True:
            digestion = ProteinDigestion(
                protein_data_path,
                self.__logs_path,
                self.__statistics_write_period,
                number_of_threads,
                self.__enzyme_name,
                self.__maximum_number_of_missed_cleavages,
                self.__minimum_peptide_length,
                self.__maximum_peptide_length,
                run_count
            )
            error_count = digestion.run(self.__database_url)

            # If not errors occure break while-loop
            if not error_count:
                break
            
            # If loop continues, there were errors. Create temporary workdir
            self.__temporary_protein_data_path.mkdir(parents=True, exist_ok=True)
            # Remove files in temporary workdir
            for file_path in self.__temporary_protein_data_path.glob('*.*'):
                if file_path.is_file():
                    file_path.unlink()

            # Copy the unprocessible proteins to workdir
            shutil.copy(
                self.__logs_path.joinpath(f"unprocessible_proteins_{run_count}.txt"),
                self.__temporary_protein_data_path.joinpath(f"unprocessible_proteins_{run_count}.txt")
            )

            # Increase the run counter
            run_count += 1
            # Set protein data path to temporary workdir
            protein_data_path = self.__temporary_protein_data_path
            # Decrease number of threads to one third
            number_of_threads = int(self.__number_of_threads / 3)
            # Make sure it is at least 1
            if number_of_threads == 0:
                number_of_threads = 1


            print(f"Digest ended with {error_count} errors. Digest the remaining proteins with {number_of_threads} threads.")

        self.__set_databse_status(self.__database_url, DatabaseStatus.METADATA_COLLECTION, False)

        # Collect peptide meta data
        collector = PeptideMetadataCollector(
            self.__logs_path,
            # Increase the number of threads by one third. The peptide updates are independent from one another, so there are no deadlocks to be expected.
            self.__statistics_write_period,
            int(self.__number_of_threads * 1.33)
        )
        collector.run(self.__database_url)

        now = datetime.datetime.utcnow()
        epoch = datetime.datetime(1970,1,1)
        last_update_timestamp = (now - epoch).total_seconds()
        self.__set_databse_status(self.__database_url, DatabaseStatus.READY, False, last_update_timestamp)

        # Cleanup by removing all files in the temporary work directory
        if self.__temporary_protein_data_path.is_dir():
            shutil.rmtree(self.__temporary_protein_data_path, ignore_errors=True)


    def __set_databse_status(self, database_url: str, status: DatabaseStatus, maintenance_mode: bool, last_update_timestamp: int = None):
        """
        Sets the database status
        @param status Indicate if set to maintenance mode
        @param last_update_timestamp UTC timestamp in seconds
        """
        database_connection = psycopg2.connect(database_url)
        with database_connection.cursor() as database_cursor:
            database_cursor.execute("SELECT values FROM maintenance_information WHERE key = %s;", (MaintenanceInformation.DATABASE_STATUS_KEY,))
            database_status_row = database_cursor.fetchone()
            if database_status_row:
                database_status_values = database_status_row[0]
                database_status_values['maintenance_mode'] = maintenance_mode
                database_status_values['status'] = status.value
                if last_update_timestamp:
                    database_status_values['last_update'] = last_update_timestamp
                database_cursor.execute("UPDATE maintenance_information SET values = %s WHERE key = %s;", (json.dumps(database_status_values), MaintenanceInformation.DATABASE_STATUS_KEY))
            else:
                database_status_values = {
                    'status': status.value,
                    'maintenance_mode': maintenance_mode 
                }
                if last_update_timestamp:
                    database_status_values['last_update'] = last_update_timestamp
                else:
                    database_status_values['last_update'] = 0
                database_cursor.execute("INSERT INTO maintenance_information (key, values) VALUES (%s, %s);", (MaintenanceInformation.DATABASE_STATUS_KEY, json.dumps(database_status_values)))
            database_connection.commit()
        database_connection.close()


    @classmethod
    def maintain_from_command_line(cls, args):
        maintenance = cls(
            args.database_url,
            pathlib.Path(args.work_dir),
            args.thread_count,
            args.statistics_write_period,
            args.enzyme_name,
            args.maximum_number_of_missed_cleavages,
            args.minimum_peptide_length,
            args.maximum_peptide_length,
        )
        maintenance.start()

    @classmethod
    def comand_line_arguments(cls, subparsers):
        parser = subparsers.add_parser('database', help="Create and update a TrypperDB.")
        parser.add_argument("--work-dir", "-w", type=str, required=True, help="Protein file (can be used multiple times)")
        parser.add_argument("--statistics-write-period", type=int, help="Seconds between writes to the statistics file (default: 900)", default=900)
        parser.add_argument("--thread-count", "-t", type=int, required=True, help="Number of concurrent digestions and inserts to the database")
        parser.add_argument("--enzyme-name", "-e", type=str, required=True, help="The name of the enzyme", choices=DigestEnzyme.get_known_enzymes())
        parser.add_argument("--maximum-number-of-missed-cleavages", "-c", type=int, help="Maximum number of missed cleavages (default: 2)", default="2")
        parser.add_argument("--minimum-peptide-length", type=int, help="Minimum peptide length (default: 5)", default="5")
        parser.add_argument("--maximum-peptide-length", type=int, help="Maximum peptide length (default: 60)", default="60")
        parser.add_argument("--database-url", "-d", type=str, required=True, help="Database url for postgres")
        parser.set_defaults(func=cls.maintain_from_command_line)
