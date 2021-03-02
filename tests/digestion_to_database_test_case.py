import pathlib
import os

from sqlalchemy import func, distinct
from sqlalchemy.sql import exists

from macpepdb.tasks.digestion import Digestion
from macpepdb.proteomics.enzymes.digest_enzyme import DigestEnzyme
from macpepdb.proteomics.file_reader.uniprot_text_reader import UniprotTextReader
from macpepdb.models.peptide import Peptide
from macpepdb.models.protein import Protein
from macpepdb.models.maintenance_information import MaintenanceInformation

from .abstract_database_test_case import AbstractDatabaseTestCase

TRYPSIN_MAX_MISSED_CLEAVAGES = 2
TRYPSIN_MIN_PEPTIDE_LENGTH = 5
TRYPSIN_MAX_PEPTIDE_LENGTH = 40

class DigestionToDatabaseTestCase(AbstractDatabaseTestCase):
    def test_digestion(self):
        test_file_path = pathlib.Path('./test_files/UP000002006.txt')

        digestion = Digestion(
            [test_file_path],
            pathlib.Path("./digestion_test.log"),
            pathlib.Path("./digestion_test.unprocessible_proteins.txt"),
            pathlib.Path("./digestion_test.statistics.csv"),
            5,
            4,
            "Trypsin",
            TRYPSIN_MAX_MISSED_CLEAVAGES,
            TRYPSIN_MIN_PEPTIDE_LENGTH,
            TRYPSIN_MAX_PEPTIDE_LENGTH
        )

        digestion.digest_to_database(os.getenv("TEST_MACPEPDB_URL"))

        EnzymeClass = DigestEnzyme.get_enzyme_by_name("Trypsin")
        trypsin = EnzymeClass(TRYPSIN_MAX_MISSED_CLEAVAGES, TRYPSIN_MIN_PEPTIDE_LENGTH, TRYPSIN_MAX_PEPTIDE_LENGTH)

        proteins = []
        peptides = set()

        with test_file_path.open("r") as test_file:
            protein_file_reader = UniprotTextReader(test_file)

            for protein in protein_file_reader:
                for peptide in trypsin.digest(protein):
                    peptides.add(peptide)
                proteins.append(protein)

        with self.database_connection:
            with self.database_connection.cursor() as database_cursor:

                # Check if protein count in database are equals to set
                database_cursor.execute(f"SELECT count(*) FROM {Protein.TABLE_NAME};")
                self.assertEqual(len(proteins), database_cursor.fetchone()[0])
                
                # Check if set count is equals db count
                database_cursor.execute(f"SELECT count(*) FROM {Peptide.TABLE_NAME};")
                self.assertEqual(len(peptides), database_cursor.fetchone()[0])
                # Check if all peptides from Set exists in database
                for peptide in peptides:
                    database_cursor.execute(f"SELECT true FROM {Peptide.TABLE_NAME} WHERE sequence = %s;", (peptide.sequence,))
                    self.assertTrue(database_cursor.fetchone()[0])
                # Check if all peptides in database are present in Set
                for peptide in Peptide.select(database_cursor, fetchall=True):
                    self.assertTrue(peptide in peptides)

                # Check if maintenance mode is false and update timestamp is greater zero
                database_status = MaintenanceInformation.select(database_cursor, MaintenanceInformation.DATABASE_STATUS_KEY)
                self.assertNotEqual(database_status, None)
                self.assertGreater(database_status.values['last_update'], 0)
                self.assertFalse(database_status.values['maintenance_mode'])

    def test_digest_with_protein(self):
        # Digest the databse
        test_file_path = pathlib.Path("./test_files/UP000002006.txt")

        digestion = Digestion(
            [test_file_path],
            pathlib.Path("./digestion_test.log"),
            pathlib.Path("./digestion_test.unprocessible_proteins.txt"),
            pathlib.Path("./digestion_test.statistics.csv"),
            5,
            4,
            "Trypsin",
            TRYPSIN_MAX_MISSED_CLEAVAGES,
            TRYPSIN_MIN_PEPTIDE_LENGTH,
            TRYPSIN_MAX_PEPTIDE_LENGTH
        )

        digestion.digest_to_database(os.getenv("TEST_MACPEPDB_URL"))

        # Digest the modified B0FIH3. Sequence, accession, taxonomy id and proteome id are updated.
        test_file_path = pathlib.Path("./test_files/B0FIH3_updated.txt")

        # Digest with different digest values to see if the digestion uses the values from the database
        digestion = Digestion(
            [test_file_path],
            pathlib.Path("./digestion_test.log"),
            pathlib.Path("./digestion_test.unprocessible_proteins.txt"),
            pathlib.Path("./digestion_test.statistics.csv"),
            5,
            4,
            "Trypsin",
            8,
            6,
            50
        )

        digestion.digest_to_database(os.getenv("TEST_MACPEPDB_URL"))

        EnzymeClass = DigestEnzyme.get_enzyme_by_name("Trypsin")
        trypsin = EnzymeClass(TRYPSIN_MAX_MISSED_CLEAVAGES, TRYPSIN_MIN_PEPTIDE_LENGTH, TRYPSIN_MAX_PEPTIDE_LENGTH)

        read_protein = None
        read_protein_peptide_sequences = {}
        with test_file_path.open("r") as test_file:
            file_reader = UniprotTextReader(test_file)

            for protein in file_reader:
                read_protein = protein
        read_protein_peptide_sequences = {peptide.sequence for peptide in trypsin.digest(read_protein)}


        with self.database_connection:
            with self.database_connection.cursor() as database_cursor:
                # Fetch peptide from database
                database_protein = Protein.select(database_cursor, ("accession = %s", [read_protein.accession]))

                # Test if peptide exists
                self.assertNotEqual(database_protein, None)

                # Test protein attributes
                self.assertEqual(read_protein.accession, database_protein.accession)
                self.assertEqual(read_protein.taxonomy_id, database_protein.taxonomy_id)
                self.assertEqual(read_protein.proteome_id, database_protein.proteome_id)
                self.assertEqual(read_protein.sequence, database_protein.sequence)

                # Test peptides
                databse_protein_peptide_sequences = {peptide.sequence for peptide in database_protein.peptides(database_cursor)}
                for sequence in read_protein_peptide_sequences:
                    self.assertIn(sequence, databse_protein_peptide_sequences)

                for sequence in databse_protein_peptide_sequences:
                    self.assertIn(sequence, read_protein_peptide_sequences)

                # Test secondary accessions
                # Convert to sets of source accession
                read_protein_secondary_accessions = set(read_protein.secondary_accessions)
                database_protein_secondary_accessions = set(database_protein.secondary_accessions)

                # Crosscheck if each accession from the on set is in the other
                for accession in read_protein_secondary_accessions:
                    self.assertIn(accession, database_protein_secondary_accessions)

                for accession in database_protein_secondary_accessions:
                    self.assertIn(accession, read_protein_secondary_accessions)