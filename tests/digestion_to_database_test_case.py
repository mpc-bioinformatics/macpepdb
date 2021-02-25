import pathlib
import os

from sqlalchemy import func, distinct
from sqlalchemy.sql import exists

from macpepdb.tasks.digestion import Digestion
from macpepdb.proteomics.enzymes.digest_enzyme import DigestEnzyme
from macpepdb.proteomics.file_reader.uniprot_text_reader import UniprotTextReader
from macpepdb.models.peptide import Peptide
from macpepdb.models.protein import Protein
from macpepdb.models.protein_merge import ProteinMerge
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
        protein_merges = []
        peptides = set()

        with test_file_path.open("r") as test_file:
            protein_file_reader = UniprotTextReader(test_file)

            # Redirects not needed, so use _
            for protein, read_protein_merges in protein_file_reader:
                for peptide in trypsin.digest(protein):
                    peptides.add(peptide)
                for merge in read_protein_merges:
                    protein_merges.append(merge)
                proteins.append(protein)

        session = self.session_factory()

        # Check if protein count in database are equals to set
        self.assertEqual(len(proteins), session.query(func.count(distinct(Protein.accession))).scalar())
        
        # Check if set count is equals db count
        self.assertEqual(len(peptides), session.query(func.count(distinct(Peptide.sequence))).scalar())
        # Check if all peptides from Set exists in database
        for peptide in peptides:
            self.assertTrue(session.query(exists().where(Peptide.sequence == peptide.sequence)).scalar())
        # Check if all peptides in database are present in Set
        db_peptides = session.query(Peptide).all()
        for peptide in db_peptides:
            self.assertTrue(peptide in peptides)

        # Check if protein merge count in database are equals to set
        self.assertEqual(len(protein_merges), session.query(func.count(ProteinMerge.source_accession)).scalar())
        # Check if all protein merges from set exists in database
        for merge in protein_merges:
            self.assertTrue(session.query(exists().where(ProteinMerge.source_accession == merge.source_accession and ProteinMerge.target_accession == merge.target_accession)).scalar())
        # Check if protein merges in db are present in the set
        db_protein_merges = session.query(ProteinMerge).all()
        for merge in db_protein_merges:
            self.assertTrue(merge in protein_merges)

        # Check if maintenance mode is false and update timestamp is greater zero
        database_status = session.query(MaintenanceInformation).filter(MaintenanceInformation.key == MaintenanceInformation.DATABASE_STATUS_KEY).one_or_none()
        self.assertNotEqual(database_status, None)
        self.assertGreater(database_status.values['last_update'], 0)
        self.assertFalse(database_status.values['maintenance_mode'])
        self.tearDown()

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
        read_protein_merges = []
        read_protein_peptide_sequences = {}
        with test_file_path.open("r") as test_file:
            file_reader = UniprotTextReader(test_file)

            # Redirects not needed, so use _
            for protein, protein_merges in file_reader:
                read_protein = protein
                read_protein_merges = protein_merges
        read_protein_peptide_sequences = {peptide.sequence for peptide in trypsin.digest(read_protein)}


        session = self.session_factory()
        # Fetch peptide from database
        database_protein = session.query(Protein).filter(Protein.accession == read_protein.accession).one_or_none()

        # Test if peptide exists
        self.assertNotEqual(database_protein, None)

        # Test protein attributes
        self.assertEqual(read_protein.accession, database_protein.accession)
        self.assertEqual(read_protein.taxonomy_id, database_protein.taxonomy_id)
        self.assertEqual(read_protein.proteome_id, database_protein.proteome_id)
        self.assertEqual(read_protein.sequence, database_protein.sequence)

        # Test peptides
        databse_protein_peptide_sequences = {peptide.sequence for peptide in database_protein.peptides.all()}
        for sequence in read_protein_peptide_sequences:
            self.assertIn(sequence, databse_protein_peptide_sequences)

        for sequence in databse_protein_peptide_sequences:
            self.assertIn(sequence, read_protein_peptide_sequences)

        # Test protein merges
        database_protein_merges = session.query(ProteinMerge).filter(ProteinMerge.target_accession == read_protein.accession).all()

        # Convert to sets of source accession
        read_protein_merges_source_accessions = {merge.source_accession for merge in read_protein_merges}
        database_protein_merges_source_accessions = {merge.source_accession for merge in database_protein_merges}

        # Crosscheck if each accession from the on set is in the other
        for accession in read_protein_merges_source_accessions:
            self.assertIn(accession, database_protein_merges_source_accessions)

        for accession in database_protein_merges_source_accessions:
            self.assertIn(accession, read_protein_merges_source_accessions)