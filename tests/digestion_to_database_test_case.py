import pathlib
import os

from sqlalchemy import func, distinct
from sqlalchemy.sql import exists

from trypperdb.tasks.digestion import Digestion
from trypperdb.proteomics.enzymes.digest_enzyme import DigestEnzyme
from trypperdb.proteomics.file_reader.file_reader import FileReader
from trypperdb.models.peptide import Peptide
from trypperdb.models.protein import Protein
from trypperdb.models.protein_merge import ProteinMerge

from .abstract_database_test_case import AbstractDatabaseTestCase

TRYPSIN_MAX_MISSED_CLEAVAGES = 2
TRYPSIN_MIN_PEPTIDE_LENGTH = 5
TRYPSIN_MAX_PEPTIDE_LENGTH = 40

class DigestionToDatabaseTestCase(AbstractDatabaseTestCase):
    def test_digestion(self):
        test_files_and_formats = [
            ("./test_files/UP000002006.fasta", "fasta"),
            ("./test_files/UP000002006.txt", "text")
        ]
        for path, format in test_files_and_formats:

            test_file_path = pathlib.Path(path)

            digestion = Digestion(
                [test_file_path],
                format,
                pathlib.Path("./digestion_test.log"),
                pathlib.Path("./digestion_test.unprocessible_proteins.fasta"),
                pathlib.Path("./digestion_test.statistics.csv"),
                5,
                4,
                "Trypsin",
                TRYPSIN_MAX_MISSED_CLEAVAGES,
                TRYPSIN_MIN_PEPTIDE_LENGTH,
                TRYPSIN_MAX_PEPTIDE_LENGTH,
                0
            )

            digestion.digest_to_database(os.getenv("TRYPPERDB_TEST_DB_URL"))

            EnzymeClass = DigestEnzyme.get_enzyme_by_name("Trypsin")
            trypsin = EnzymeClass(TRYPSIN_MAX_MISSED_CLEAVAGES, TRYPSIN_MIN_PEPTIDE_LENGTH, TRYPSIN_MAX_PEPTIDE_LENGTH)

            proteins = []
            protein_merges = []
            peptides = set()

            with test_file_path.open("r") as test_file:
                FileReaderClass = FileReader.get_reader_by_file_format(format)
                file_reader = FileReaderClass(test_file)

                # Redirects not needed, so use _
                for protein, read_protein_merges in file_reader:
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
