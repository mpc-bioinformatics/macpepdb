import unittest
import pathlib
from trypperdb.proteomics.file_reader.fasta_reader import FastaReader

class FastaReaderTestCase(unittest.TestCase):
    def test_reading(self):
        fasta_path = pathlib.Path("./test_files/UP000002006.fasta")
        with fasta_path.open("r") as fasta_file:
            fasta = FastaReader(fasta_file)
            proteins = []
            for protein in fasta:
                proteins.append(protein)
        
        with fasta_path.open() as fasta_file:
            fasta_plain_content = fasta_file.read()

        number_of_headers = fasta_plain_content.count(">")

        # Make sure all proteins are read
        self.assertEqual(number_of_headers, len(proteins))