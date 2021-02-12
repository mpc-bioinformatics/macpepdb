import unittest
import pathlib

from macpep_db.proteomics.file_reader.uniprot_text_reader import UniprotTextReader

class UniprotTextReaderTestCase(unittest.TestCase):
    def test_reading(self):
        test_file_path = pathlib.Path("./test_files/UP000002006.txt")

        proteins = []

        with test_file_path.open("r") as test_file:
            reader = UniprotTextReader(test_file)

            for protein in reader:
                proteins.append(protein)
        
        with test_file_path.open() as test_file:
            test_file_plain_content = test_file.read()

        id_lines = test_file_plain_content.count("ID   ")

        # Make sure all proteins are read
        self.assertEqual(id_lines, len(proteins))