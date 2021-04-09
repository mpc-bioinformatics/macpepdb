import unittest
import pathlib
import re

from macpepdb.proteomics.file_reader.uniprot_text_reader import UniprotTextReader

class UniprotTextReaderTestCase(unittest.TestCase):
    def test_reading(self):
        test_file_path = pathlib.Path("./test_files/proteins.txt")

        proteins = []

        with test_file_path.open("r") as test_file:
            reader = UniprotTextReader(test_file)

            for protein in reader:
                proteins.append(protein)
        
        with test_file_path.open() as test_file:
            test_file_plain_content = test_file.read()

        id_line_regex = re.compile(r"^ID\W{3}", re.MULTILINE)
        id_line_matches = id_line_regex.findall(test_file_plain_content)

        # Make sure all proteins are read
        self.assertEqual(len(id_line_matches), len(proteins))