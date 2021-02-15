import unittest
import pathlib

from macpepdb.proteomics.modification import Modification, ModificationPosition
from macpepdb.proteomics.amino_acid import AminoAcid

class ModificationTestCase(unittest.TestCase):
    def test_csv_read(self):
        csv_file_path = pathlib.Path("./test_files/modifications.csv")

        modifications = Modification.read_from_csv_file(csv_file_path)

        self.assertEqual(2, len(modifications))

    def test_modification_positions(self):
        # Assert KeyError
        self.assertRaises(KeyError, ModificationPosition.from_string, "x_terminus")
        # Do not assert any problems here
        for position in ModificationPosition:
            position = ModificationPosition.from_string(str(position))