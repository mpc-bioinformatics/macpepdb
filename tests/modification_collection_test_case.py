# std imports
import unittest
import tempfile
import pathlib

# internal imports
from macpepdb.proteomics.modification_collection import ModificationCollection, ModificationLimitError, InvalidModificationCombinationError

CORRECT = """
"accession","name","aa","delta","static","position"
"unimod:35","Oxidation of M","M",15.994915,"variable","anywhere"
"unimod:4","Carbamidomethyl of C","C",57.021464,"static","anywhere"
"""

TOO_MANY_VARIABLE_MODIFICATIONS = """
"accession","name","aa","delta","static","position"
"fakemod:1","Fake of P","P",1,"variable","anywhere"
"fakemod:2","Fake of C","C",2,"variable","c_terminus"
"fakemod:3","Fake of V","V",3,"variable","n_terminus"
"fakemod:4","Fake of K","K",4,"variable","anywhere"
"fakemod:1","Fake of Q","Q",10,"static","n_terminus"
"fakemod:5","Fake of R","R",5,"variable","anywhere"
"fakemod:6","Fake of E","F",6,"variable","c_terminus"
"fakemod:7","Fake of O","O",7,"variable","anywhere"
"fakemod:8","Fake of S","S",8,"variable","anywhere"
"fakemod:9","Fake of U","U",9,"variable","anywhere"
"fakemod:10","Fake of W","W",10,"variable","n_terminus"
"""

TOO_MANY_N_TERMINUS_MODIFICATIONS = """
"accession","name","aa","delta","static","position"
"unimod:35","Oxidation of M","M",15.994915,"variable","anywhere"
"fakemod:1","Fake of C","C",10,"static","n_terminus"
"fakemod:2","Fake of P","P",20,"static","n_terminus"
"""

TOO_MANY_C_TERMINUS_MODIFICATIONS = """
"accession","name","aa","delta","static","position"
"unimod:35","Oxidation of M","M",15.994915,"variable","anywhere"
"fakemod:1","Fake of C","C",10,"static","c_terminus"
"fakemod:2","Fake of P","P",20,"static","c_terminus"
"""

STATIC_AND_VARIABLE_MODIFICATION_FOR_SAME_AMINO_ACID = """
"accession","name","aa","delta","static","position"
"unimod:35","Oxidation of M","M",15.994915,"variable","anywhere"
"fakemod:1","Static fake of C","C",10,"static","anywhere"
"fakemod:2","Variable fake of C","C",20,"variable","anywhere"
"""


class ModificationCollectionTestCase(unittest.TestCase):
    def test_creation(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            modification_file_path = pathlib.Path(temp_dir).joinpath("modifications.csv")

            with modification_file_path.open("w") as modification_file:

                modification_file.write(CORRECT.strip())
                modification_file.flush()

            ModificationCollection.read_from_csv_file(modification_file_path)


            for csv_content in [ TOO_MANY_VARIABLE_MODIFICATIONS, TOO_MANY_N_TERMINUS_MODIFICATIONS, TOO_MANY_C_TERMINUS_MODIFICATIONS ]:
                with modification_file_path.open("w") as modification_file:
                    modification_file.write(csv_content.strip())
                    modification_file.flush()
                self.assertRaises(ModificationLimitError, ModificationCollection.read_from_csv_file, modification_file_path)
        
            with modification_file_path.open("w") as modification_file:
                modification_file.write(STATIC_AND_VARIABLE_MODIFICATION_FOR_SAME_AMINO_ACID.strip())
                modification_file.flush()
            self.assertRaises(InvalidModificationCombinationError, ModificationCollection.read_from_csv_file, modification_file_path)