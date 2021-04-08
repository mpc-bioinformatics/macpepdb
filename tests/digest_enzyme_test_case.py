import unittest

from macpepdb.proteomics.enzymes.digest_enzyme import DigestEnzyme

class DigestEnzymeTestCase(unittest.TestCase):
    def test_differentiate_ambigous_sequences(self):
        AMBIGOUS_SEQUENCE = "MDQZTLABBQQILASLZPSR"
        EXPECTED_DIFFERENTIATED_SEQEUNCES = [
            "MDQETLANDQQILASLQPSR",
            "MDQETLANDQQILASLEPSR",
            "MDQQTLADNQQILASLQPSR",
            "MDQQTLANDQQILASLQPSR",
            "MDQQTLANNQQILASLEPSR",
            "MDQETLADNQQILASLEPSR",
            "MDQETLANNQQILASLQPSR",
            "MDQQTLANDQQILASLEPSR",
            "MDQETLADNQQILASLQPSR",
            "MDQQTLADDQQILASLQPSR",
            "MDQQTLADDQQILASLEPSR",
            "MDQQTLANNQQILASLQPSR",
            "MDQQTLADNQQILASLEPSR",
            "MDQETLANNQQILASLEPSR",
            "MDQETLADDQQILASLEPSR",
            "MDQETLADDQQILASLQPSR"
        ]

        differentiated_sequences = DigestEnzyme.differentiate_ambigous_sequences(AMBIGOUS_SEQUENCE)
        for sequence in differentiated_sequences:
            self.assertIn(sequence, EXPECTED_DIFFERENTIATED_SEQEUNCES)


    def test_is_sequence_containing_replaceable_ambigous_amino_acids(self):
        AMBIGOUS_SEQUENCE = "MDQZTLABBQQILASLZPSR"
        AMBIGOUS_SEQUENCE_WITH_B = "MDQTLABBQQILASLPSR"
        AMBIGOUS_SEQUENCE_WITH_Z = "MDQZTLAQQILASLZPSR"
        UNAMBIGOUS_SEQUENCE = "MDQTLAQQILASLPSR"

        self.assertTrue(DigestEnzyme.is_sequence_containing_replaceable_ambigous_amino_acids(AMBIGOUS_SEQUENCE_WITH_B))
        self.assertTrue(DigestEnzyme.is_sequence_containing_replaceable_ambigous_amino_acids(AMBIGOUS_SEQUENCE_WITH_Z))
        self.assertTrue(DigestEnzyme.is_sequence_containing_replaceable_ambigous_amino_acids(AMBIGOUS_SEQUENCE))
        self.assertFalse(DigestEnzyme.is_sequence_containing_replaceable_ambigous_amino_acids(UNAMBIGOUS_SEQUENCE));

