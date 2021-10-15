# std imports
import unittest

# internal imports
from macpepdb.models.peptide import Peptide
from macpepdb.peptide_mass_validator import PeptideMassValidator
from macpepdb.proteomics.amino_acid import AminoAcid
from macpepdb.proteomics.mass.convert import to_int as mass_to_int
from macpepdb.proteomics.mass.precursor_range import PrecursorRange
from macpepdb.proteomics.modification import Modification, ModificationPosition
from macpepdb.proteomics.modification_collection import ModificationCollection

# Leptin peptide with 2 missed cleavages with a weigh of 5818137657950 (5818.137657950 Da)
LEPTIN_PEPTIDE_SEQUENCE = "DJJHJJAASKSCPJPQVRAJESJESJGVVJEASJYSTEVVAJSRJQGSJQDMJR"

class PeptideMassValidatorTestCase(unittest.TestCase):
    def test_validation(self):
        static_carbamidomethylation_of_c = Modification('unimod:4', 'carbamidomethylation of cysteine', AminoAcid.get_by_one_letter_code('C'), mass_to_int(57.021464), True, ModificationPosition.ANYWHERE)
        variable_oxidation_of_m = Modification('unimod:35', 'oxidation of methionine', AminoAcid.get_by_one_letter_code('M'), mass_to_int(15.994915), False, ModificationPosition.ANYWHERE)
        static_custom_modification_of_n_terminal_d = Modification('custom:1', 'custom of aspartic acid', AminoAcid.get_by_one_letter_code('D'), mass_to_int(10.01541), True, ModificationPosition.N_TERMINUS)
        variable_custom_modification_of_n_terminal_d = Modification('custom:2', 'custom of aspartic acid', AminoAcid.get_by_one_letter_code('D'), mass_to_int(10.01541), False, ModificationPosition.N_TERMINUS)
        static_custom_modification_of_c_terminal_r = Modification('custom:3', 'custom of arginine', AminoAcid.get_by_one_letter_code('R'), mass_to_int(6.153215), True, ModificationPosition.C_TERMINUS)
        variable_custom_modification_of_c_terminal_r = Modification('custom:4', 'custom of arginine', AminoAcid.get_by_one_letter_code('R'), mass_to_int(6.153215), False, ModificationPosition.C_TERMINUS)

        peptide = Peptide(LEPTIN_PEPTIDE_SEQUENCE, 2)

        # Static carbamidomethylation of C
        expected_peptide_mass = peptide.mass + peptide.c_count * static_carbamidomethylation_of_c.delta
        modification_collection = ModificationCollection([static_carbamidomethylation_of_c])
        precursor_range = PrecursorRange(expected_peptide_mass, 0, 0)
        validator = PeptideMassValidator(modification_collection, 0, precursor_range)
        self.assertTrue(validator.validate(peptide))

        # This should als match with allowed variable modification (where actually none is applied)
        # Static carbamidomethylation of C
        # Variable oxidation of M (not considered in expected_mass)
        modification_collection = ModificationCollection([static_carbamidomethylation_of_c, variable_oxidation_of_m])
        validator = PeptideMassValidator(modification_collection, 3, precursor_range)
        self.assertTrue(validator.validate(peptide))


        # Static carbamidomethylation of C
        # 1 variable oxidation of M
        expected_peptide_mass = peptide.mass \
            + peptide.c_count * static_carbamidomethylation_of_c.delta \
            + 1 * variable_oxidation_of_m.delta
        modification_collection = ModificationCollection([static_carbamidomethylation_of_c, variable_oxidation_of_m])
        precursor_range = PrecursorRange(expected_peptide_mass, 0, 0)
        validator = PeptideMassValidator(modification_collection, 3, precursor_range)
        self.assertTrue(validator.validate(peptide))

        # This should not match if no variable modifiations are allowed
        # Static carbamidomethylation of C
        # Variable oxidation of M (considered in expected_mass but no variable modification allowed in validation)
        validator.set_maximum_number_of_variable_modifications(0)
        self.assertFalse(validator.validate(peptide))

        # Lets replace two Js with Ms and test 3 applied variable oxidations of M
        # Static carbamidomethylation of C
        # 3 Variable oxidation of M
        peptide = Peptide(LEPTIN_PEPTIDE_SEQUENCE.replace('J', 'M', 2), 2)
        expected_peptide_mass = peptide.mass \
            + peptide.c_count * static_carbamidomethylation_of_c.delta \
            + 3 * variable_oxidation_of_m.delta
        modification_collection = ModificationCollection([static_carbamidomethylation_of_c, variable_oxidation_of_m])
        precursor_range = PrecursorRange(expected_peptide_mass, 0, 0)
        validator = PeptideMassValidator(modification_collection, 3, precursor_range)
        self.assertTrue(validator.validate(peptide))

        # This should fail with only 2 allowed variable modifications
        validator.set_maximum_number_of_variable_modifications(2)
        self.assertFalse(validator.validate(peptide))


        # Test variable n-terminal
        # Variable n-terminal modification of D
        # Static carbamidomethylation of C
        # 2 variable oxidation of M
        expected_peptide_mass = peptide.mass \
            + variable_custom_modification_of_n_terminal_d.delta \
            + peptide.c_count * static_carbamidomethylation_of_c.delta \
            + 2 * variable_oxidation_of_m.delta
        modification_collection = ModificationCollection([static_carbamidomethylation_of_c, variable_oxidation_of_m, variable_custom_modification_of_n_terminal_d])
        precursor_range = PrecursorRange(expected_peptide_mass, 0, 0)
        validator = PeptideMassValidator(modification_collection, 3, precursor_range)
        self.assertTrue(validator.validate(peptide))

        # This should fail with only 2 allowed variable modifications
        validator.set_maximum_number_of_variable_modifications(2)
        self.assertFalse(validator.validate(peptide))


        # Test static n-terminal modification
        # Static n-terminal modification of D
        # Static carbamidomethylation of C
        # 2 variable oxidation of M
        expected_peptide_mass = peptide.mass \
            + static_custom_modification_of_n_terminal_d.delta \
            + peptide.c_count * static_carbamidomethylation_of_c.delta \
            + 2 * variable_oxidation_of_m.delta
        modification_collection = ModificationCollection([static_carbamidomethylation_of_c, variable_oxidation_of_m, static_custom_modification_of_n_terminal_d])
        precursor_range = PrecursorRange(expected_peptide_mass, 0, 0)
        validator = PeptideMassValidator(modification_collection, 3, precursor_range)
        self.assertTrue(validator.validate(peptide))

        # Test variable n-terminal
        # Variable c-terminal modification of R
        # Static carbamidomethylation of C
        # 2 variable oxidation of M
        expected_peptide_mass = peptide.mass \
            + variable_custom_modification_of_c_terminal_r.delta \
            + peptide.c_count * static_carbamidomethylation_of_c.delta \
            + 2 * variable_oxidation_of_m.delta
        modification_collection = ModificationCollection([static_carbamidomethylation_of_c, variable_oxidation_of_m, variable_custom_modification_of_c_terminal_r])
        precursor_range = PrecursorRange(expected_peptide_mass, 0, 0)
        validator = PeptideMassValidator(modification_collection, 3, precursor_range)
        self.assertTrue(validator.validate(peptide))

        # This should fail with only 2 allowed variable modifications
        validator.set_maximum_number_of_variable_modifications(2)
        self.assertFalse(validator.validate(peptide))


        # Test static n-terminal modification
        # Static c-terminal modification of R
        # Static carbamidomethylation of C
        # 2 variable oxidation of M
        expected_peptide_mass = peptide.mass \
            + static_custom_modification_of_c_terminal_r.delta \
            + peptide.c_count * static_carbamidomethylation_of_c.delta \
            + 2 * variable_oxidation_of_m.delta
        modification_collection = ModificationCollection([static_carbamidomethylation_of_c, variable_oxidation_of_m, static_custom_modification_of_c_terminal_r])
        precursor_range = PrecursorRange(expected_peptide_mass, 0, 0)
        validator = PeptideMassValidator(modification_collection, 3, precursor_range)
        self.assertTrue(validator.validate(peptide))