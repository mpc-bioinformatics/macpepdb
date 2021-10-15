# std imports
import unittest

# internal imports
from macpepdb.models.peptide import Peptide
from macpepdb.proteomics.amino_acid import AminoAcid
from macpepdb.proteomics.mass.convert import to_float as mass_to_float
from macpepdb.proteomics.neutral_loss import H2O

# Fictional sequence 
FICTIONAL_SEQUENCE = "AZVPIRKVJQODDEEWWTBKTLIKTUOIVCTRINDISHTJQSCCVSSKQRVTGJLDFYIPGLOHPLLBSLSKMCDWWOQTLA"

class PeptideTestCase(unittest.TestCase):
    def test_termini(self):
        fictional_peptide = Peptide(FICTIONAL_SEQUENCE, 0)
        self.assertEqual("A", fictional_peptide.n_terminus)
        self.assertEqual("A", fictional_peptide.c_terminus)

    def test_amino_acid_counts(self):
        fictional_peptide = Peptide(FICTIONAL_SEQUENCE, 0)
        self.assertEqual(fictional_peptide.a_count, 2)
        self.assertEqual(fictional_peptide.b_count, 2)
        self.assertEqual(fictional_peptide.c_count, 4)
        self.assertEqual(fictional_peptide.d_count, 5)
        self.assertEqual(fictional_peptide.e_count, 2)
        self.assertEqual(fictional_peptide.f_count, 1)
        self.assertEqual(fictional_peptide.g_count, 2)
        self.assertEqual(fictional_peptide.h_count, 2)
        self.assertEqual(fictional_peptide.i_count, 6)
        self.assertEqual(fictional_peptide.j_count, 3)
        self.assertEqual(fictional_peptide.k_count, 5)
        self.assertEqual(fictional_peptide.l_count, 7)
        self.assertEqual(fictional_peptide.m_count, 1)
        self.assertEqual(fictional_peptide.n_count, 1)
        self.assertEqual(fictional_peptide.o_count, 4)
        self.assertEqual(fictional_peptide.p_count, 3)
        self.assertEqual(fictional_peptide.q_count, 4)
        self.assertEqual(fictional_peptide.r_count, 3)
        self.assertEqual(fictional_peptide.s_count, 6)
        self.assertEqual(fictional_peptide.t_count, 7)
        self.assertEqual(fictional_peptide.u_count, 1)
        self.assertEqual(fictional_peptide.v_count, 5)
        self.assertEqual(fictional_peptide.w_count, 4)
        self.assertEqual(fictional_peptide.y_count, 1)
        self.assertEqual(fictional_peptide.z_count, 1)

        amino_acid_sum = fictional_peptide.a_count \
            + fictional_peptide.b_count \
            + fictional_peptide.c_count \
            + fictional_peptide.d_count \
            + fictional_peptide.e_count \
            + fictional_peptide.f_count \
            + fictional_peptide.g_count \
            + fictional_peptide.h_count \
            + fictional_peptide.i_count \
            + fictional_peptide.j_count \
            + fictional_peptide.k_count \
            + fictional_peptide.l_count \
            + fictional_peptide.m_count \
            + fictional_peptide.n_count \
            + fictional_peptide.o_count \
            + fictional_peptide.p_count \
            + fictional_peptide.q_count \
            + fictional_peptide.r_count \
            + fictional_peptide.s_count \
            + fictional_peptide.t_count \
            + fictional_peptide.u_count \
            + fictional_peptide.v_count \
            + fictional_peptide.w_count \
            + fictional_peptide.y_count \
            + fictional_peptide.z_count

        self.assertEqual(len(fictional_peptide.sequence), amino_acid_sum)

    def test_mass_calculation(self):
        fictional_peptide = Peptide(FICTIONAL_SEQUENCE, 0)

        # If the amino acid cound test passes, we can use the counts to calculate the mass manally.
        # Actually there is no external tool which supports all of our known amino acids, so we can double check the weigth.
        mass = fictional_peptide.a_count * AminoAcid.get_by_one_letter_code('A').mono_mass \
            + fictional_peptide.b_count * AminoAcid.get_by_one_letter_code('B').mono_mass \
            + fictional_peptide.c_count * AminoAcid.get_by_one_letter_code('C').mono_mass \
            + fictional_peptide.d_count * AminoAcid.get_by_one_letter_code('D').mono_mass \
            + fictional_peptide.e_count * AminoAcid.get_by_one_letter_code('E').mono_mass \
            + fictional_peptide.f_count * AminoAcid.get_by_one_letter_code('F').mono_mass \
            + fictional_peptide.g_count * AminoAcid.get_by_one_letter_code('G').mono_mass \
            + fictional_peptide.h_count * AminoAcid.get_by_one_letter_code('H').mono_mass \
            + fictional_peptide.i_count * AminoAcid.get_by_one_letter_code('I').mono_mass \
            + fictional_peptide.j_count * AminoAcid.get_by_one_letter_code('J').mono_mass \
            + fictional_peptide.k_count * AminoAcid.get_by_one_letter_code('K').mono_mass \
            + fictional_peptide.l_count * AminoAcid.get_by_one_letter_code('L').mono_mass \
            + fictional_peptide.m_count * AminoAcid.get_by_one_letter_code('M').mono_mass \
            + fictional_peptide.n_count * AminoAcid.get_by_one_letter_code('N').mono_mass \
            + fictional_peptide.o_count * AminoAcid.get_by_one_letter_code('O').mono_mass \
            + fictional_peptide.p_count * AminoAcid.get_by_one_letter_code('P').mono_mass \
            + fictional_peptide.q_count * AminoAcid.get_by_one_letter_code('Q').mono_mass \
            + fictional_peptide.r_count * AminoAcid.get_by_one_letter_code('R').mono_mass \
            + fictional_peptide.s_count * AminoAcid.get_by_one_letter_code('S').mono_mass \
            + fictional_peptide.t_count * AminoAcid.get_by_one_letter_code('T').mono_mass \
            + fictional_peptide.u_count * AminoAcid.get_by_one_letter_code('U').mono_mass \
            + fictional_peptide.v_count * AminoAcid.get_by_one_letter_code('V').mono_mass \
            + fictional_peptide.w_count * AminoAcid.get_by_one_letter_code('W').mono_mass \
            + fictional_peptide.y_count * AminoAcid.get_by_one_letter_code('Y').mono_mass \
            + fictional_peptide.z_count * AminoAcid.get_by_one_letter_code('Z').mono_mass \
            + H2O.mono_mass

        self.assertEqual(mass, fictional_peptide.mass)

