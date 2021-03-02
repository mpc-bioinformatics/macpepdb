from macpepdb.models.peptide import Peptide
from macpepdb.proteomics.mass.convert import to_float as mass_to_float

# Part (20 - 60) of Leptin (UniProt  accession: Q257X2)
LEPTIN_SEQUENCE = "AVPIRKVQDDTKTLIKTIVTRINDISHTQSVSSKQRVTGLDFIPGLHPLLSLSKMDQTLA"
# Reference calculated with https://web.expasy.org/compute_pi/
LEPTIN_SEQUENCE_WEIGHT = 6622.66

class PeptideTestCase:
    def test_termini(self):
        leptin = Peptide(LEPTIN_SEQUENCE, 0)
        self.assertEqual("A", leptin.n_terminus)
        self.assertEqual("A", leptin.c_terminus)

    def test_amino_acid_counts(self):
        leptin = Peptide(LEPTIN_SEQUENCE, 0)
        self.assertEqual(leptin.a_count, 2)
        self.assertEqual(leptin.c_count, 0)
        self.assertEqual(leptin.d_count, 5)
        self.assertEqual(leptin.e_count, 0)
        self.assertEqual(leptin.f_count, 1)
        self.assertEqual(leptin.g_count, 2)
        self.assertEqual(leptin.h_count, 2)
        self.assertEqual(leptin.i_count, 6)
        self.assertEqual(leptin.k_count, 5)
        self.assertEqual(leptin.l_count, 7)
        self.assertEqual(leptin.m_count, 1)
        self.assertEqual(leptin.n_count, 1)
        self.assertEqual(leptin.o_count, 0)
        self.assertEqual(leptin.p_count, 3)
        self.assertEqual(leptin.q_count, 4)
        self.assertEqual(leptin.r_count, 3)
        self.assertEqual(leptin.s_count, 6)
        self.assertEqual(leptin.t_count, 7)
        self.assertEqual(leptin.u_count, 0)
        self.assertEqual(leptin.v_count, 5)
        self.assertEqual(leptin.w_count, 0)
        self.assertEqual(leptin.y_count, 0)

        amino_acid_sum = leptin.a_count \
            + leptin.c_count \
            + leptin.d_count \
            + leptin.e_count \
            + leptin.f_count \
            + leptin.g_count \
            + leptin.h_count \
            + leptin.i_count \
            + leptin.k_count \
            + leptin.l_count \
            + leptin.m_count \
            + leptin.n_count \
            + leptin.o_count \
            + leptin.p_count \
            + leptin.q_count \
            + leptin.r_count \
            + leptin.s_count \
            + leptin.t_count \
            + leptin.u_count \
            + leptin.v_count \
            + leptin.w_count \
            + leptin.y_count

        self.assertEqual(len(leptin.sequence), amino_acid_sum)

    def test_weight_calculation(self):
        peptide = Peptide(LEPTIN_SEQUENCE, 0)
        self.assertEqual(LEPTIN_SEQUENCE_WEIGHT, round(mass_to_float(peptide.weight), 2))

