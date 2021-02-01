from macpepdb.models.peptide import Peptide
from macpepdb.proteomics.mass.convert import to_float as mass_to_float

from .abstract_database_test_case import AbstractDatabaseTestCase

# Part (20 - 60) of Leptin (UniProt  accession: Q257X2)
LEPTIN_SEQUENCE = "AVPIRKVQDDTKTLIKTIVTRINDISHTQSVSSKQRVTGLDFIPGLHPLLSLSKMDQTLA"
# Reference calculated with https://web.expasy.org/compute_pi/
LEPTIN_SEQUENCE_WEIGHT = 6622.66

class PeptideTestCase(AbstractDatabaseTestCase):
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
    
    def test_create_read_update_delete_cycle(self):
        leptin = Peptide(LEPTIN_SEQUENCE, 0)

        ## create
        # start db session
        session = self.session_factory()
        session.add(leptin)
        session.commit()
        # check if id is now an integer (autoincrement from db)
        self.assertTrue(isinstance(leptin.id, int))
        # save autoincremented id
        LEPTIN_ID = leptin.id
        # close session and set Leptin to None, so all connections to db are lost
        session.close()
        leptin = None

        ## read
        # start new db session
        session = self.session_factory()
        # get Leptin-peptide by sequence
        leptin = session.query(Peptide).filter(Peptide.sequence == LEPTIN_SEQUENCE).one()
        self.assertEqual(LEPTIN_ID, leptin.id)
        session.close()

        ## update
        # start new db session
        session = self.session_factory()
        # update sequence
        NEW_SEQUENCE = "THESEQUENCEHASCHANGEDEFORTHETEST"
        leptin.sequence = NEW_SEQUENCE
        session.add(leptin)
        session.commit()
        # close session to make sure nothing is cached
        session.close()
        # start new session
        session = self.session_factory()
        # query sequence of Leptin
        leptin_sequence = session.query(Peptide.sequence).filter(Peptide.id == LEPTIN_ID).scalar()
        self.assertEqual(NEW_SEQUENCE, leptin_sequence)
        session.close()


        ## delete
        # start new session
        session = self.session_factory()
        session.delete(leptin)
        session.commit()
        self.assertEqual(0, session.query(Peptide).count())
        session.close()

    def test_weight_calculation(self):
        peptide = Peptide(LEPTIN_SEQUENCE, 0)
        self.assertEqual(LEPTIN_SEQUENCE_WEIGHT, round(mass_to_float(peptide.weight), 2))

