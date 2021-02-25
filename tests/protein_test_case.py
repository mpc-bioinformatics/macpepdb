from macpepdb.models.protein import Protein

from .abstract_database_test_case import AbstractDatabaseTestCase

class ProteinTestCase(AbstractDatabaseTestCase):
    def test_create_read_update_delete_cycle(self):
        # Using Leptin (UniProt  accession: Q257X2)
        leptin = Protein('Q257X2', 'LEP_CAPHI', 'Leptin', 'MRCGPLYRFLWLWPYLSYVEAVPIRKVQDDTKTLIKTIVTRINDISHTQSVSSKQRVTGLDFIPGLHPLLSLSKMDQTLAIYQQILASLPSRNVIQISNDLENLRDLLHLLAASKSCPLPQVRALESLESLGVVLEASLYSTEVVALSRLQGSLQDMLRQLDLSPGC', 9925, 'UP000291000', True)

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
        # get Leptin by accession
        leptin = session.query(Protein).filter(Protein.accession == "Q257X2").one()
        self.assertEqual(LEPTIN_ID, leptin.id)
        session.close()

        ## update
        # start new db session
        session = self.session_factory()
        # update accession
        leptin.accession = "CHANGED"
        session.add(leptin)
        session.commit()
        # close session to make sure nothing is cached
        session.close()
        # start new session
        session = self.session_factory()
        # query accession of Leptin
        leptin_accession = session.query(Protein.accession).filter(Protein.id == LEPTIN_ID).scalar()
        self.assertEqual("CHANGED", leptin_accession)
        session.close()


        ## delete
        # start new session
        session = self.session_factory()
        session.delete(leptin)
        session.commit()
        self.assertEqual(0, session.query(Protein).count())
        session.close()

