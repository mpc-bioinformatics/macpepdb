from macpepdb.models.protein import Protein
from macpepdb.models.peptide import Peptide
from macpepdb.proteomics.enzymes.trypsin import Trypsin

from .abstract_database_test_case import AbstractDatabaseTestCase

class ProteinPeptideAssociationTestCase(AbstractDatabaseTestCase):
    def test_create_read_update_delete_cycle(self):
        # Using Leptin (UniProt  accession: Q257X2)
        leptin = Protein('Q257X2', ['TESTACC'], 'LEP_CAPHI', 'Leptin', 'MRCGPLYRFLWLWPYLSYVEAVPIRKVQDDTKTLIKTIVTRINDISHTQSVSSKQRVTGLDFIPGLHPLLSLSKMDQTLAIYQQILASLPSRNVIQISNDLENLRDLLHLLAASKSCPLPQVRALESLESLGVVLEASLYSTEVVALSRLQGSLQDMLRQLDLSPGC', 9925, 'UP000291000', True)

        ## Create
        # Start db session
        session = self.session_factory()
        session.add(leptin)
        try:
            session.commit()
        except:
            print("leptin already exists")
        # Check if id is now an integer (autoincrement from db)
        self.assertTrue(isinstance(leptin.id, int))
        # Digest protein and save peptides with association to protein
        trypsin = Trypsin(3, 0, 60)
        leptin.peptides = trypsin.digest(leptin)
        PEPTIDE_COUNT = leptin.peptides.count()
        session.commit()
        for peptide in leptin.peptides:
            # Check if id is now an integer (autoincrement from db)
            self.assertTrue(isinstance(peptide.id, int))
        # Save autoincremented id
        LEPTIN_ID = leptin.id
        # Close session and set Leptin to None, so all connections to db are lost
        session.close()
        leptin = None

        ## Read
        session = self.session_factory()
        # Get Leptin by accession
        leptin = session.query(Protein).filter(Protein.accession == "Q257X2").one()
        leptin_petides = leptin.peptides.all()
        self.assertEqual(LEPTIN_ID, leptin.id)
        self.assertEqual(PEPTIDE_COUNT, len(leptin_petides))
        session.close()

        ## Update
        # Not implemented yes

        ## Delete
        # Start new session
        session = self.session_factory()
        # Bound letpin to the new session
        session.add(leptin)
        # Remove association between leptin and peptides
        for peptide in leptin.peptides.all():
            session.delete(peptide)
        session.delete(leptin)
        session.commit()
        self.assertEqual(0, session.query(Protein).count())
        self.assertEqual(0, session.query(Peptide).count())
        session.close()

