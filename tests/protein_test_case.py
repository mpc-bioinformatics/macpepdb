from macpepdb.models.protein import Protein
from macpepdb.models.protein_peptide_association import ProteinPeptideAssociation
from macpepdb.models.peptide import Peptide
from macpepdb.proteomics.enzymes.trypsin import Trypsin

from .abstract_database_test_case import AbstractDatabaseTestCase

class ProteinTestCase(AbstractDatabaseTestCase):
    def test_lifecycle(self):
        trypsin = Trypsin(2, 6, 50)
        # Using Leptin (UniProt  accession: Q257X2)
        leptin = Protein('Q257X2', ['TESTACC'], 'LEP_CAPHI', 'Leptin', 'MRCGPLYRFLWLWPYLSYVEAVPIRKVQDDTKTLIKTIVTRINDISHTQSVSSKQRVTGLDFIPGLHPLLSLSKMDQTLAIYQQILASLPSRNVIQISNDLENLRDLLHLLAASKSCPLPQVRALESLESLGVVLEASLYSTEVVALSRLQGSLQDMLRQLDLSPGC', 9925, 'UP000291000', True)
        leptin_peptides = trypsin.digest(leptin)
        # Letpin with a new accession, old accession moved to secondary accessions, and new sequence where the first leucine is replaced by an isoleucine which creates a new peptide.
        updated_leptin = Protein('Q257X2V2', ['Q257X2', 'TESTACC'], 'LEP_CAPHI', 'Leptin', 'MRCGPIYRFLWLWPYLSYVEAVPIRKVQDDTKTLIKTIVTRINDISHTQSVSSKQRVTGLDFIPGLHPLLSLSKMDQTLAIYQQILASLPSRNVIQISNDLENLRDLLHLLAASKSCPLPQVRALESLESLGVVLEASLYSTEVVALSRLQGSLQDMLRQLDLSPGC', 9925, 'UP000291000', True)
        updated_leptin_peptides = {peptide.sequence: peptide for peptide in trypsin.digest(updated_leptin)}

        inserted_leptin_peptide_count = 0    
        ## Create
        # Start db session
        with self.database_connection:
            with self.database_connection.cursor() as database_cursor:
                inserted_leptin_peptide_count = Protein.create(database_cursor, leptin, trypsin)
                self.database_connection.commit()
                database_cursor.execute(f"SELECT id FROM {Protein.TABLE_NAME} WHERE accession = %s;", (leptin.accession,))
                leptin.id = database_cursor.fetchone()[0]
                # Check if id is now an integer (autoincrement from db)
                self.assertTrue(isinstance(leptin.id, int))
                # Check if alle peptides were inserted
                self.assertEqual(inserted_leptin_peptide_count, len(leptin_peptides))

                # Database should contain exactly the amout of leptin peptides
                database_cursor.execute(f"SELECT count(*) FROM {Peptide.TABLE_NAME};")
                self.assertEqual(len(leptin_peptides), database_cursor.fetchone()[0])

                # Database should contain also exactly one association per leptin peptide
                database_cursor.execute(f"SELECT count(*) FROM {ProteinPeptideAssociation.TABLE_NAME};")
                self.assertEqual(len(leptin_peptides), database_cursor.fetchone()[0])

                for peptide in leptin_peptides:
                    database_cursor.execute(f"SELECT true FROM {Peptide.TABLE_NAME} WHERE sequence = %s;", (peptide.sequence,))
                    self.assertTrue(database_cursor.fetchone()[0])

        ## Read
        with self.database_connection:
            with self.database_connection.cursor() as database_cursor:
                # Get Leptin by accession
                database_leptin = Protein.select(
                    database_cursor,
                    ("accession = %s", [leptin.accession])
                )
                database_leptin_petides = database_leptin.peptides(database_cursor)
                self.assertEqual(database_leptin.id, leptin.id)
                self.assertEqual(len(database_leptin_petides), len(leptin_peptides))


        ## Update
        with self.database_connection:
            with self.database_connection.cursor() as database_cursor:
                database_leptin = Protein.select(
                    database_cursor,
                    ("accession = %s", [leptin.accession])
                )
                database_leptin.update(database_cursor, updated_leptin, trypsin)
                self.database_connection.commit()
                database_cursor.execute(f"SELECT COUNT(*) FROM {Protein.TABLE_NAME};")
                protein_count = database_cursor.fetchone()[0]
                # There should still be only one protein (updated letpin)
                self.assertEqual(protein_count, 1)
                updated_database_leptin = Protein.select(database_cursor, ("accession = %s", [updated_leptin.accession]))
                # Check the updated attributes
                self.assertEqual(updated_database_leptin.accession, updated_leptin.accession)
                self.assertEqual(updated_database_leptin.secondary_accessions, updated_leptin.secondary_accessions)
                self.assertEqual(updated_database_leptin.sequence, updated_leptin.sequence)
                self.assertNotEqual(updated_database_leptin, None)
                # Fetch peptides
                updated_database_leptin_peptides = {peptide.sequence: peptide for peptide in updated_database_leptin.peptides(database_cursor)}
                self.assertEqual(len(updated_database_leptin_peptides), len(updated_leptin_peptides))
                # Cross check if only the updated leptin peptides are returned
                for sequence in updated_database_leptin_peptides.keys():
                    self.assertIn(sequence, updated_leptin_peptides)



        ## Delete
        with self.database_connection:
            with self.database_connection.cursor() as database_cursor:
                database_leptin = Protein.select(
                    database_cursor, 
                    ("accession = %s", [updated_leptin.accession])
                )
                Protein.delete(database_cursor, database_leptin)
                self.database_connection.commit()

                database_cursor.execute(f"SELECT COUNT(*) FROM {Protein.TABLE_NAME};")
                protein_count = database_cursor.fetchone()[0]
                self.assertEqual(0, protein_count)
                database_cursor.execute(f"SELECT COUNT(*) FROM {ProteinPeptideAssociation.TABLE_NAME};")
                peptide_count = database_cursor.fetchone()[0]
                self.assertEqual(0, peptide_count)
                # Peptides will not be deleted on protein delete,  because thay may be referenced by another protein

