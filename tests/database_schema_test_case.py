import re

from macpep_db.models.peptide import Peptide
from macpep_db.proteomics.amino_acid import AminoAcid
from macpep_db.proteomics.neutral_loss import NeutralLoss
from macpep_db.tasks.statistics import Statistics

from .abstract_database_test_case import AbstractDatabaseTestCase

PEPTIDE_PARTITION_COUNT = 100

class DatabaseSchemaTestCase(AbstractDatabaseTestCase):
    def create_and_delete_peptide_with_weight(self, peptide_class, weight, session):
            peptide = peptide_class("ABC", 0)
            peptide.weight = weight
            session.add(peptide)
            session.commit()
            self.assertIsInstance(peptide.id, int)
            session.delete(peptide)
            session.commit()

    def check_partition_boundaries_sequence(self, boundaries):
        # Beginning with 0
        self.assertEqual(boundaries[0][1], 0)
        # Each upper boundary of the current partition must match lower boundary of the next partition
        for idx in range(len(boundaries) - 1):
            self.assertEqual(boundaries[idx][2], boundaries[idx + 1][1])
        # Last partition needs a upper limit of 60 times Tryptophan + 1 Water + 1 for the explicit boundary
        upper_boundary = AminoAcid.get_haviest().mono_mass * 60 + NeutralLoss.get_by_name("H2O").mono_mass  + 1
        self.assertEqual(boundaries[len(boundaries) - 1][2], upper_boundary)


    def test_partition_ranges(self):
        session = self.session_factory()

        boundaries = Statistics.get_partition_boundaries(session)
        self.assertEqual(len(boundaries), PEPTIDE_PARTITION_COUNT)
        self.check_partition_boundaries_sequence(boundaries)
        for idx, partition in enumerate(boundaries):
            # Check if a peptide with a weight of the upper partition boundary minus 1 fits into the database (cause upper boundary is exclusive) 
            self.create_and_delete_peptide_with_weight(Peptide, partition[2] - 1, session)
            if idx < len(boundaries) - 1:
                # Check if a peptide with a weight of the upper partition boundary fits into the database 
                self.create_and_delete_peptide_with_weight(Peptide, partition[2], session)