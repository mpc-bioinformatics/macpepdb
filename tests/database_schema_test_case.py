from macpepdb.models.peptide import Peptide
from macpepdb.proteomics.amino_acid import AminoAcid
from macpepdb.proteomics.neutral_loss import NeutralLoss
from macpepdb.tasks.statistics import Statistics

from .abstract_database_test_case import AbstractDatabaseTestCase

PEPTIDE_PARTITION_COUNT = 100

class DatabaseSchemaTestCase(AbstractDatabaseTestCase):
    def test_partition_ranges(self):
        with self.database_connection:
            with self.database_connection.cursor() as database_cursor:
                boundaries = Statistics.get_partition_boundaries(database_cursor)
                self.assertEqual(len(boundaries), PEPTIDE_PARTITION_COUNT)
                # Beginning with 0
                self.assertEqual(boundaries[0][1], 0)
                # Each upper boundary of the current partition must match lower boundary of the next partition
                for idx in range(len(boundaries) - 1):
                    self.assertEqual(boundaries[idx][2], boundaries[idx + 1][1])
                # Last partition needs a upper limit of 60 times Tryptophan + 1 Water + 1 for the explicit boundary
                upper_boundary = AminoAcid.get_haviest().mono_mass * 60 + NeutralLoss.get_by_name("H2O").mono_mass  + 1
                self.assertEqual(boundaries[len(boundaries) - 1][2], upper_boundary)