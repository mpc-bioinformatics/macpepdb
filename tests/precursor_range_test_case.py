import unittest

from macpepdb.proteomics.mass.precursor_range import PrecursorRange

PRECURSOR = 1325887444084
PRECURSOR_TOLERANCE = 5
LOWER_LIMIT = 1325880814647
UPPER_LIMIT = 1325894073521

class PrecursorRangeTestCase(unittest.TestCase):
    def test_mass_calculation(self):
        precursor_range = PrecursorRange(PRECURSOR, 0, 0)
        self.assertEqual(PRECURSOR, precursor_range.precursor)
        self.assertEqual(PRECURSOR, precursor_range.lower_limit)
        self.assertEqual(PRECURSOR, precursor_range.upper_limit)

        precursor_range = PrecursorRange(PRECURSOR, 0, PRECURSOR_TOLERANCE)
        self.assertEqual(PRECURSOR, precursor_range.precursor)
        self.assertEqual(PRECURSOR, precursor_range.lower_limit)
        self.assertEqual(UPPER_LIMIT, precursor_range.upper_limit)
        self.assertGreaterEqual(precursor_range.upper_limit, PRECURSOR)

        precursor_range = PrecursorRange(PRECURSOR, PRECURSOR_TOLERANCE, 0)
        self.assertEqual(PRECURSOR, precursor_range.precursor)
        self.assertEqual(LOWER_LIMIT, precursor_range.lower_limit)
        self.assertEqual(PRECURSOR, precursor_range.upper_limit)
        self.assertLess(precursor_range.lower_limit, PRECURSOR)

        precursor_range = PrecursorRange(PRECURSOR, PRECURSOR_TOLERANCE, PRECURSOR_TOLERANCE)
        self.assertEqual(PRECURSOR, precursor_range.precursor)
        self.assertEqual(LOWER_LIMIT, precursor_range.lower_limit)
        self.assertEqual(UPPER_LIMIT, precursor_range.upper_limit)
        self.assertLess(precursor_range.lower_limit, PRECURSOR)
        self.assertGreater(precursor_range.upper_limit, PRECURSOR)


    def test_operators(self):
        precursor_range = PrecursorRange(PRECURSOR, PRECURSOR_TOLERANCE, PRECURSOR_TOLERANCE)

        self.assertIn(LOWER_LIMIT, precursor_range)
        self.assertIn(PRECURSOR, precursor_range)
        self.assertIn(UPPER_LIMIT, precursor_range)
        self.assertNotIn(LOWER_LIMIT-1, precursor_range)
        self.assertNotIn(UPPER_LIMIT+1, precursor_range)
        
        self.assertLess(LOWER_LIMIT-1, precursor_range)
        self.assertFalse(LOWER_LIMIT < precursor_range)

        self.assertLessEqual(UPPER_LIMIT, precursor_range)
        self.assertLessEqual(LOWER_LIMIT, precursor_range)
        self.assertLessEqual(LOWER_LIMIT-1, precursor_range)
        self.assertFalse(UPPER_LIMIT+1 <= precursor_range)

        self.assertGreater(UPPER_LIMIT + 1, precursor_range)
        self.assertFalse(UPPER_LIMIT > precursor_range)

        self.assertGreaterEqual(UPPER_LIMIT, precursor_range)
        self.assertGreaterEqual(LOWER_LIMIT, precursor_range)
        self.assertGreaterEqual(UPPER_LIMIT+1, precursor_range)
        self.assertFalse(LOWER_LIMIT-1 >= precursor_range)

        self.assertEqual(LOWER_LIMIT, precursor_range)
        self.assertEqual(PRECURSOR, precursor_range)
        self.assertEqual(UPPER_LIMIT, precursor_range)
        self.assertNotEqual(LOWER_LIMIT-1, precursor_range)
        self.assertNotEqual(UPPER_LIMIT+1, precursor_range)