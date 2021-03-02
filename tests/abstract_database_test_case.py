import unittest
import os
import time

import psycopg2

# Abstrct class for database dependent test cases.
class AbstractDatabaseTestCase(unittest.TestCase):
    def setUp(self):
        """
        Creates engine and session_factory before each test.
        """
        self.database_connection = psycopg2.connect(os.getenv("TEST_MACPEPDB_URL"))

    def tearDown(self):
        """
        Removes data from all tables after each tests.
        """
        for table in ['proteins_peptides', 'peptides', 'proteins', 'taxonomies', 'taxonomy_merges', 'maintenance_information']:
            with self.database_connection:
                with self.database_connection.cursor() as database_cursor:
                    database_cursor.execute(f"delete from {table};")
        self.database_connection.close()