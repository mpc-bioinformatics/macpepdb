from __future__ import annotations#
from psycopg2.extras import execute_values

class TaxonomyMerge:
    TABLE_NAME = "taxonomy_merges"

    def __init__(self, source_id: int, target_id: int):
        self.source_id = source_id
        self.target_id = target_id

    # This method is implemented to use bot accessions as hash when a protein merge is stored in a hashable collection (Set, Dictionary, ...)
    def __hash__(self):
        return hash(f"{self.source_id},{self.target_id}")
    
    # According to the Python documentation this should be implemented if __hash__() is implemented.
    def __eq__(self, other):
        return self.source_id == other.source_id and self.source_id == other.source_id

    @staticmethod
    def insert(database_cursor, taxonomy_merge: TaxonomyMerge):
        """
        @param database_cursor Database cursor
        @param taxonomy_merges TaxonomyMerge to insert
        """
        INSERT_QUERY = f"INSERT INTO {TaxonomyMerge.TABLE_NAME} (source_id, target_id) VALUES (%s, %s);"
        database_cursor.execute(
            INSERT_QUERY,
            (taxonomy_merge.source_id, taxonomy_merge.target_id)
        )

    @classmethod
    def bulk_insert(cls, database_cursor, taxonomy_merges: list):
        """
        @param database_cursor Database cursor with open transaction.
        @param taxonomy_merges TaxonomyMerges for bulk insert.
        """
        BULK_INSERT_QUERY = (
            f"INSERT INTO {cls.TABLE_NAME} (source_id, target_id) "
            "VALUES %s ON CONFLICT DO NOTHING;"
        )
        # Bulk insert the new peptides
        execute_values(
            database_cursor,
            BULK_INSERT_QUERY,
            [
                (
                    taxonomy_merge.source_id, 
                    taxonomy_merge.target_id
                ) for taxonomy_merge in taxonomy_merges
            ]
        )
    
    @classmethod
    def select(cls, database_cursor, select_conditions: tuple = ("", []), fetchall: bool = False):
        """
        @param database_cursor
        @param select_conditions A tupel with the where statement (without WHERE) and a list of parameters, e.g. ("source_id = %s", [1])
        @param fetchall Indicates if multiple rows should be fetched
        @return TaxonomyMerge or list of taxonomy_merges
        """
        select_query = f"SELECT source_id, target_id FROM {cls.TABLE_NAME}"
        if len(select_conditions) == 2 and len(select_conditions[0]):
            select_query += f" WHERE {select_conditions[0]}"
        select_query += ";"
        database_cursor.execute(select_query, select_conditions[1])
        
        if fetchall:
            return [cls(row[0], row[1]) for row in database_cursor.fetchall()]
        else:
            row = database_cursor.fetchone()
            if row:
                return cls(row[0], row[1])
            else:
                return None
