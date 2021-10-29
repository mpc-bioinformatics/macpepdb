# std imports
from __future__ import annotations

# external imports
from psycopg2.extras import execute_values

class TaxonomyMerge:
    """
    Defines merged taxonomy.

    Parameters
    ----------
    source_id : int
        Old ID
    target_id : int
        New ID
    """

    TABLE_NAME = "taxonomy_merges"

    def __init__(self, source_id: int, target_id: int):
        self.source_id = source_id
        self.target_id = target_id

    def __hash__(self):
        """
        Implements the ability to use as key in dictionaries and sets.
        """
        return hash(f"{self.source_id},{self.target_id}")
    
    def __eq__(self, other):
        """
        Implements the equals operator.
        According to the Python documentation this should be implemented if __hash__() is implemented.
        """
        return self.source_id == other.source_id and self.source_id == other.source_id

    @staticmethod
    def insert(database_cursor, taxonomy_merge: TaxonomyMerge):
        """
        Insers a taxonomy merge into the database.

        Parameters
        ----------
        database_cursor
            Database cursor
        taxonomy_merge : TaxonomyMerge
            TaxonomyMerge to insert
        """
        INSERT_QUERY = f"INSERT INTO {TaxonomyMerge.TABLE_NAME} (source_id, target_id) VALUES (%s, %s);"
        database_cursor.execute(
            INSERT_QUERY,
            (taxonomy_merge.source_id, taxonomy_merge.target_id)
        )

    @classmethod
    def bulk_insert(cls, database_cursor, taxonomy_merges: list) -> int:
        """
        Efficiently inserts multiple taxonomy merges.

        Argument
        ========
        database_cursor
            Database cursor with open transaction.
        taxonomy_merges : List[TaxonomyMerge]
            TaxonomyMerges for bulk insert.
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
            ],
            page_size=len(taxonomy_merges)
        )
        # rowcount is only accurate, because the page size is as high as the number of inserted data. If the page size would be smaller rowcount would only return the rowcount of the last processed page.
        return database_cursor.rowcount
    
    @classmethod
    def select(cls, database_cursor, select_conditions: tuple = ("", []), fetchall: bool = False):
        """
        Selects one or many taxonomy merges

        Parameters
        ----------
        database_cursor
        select_conditions : Tuble[str, List[Any]]
            A tupel with the where statement (without WHERE) and a list of parameters, e.g. ("source_id = %s", [1])
        fetchall : bool
            Indicates if multiple rows should be fetched

        Returns
        -------
        TaxonomyMerge or list of taxonomy_merges
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
