from __future__ import annotations

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