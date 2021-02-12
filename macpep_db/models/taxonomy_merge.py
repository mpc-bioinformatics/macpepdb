from .database_record import DatabaseRecord

from sqlalchemy import Column, Integer, Text, ForeignKey
from sqlalchemy.orm import relationship

class TaxonomyMerge(DatabaseRecord):
    __tablename__ = "taxonomy_merges"

    source_id = Column(Integer, primary_key=True)
    target_id = Column(Integer, primary_key=True)

    def __init__(self, source_id: int, target_id: int):
        self.source_id = source_id
        self.target_id = target_id

    # This method is implemented to use bot accessions as hash when a protein merge is stored in a hashable collection (Set, Dictionary, ...)
    def __hash__(self):
        return hash(f"{self.source_id},{self.target_id}")
    
    # According to the Python documentation this should be implemented if __hash__() is implemented.
    def __eq__(self, other):
        return self.source_id == other.source_id and self.source_id == other.source_id