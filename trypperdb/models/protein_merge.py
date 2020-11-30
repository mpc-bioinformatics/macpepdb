from .database_record import DatabaseRecord

from sqlalchemy import Column, VARCHAR

class ProteinMerge(DatabaseRecord):
    __tablename__ = "protein_merges"

    source_accession = Column(VARCHAR(10), primary_key=True)
    target_accession = Column(VARCHAR(10), primary_key=True)

    def __init__(self, source_accession: str, target_accession: str):
        self.source_accession = source_accession
        self.target_accession = target_accession

    # This method is implemented to use bot accessions as hash when a protein merge is stored in a hashable collection (Set, Dictionary, ...)
    def __hash__(self):
        return hash(self.source_accession + self.target_accession)
    
    # According to the Python documentation this should be implemented if __hash__() is implemented.
    def __eq__(self, other):
        return self.source_accession == other.source_accession and self.target_accession == other.target_accession