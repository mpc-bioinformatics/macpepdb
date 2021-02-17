import re
from sqlalchemy import Column, BigInteger, String, VARCHAR, Boolean, Integer
from sqlalchemy.orm import relationship

from .database_record import DatabaseRecord
from .associacions import proteins_peptides

class Protein(DatabaseRecord):
    __tablename__ = "proteins"

    id = Column(BigInteger, primary_key = True)
    # Attentions: Do not change the accession if the protein is stored in a hashabale collection, because accession is used as hash-key,
    # therefore a change of the accession will result in undetectability of this protein within the collection.
    accession = Column(VARCHAR(10))
    entry_name = Column(VARCHAR(16))
    name = Column(String)
    sequence = Column(String)
    # https://www.uniprot.org/help/taxonomic_identifier
    taxonomy_id = Column(Integer)
    # https://www.uniprot.org/help/proteome_id
    proteome_id = Column(VARCHAR(11))
    is_reviewed = Column(Boolean)

    # It is sufficient to join Peptide.id with proteins_peptides.c.peptide_id without considering the weight, because weight is contained in the primary for partitioning only. This reduces time consumption for the join by over 50 %
    peptides = relationship('Peptide', secondary=proteins_peptides, secondaryjoin="Peptide.id == proteins_peptides.c.peptide_id", back_populates='proteins', lazy='dynamic')


    def __init__(self, accession: str, entry_name: str, name: str, sequence: str, taxonomy_id: int, proteome_id: str, is_reviewed: bool):
        self.accession = accession
        self.entry_name = entry_name
        self.name = name
        self.sequence = sequence
        self.taxonomy_id = taxonomy_id
        self.proteome_id = proteome_id
        self.is_reviewed = is_reviewed

    def to_fasta_entry(self, protein_merges: list) -> str:
        review_type = "sp" if self.is_reviewed else "tr"
        non_fasta_attributes = []
        if self.proteome_id:
            non_fasta_attributes.append(f"TDBPI={self.proteome_id}")
        if len(protein_merges):
            accessions = ",".join([merge.source_accession for merge in protein_merges])
            non_fasta_attributes.append(f"TDBPM={accessions}")
        if len(non_fasta_attributes):
            non_fasta_attributes.insert(0, "")
            non_fasta_attributes = " ".join(non_fasta_attributes)
        else:
            non_fasta_attributes = ""
        return f">{review_type}|{self.accession}|{self.entry_name} {self.name} OX={self.taxonomy_id}{non_fasta_attributes}\n{self.sequence}"

    # This method is implemented to make sure only the accession is used as hash when a protein is stored in a hashable collection (Set, Dictionary, ...)
    def __hash__(self):
        return hash(self.accession)
    
    # According to the Python documentation this should be implemented if __hash__() is implemented.
    def __eq__(self, other):
        if not isinstance(other, Protein):
            return False
        return self.accession == other.accession

    def to_dict(self):
        return {
            "id": self.id,
            "accession": self.accession,
            "entry_name": self.entry_name,
            "name": self.name,
            "sequence": self.sequence,
            "taxonomy_id": self.taxonomy_id,
            "proteome_id": self.proteome_id,
            "is_reviewed": self.is_reviewed
        }

# Prevent circular import problem by SQLAlchemy relation() by put this import after the Protein definition
from .peptide import Peptide