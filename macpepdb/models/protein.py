import re
from sqlalchemy import Column, BigInteger, String, VARCHAR, Boolean, Integer
from sqlalchemy.dialects.postgresql import ARRAY as psql_array
from sqlalchemy.orm import relationship

from .database_record import DatabaseRecord
from .associacions import proteins_peptides

class Protein(DatabaseRecord):
    EMBL_AMINO_ACID_GROUPS_PER_LINE = 6
    EMBL_AMINO_ACID_GROUP_LEN = 10
    EMBL_ACCESSIONS_PER_LINE = 8

    __tablename__ = "proteins"

    id = Column(BigInteger, primary_key = True)
    # Attentions: Do not change the accession if the protein is stored in a hashabale collection, because accession is used as hash-key,
    # therefore a change of the accession will result in undetectability of this protein within the collection.
    accession = Column(VARCHAR(10))
    secondary_accessions = Column(psql_array(VARCHAR(10)), nullable=False)
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


    def __init__(self, accession: str, secondary_accessions: list, entry_name: str, name: str, sequence: str, taxonomy_id: int, proteome_id: str, is_reviewed: bool):
        self.accession = accession
        self.secondary_accessions = secondary_accessions
        self.entry_name = entry_name
        self.name = name
        self.sequence = sequence
        self.taxonomy_id = taxonomy_id
        self.proteome_id = proteome_id
        self.is_reviewed = is_reviewed

    def to_embl_entry(self) -> str:
        embl_entry = f"ID   {self.entry_name}    {'Reviewed' if self.is_reviewed else 'Unreviewed'};    {len(self.sequence)}\n"

        embl_accessions = [self.accession] + self.secondary_accessions
        embl_accessions_start = 0
        while embl_accessions_start < len(embl_accessions):
            # Add only 1 whitespace after AC, because each accession will be prepended by one whitespace
            embl_entry += "AC  "
            for accession in embl_accessions[embl_accessions_start:embl_accessions_start+Protein.EMBL_ACCESSIONS_PER_LINE]:
                embl_entry += f" {accession};"
            embl_entry += "\n"
            embl_accessions_start += Protein.EMBL_ACCESSIONS_PER_LINE

        embl_entry += f"OX   NCBI_TaxID={self.taxonomy_id};\n"
        embl_entry += f"DR   Proteomes; {self.proteome_id};\n"
        embl_entry += f"DE   RecName: Full={self.name};\n"
        embl_entry += f"SQ   SEQUENCE\n"

        sequence_chunk_size = Protein.EMBL_AMINO_ACID_GROUP_LEN * Protein.EMBL_AMINO_ACID_GROUPS_PER_LINE
        seq_group_start = 0
        while seq_group_start < len(self.sequence):
            embl_entry += ' ' * 5
            for idx, amino_acid in enumerate(self.sequence[seq_group_start:seq_group_start+sequence_chunk_size]):
                if (idx + 1) % Protein.EMBL_AMINO_ACID_GROUP_LEN:
                    # If this is not the last  amino acid of this group add only the amino acid
                    embl_entry += amino_acid
                else:
                    # If this is the last amino acid of this group add the amino acid and a whitespace
                    embl_entry += f"{amino_acid} "
            embl_entry += "\n"
            seq_group_start += sequence_chunk_size

        embl_entry += "//"
            
        return embl_entry

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