from sqlalchemy.orm import relationship

from .database_record import DatabaseRecord
from .peptide_base import PeptideBase
from .associacions import proteins_peptides

class Peptide(PeptideBase, DatabaseRecord):
    __tablename__ = "peptides"
    # It is sufficient to join Peptide.id with proteins_peptides.c.peptide_id without considering the weight, because weight is contained in the primary for partitioning only. This reduces time consumption for the join by over 50 %
    proteins = relationship('Protein', secondary=proteins_peptides, primaryjoin="Peptide.id == proteins_peptides.c.peptide_id", back_populates='peptides', lazy='dynamic')

    def __init__(self, sequence: str, number_of_missed_cleavages: int):
        PeptideBase.__init__(self, sequence, number_of_missed_cleavages)

# Prevent circular import problem by SQLAlchemy relation() by put this import after the Peptide definition
from .protein import Protein