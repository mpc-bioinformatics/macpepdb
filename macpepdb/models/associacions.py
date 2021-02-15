from sqlalchemy import Table, Column, ForeignKey, BigInteger, ForeignKeyConstraint
from .database_record import DatabaseRecord

proteins_peptides = Table('proteins_peptides', DatabaseRecord.metadata,
    Column('protein_id', ForeignKey('proteins.id')),
    Column('peptide_id', BigInteger),
    Column('peptide_weight', BigInteger),
    ForeignKeyConstraint(['peptide_id', 'peptide_weight'], ['peptides.id', 'peptides.weight'])
)