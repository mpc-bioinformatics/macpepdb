"""create protein peptide association

Revision ID: e7b8ea6c30af
Revises: 09b1d715d8c6
Create Date: 2020-09-10 11:48:56.286670+00:00

"""
from alembic import op
import sqlalchemy as sa

# revision identifiers, used by Alembic.
revision = 'e7b8ea6c30af'
down_revision = '09b1d715d8c6'
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        'proteins_peptides',
        sa.Column('protein_accession', sa.VARCHAR(10), sa.ForeignKey('proteins.accession', onupdate='CASCADE')),
        sa.Column('peptide_weight', sa.BigInteger),
        sa.Column('peptide_sequence', sa.VARCHAR(60)),
        sa.ForeignKeyConstraint(['peptide_weight', 'peptide_sequence'], ['peptides.weight', 'peptides.sequence']),
        sa.PrimaryKeyConstraint('protein_accession', 'peptide_weight', 'peptide_sequence')
    )

    op.create_index('proteins_peptides_peptide_seqeunce_idx', 'proteins_peptides', ['peptide_sequence'])

def downgrade():
    op.drop_table('proteins_peptides')
