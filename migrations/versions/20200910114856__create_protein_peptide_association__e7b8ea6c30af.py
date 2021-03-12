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
        sa.Column('protein_id', sa.BigInteger, sa.ForeignKey('proteins.id')),
        sa.Column('peptide_id', sa.BigInteger),
        sa.Column('peptide_weight', sa.BigInteger),
        sa.ForeignKeyConstraint(['peptide_id', 'peptide_weight'], ['peptides.id', 'peptides.weight'])
    )
    op.create_index('proteins_peptides_protein_id_idx', 'proteins_peptides', ['protein_id'])
    op.create_index('proteins_peptides_peptide_id_peptide_weight_idx', 'proteins_peptides', ['peptide_id', 'peptide_weight'])


def downgrade():
    op.drop_table('proteins_peptides')
