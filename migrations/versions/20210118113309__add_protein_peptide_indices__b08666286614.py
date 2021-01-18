"""add protein-peptide indices

Revision ID: b08666286614
Revises: d74dd05f837b
Create Date: 2021-01-18 11:33:09.424540+00:00

"""
from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision = 'b08666286614'
down_revision = 'd74dd05f837b'
branch_labels = None
depends_on = None


def upgrade():
    op.create_index('proteins_peptides_protein_id_idx', 'proteins_peptides', ['protein_id'])
    op.create_index('proteins_peptides_peptide_id_peptide_weight_idx', 'proteins_peptides', ['peptide_id', 'peptide_weight'])


def downgrade():
    op.drop_index('proteins_peptides_protein_id_idx', 'proteins_peptides')
    op.drop_index('proteins_peptides_peptide_id_peptide_weight_idx', 'proteins_peptides')
