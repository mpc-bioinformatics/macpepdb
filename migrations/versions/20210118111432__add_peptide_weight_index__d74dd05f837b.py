"""add peptide weight index

Revision ID: d74dd05f837b
Revises: c0986aa71a43
Create Date: 2021-01-18 11:14:32.273967+00:00

"""
from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision = 'd74dd05f837b'
down_revision = 'c0986aa71a43'
branch_labels = None
depends_on = None


def upgrade():
    op.create_index('peptide_weight_idx', 'peptides', ['weight'])


def downgrade():
    op.drop_index('peptide_weight_idx', 'peptides')
