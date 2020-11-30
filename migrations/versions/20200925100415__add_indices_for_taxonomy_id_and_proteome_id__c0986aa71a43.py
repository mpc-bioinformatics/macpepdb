"""add indices for taxonomy_id and proteome_id

Revision ID: c0986aa71a43
Revises: a450785500a9
Create Date: 2020-09-25 10:04:15.336571+00:00

"""
from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision = 'c0986aa71a43'
down_revision = 'a450785500a9'
branch_labels = None
depends_on = None


def upgrade():
    op.create_index('protein_taxonomy_id_idx', 'proteins', ['taxonomy_id'])
    op.create_index('protein_proteome_id_idx', 'proteins', ['proteome_id'])


def downgrade():
    op.drop_index('protein_taxonomy_id_idx', 'proteins')
    op.drop_index('protein_proteome_id_idx', 'proteins')
