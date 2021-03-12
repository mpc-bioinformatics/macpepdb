"""create taxonomies and taxonomy_merges

Revision ID: 49f8a0465f91
Revises: e7b8ea6c30af
Create Date: 2020-09-10 11:50:27.073067+00:00

"""
from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision = '49f8a0465f91'
down_revision = 'e7b8ea6c30af'
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        'taxonomies',
        sa.Column('id', sa.Integer, primary_key=True),
        sa.Column('parent_id', sa.Integer),
        sa.Column('name', sa.Text, nullable=False),
        sa.Column('rank', sa.SmallInteger, nullable=True)
    )
    op.create_index('taxonomies_parent_id_idx', 'taxonomies', ['parent_id'])
    op.create_index('taxonomies_name_idx', 'taxonomies', ['name'])
    op.create_index('taxonomies_rank_idx', 'taxonomies', ['rank'])

    op.create_table(
        'taxonomy_merges',
        sa.Column('source_id', sa.Integer),
        sa.Column('target_id', sa.Integer),
        sa.PrimaryKeyConstraint('source_id', 'target_id')
    )
    op.create_index('taxonomy_merges_target_id_idx', 'taxonomy_merges', ['target_id'])


def downgrade():
    op.drop_table('taxonomies')
    op.drop_table('taxonomy_merges')