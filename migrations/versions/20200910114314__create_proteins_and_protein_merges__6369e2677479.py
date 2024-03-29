"""create proteins and protein_merges

Revision ID: 6369e2677479
Revises: 
Create Date: 2020-09-10 11:43:14.099297+00:00

"""
from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision = '6369e2677479'
down_revision = None
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        'proteins',
        sa.Column('accession', sa.VARCHAR(10), nullable=False, primary_key=True),
        sa.Column('secondary_accessions', sa.dialects.postgresql.ARRAY(sa.VARCHAR(10)), nullable=False),
        sa.Column('entry_name', sa.VARCHAR(16), nullable=False),
        sa.Column('name', sa.Text),
        sa.Column('sequence', sa.Text, nullable=False),
        sa.Column('taxonomy_id', sa.Integer),
        sa.Column('proteome_id', sa.VARCHAR(11)),
        sa.Column('is_reviewed', sa.Boolean, nullable=False),
        sa.Column('updated_at', sa.BigInteger, nullable=False, server_default="0")
    )
    op.create_index('protein_secondary_accessions_idx', 'proteins', ['secondary_accessions'], postgresql_using='gin')
    op.create_index('protein_taxonomy_id_idx', 'proteins', ['taxonomy_id'])
    op.create_index('protein_proteome_id_idx', 'proteins', ['proteome_id'])

    # Make 'proteins' a reference table (copy to all nodes)
    connection = op.get_bind()
    connection.execute(sa.text("SELECT create_reference_table('proteins');"))

def downgrade():
    op.drop_table('proteins')
