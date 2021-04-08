"""create peptides

Revision ID: 09b1d715d8c6
Revises: 6369e2677479
Create Date: 2020-09-10 11:43:55.054433+00:00

"""
from alembic import op
import sqlalchemy as sa
import pathlib
import imp


# revision identifiers, used by Alembic.
revision = '09b1d715d8c6'
down_revision = '6369e2677479'
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        'peptides',
        sa.Column('weight', sa.BigInteger, nullable=False),
        sa.Column('sequence', sa.VARCHAR(60), nullable=False),
        sa.Column('length', sa.SmallInteger, nullable=False),
        sa.Column('number_of_missed_cleavages', sa.SmallInteger, nullable=False),
        sa.Column('a_count', sa.SmallInteger, default=0, nullable=False),
        sa.Column('b_count', sa.SmallInteger, default=0, nullable=False),
        sa.Column('c_count', sa.SmallInteger, default=0, nullable=False),
        sa.Column('d_count', sa.SmallInteger, default=0, nullable=False),
        sa.Column('e_count', sa.SmallInteger, default=0, nullable=False),
        sa.Column('f_count', sa.SmallInteger, default=0, nullable=False),
        sa.Column('g_count', sa.SmallInteger, default=0, nullable=False),
        sa.Column('h_count', sa.SmallInteger, default=0, nullable=False),
        sa.Column('i_count', sa.SmallInteger, default=0, nullable=False),
        sa.Column('j_count', sa.SmallInteger, default=0, nullable=False),
        sa.Column('k_count', sa.SmallInteger, default=0, nullable=False),
        sa.Column('l_count', sa.SmallInteger, default=0, nullable=False),
        sa.Column('m_count', sa.SmallInteger, default=0, nullable=False),
        sa.Column('n_count', sa.SmallInteger, default=0, nullable=False),
        sa.Column('o_count', sa.SmallInteger, default=0, nullable=False),
        sa.Column('p_count', sa.SmallInteger, default=0, nullable=False),
        sa.Column('q_count', sa.SmallInteger, default=0, nullable=False),
        sa.Column('r_count', sa.SmallInteger, default=0, nullable=False),
        sa.Column('s_count', sa.SmallInteger, default=0, nullable=False),
        sa.Column('t_count', sa.SmallInteger, default=0, nullable=False),
        sa.Column('u_count', sa.SmallInteger, default=0, nullable=False),
        sa.Column('v_count', sa.SmallInteger, default=0, nullable=False),
        sa.Column('w_count', sa.SmallInteger, default=0, nullable=False),
        sa.Column('y_count', sa.SmallInteger, default=0, nullable=False),
        sa.Column('z_count', sa.SmallInteger, default=0, nullable=False),
        sa.Column('n_terminus', sa.CHAR(1), nullable=False),
        sa.Column('c_terminus', sa.CHAR(1), nullable=False),
        # Flag to mark peptide for metadata update
        sa.Column('is_metadata_up_to_date', sa.Boolean, server_default=sa.text('false'), nullable=False),
        # Peptide metadata, collected/duplicated from proteins-table in a second step.
        sa.Column('is_swiss_prot', sa.Boolean, server_default=sa.text('false'), nullable=False),
        sa.Column('is_trembl', sa.Boolean, server_default=sa.text('false'), nullable=False),
        sa.Column('taxonomy_ids', sa.dialects.postgresql.ARRAY(sa.Integer), server_default='{}', nullable=False),
        sa.Column('unique_taxonomy_ids', sa.dialects.postgresql.ARRAY(sa.Integer), server_default='{}', nullable=False),
        sa.Column('proteome_ids', sa.dialects.postgresql.ARRAY(sa.VARCHAR(11)), server_default='{}', nullable=False),
        sa.PrimaryKeyConstraint('weight', 'sequence'),
        postgresql_partition_by='RANGE (weight)'
    )

    # The migration folder is not a Python module, so we need to import the partition_boundaries file directly
    partition_boundaries = imp.load_source('partition_boundaries', str(pathlib.Path(__file__).parent.parent.joinpath('constants').joinpath('partition_boundaries.py')))

    connection = op.get_bind()
    for idx, (lower, upper) in enumerate(partition_boundaries.PEPTIDE_WEIGHTS):
        connection.execute(f"CREATE TABLE peptides_{str(idx).zfill(3)} PARTITION OF peptides FOR VALUES FROM ('{lower}') TO ('{upper}');")

    op.create_index('peptide_seqeunce_idx', 'peptides', ['sequence'])

def downgrade():
    op.drop_table('peptides')
