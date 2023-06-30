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
        sa.Column('partition', sa.SmallInteger, nullable=False),
        sa.Column('mass', sa.BigInteger, nullable=False),
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
        sa.Column('n_terminus', sa.SmallInteger, nullable=False),
        sa.Column('c_terminus', sa.SmallInteger, nullable=False),
        # Flag to mark peptide for metadata update
        sa.Column('is_metadata_up_to_date', sa.Boolean, server_default=sa.text('false'), nullable=False),
        sa.PrimaryKeyConstraint('partition', 'mass', 'sequence')
    )

    # Distribute 'peptides' by mass
    connection = op.get_bind()
    connection.execute(sa.text("SELECT create_distributed_table('peptides', 'partition');"))

    op.create_index("peptides_mass_aa_count_idx", "peptides", [
        "partition",
        "mass",
        "m_count",
        "c_count",
        "s_count",
        "t_count",
        "y_count",
        "n_terminus",
        "c_terminus",
        "r_count", 
        "k_count",
        "a_count",
        "b_count",
        "d_count",
        "e_count",
        "f_count",
        "g_count",
        "h_count",
        "i_count",
        "j_count",
        "l_count",
        "n_count",
        "o_count",
        "p_count",
        "q_count",
        "u_count",
        "v_count",
        "w_count",
        "z_count"
    ])

    

def downgrade():
    op.drop_table('peptides')
