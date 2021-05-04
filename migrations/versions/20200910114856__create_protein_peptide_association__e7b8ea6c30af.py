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
        # Column 'mass' contains the peptide mass. Because the colocation feature of Citus needs the same column name in both distributed tables, this column is called 'weight' instead of 'peptide_weight'.
        sa.Column('mass', sa.BigInteger),
        sa.Column('peptide_sequence', sa.VARCHAR(60)),
        sa.ForeignKeyConstraint(['mass', 'peptide_sequence'], ['peptides.mass', 'peptides.sequence']),
        sa.PrimaryKeyConstraint('protein_accession', 'mass', 'peptide_sequence')
    )

    op.create_index('proteins_peptides_peptide_mass_seqeunce_idx', 'proteins_peptides', ['mass', 'peptide_sequence'])

    connection = op.get_bind()
    connection.execute("SELECT create_distributed_table('proteins_peptides', 'mass', colocate_with => 'peptides');",)


def downgrade():
    op.drop_table('proteins_peptides')
