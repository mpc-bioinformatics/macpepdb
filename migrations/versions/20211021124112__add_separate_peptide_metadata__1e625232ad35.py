"""add separate peptide metadata

Revision ID: 1e625232ad35
Revises: 54dc2b03bfd2
Create Date: 2021-10-21 12:41:12.161707+00:00

"""
from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision = '1e625232ad35'
down_revision = '54dc2b03bfd2'
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        "peptide_metadata",
        sa.Column("partition", sa.SmallInteger, nullable=False),
        sa.Column("mass", sa.BigInteger, nullable=False),
        sa.Column("sequence", sa.VARCHAR(60), nullable=False),
        sa.Column("is_swiss_prot", sa.Boolean, server_default=sa.text("false"), nullable=False),
        sa.Column("is_trembl", sa.Boolean, server_default=sa.text("false"), nullable=False),
        sa.Column("taxonomy_ids", sa.dialects.postgresql.ARRAY(sa.Integer), server_default="{}", nullable=False),
        sa.Column("unique_taxonomy_ids", sa.dialects.postgresql.ARRAY(sa.Integer), server_default="{}", nullable=False),
        sa.Column("proteome_ids", sa.dialects.postgresql.ARRAY(sa.VARCHAR(11)), server_default="{}", nullable=False),
        sa.PrimaryKeyConstraint("partition", "mass", "sequence")
    )

    connection = op.get_bind()
    connection.execute("SELECT create_distributed_table('peptide_metadata', 'partition', colocate_with => 'peptides');")


def downgrade():
    op.drop_table("peptide_metadata")
