"""add maintenance information table

Revision ID: 54dc2b03bfd2
Revises: b08666286614
Create Date: 2021-02-17 15:37:20.411377+00:00

"""
from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision = '54dc2b03bfd2'
down_revision = '49f8a0465f91'
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        'maintenance_information',
        sa.Column('key', sa.VARCHAR(256), primary_key=True),
        sa.Column('values', sa.dialects.postgresql.JSONB, nullable=False)
    )


def downgrade():
    op.drop_table('maintenance_information')
