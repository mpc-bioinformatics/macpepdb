"""enable nullable for taxonomies.rank

Revision ID: a450785500a9
Revises: 49f8a0465f91
Create Date: 2020-09-24 14:59:20.462315+00:00

"""
from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision = 'a450785500a9'
down_revision = '49f8a0465f91'
branch_labels = None
depends_on = None


def upgrade():
    op.alter_column('taxonomies', 'rank', nullable=True)


def downgrade():
    op.alter_column('taxonomies', 'rank', nullable=False)
