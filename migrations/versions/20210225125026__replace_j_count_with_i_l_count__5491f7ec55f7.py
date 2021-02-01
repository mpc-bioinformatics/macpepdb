"""replace j_count with i/l_count

Revision ID: 5491f7ec55f7
Revises: 54dc2b03bfd2
Create Date: 2021-01-29 15:00:26.858274+00:00

"""
from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision = '5491f7ec55f7'
down_revision = '54dc2b03bfd2'
branch_labels = None
depends_on = None


def upgrade():
    op.drop_column('peptides', 'j_count')
    op.add_column('peptides', sa.Column('i_count', sa.SmallInteger, default=0, nullable=False))
    op.add_column('peptides', sa.Column('l_count', sa.SmallInteger, default=0, nullable=False))


def downgrade():
    op.add_column('peptides', sa.Column('j_count', sa.SmallInteger, default=0, nullable=False))
    op.drop_column('peptides', 'i_count')
    op.drop_column('peptides', 'l_count')
