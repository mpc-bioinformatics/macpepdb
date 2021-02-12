from sqlalchemy.ext.declarative import declarative_base

Base = declarative_base()

class DatabaseRecord(Base):
    __abstract__ = True
