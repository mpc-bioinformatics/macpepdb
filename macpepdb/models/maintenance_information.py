from sqlalchemy import VARCHAR, Column
from sqlalchemy.dialects.postgresql import JSONB

from .database_record import DatabaseRecord

class MaintenanceInformation(DatabaseRecord):
    __tablename__ = "maintenance_information"

    key = Column(VARCHAR(256), primary_key = True)
    values = Column(JSONB)

    def __init__(self, key: str, values: dict):
        self.key = key
        self.values = values