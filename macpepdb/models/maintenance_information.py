from sqlalchemy import VARCHAR, Column
from sqlalchemy.dialects.postgresql import JSONB
from sqlalchemy.orm import Session

from .database_record import DatabaseRecord

class MaintenanceInformation(DatabaseRecord):
    DATABASE_STATUS_KEY = 'database_status'
    DIGESTION_PARAMTERS_KEY = 'digestion_parameters'

    __tablename__ = "maintenance_information"

    key = Column(VARCHAR(256), primary_key = True)
    values = Column(JSONB)

    def __init__(self, key: str, values: dict):
        self.key = key
        self.values = values