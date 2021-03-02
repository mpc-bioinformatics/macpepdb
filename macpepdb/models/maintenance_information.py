from __future__ import annotations
import json

class MaintenanceInformation:
    DATABASE_STATUS_KEY = 'database_status'
    DIGESTION_PARAMTERS_KEY = 'digestion_parameters'

    TABLE_NAME = "maintenance_information"

    def __init__(self, key: str, values: dict):
        self.key = key
        self.values = values

    @staticmethod
    def select(database_cursor, key: str) -> MaintenanceInformation:
        """
        @param database_cursor Database cursor with open transaction.
        @param key Key for identification
        @return MaintenanceInformation If the SELECT returns something, otherwise None
        """
        SELECT_QUERY = f"SELECT values FROM {MaintenanceInformation.TABLE_NAME} WHERE key = %s;"
        database_cursor.execute(
            SELECT_QUERY,
            (key,)
        )
        row = database_cursor.fetchone()
        if row:
            return MaintenanceInformation(key, row[0])
        else:
            return None

    @staticmethod
    def insert(database_cursor, maintenance_information: MaintenanceInformation):
        """
        @param database_cursor Database cursor with open transaction.
        @param mainenance_infromation MaintenanceInformation to insert
        """
        INSERT_QUERY = f"INSERT INTO {TMaintenanceInformation.TABLE_NAME} (key, values) VALUES (%s, %s);"
        database_cursor.execute(
            INSERT_QUERY,
            (maintenance_information.key, maintenance_information.value)
        )

    @staticmethod
    def update(database_cursor, maintenance_information: MaintenanceInformation):
        """
        @param database_cursor Database cursor with open transaction.
        @param mainenance_infromation MaintenanceInformation to insert
        """
        UPDATE_QUERY = f"UPDATE {MaintenanceInformation.TABLE_NAME} SET values = %s WHERE key = %s;"
        database_cursor.execute(
            UPDATE_QUERY,
            (json.dumps(database_status_values), MaintenanceInformation.DATABASE_STATUS_KEY)
        )