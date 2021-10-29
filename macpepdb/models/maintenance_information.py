# std imports
from __future__ import annotations
import json

class MaintenanceInformation:
    """
    Holds information of the database, e.g. which protein was used for digestion, min/max length of peptides, ...

    Parameters
    ----------
    key : str
        Key of the maintenance information
    value : dict
        Actual information
    """

    DATABASE_STATUS_KEY = 'database_status'
    DIGESTION_PARAMTERS_KEY = 'digestion_parameters'

    TABLE_NAME = "maintenance_information"

    def __init__(self, key: str, values: dict):
        self.key = key
        self.values = values

    @staticmethod
    def select(database_cursor, key: str) -> MaintenanceInformation:
        """
        Selectes maintenance information from the database.
        Parameters
        ----------
        database_cursor
            Database cursor with open transaction.
        key : str
            Key for identification
        
        Returns
        
        -------
        Returns MaintenanceInformation if the SELECT returns something, otherwise None
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
        Inserts maintenance information into the database.

        Parameters
        ----------
        database_cursor :
            Database cursor with open transaction.
        mainenance_infromation : MaintenanceInformation
            MaintenanceInformation to insert
        """
        INSERT_QUERY = f"INSERT INTO {MaintenanceInformation.TABLE_NAME} (key, values) VALUES (%s, %s);"
        database_cursor.execute(
            INSERT_QUERY,
            (maintenance_information.key, maintenance_information.value)
        )

    @staticmethod
    def update(database_cursor, maintenance_information: MaintenanceInformation):
        """
        Updates maintenance information in the database.

        Parameters
        ----------
        database_cursor : 
            Database cursor with open transaction.
        mainenance_infromation : MaintenanceInformation
            MaintenanceInformation to insert
        """
        UPDATE_QUERY = f"UPDATE {MaintenanceInformation.TABLE_NAME} SET values = %s WHERE key = %s;"
        database_cursor.execute(
            UPDATE_QUERY,
            (json.dumps(maintenance_information.values), MaintenanceInformation.DATABASE_STATUS_KEY)
        )