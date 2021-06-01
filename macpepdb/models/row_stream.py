import uuid

from macpepdb.models.column_conditions_list import ColumnConditionsList
from psycopg2.extensions import connection as DatabaseConnection
from typing import List, Any

class RowStream:
    def __init__(self, database_connection: DatabaseConnection, sql_query: str, query_values: List[Any], column_names: List[str], column_conditions_list: ColumnConditionsList, chunk_size: int = 10000):
        self.__database_connection = database_connection
        self.__sql_query = sql_query
        self.__query_values = query_values
        self.__column_names = column_names
        self.__column_conditions_list = column_conditions_list
        self.__chunk_size = chunk_size
        self.__database_cursor = None
        self.__cursor_uuid = None

    def __check_column_conditions(self, row: tuple) -> bool:
        # Checking columns for conditions
        for column_idx, column_name in enumerate(self.__column_names):
            if not self.__column_conditions_list.check_column_value(column_name, row[column_idx]):
                return False
        return True

    def __close_iterator(self):
        if self.__database_cursor and not self.__database_cursor.closed:
            self.__database_cursor.close()

    
    def __del__(self):
        self.__close_iterator()

    def __iter__(self):
        # Close old iterator (maybe from interrupted loop)
        self.__close_iterator()

        # Create a named cursor with a uuid based on the hostname and the timestamp.
        self.__cursor_uuid = str(uuid.uuid1())
        self.__database_cursor = self.__database_connection.cursor(name=f"row_stream_{self.__cursor_uuid}")
        # Set itersize
        self.__database_cursor.itersize = self.__chunk_size
        # Execute query
        self.__database_cursor.execute(self.__sql_query, self.__query_values)
        return self

    def __next__(self):
        while True:
            try:
                row = next(self.__database_cursor)
                if self.__check_column_conditions(row):
                    return row
            except StopIteration:
                self.__close_iterator()
                raise StopIteration