from macpepdb.models.column_conditions_list import ColumnConditionsList
from psycopg2.extensions import cursor as DatabaseCursor
from typing import List, Any

class RowStream:
    def __init__(self, database_cursor: DatabaseCursor, sql_query: str, query_values: List[Any], column_names: List[str], column_conditions_list: ColumnConditionsList, chunk_size: int = 10000):
        self.__database_cursor = database_cursor
        self.__sql_query = sql_query
        self.__query_values = query_values
        self.__column_names = column_names
        self.__column_conditions_list = column_conditions_list
        self.__chunk_size = chunk_size
        self.__reset_row_chunk_and_row_index()

    def __reset_row_chunk_and_row_index(self):
        self.__row_chunk = []
        self.__row_chunk_idx = -1

    def __check_column_conditions(self, row: tuple) -> bool:
        # Checking columns for conditions
        for column_idx, column_name in enumerate(self.__column_names):
            if not self.__column_conditions_list.check_column_value(column_name, row[column_idx]):
                return False
        return True

    def __iter__(self):
        self.__database_cursor.execute(self.__sql_query, self.__query_values)
        # Reset chunk and chunk indexs
        self.__reset_row_chunk_and_row_index()
        return self

    def __next__(self):
        while True:
            self.__row_chunk_idx += 1
            if self.__row_chunk_idx < len(self.__row_chunk):
                if self.__check_column_conditions(self.__row_chunk[self.__row_chunk_idx]):
                    return self.__row_chunk[self.__row_chunk_idx]
            else:
                self.__reset_row_chunk_and_row_index()
                self.__row_chunk = self.__database_cursor.fetchmany(self.__chunk_size)
                # If chunk is empty stop iteration
                if not len(self.__row_chunk):
                    raise StopIteration