# std imports
import itertools
from dataclasses import dataclass
from typing import List, Any, Tuple, Type

# internal imports
from macpepdb.database.indexes.abstract_index import AbstractIndex
from macpepdb.database.query_helpers.column_condition import ColumnCondition

@dataclass
class DatabaseIndexWhereClauseBuilder:
    """
    Creates the where clause which utilize the given database index.
    """

    __slots__ = ["__database_index", "__column_conditions", "__highest_used_column_condition_idx"]

    __database_index: Type[AbstractIndex]
    __column_conditions: List[ColumnCondition]
    __highest_used_column_condition_idx: int

    def __init__(self, database_index: AbstractIndex):
        self.__database_index = database_index
        self.__column_conditions = [ColumnCondition.from_column_defintion(column_definition) for column_definition in self.__database_index.COLUMN_CONDITIONS_TEMPLATE]
        self.__highest_used_column_condition_idx = -1

    def set_condition(self, column_condition: ColumnCondition) -> bool:
        """
        Sets a new column condition.

        Arguments
        =========
        column_condition : ColumnCondition
            Column condition

        Return
        ======
        True/False indicating if the condition is set.
        """
        try:
            column_index = self.__database_index.get_column_index(column_condition.column_name)
            self.__highest_used_column_condition_idx = max(self.__highest_used_column_condition_idx, column_index)
            self.__column_conditions[column_index] = column_condition
            return True
        except ValueError:
            return False

    def to_sql(self) -> Tuple[str, List[Any]]:
        if self.__highest_used_column_condition_idx < 0:
            self.__highest_used_column_condition_idx = len(self.__column_conditions) - 1
        return (
            " AND ".join(
                [column_condition.get_sql_definition() for column_condition_index, column_condition in enumerate(self.__column_conditions) if column_condition_index <= self.__highest_used_column_condition_idx]
            ),
            itertools.chain(*[column_condition.values for column_condition_index, column_condition in enumerate(self.__column_conditions) if column_condition_index <= self.__highest_used_column_condition_idx])
        )
