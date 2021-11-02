# std imports
import itertools
from dataclasses import dataclass
from typing import List, Type

# internal imports
from macpepdb.database.indexes.abstract_index import AbstractIndex
from macpepdb.database.query_helpers.column_condition import ColumnCondition
from macpepdb.database.query_helpers.where_condition import WhereCondition

@dataclass
class DatabaseIndexWhereClauseBuilder:
    """
    Creates the where clause which utilize the given database index.

    Parameter
    ---------
    database_index : AbstractIndex
        Database index to use for the where clause
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

        Parameters
        ----------
        column_condition : ColumnCondition
            Column condition

        Returns
        -------
        True/False indicating if the condition is set.
        """
        try:
            column_index = self.__database_index.get_column_index(column_condition.column_name)
            self.__highest_used_column_condition_idx = max(self.__highest_used_column_condition_idx, column_index)
            self.__column_conditions[column_index] = column_condition
            return True
        except ValueError:
            return False

    def to_where_condition(self) -> WhereCondition:
        """
        Returns WhereCondition which can be used create a SQL-query which utilizes the PTM index.

        Returns
        -------
        WhereCondition
        """
        if self.__highest_used_column_condition_idx < 0:
            self.__highest_used_column_condition_idx = len(self.__column_conditions) - 1
        return WhereCondition (
            # Concatenate all column conditions with `AND`
            " AND ".join(
                [column_condition.get_sql_definition() for column_condition_index, column_condition in enumerate(self.__column_conditions) if column_condition_index <= self.__highest_used_column_condition_idx]
            ),
            # Concatenate column conditions
            list(itertools.chain(*[column_condition.values for column_condition_index, column_condition in enumerate(self.__column_conditions) if column_condition_index <= self.__highest_used_column_condition_idx]))
        )
