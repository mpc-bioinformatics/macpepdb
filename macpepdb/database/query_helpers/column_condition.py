# std imports
from __future__ import annotations
from dataclasses import dataclass
from typing import Any, Tuple

# internal imports
from macpepdb.database.indexes.column_definition import ColumnDefinition


@dataclass
class ColumnCondition(ColumnDefinition):
    """
    Subclass of `macpepdb.database.indexes.column_definition.ColumnDefinition` but with settable operator and values attributes.

    Parameters
    ----------
    column_name : str
        Name of the targeted column
    operator : str
        Comparison operator with placeholder for values, e.g. `= %s` or `BETWEEN %s AND %s`
    values : Tuple[Any]
        Values for the operator
    """
    
    def __init__(self, column_name: str, operator: str, values: Any):
        ColumnDefinition.__init__(self, column_name, operator, values)

    @ColumnDefinition.operator.setter
    def operator(self, operator: str):
        self.__operator = operator

    @ColumnDefinition.values.setter
    def values(self, values: Tuple[Any]):
        self.__values = values

    @classmethod
    def from_column_defintion(cls, column_definition: ColumnDefinition) -> ColumnCondition:
        """
        Uses a ColumnDefinition to create a ColumnCondition

        Parameters
        ----------
        column_definition : ColumnDefinition
            Column defintion

        Returns
        -------
        Returns a column condition
        """
        return ColumnCondition(
            column_definition.column_name,
            column_definition.operator,
            column_definition.values
        )