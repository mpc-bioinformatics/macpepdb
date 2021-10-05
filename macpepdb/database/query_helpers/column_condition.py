from __future__ import annotations

from dataclasses import dataclass
from macpepdb.database.indexes.column_definition import ColumnDefinition
from typing import Any, Tuple

@dataclass
class ColumnCondition(ColumnDefinition):
    """
    Subclass of macpepdb.database.indexes.column_definition.ColumnDefinition but with settable operator and values attributes.
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

        Arguments
        =========
        column_definition : ColumnDefinition
            Column defintion

        Return
        ======
        Returns a column condition
        """
        return ColumnCondition(
            column_definition.column_name,
            column_definition.operator,
            column_definition.values
        )