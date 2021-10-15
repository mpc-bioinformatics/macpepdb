# std imports
from dataclasses import dataclass
from typing import Any, Tuple

@dataclass
class ColumnDefinition:
    """
    Defines the column name, default operator and values to build a query which uses a specific index.
    """
    __slots__ = ["__column_name", "__operator", "__values"]
    __column_name: str
    __operator: str
    __values: Any

    def __init__(self, column_name: str, operator: str, values: Tuple[Any]):
        self.__column_name = column_name
        self.__operator = operator
        self.__values = values

    @property
    def column_name(self) -> str:
        return self.__column_name

    @property
    def operator(self) -> str:
        return self.__operator
    
    @property
    def values(self) -> Tuple[Any]:
        return self.__values

    def get_sql_definition(self) -> str:
        """
        Put column name and operation together

        Return
        ======
        E.g. "mass BETWEEN %s AND %s" when column name was "mass" and operator was "BETWEEN %s AND %s".
        """
        return f"{self.__column_name} {self.__operator}"
