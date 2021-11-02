# std imports
from __future__ import annotations
from dataclasses import dataclass
from typing import List, Any

@dataclass
class WhereCondition:
    """
    Instances of this class keep a SQL-WHERE conditions, with placeholders for the actual values and a list
    of these values.

    Paramters
    ---------
    condition : str
        SQL-definition (without `WHERE`-keyword) with placeholders for the values, e.g. mass BETWEEN %s AND %s AND m_count >= %s
    values : List[ANY]
        List of values, used in the condition

    Raises
    ------
    ValueError
        If number of placeholder and values are different.
    """
    __slots__ = ["__condition", "__values"]

    __condition: str
    __values: List[Any]

    def __init__(self, condition: str, condtition_values: List[Any]):
        if condition.count("%s") != len(condtition_values):
            raise ValueError(f"condition has {condition.count('%s')} placeholders but {len(condtition_values)} values were given")
        self.__condition = condition
        self.__values = condtition_values

    @property
    def condition(self) -> str:
        """
        Returns
        -------
        SQL-condition
        """
        return self.__condition


    @property
    def values(self) -> List[Any]:
        """
        Returns
        -------
        Values uses in the conditions.
        """
        return self.__values

    def concatenate(self, other_where_condition: WhereCondition, logical_sql_operator: str):
        """
        Concatenates another WhereCondition with self.

        Parameters
        ----------
        other_where_condition : WhereCondition
            Other where condition
        logical_sql_operator : str
            Logical SQL-operator for the concatentation, e.g. `AND` or `OR`.
            This is only added if the condition is not and empty string
        """

        if len(self.__condition) > 0:
            self.__condition = f"{self.__condition} {logical_sql_operator} {other_where_condition.condition}"
        else:
            self.__condition = other_where_condition.condition
        self.__values += other_where_condition.values
