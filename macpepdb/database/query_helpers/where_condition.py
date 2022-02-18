# std imports
from __future__ import annotations
from dataclasses import dataclass
from typing import List, Any, Optional

@dataclass
class WhereCondition:
    """
    Instances of this class keep a SQL-WHERE conditions, with placeholders for the actual values and a list
    of these values.

    Paramters
    ---------
    condition : List[str]
        SQL-condition (without `WHERE`-keyword) with placeholders for the values.
        Each element should be a column condition or a logical operator between two columns conditions,
        e.g. ["mass BETWEEN %s AND %s", "AND", "m_count >= %s"]
    values : List[ANY]
        List of values, used in the condition

    Raises
    ------
    ValueError
        If number of placeholder and values are different.
    """

    __slots__ = [
        "__condition",
        "__values"
    ]

    __condition: List[str]
    __values: List[Any]

    def __init__(self, condition: List[str], condtition_values: List[Any]):
        placeholder_count = sum([part.count("%s") for part in condition])
        if placeholder_count != len(condtition_values):
            raise ValueError(f"condition has {placeholder_count} placeholders but {len(condtition_values)} values were given")
        self.__condition = condition
        self.__values = condtition_values

    @property
    def condition(self) -> List[str]:
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

    def add(self, condition: str, value: Any):
        self.__condition.append(condition)
        self.__values.append(value)

    def get_condition_str(self, table: Optional[str] = None) -> str:
        """
        Returns the condition as string.

        Parameters
        ----------
        table : Optional[str], optional
            Table name which is added to the column names, by default None.
            E.g. if conditions is `["mass BETWEEN %s AND %s", "AND", "m_count >= %s"]` and table is `foo`
            output become: `foo.mass BETWEEN %s AND %s AND foo.m_count >= %s`

        Returns
        -------
        Condition string for SQL-query
        """
        if table is None:
            return " ".join(self.__condition)
        else:
            return " ".join([f"{table}.{part}" if part_idx % 2 == 0 and part[0].isalnum() else part for part_idx, part in enumerate(self.__condition)])

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
            self.__condition += [logical_sql_operator] + other_where_condition.condition
        else:
            self.__condition = other_where_condition.condition
        self.__values += other_where_condition.values
