from typing import List, Tuple, Any, Callable

class ColumnConditionsList:
    def __init__(self):
        self.__column_names = []
        self.__sql_operatos = []
        self.__func_operators = []
        self.__values = []
        self.__flatten_value_for_sql = []
        # Create a index for faster lookup
        self.__column_index = {}

    @property
    def column_names(self) -> List[str]:
        """Returns column names
        """
        return self.__column_names

    def append(self, column_name: str, sql_operator: str, func_operator: Callable[[Any, Any], bool], value: Any, flatten_value_for_sql: bool = False):
        """Appends a new conditions.
        Throws ValueError if the given column already has a condition 

        Arguments:
        column_name -- Column name, e.g. 'mass'
        sql_operator -- SQL-representation of condition operation, e.g. 'BETWEEN %s AND %s'
        func_operator -- Function which has two arguments (first argument is the value from the database, second argument is the value to check agains) and returns a bool. e.g. `lambda mass, tolerance: tolerance[0] <= mass <= tolerance[1]`
        value -- Value to check against, e.g. [200, 300]
        flatten_value_for_sql -- If value is a list or tupel and this is True the elements of value are appended to the value list returned by to_sql(). If this is False, value is added as one element to the value list of to_sql()
        """
        if column_name in self.__column_index:
            raise ValueError(f"'{column_name}' already included.")
        self.__column_names.append(column_name)
        self.__sql_operatos.append(sql_operator)
        self.__func_operators.append(func_operator)
        self.__values.append(value)
        self.__flatten_value_for_sql.append(flatten_value_for_sql)
        self.__column_index[column_name] = len(self.__column_names) - 1

    def contains_column_name(self, column_name: str) -> bool:
        return column_name in self.__column_index

    def __len__(self):
        return len(self.__column_names) 

    def to_sql(self) -> Tuple[str, List[Any]]:
        """Returns a tuple where the first element is a SQL-where-string with placeholders for the actual values
        and the second element is an array with values. Both can be used to create the actual query with psycopg.cursor.mogrify/create.
        """
        where_string = []
        where_values = []
        for idx, column_name in enumerate(self.column_names):
            where_string.append(f"{column_name} {self.__sql_operatos[idx]}")
            # Append value if it should not be flattened or if the value is not a tuple or a list
            if not self.__flatten_value_for_sql[idx] or (not isinstance(self.__values[idx], tuple) and not isinstance(self.__values[idx], list)):
                where_values.append(self.__values[idx])
            else:
                where_values += self.__values[idx]
        return (
            " AND ".join(where_string),
            where_values
        )

    def check_column_value(self, column_name: str, column_value: Any) -> bool:
        """Applies the operator on the given column value and the value given for this  
        Throws a KeyError if the column has no condition.

        Arguments:
        column_name -- Name of the column.
        column_value -- Value on the left side of the operator.
        """
        column_idx = self.__column_index[column_name]
        return self.__func_operators[column_idx](column_value, self.__values[column_idx])

