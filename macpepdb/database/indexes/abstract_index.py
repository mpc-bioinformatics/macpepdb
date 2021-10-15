# std imports
from dataclasses import dataclass
from typing import ClassVar, List


# inner imports
from macpepdb.database.indexes.column_definition import ColumnDefinition


@dataclass
class AbstractIndex:
    """
    This class and its subclasses should provide a order for columns used in the represented database index.
    This class can than be used to order columns in a where clause so the acutal query uses the index.
    """
    COLUMN_CONDITIONS_TEMPLATE: ClassVar[List[ColumnDefinition]] = []

    @classmethod
    def get_column_index(cls, column_name: str) -> int:
        """
        Returns the index of the column within the index.

        Arguments
        =========
        column_name : str
            Name of the column to search for.

        Return
        ======
        Returns the column index.

        Raise
        =====
        Raises ValueError when the column was not found
        """
        for column_index, column_condition in enumerate(cls.COLUMN_CONDITIONS_TEMPLATE):
            if column_condition.column_name == column_name:
                return column_index
        raise ValueError
