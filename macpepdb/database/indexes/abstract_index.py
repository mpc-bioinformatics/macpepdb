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
    """List of column definitions for the index.
    """

    @classmethod
    def get_column_index(cls, column_name: str) -> int:
        """
        Returns the index of the column within the index.

        Parameters
        ----------
        column_name : str
            Name of the column to search for.

        Returns
        -------
        Returns the column index.

        Raises
        ------
        ValueError
            If the column was not found
        """
        for column_index, column_condition in enumerate(cls.COLUMN_CONDITIONS_TEMPLATE):
            if column_condition.column_name == column_name:
                return column_index
        raise ValueError
