# # internal imports
# from macpepdb.database.indexes.abstract_index import AbstractIndex
# from macpepdb.database.indexes.column_definition import ColumnDefinition

# class PeptidesPrimaryKeyIndex(AbstractIndex):
#     """
#     Defines the index for the peptide primary key index.
#     """
    
#     COLUMN_CONDITIONS_TEMPLATE = [
#         ColumnDefinition("partition", "> %s", (-1,)),
#         ColumnDefinition("mass", "> %s", (-1,)),
#         # ColumnDefinition("sequence", "> %s", (-1,))
#     ]