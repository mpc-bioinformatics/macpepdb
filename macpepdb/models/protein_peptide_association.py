from psycopg2.extras import execute_values

class ProteinPeptideAssociation:
    TABLE_NAME = 'proteins_peptides'

    def __init__(self, protein_id, peptide_id):
        self.__protein_id = protein_id
        self.__peptide_id = peptide_id

    @property
    def protein_id(self):
        return self.__protein_id
    
    @property
    def peptide_id(self):
        return self.__peptide_id

    @staticmethod
    def bulk_insert(database_cursor, protein_peptide_associations: list):
        """
        @param database_cursor Database cursor with open transaction.
        @param protein_peptide_associations List of ProteinPeptideAssciation.
        """
        BULK_INSERT_QUERY = f"INSERT INTO {ProteinPeptideAssociation.TABLE_NAME} (protein_id, peptide_id, peptide_weight) VALUES %s;"
        # The peptide weigth is indeed part of the foreign key, but we reduce the JOIN to protein_id/peptide_id to be faster. So we set it to None as SQLAlchemy has done it before.
        execute_values(
            database_cursor,
            BULK_INSERT_QUERY,
            [(association.protein_id, association.peptide_id, None) for association in protein_peptide_associations ]
        )
    
    @staticmethod
    def delete(database_cursor, where_conditions: list):
        """
        Deletes ProteinPeptideAssociation from the database. Does nothing if no where_conditions were given.
        @param database_cursor
        @param where_conditions List of tupel, where each's tupel first element is the condition and the second element is the value, e.g. ("peptide = %s", "Q257X2)
        """
        if len(where_conditions):
            delete_query = f"DELETE FROM {ProteinPeptideAssociation.TABLE_NAME} WHERE "
            delete_query += " AND ".join([condition[0] for condition in where_conditions])
            delete_query += ";"
            database_cursor.execute(delete_query, [condition[1] for condition in where_conditions])
