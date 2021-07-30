from psycopg2.extras import execute_values

class ProteinPeptideAssociation:
    TABLE_NAME = 'proteins_peptides'

    def __init__(self, protein, peptide):
        self.__protein_accession = protein.accession
        self.__peptide_sequence = peptide.sequence
        self.__peptide_mass = peptide.mass

    @property
    def protein_accession(self):
        return self.__protein_accession
    
    @property
    def peptide_sequence(self):
        return self.__peptide_sequence

    @property
    def peptide_mass(self):
        return self.__peptide_mass

    @staticmethod
    def bulk_insert(database_cursor, protein_peptide_associations: list) -> int:
        """
        @param database_cursor Database cursor with open transaction.
        @param protein_peptide_associations List of ProteinPeptideAssciation.
        """
        BULK_INSERT_QUERY = f"INSERT INTO {ProteinPeptideAssociation.TABLE_NAME} (protein_accession, peptide_sequence, peptide_mass) VALUES %s ON CONFLICT DO NOTHING;"
        execute_values(
            database_cursor,
            BULK_INSERT_QUERY,
            [(association.protein_accession, association.peptide_sequence, association.peptide_mass) for association in protein_peptide_associations ],
            page_size=len(protein_peptide_associations)
        )
        # rowcount is only accurate, because the page size is as high as the number of inserted data. If the page size would be smaller rowcount would only return the rowcount of the last processed page.
        return database_cursor.rowcount
    
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
