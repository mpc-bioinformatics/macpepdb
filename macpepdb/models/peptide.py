from .peptide_base import PeptideBase
from .protein_peptide_association import ProteinPeptideAssociation
from . import protein

class Peptide(PeptideBase):
    TABLE_NAME = 'peptides'
    
    def __init__(self, sequence: str, number_of_missed_cleavages: int):
        PeptideBase.__init__(self, sequence, number_of_missed_cleavages)

    def proteins(self, database_cursor):
        PROTEIN_QUERY = (
            f"SELECT id, accession, secondary_accessions, entry_name, name, sequence, taxonomy_id, proteome_id, is_reviewed FROM {protein.Protein.TABLE_NAME} "
            f"WHERE id = ANY(SELECT protein_id FROM {ProteinPeptideAssociation.TABLE_NAME} WHERE peptide_sequence = %s);"
        )
        database_cursor.execute(
            PROTEIN_QUERY,
            (self.__sequence)
        )
        return [protein.Protein.from_sql_row(row) for row in database_cursor.fetchall()]

    @staticmethod
    def flag_for_metadata_update(database_cursor, peptides: list):
        database_cursor.execute(f"UPDATE {Peptide.TABLE_NAME} SET is_metadata_up_to_date = %s WHERE sequence = ANY(%s);", (False, [peptide.sequence for peptide in peptides]))

    @staticmethod
    def update_metadata(database_cursor, peptide_sequence: str):
        update_values = Peptide.generate_metadata_update_values(database_cursor, peptide_sequence)
        database_cursor.execute(
            f"UPDATE {Peptide.TABLE_NAME} SET is_metadata_up_to_date = true, is_swiss_prot = %s, is_trembl = %s, taxonomy_ids = %s, unique_taxonomy_ids = %s, proteome_ids = %s WHERE sequence = %s;",
            (
                update_values[0],
                update_values[1],
                update_values[2],
                update_values[3],
                update_values[4],
                update_values[5]
            )
        )

    @staticmethod
    def generate_metadata_update_values(database_cursor, peptide_sequence: str) -> tuple():
        """
        Generates a tuple with the updated metadata for peptides.

        @param database_cursor
        @param peptide_sequence
        @return tuple Form (is_swiss_prot, is_trembl, taxonomy_ids, unique_taxonomy_ids, proteome_ids, sequence)
        """
        review_statuses = []
        proteome_ids = set()
        # Key is a taxonomy id, value is a counter which indicates how often the taxonomy among the referenced proteins
        taxonomy_id_count_map = {} 
        database_cursor.execute(f"SELECT is_reviewed, taxonomy_id, proteome_id FROM {protein.Protein.TABLE_NAME} WHERE accession = ANY(SELECT protein_accession FROM {ProteinPeptideAssociation.TABLE_NAME} WHERE peptide_sequence = %s);", (peptide_sequence,))
        for row in database_cursor.fetchall():
            review_statuses.append(row[0])
            proteome_ids.add(row[2])
            if row[1] in taxonomy_id_count_map:
                taxonomy_id_count_map[row[1]] += 1
            else:
                taxonomy_id_count_map[row[1]] = 1
        unique_taxonomy_ids = [taxonomy_id for taxonomy_id, taxonomy_counter in taxonomy_id_count_map.items() if taxonomy_counter == 1]
        return (
            # is_swiss_prot when at least one status is true
            any(review_statuses),
            # is_trembl when not all are true
            not all(review_statuses),
            list(taxonomy_id_count_map.keys()),
            unique_taxonomy_ids,
            list(proteome_ids),
            peptide_sequence
        )