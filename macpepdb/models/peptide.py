# std imports
from __future__ import annotations

# internal imports
from macpepdb.models.peptide_base import PeptideBase
from macpepdb.models.protein_peptide_association import ProteinPeptideAssociation
from macpepdb.models import protein

class Metadata:
    def __init__(self, is_swiss_prot: bool, is_trembl: bool, taxonomy_ids: list, unique_taxonomy_ids: list, proteome_ids: list):
        self.__is_swiss_prot = is_swiss_prot
        self.__is_trembl = is_trembl
        self.__taxonomy_ids = taxonomy_ids
        self.__unique_taxonomy_ids = unique_taxonomy_ids
        self.__proteome_ids = proteome_ids

    @property
    def is_swiss_prot(self):
        return self.__is_swiss_prot
    
    @property
    def is_trembl(self):
        return self.__is_trembl
    
    @property
    def taxonomy_ids(self):
        return self.__taxonomy_ids
    
    @property
    def unique_taxonomy_ids(self):
        return self.__unique_taxonomy_ids
    
    @property
    def proteome_ids(self):
        return self.__proteome_ids
    

class Peptide(PeptideBase):
    TABLE_NAME = 'peptides'
    
    def __init__(self, sequence: str, number_of_missed_cleavages: int, metadata: Metadata = None):
        PeptideBase.__init__(self, sequence, number_of_missed_cleavages)
        self.__metadata = metadata

    @property
    def metadata(self):
        return self.__metadata

    @classmethod
    def select(cls, database_cursor, select_conditions: tuple = ("", []), fetchall: bool = False, include_metadata: bool = False):
        """
        @param database_cursor
        @param select_conditions A tupel with the where statement (without WHERE) and a list of parameters, e.g. ("accession = %s AND taxonomy_id = %s",["Q257X2", 6909])
        @param fetchall Indicates if multiple rows should be fetched
        @param include_metadata Indicates if peptides is returned with metadata (is_swiss_prot, is_trembl, taxonomy_ids, unique_taxonomy_ids, proteome_ids)
        @return Petide or list of peptides
        """
        if not include_metadata:
            return super().select(database_cursor, select_conditions, fetchall)
        else:
            select_query = f"SELECT sequence, number_of_missed_cleavages, is_swiss_prot, is_trembl, taxonomy_ids, unique_taxonomy_ids, proteome_ids FROM {cls.TABLE_NAME}"
            if len(select_conditions) == 2 and len(select_conditions[0]):
                select_query += f" WHERE {select_conditions[0]}"
            select_query += ";"
            database_cursor.execute(select_query, select_conditions[1])
            if fetchall:
                return [cls(row[0], row[1], Metadata(row[2], row[3], row[4], row[5], row[6])) for row in database_cursor.fetchall()]
            else:
                row = database_cursor.fetchone()
                if row:
                    return cls(row[0], row[1], Metadata(row[2], row[3], row[4], row[5], row[6]))
                else:
                    return None


    def proteins(self, database_cursor):
        PROTEIN_QUERY = (
            f"SELECT accession, secondary_accessions, entry_name, name, sequence, taxonomy_id, proteome_id, is_reviewed FROM {protein.Protein.TABLE_NAME} "
            f"WHERE accession = ANY(SELECT protein_accession FROM {ProteinPeptideAssociation.TABLE_NAME} WHERE partition = %s peptide_mass = %s AND peptide_sequence = %s);"
        )
        database_cursor.execute(
            PROTEIN_QUERY,
            (self.partition, self.mass, self.sequence)
        )
        return [protein.Protein.from_sql_row(row) for row in database_cursor.fetchall()]

    @staticmethod
    def flag_for_metadata_update(database_cursor, peptides: list):
        database_cursor.execute(f"UPDATE {Peptide.TABLE_NAME} SET is_metadata_up_to_date = %s WHERE sequence = ANY(%s);", (False, [peptide.sequence for peptide in peptides]))

    @staticmethod
    def update_metadata(database_cursor, peptide: Peptide):
        update_values = Peptide.generate_metadata_update_values(database_cursor, peptide)
        database_cursor.execute(
            f"UPDATE {Peptide.TABLE_NAME} SET is_metadata_up_to_date = true, is_swiss_prot = %s, is_trembl = %s, taxonomy_ids = %s, unique_taxonomy_ids = %s, proteome_ids = %s WHERE mass = %s AND sequence = %s;",
            (
                update_values[0],
                update_values[1],
                update_values[2],
                update_values[3],
                update_values[4],
                update_values[5],
                update_values[6],
                update_values[7]
            )
        )

    @staticmethod
    def generate_metadata_update_values(database_cursor, peptide: Peptide) -> tuple():
        """
        Generates a tuple with the updated metadata for peptides.

        @param database_cursor
        @param peptide
        @return tuple Form (is_swiss_prot, is_trembl, taxonomy_ids, unique_taxonomy_ids, proteome_ids, mass, sequence)
        """
        review_statuses = []
        proteome_ids = set()
        # Key is a taxonomy id, value is a counter which indicates how often the taxonomy among the referenced proteins
        taxonomy_id_count_map = {} 
        database_cursor.execute(f"SELECT is_reviewed, taxonomy_id, proteome_id FROM {protein.Protein.TABLE_NAME} WHERE accession = ANY(SELECT protein_accession FROM {ProteinPeptideAssociation.TABLE_NAME} WHERE partition = %s AND peptide_mass = %s AND peptide_sequence = %s);", (peptide.partition, peptide.mass, peptide.sequence))
        for row in database_cursor.fetchall():
            review_statuses.append(row[0])
            # Some proteins do not seeem to have an proteome ID
            if row[2] != None:
                proteome_ids.add(row[2])
            if not row[1] in taxonomy_id_count_map:
                taxonomy_id_count_map[row[1]] = 0
            taxonomy_id_count_map[row[1]] += 1
        unique_taxonomy_ids = [taxonomy_id for taxonomy_id, taxonomy_counter in taxonomy_id_count_map.items() if taxonomy_counter == 1]
        return (
            # is_swiss_prot when at least one status is true
            any(review_statuses),
            # is_trembl when not all are true
            not all(review_statuses),
            list(taxonomy_id_count_map.keys()),
            unique_taxonomy_ids,
            list(proteome_ids),
            peptide.partition,
            peptide.mass,
            peptide.sequence
        )