# std imports
from __future__ import annotations
from typing import Iterator, Optional, List, Union


# internal imports
from macpepdb.database.query_helpers.where_condition import WhereCondition
from macpepdb.models.peptide_base import PeptideBase
from macpepdb.models.protein_peptide_association import ProteinPeptideAssociation
from macpepdb.models import protein
# Do no import PeptideMetadata directly to prevent circula import
from macpepdb.models import peptide_metadata as metadata_module

class Peptide(PeptideBase):
    """
    Extend the peptide base by protein associations.

    Parameters
    ----------
    sequence : str
        Amino acid sequence
    number_of_missed_cleavages : int
        Number of missed cleavages
    metadata : metadata_module.PeptideMetadata
        Peptide metadata (optional)
    """

    TABLE_NAME = 'peptides'
    """Database table name
    """
    
    def __init__(self, sequence: str, number_of_missed_cleavages: int, metadata: metadata_module.PeptideMetadata = None):
        PeptideBase.__init__(self, sequence, number_of_missed_cleavages)
        self.__metadata = metadata

    @property
    def metadata(self) -> Optional[metadata_module.PeptideMetadata]:
        """
        Returns
        -------
        Peptide metadata if present, otherwise None.
        """
        return self.__metadata

    @classmethod
    # pylint: disable=arguments-differ
    def select(cls, database_cursor, where_condition: Optional[WhereCondition] = None,
        order_by: Optional[str] = None, fetchall: bool = False, stream: bool = False, include_metadata: bool = False) -> Optional[Union[PeptideBase, List[PeptideBase], Iterator[PeptideBase]]]:
        """
        Selects peptides.
        
        database_cursor
            Active database cursor
        where_condition : WhereCondition
            Where condition (optional)
        order_by : str
            Order by instruction, e.g `mass DESC` (optional)
        fetchall : bool
            Indicates if multiple rows should be fetched
        stream : bool
            If true, a generator is returned which yields all matching PeptideBase records
        include_metadata : bool
            Indicates if peptides is returned with metadata (is_swiss_prot, is_trembl, taxonomy_ids, unique_taxonomy_ids, proteome_ids)
        
        Returns
        -------
        None, Petide, list of peptides or generator which yield peptides
        """
        if not include_metadata:
            return super().select(database_cursor, where_condition, order_by, fetchall, stream)
        else:
            select_query = (
                f"SELECT peps.partition, peps.mass, peps.sequence, peps.number_of_missed_cleavages, meta.is_swiss_prot, meta.is_trembl, meta.taxonomy_ids, meta.unique_taxonomy_ids, meta.proteome_ids FROM {cls.TABLE_NAME} as peps "
                f"INNER JOIN {metadata_module.PeptideMetadata.TABLE_NAME} as meta ON meta.partition = peps.partition AND meta.mass = peps.mass AND meta.sequence = peps.sequence"
            )
            select_values = ()
            if where_condition is not None:
                select_query += f" WHERE {where_condition.get_condition_str(table='peps')}"
                select_values = where_condition.values
            if order_by is not None:
                select_query += f" ORDER BY peps.{order_by}"
            select_query += ";"
            database_cursor.execute(select_query, select_values)
            if not stream:
                if fetchall:
                    return [cls(row[2], row[3], metadata_module.PeptideMetadata(row[4], row[5], row[6], row[7], row[8])) for row in database_cursor.fetchall()]
                else:
                    row = database_cursor.fetchone()
                    if row:
                        return cls(row[2], row[3], metadata_module.PeptideMetadata(row[4], row[5], row[6], row[7], row[8]))
                    else:
                        return None
            else:
                def gen():
                    for row in database_cursor:
                        yield cls(row[2], row[3], metadata_module.PeptideMetadata(row[4], row[5], row[6], row[7], row[8]))
                return gen()

    def proteins(self, database_cursor):
        """
        Selects the proteins which contain this peptide.

        Returns
        -------
        List of proteins.
        """
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
        """
        Marks peptides for a metadata update.

        Parameters
        ----------
        databas_cursor
            Active database cursor
        peptides : List[Peptide]
            List of peptides to update
        """
        database_cursor.execute(f"UPDATE {Peptide.TABLE_NAME} SET is_metadata_up_to_date = %s WHERE sequence = ANY(%s);", (False, [peptide.sequence for peptide in peptides]))

    def fetch_metadata_from_proteins(self, database_cursor):
        """
        Fetches and set metadata from proteins if metadata is currently.

        Parameters
        ----------
        database_cursor
            Active database cursor
        """
        if self.metadata is None:
            review_statuses = []
            proteome_ids = set()
            # Key is a taxonomy id, value is a counter which indicates how often the taxonomy among the referenced proteins
            taxonomy_id_count_map = {} 
            database_cursor.execute(
                f"SELECT is_reviewed, taxonomy_id, proteome_id FROM {protein.Protein.TABLE_NAME} WHERE accession = ANY(SELECT protein_accession FROM {ProteinPeptideAssociation.TABLE_NAME} WHERE partition = %s AND peptide_mass = %s AND peptide_sequence = %s);", 
                (self.partition, self.mass, self.sequence)
            )
            for row in database_cursor.fetchall():
                review_statuses.append(row[0])
                # Some proteins do not seeem to have an proteome ID
                if row[2] is not None:
                    proteome_ids.add(row[2])
                if not row[1] in taxonomy_id_count_map:
                    taxonomy_id_count_map[row[1]] = 0
                taxonomy_id_count_map[row[1]] += 1
            unique_taxonomy_ids = [taxonomy_id for taxonomy_id, taxonomy_counter in taxonomy_id_count_map.items() if taxonomy_counter == 1]
            self.__metadata = metadata_module.PeptideMetadata(
                # is_swiss_prot when at least one status is true
                any(review_statuses),
                # is_trembl when not all are true
                not all(review_statuses),
                list(taxonomy_id_count_map.keys()),
                unique_taxonomy_ids,
                list(proteome_ids)
            )