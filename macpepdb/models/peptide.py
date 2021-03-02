from .peptide_base import PeptideBase
from .protein_peptide_association import ProteinPeptideAssociation
from . import protein

class Peptide(PeptideBase):
    TABLE_NAME = 'peptides'
    
    def __init__(self, sequence: str, number_of_missed_cleavages: int, id = None):
        PeptideBase.__init__(self, sequence, number_of_missed_cleavages, id)

    def proteins(self, database_cursor):
        if not self.__id:
            return []
        PROTEIN_QUERY = (
            f"SELECT id, accession, secondary_accessions, entry_name, name, sequence, taxonomy_id, proteome_id, is_reviewed FROM {protein.Protein.TABLE_NAME} "
            f"WHERE id = ANY(SELECT protein_id FROM {ProteinPeptideAssociation.TABLE_NAME} WHERE peptide_id = %s);"
        )
        database_cursor.execute(
            PROTEIN_QUERY,
            (self.__id)
        )
        return [protein.Protein.from_sql_row(row) for row in database_cursor.fetchall()]
