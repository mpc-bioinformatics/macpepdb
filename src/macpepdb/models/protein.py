from __future__ import annotations
import re
from datetime import datetime
from typing import ByteString, List, Tuple, Iterator

# external imports
from psycopg2.extras import execute_values
from macpepdb.database.query_helpers.where_condition import WhereCondition


# internal imports
from macpepdb.models.protein_peptide_association import ProteinPeptideAssociation
from macpepdb.models import peptide as peptide_module

class Protein:
    """
    Defines a protein.

    Parameters
    ----------
    accession : str
        Primary accession
    secondary_accessions : List[str]
        Secondary accessions
    entry_name : str
        Entry name
    name : str
        Name
    sequence : str
        Amino acid sequence
    taxonomy_id : ID
        Taxonomy ID
    proteome_id : str
        Proteome ID
    is_reviewed : bool
        Review status
    """
    EMBL_AMINO_ACID_GROUPS_PER_LINE = 6
    EMBL_AMINO_ACID_GROUP_LEN = 10
    EMBL_ACCESSIONS_PER_LINE = 8
    # Lookup for month name by number. So no locale change is necessary
    DT_MONTH_LOOKUP_TABLE = {
        1:  "JAN",
        2:  "FEB",
        3:  "MAR",
        4:  "APR",
        5:  "MAY",
        6:  "JUN",
        7:  "JUL",
        8:  "AUG",
        9:  "SEP",
        10: "OCT",
        11: "NOV",
        12: "DEC"
    }

    TABLE_NAME = 'proteins'

    def __init__(self, accession: str, secondary_accessions: list, entry_name: str, name: str, sequence: str, taxonomy_id: int, proteome_id: str, is_reviewed: bool, updated_at: int):
        self.accession = accession
        self.secondary_accessions = secondary_accessions
        self.entry_name = entry_name
        self.name = name
        self.sequence = sequence
        self.taxonomy_id = taxonomy_id
        self.proteome_id = proteome_id
        self.is_reviewed = is_reviewed
        self.updated_at = updated_at

    def to_embl_entry(self) -> str:
        """
        Creates an EMBL entry of the protein.

        Returns
        -------
        EMBL entry
        """
        embl_entry = f"ID   {self.entry_name}    {'Reviewed' if self.is_reviewed else 'Unreviewed'};    {len(self.sequence)}\n"

        embl_accessions = [self.accession] + self.secondary_accessions
        embl_accessions_start = 0
        while embl_accessions_start < len(embl_accessions):
            # Add only 1 whitespace after AC, because each accession will be prepended by one whitespace
            embl_entry += "AC  "
            for accession in embl_accessions[embl_accessions_start:embl_accessions_start+Protein.EMBL_ACCESSIONS_PER_LINE]:
                embl_entry += f" {accession};"
            embl_entry += "\n"
            embl_accessions_start += Protein.EMBL_ACCESSIONS_PER_LINE
        

        last_update = datetime.utcfromtimestamp(self.updated_at)
        dt_day = str(last_update.day).zfill(2)
        embl_entry += f"DT   {dt_day}-{self.__class__.DT_MONTH_LOOKUP_TABLE.get(last_update.month, 'JAN')}-{last_update.year}\n"

        embl_entry += f"OX   NCBI_TaxID={self.taxonomy_id};\n"
        embl_entry += f"DR   Proteomes; {self.proteome_id};\n"
        embl_entry += f"DE   RecName: Full={self.name};\n"
        embl_entry += f"SQ   SEQUENCE\n"

        sequence_chunk_size = Protein.EMBL_AMINO_ACID_GROUP_LEN * Protein.EMBL_AMINO_ACID_GROUPS_PER_LINE
        seq_group_start = 0
        while seq_group_start < len(self.sequence):
            embl_entry += ' ' * 5
            for idx, amino_acid in enumerate(self.sequence[seq_group_start:seq_group_start+sequence_chunk_size]):
                if (idx + 1) % Protein.EMBL_AMINO_ACID_GROUP_LEN:
                    # If this is not the last  amino acid of this group add only the amino acid
                    embl_entry += amino_acid
                else:
                    # If this is the last amino acid of this group add the amino acid and a whitespace
                    embl_entry += f"{amino_acid} "
            embl_entry += "\n"
            seq_group_start += sequence_chunk_size

        embl_entry += "//"
            
        return embl_entry

    def __hash__(self):
        """
        Implements the ability to use proteins as key in sets and dictionaries.
        """
        return hash(self.accession)
    
    def __eq__(self, other):
        """
        Implements the equals operatis.
        According to the Python documentation this should be implemented if __hash__() is implemented.
        """
        if not isinstance(other, Protein):
            return False
        return self.accession == other.accession

    def peptides(self, database_cursor, order_by = None, order_descending: bool = False, offset: int = None, limit: int = None):
        """
        Selects the associated peptides of this protein.

        Parameters
        ----------
        database_cursor
            Active database cursor
        order_by : str
            adds an order column to the query
        order_descending : bool
            Descending if True, otherwise False. Applies only if order_by is not None
        offset : int
            Adds an offset to the query.
        limit : int
            Adds a limit to the query.

        Returns
        -------
        List of peptides
        """
        referenced_peptides_query = (
            f"SELECT sequence, number_of_missed_cleavages "
            f"FROM {peptide_module.Peptide.TABLE_NAME} "
            f"WHERE (partition, mass, sequence) IN (SELECT partition, peptide_mass, peptide_sequence FROM {ProteinPeptideAssociation.TABLE_NAME} WHERE protein_accession = %s)"
        )
        if order_by:
            order_type = "ASC" if not order_descending else "DESC"
            referenced_peptides_query += f" ORDER BY {order_by} {order_type}"
        if offset:
            referenced_peptides_query += f" OFFSET {offset}"
        if limit:
            referenced_peptides_query += f" LIMIT {limit}"
        referenced_peptides_query += ";"
        database_cursor.execute(
            referenced_peptides_query,
            (self.accession,)
        )
        return [
            peptide_module.Peptide(
                row[0],
                row[1]
            ) for row in database_cursor.fetchall()
        ]

    @staticmethod
    def select(database_cursor, where_condition: WhereCondition = None, fetchall: bool = False):
        """
        Selects one or many proteins.

        Parameters
        ----------
        database_cursor
            Active database cursor
        where_condition : WhereCondition
            Where condition
        fetchall : bool
            Indicates if multiple rows should be fetched

        Returns
        -------
        Protein or list of proteins
        """
        select_query = f"SELECT accession, secondary_accessions, entry_name, name, sequence, taxonomy_id, proteome_id, is_reviewed, updated_at FROM {Protein.TABLE_NAME}"
        select_values = ()
        if where_condition is not None:
            select_query += f" WHERE {where_condition.get_condition_str()}"
            select_values = where_condition.values
        select_query += ";"
        database_cursor.execute(select_query, select_values)
        if fetchall:
            return [Protein.from_sql_row(row) for row in database_cursor.fetchall()]
        else:
            row = database_cursor.fetchone()
            if row:
                return Protein.from_sql_row(row)
            else:
                return None

    @staticmethod
    def from_sql_row(sql_row: Tuple[str, List[str], str, str, str, int, str, bool, int]) -> Protein:
        """
        Parameters
        ----------
        sql_row : tuple
            Contains the protein columns in the following order: accession, secondary_accessions, entry_name, name, sequence, taxonomy_id, proteome_id, is_reviewed, updated_at

        Returns
        -------
        Protein
        """
        return Protein(
            sql_row[0],
            sql_row[1],
            sql_row[2],
            sql_row[3],
            sql_row[4],
            sql_row[5],
            sql_row[6],
            sql_row[7],
            sql_row[8]
        )

    @staticmethod
    def insert(database_cursor, protein):
        """
        Inserts a protein.

        Parameters
        ----------
        database_cursor
            Active database cursor
        protein : Protein
            Protein to insert
        """
        INSERT_QUERY = f"INSERT INTO {Protein.TABLE_NAME} (accession, secondary_accessions, entry_name, name, sequence, taxonomy_id, proteome_id, is_reviewed, updated_at) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s);"
        database_cursor.execute(
            INSERT_QUERY,
            (
                protein.accession,
                protein.secondary_accessions,
                protein.entry_name,
                protein.name,
                protein.sequence,
                protein.taxonomy_id,
                protein.proteome_id,
                protein.is_reviewed,
                protein.updated_at
            )
        )
        return database_cursor.rowcount

    @staticmethod
    def delete(database_cursor, protein):
        """
        Deletes a peptides

        Parameters
        ----------
        database_cursor
            Active database cursor
        protein : Protein
            Protein to delete
        """
        DELETE_QUERY = f"DELETE FROM {Protein.TABLE_NAME} WHERE accession = %s;"
        ProteinPeptideAssociation.delete(
            database_cursor,
            [
                ("protein_accession = %s", protein.accession)
            ]
        )
        database_cursor.execute(DELETE_QUERY, (protein.accession,))

    @staticmethod
    def create(database_cursor, protein, enzyme) -> int:
        """
        Creates a new protein, by storing insert it and its peptides to the database. Make sure the protein does not already exists

        Parameters
        ----------
        database_cursor
            Database cursor with open transaction.
        protein : Protein
            Protein to digest
        enzyme : DigestEnzym
            Digest enzym

        Returns
        -------
        Number of newly inserted peptides
        """
        # Create protein and get protein ID
        Protein.insert(database_cursor, protein)

        inserted_peptide_count = 0

        # Digest protein and create sequence => peptide map
        new_peptides = {peptide.sequence: peptide for peptide in enzyme.digest(protein)}

        # Some proteins may be to short or have to few cleavage sides to produce peptides for the allowed length. If not peptides where returned, we can omit the peptide handling.
        if len(new_peptides):
            protein_peptide_associations = []
            peptides_for_metadata_update = []

            stored_peptides = Protein.__select_existing_peptides_with_metadata_status(database_cursor, new_peptides.values())

            # Remove already existing peptides form new peptides and associate already stored peptides
            for peptide, is_metadata_up_to_date in stored_peptides:
                new_peptides.pop(peptide.sequence, None)
                protein_peptide_associations.append(ProteinPeptideAssociation(protein, peptide))
                if is_metadata_up_to_date:
                    peptides_for_metadata_update.append(peptide)

            # Convert sequence => peptide map to peptide list
            new_peptides = list(new_peptides.values())

            if len(new_peptides):
                inserted_peptide_count = peptide_module.Peptide.bulk_insert(database_cursor, new_peptides)

                for peptide in new_peptides:
                    protein_peptide_associations.append(ProteinPeptideAssociation(protein, peptide))

            if len(protein_peptide_associations):
                ProteinPeptideAssociation.bulk_insert(database_cursor, protein_peptide_associations)

            if len(peptides_for_metadata_update):
                peptide_module.Peptide.flag_for_metadata_update(database_cursor, peptides_for_metadata_update)

        return inserted_peptide_count

    def update(self, database_cursor, updated_protein: Protein, enzyme: digest_enzyme.DigestEnzyme) -> int:
        """
        Updates the protein with the updated_protein if the updated_at timestamp of the given protein is higher than the updated_at timestamp from the current protein.
        
        Parameters
        ----------
        database_cursor : DatabaseCursor
            Database cursor with open transaction.
        updated_protein : Protein
            Protein with updates from file
        stored_protein : Protein
            Protein from database
        enzyme : Digest
            Digest enzym

        Return
        ------
        Number of newly inserted peptides
        """
        # This protein has the more current update date skip the update
        if self.updated_at >= updated_protein.updated_at:
            return 0


        inserted_peptide_count = 0

        update_columns = []
        update_values = []

        if updated_protein.accession != self.accession:
            update_columns.append('accession = %s')
            update_values.append(updated_protein.accession)
        if updated_protein.secondary_accessions != self.secondary_accessions:
            update_columns.append('secondary_accessions = %s')
            update_values.append(updated_protein.secondary_accessions)
        if updated_protein.taxonomy_id != self.taxonomy_id:
            update_columns.append('taxonomy_id = %s')
            update_values.append(updated_protein.taxonomy_id)
        if updated_protein.proteome_id != self.proteome_id:
            update_columns.append('proteome_id = %s')
            update_values.append(updated_protein.proteome_id)
        if updated_protein.sequence != self.sequence:
            update_columns.append('sequence = %s')
            update_values.append(updated_protein.sequence)
        if updated_protein.updated_at != self.updated_at:
            update_columns.append('updated_at = %s')
            update_values.append(updated_protein.updated_at)

            # Create sequence => peptide map
            new_peptides = {peptide.sequence: peptide for peptide in enzyme.digest(updated_protein)}

            peptides_for_metadata_update = []

            ### Remove already referenced peptides from new_peptides and dereference peptides which no longer part of the protein
            referenced_peptides_query = (
                f"SELECT {peptide_module.Peptide.TABLE_NAME}.sequence, {peptide_module.Peptide.TABLE_NAME}.number_of_missed_cleavages, {peptide_module.Peptide.TABLE_NAME}.is_metadata_up_to_date "
                f"FROM {peptide_module.Peptide.TABLE_NAME} "
                f"WHERE (partition, mass, sequence) IN (SELECT partition, peptide_mass, peptide_sequence FROM {ProteinPeptideAssociation.TABLE_NAME} WHERE protein_accession = %s);"
            )
            database_cursor.execute(referenced_peptides_query, (self.accession,))
            currently_referenced_peptides = [(peptide_module.Peptide(row[0], row[1]), row[2]) for row in database_cursor.fetchall()]
            
            peptides_to_unreference = []
            # Check if the referenced peptides are in the map of new peptides
            for peptide, is_metadata_up_to_date in currently_referenced_peptides:
                if not peptide.sequence in new_peptides:
                    # If this peptide is no longer in the peptide list of the protein, delete the reference
                    peptides_to_unreference.append(peptide)
                    if is_metadata_up_to_date:
                        peptides_for_metadata_update.append(peptide)
                else:
                    # Remove it from new peptides, because it already exists and is associated with this protein
                    new_peptides.pop(peptide.sequence, None)
            if len(peptides_to_unreference):
                database_cursor.execute("DELETE FROM proteins_peptides WHERE protein_accession = %s AND peptide_sequence  = ANY(%s);", (self.accession, [peptide.sequence for peptide in peptides_to_unreference]))
            ### At this point new_peptides contain new and not referenced peptides

            # Check if there are peptides left
            if len(new_peptides):
                ### Remove already existing peptides from new peptides and create association value
                stored_peptides = Protein.__select_existing_peptides_with_metadata_status(database_cursor, new_peptides.values())

                protein_peptide_associations = []

                # Remove existing peptides from new_peptides and create association value
                for peptide, is_metadata_up_to_date in stored_peptides:
                    # Remove the peptide from new peptides and create association values
                    if peptide.sequence in new_peptides:
                        new_peptides.pop(peptide.sequence, None)
                        protein_peptide_associations.append(ProteinPeptideAssociation(self, peptide))
                        if is_metadata_up_to_date:
                            peptides_for_metadata_update.append(peptide)
                
                if len(new_peptides):
                    # Insert new peptides
                    inserted_peptide_count = peptide_module.Peptide.bulk_insert(database_cursor, new_peptides.values())

                # Create association values for new peptides
                for peptide in new_peptides.values():
                    protein_peptide_associations.append(ProteinPeptideAssociation(self, peptide))

                if len(protein_peptide_associations):
                    # Bulk insert new peptides
                    ProteinPeptideAssociation.bulk_insert(database_cursor, protein_peptide_associations)

            if len(peptides_for_metadata_update):
                peptide_module.Peptide.flag_for_metadata_update(database_cursor, peptides_for_metadata_update)

        # Update protein
        if len(update_columns):
            update_columns = ", ".join(update_columns)
            update_values.append(self.accession)
            database_cursor.execute(f"UPDATE proteins SET {update_columns} WHERE accession = %s", update_values)

        return inserted_peptide_count

    @staticmethod
    def __select_existing_peptides_with_metadata_status(database_cursor, peptides: list) -> list:
        """
        Updates the stored_protein with the updated_protein by comparing its attributes.

        Parameters
        ----------
        database_cursor
            Database cursor.
        peptides : List[Peptide]
            List of peptides

        Returns
        -------
        List of tuple, each tuple contains a peptide and the peptides metadata status
        """

        # %s after VLAUES is substitutet by "(mass, sequence), (mass, sequence)"
        EXISTING_PEPTIDE_QUERY = (
            f"SELECT {peptide_module.Peptide.TABLE_NAME}.sequence, {peptide_module.Peptide.TABLE_NAME}.number_of_missed_cleavages, {peptide_module.Peptide.TABLE_NAME}.is_metadata_up_to_date "
            f"FROM {peptide_module.Peptide.TABLE_NAME} "
            f"WHERE ({peptide_module.Peptide.TABLE_NAME}.partition, {peptide_module.Peptide.TABLE_NAME}.mass, {peptide_module.Peptide.TABLE_NAME}.sequence) IN (VALUES %s);"
        )
        return [
            (
                peptide_module.Peptide(row[0], row[1]), 
                row[2]
            ) for row in execute_values(
                database_cursor,
                EXISTING_PEPTIDE_QUERY,
                [
                    (
                        peptide.partition,
                        peptide.mass,
                        peptide.sequence
                    ) for peptide in peptides
                ],
                template="(%s, %s, %s)",
                page_size=len(peptides),
                fetch=True
            )
        ]

    def to_json(self) -> Iterator[ByteString]:
        """
        Generator which yields the protein as a json formatted string

        Yields
        ------
        Iterator[ByteString]
            JSON formatted string.
        """
        yield b"{\"accession\":\""
        yield self.accession.encode("utf-8")
        yield b"\",\"entry_name\":\""
        yield self.accession.encode("utf-8")
        yield b"\",\"name\":\""
        yield self.name.replace("\\", "\\\\").encode("utf-8")   # Proteins with "\" in name exists (A0A084VXA8), need to escape it.
        yield b"\",\"sequence\":\""
        yield self.sequence.encode("utf-8")
        yield b"\",\"taxonomy_id\":"
        yield str(self.taxonomy_id).encode("utf-8")
        yield b",\"proteome_id\":"
        yield f"\"{self.proteome_id}\"".encode("utf-8") if self.proteome_id is not None else b"null"
        yield b",\"is_reviewed\":"
        yield b"true" if self.is_reviewed else b"false"
        yield b"}"