import re

from .protein_peptide_association import ProteinPeptideAssociation
from ..proteomics.enzymes import digest_enzyme
from . import peptide as peptide_module

class Protein:
    EMBL_AMINO_ACID_GROUPS_PER_LINE = 6
    EMBL_AMINO_ACID_GROUP_LEN = 10
    EMBL_ACCESSIONS_PER_LINE = 8

    TABLE_NAME = 'proteins'

    def __init__(self, accession: str, secondary_accessions: list, entry_name: str, name: str, sequence: str, taxonomy_id: int, proteome_id: str, is_reviewed: bool, id = None):
        self.accession = accession
        self.secondary_accessions = secondary_accessions
        self.entry_name = entry_name
        self.name = name
        self.sequence = sequence
        self.taxonomy_id = taxonomy_id
        self.proteome_id = proteome_id
        self.is_reviewed = is_reviewed
        self.__id = id

    @property
    def id(self):
        return self.__id

    @id.setter
    def id(self, value):
        if self.__id == None:
            self.__id = value

    def to_embl_entry(self) -> str:
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

    # This method is implemented to make sure only the accession is used as hash when a protein is stored in a hashable collection (Set, Dictionary, ...)
    def __hash__(self):
        return hash(self.accession)
    
    # According to the Python documentation this should be implemented if __hash__() is implemented.
    def __eq__(self, other):
        if not isinstance(other, Protein):
            return False
        return self.accession == other.accession

    def to_dict(self):
        return {
            "id": self.id,
            "accession": self.accession,
            "entry_name": self.entry_name,
            "name": self.name,
            "sequence": self.sequence,
            "taxonomy_id": self.taxonomy_id,
            "proteome_id": self.proteome_id,
            "is_reviewed": self.is_reviewed
        }

    def peptides(self, database_cursor):
        if not self.__id:
            return []
        PROTEIN_QUERY = (
            f"SELECT {peptide_module.Peptide.TABLE_NAME}.id, {peptide_module.Peptide.TABLE_NAME}.sequence, {peptide_module.Peptide.TABLE_NAME}.number_of_missed_cleavages " f"FROM {peptide_module.Peptide.TABLE_NAME}, {ProteinPeptideAssociation.TABLE_NAME} "
            f"WHERE {ProteinPeptideAssociation.TABLE_NAME}.protein_id = %s AND {ProteinPeptideAssociation.TABLE_NAME}.peptide_id = {peptide_module.Peptide.TABLE_NAME}.id"
        )
        database_cursor.execute(
            PROTEIN_QUERY,
            (self.__id,)
        )
        return [
            peptide_module.Peptide(
                row[1],
                row[2],
                row[0]
            ) for row in database_cursor.fetchall()
        ]

    @staticmethod
    def select(database_cursor, select_conditions: tuple = ("", []), fetchall: bool = False):
        """
        @param database_cursor
        @param select_conditions A tupel with the where statement (without WHERE) and a list of parameters, e.g. ("accession = %s AND taxonomy_id = %s",["Q257X2", 6909])
        @param fetchall Indicates if multiple rows should be fetched
        @return Protein or list of proteins
        """
        select_query = f"SELECT id, accession, secondary_accessions, entry_name, name, sequence, taxonomy_id, proteome_id, is_reviewed FROM {Protein.TABLE_NAME}"
        if len(select_conditions) == 2 and len(select_conditions[0]):
            select_query += f" WHERE {select_conditions[0]}"
        select_query += ";"
        database_cursor.execute(select_query, select_conditions[1])
        if fetchall:
            return [Protein.from_sql_row(row) for row in database_cursor.fetchall()]
        else:
            row = database_cursor.fetchone()
            if row:
                return Protein.from_sql_row(row)
            else:
                return None

    @staticmethod
    def from_sql_row(sql_row):
        """
        @param sql_row Contains the protein columns in the following order: id, accession, secondary_accessions, entry_name, name, sequence, taxonomy_id, proteome_id, is_reviewed
        @return Protein
        """
        return Protein(
            sql_row[1],
            sql_row[2],
            sql_row[3],
            sql_row[4],
            sql_row[5],
            sql_row[6],
            sql_row[7],
            sql_row[8],
            sql_row[0]
        )

    @staticmethod
    def insert(database_cursor, protein) -> int:
        INSERT_QUERY = f"INSERT INTO {Protein.TABLE_NAME} (accession, secondary_accessions, entry_name, name, sequence, taxonomy_id, proteome_id, is_reviewed) VALUES (%s, %s, %s, %s, %s, %s, %s, %s) RETURNING id;"
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
                protein.is_reviewed
            )
        )
        return database_cursor.fetchone()[0]

    @staticmethod
    def delete(database_cursor, protein):
        if not protein.id:
            return
        DELETE_QUERY = f"DELETE FROM {Protein.TABLE_NAME} WHERE id = %s;"
        ProteinPeptideAssociation.delete(
            database_cursor,
            [
                ("protein_id = %s", protein.id)
            ]
        )
        database_cursor.execute(DELETE_QUERY, (protein.id,))

    @staticmethod
    def create(database_cursor, protein, enzyme) -> int:
        """
        Creates a new protein, by storing insert it and its peptides to the database. Make sure the protein does not already exists
        @param database_cursor Database cursor with open transaction.
        @param protein Protein to digest
        @param enzyme Digest enzym
        @return int Number of newly inserted peptides
        """
        # Create protein and get protein ID
        protein.id = Protein.insert(database_cursor, protein)
        

        # Digest protein and create sequence => peptide map
        new_peptides = {peptide.sequence: peptide for peptide in enzyme.digest(protein)}

        # Some proteins may be to short or have to few cleavage sides to produce peptides for the allowed length. If not peptides where returned, we can omit the peptide handling.
        if len(new_peptides):
            protein_peptide_associations = []

            stored_peptides = peptide_module.Peptide.select(database_cursor, ("sequence = ANY(%s)", [list(new_peptides.keys())]), True)

            # Remove already existing peptides form new peptides and associate already stored peptides
            for peptide in stored_peptides:
                new_peptides.pop(peptide.sequence, None)
                protein_peptide_associations.append(ProteinPeptideAssociation(protein.id, peptide.id))

            # Convert sequence => peptide map to peptide list
            new_peptides = list(new_peptides.values())

            new_peptide_ids = peptide_module.Peptide.bulk_insert(database_cursor, new_peptides)

            for new_peptide_id in new_peptide_ids:
                protein_peptide_associations.append(ProteinPeptideAssociation(protein.id, new_peptide_id))

            ProteinPeptideAssociation.bulk_insert(database_cursor, protein_peptide_associations)

        return len(new_peptides)

    def update(self, database_cursor, updated_protein: 'Protein', enzyme) -> int:
        """
        Updates the stored_protein with the updated_protein by comparing its attributes.
        @param database_cursor Database cursor with open transaction.
        @param updated_protein Protein with updates from file
        @param stored_protein Protein from database
        @param enzyme Digest enzym
        @return int First element indicates if the session needs to commit, seconds is the number of newly inserted peptides
        """
        number_of_new_peptides = 0

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

            # Create sequence => peptide map
            new_peptides = {peptide.sequence: peptide for peptide in enzyme.digest(updated_protein)}

            ### Remove already referenced peptides from new_peptides and dereference peptides which no longer part of the protein
            currently_referenced_peptides = self.peptides(database_cursor)
            
            peptides_ids_to_unreference = []
            # Check if the referenced peptides are in the map of new peptides
            for peptide in currently_referenced_peptides:
                if not peptide.sequence in new_peptides:
                    # If this peptide is no longer in the peptide list of the protein, delete the reference
                    peptides_ids_to_unreference.append(peptide.id)
                else:
                    # Remove it from new peptides, because it already exists and is associated with this protein
                    new_peptides.pop(peptide.sequence, None)
            if len(peptides_ids_to_unreference):
                database_cursor.execute("DELETE FROM proteins_peptides WHERE protein_id = %s AND peptide_id IN %s;", (self.id, tuple(peptides_ids_to_unreference)))
            ### At this point new_peptides contain new and not referenced peptides

            # Check if there are peptides left
            if len(new_peptides):
                ### Remove already existing peptides from new peptides and create association value
                stored_peptides = peptide_module.Peptide.select(database_cursor, ("sequence = ANY(%s)", [list(new_peptides.keys())]), True)

                protein_peptide_associations = []

                # Remove existing peptides from new_peptides and create association value
                for peptide in stored_peptides:
                    # Remove the peptide from new peptides and create association values
                    if peptide.sequence in new_peptides:
                        new_peptides.pop(peptide.sequence, None)
                        protein_peptide_associations.append(ProteinPeptideAssociation(self.id, peptide.id))

                ### At this point new_peptides should only contain peptides which not exist in the database
                number_of_new_peptides = len(new_peptides)

                # Insert new peptides
                new_peptide_ids = peptide_module.Peptide.bulk_insert(database_cursor, new_peptides.values())

                # Create association values for new peptides
                for new_peptide_id in new_peptide_ids:
                    protein_peptide_associations.append(ProteinPeptideAssociation(self.id, new_peptide_id))

                # Bulk insert new peptides
                ProteinPeptideAssociation.bulk_insert(database_cursor, protein_peptide_associations)

        # Update protein
        if len(update_columns):
            update_columns = ", ".join(update_columns)
            update_values.append(self.id)
            database_cursor.execute(f"UPDATE proteins SET {update_columns} WHERE id = %s", update_values)

        return number_of_new_peptides