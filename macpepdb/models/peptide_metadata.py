"""
Contains all necessary classes and constants for storing and retrieving peptide metadata.
"""

# std imports
from __future__ import annotations
from dataclasses import dataclass
from typing import List, ClassVar, Optional

# external imports
from psycopg2.extras import execute_batch, execute_values

# internal imports
from macpepdb.models import peptide as peptide_module # Do not import Peptide to prevent circular imports

@dataclass
class PeptideMetadata:
    """
    Holds metadata of peptides.

    Parameters
    ----------
    partition : int
        Database partition
    mass : int
        Mass in internal integer representation
    sequence : str
        Amino acid sequence
    is_swiss_prot : bool
        If this peptide is swissprot
    is_trembl : bool
        If this peptide is in trembl
    taxonomy_ids : List[int]
        IDs of taxonomies where this peptide is part of
    unique_taxonomy_ids : List[int]
        IDs of taxonomies where this peptide is part of and is only contained by one protein.
    proteome_ids : List [str]
        IDs of proteomes where this peptide is part of
    """

    TABLE_NAME: ClassVar[str] = "peptide_metadata"
    """Table name
    """

    SELECT_BY_PEPTIDE_QUERY: ClassVar[str] = (
        "SELECT is_swiss_prot, is_trembl, taxonomy_ids, unique_taxonomy_ids, proteome_ids "
        f"FROM {TABLE_NAME} WHERE partition = %s AND mass = %s AND sequence = %s;"
    )
    """Query for selecting metadata by peptide partition, mass and sequence.
    """

    BULK_INSERT_QUERY: ClassVar[str] = (
        f"INSERT INTO {TABLE_NAME} (partition, mass, sequence, is_swiss_prot, is_trembl, taxonomy_ids, unique_taxonomy_ids, proteome_ids) "
        "VALUES (%s, %s, %s, %s, %s, %s, %s, %s);"
    )
    """Query for bulk insert multiple peptide metadata records.
    """

    __slots__ = [
        "__is_swiss_prot",
        "__is_trembl",
        "__taxonomy_ids",
        "__unique_taxonomy_ids",
        "__proteome_ids"
    ]

    __is_swiss_prot: bool
    __is_trembl: bool
    __taxonomy_ids: List[int]
    __unique_taxonomy_ids: List[int]
    __proteome_ids: List[str]


    def __init__(self, is_swiss_prot: bool, is_trembl: bool, taxonomy_ids: List[int], unique_taxonomy_ids: List[int], proteome_ids: List [str]):
        self.__is_swiss_prot = is_swiss_prot
        self.__is_trembl = is_trembl
        self.__taxonomy_ids = taxonomy_ids
        self.__unique_taxonomy_ids = unique_taxonomy_ids
        self.__proteome_ids = proteome_ids

    @property
    def is_swiss_prot(self) -> bool:
        """
        Returns
        -------
        True if this peptide is swissprot
        """
        return self.__is_swiss_prot
        
    @property
    def is_trembl(self) -> bool:
        """
        Returns
        -------
        True if this peptide is in trembl
        """
        return self.__is_trembl

    @property
    def taxonomy_ids(self) -> List[int]:
        """
        Returns
        -------
        IDs of taxonomies where this peptide is part of
        """
        return self.__taxonomy_ids

    @property
    def unique_taxonomy_ids(self) -> List[int]:
        """
        Returns
        -------
        IDs of taxonomies where this peptide is part of and is only contained by one protein.
        """
        return self.__unique_taxonomy_ids

    @property
    def proteome_ids(self) -> List[str]:
        """
        Returns
        -------
        IDs of proteomes where this peptide is part of
        """
        return self.__proteome_ids

    @classmethod
    def select(cls, database_cursor, peptide: peptide_module.Peptide) -> Optional[PeptideMetadata]:
        """
        Selects peptide metadata.

        Parameters
        ----------
        database_cursor
            Database cursor
        peptide : Peptide
            Peptide for which the metadata is to be searched for

        Returns
        -------
        None or PeptideMetadata
        """
        database_cursor.execute(
            cls.SELECT_BY_PEPTIDE_QUERY,
            (peptide.partition, peptide.mass, peptide.sequence)
        )
        record = database_cursor.fetchone()
        if record:
            return cls(
                record[0],
                record[1],
                record[2],
                record[3],
                record[4]
            )
        return None

    @classmethod
    def bulk_insert(cls, database_cursor, peptides: List[peptide_module.Peptide], statement: str = BULK_INSERT_QUERY, is_prepared_statement: bool = False) -> int:
        """
        Insert multiple peptide metadata.

        Parameters
        ----------
        database_cursor
            Active database cursor
        peptides : List[peptide_module.Peptide]
            [description]
        statement : str
            Insert statement, make sure the statements covers the following parameters if you override it:
                1. partition,
                2. mass,
                3. sequence,
                4. is_swiss_prot,
                5. is_trembl,
                6. taxonomy_ids,
                7. unique_taxonomy_ids,
                8. proteome_ids
        is_prepared_statement : bool
            If False (default) psycopg2.extras.execute_values() is used, otherwise psycopg2.extras.execute_batch

        Returns
        -------
        Number of inserted metadata
        """
        execute_method = execute_values if not is_prepared_statement else execute_batch

        batch_parameters = [
            (
                peptide.partition,
                peptide.mass,
                peptide.sequence,
                peptide.metadata.is_swiss_prot,
                peptide.metadata.is_trembl,
                peptide.metadata.taxonomy_ids,
                peptide.metadata.unique_taxonomy_ids,
                peptide.metadata.proteome_ids
            )
            for peptide in peptides if peptide.metadata is not None
        ]

        execute_method(
            database_cursor,
            statement,
            batch_parameters,
            page_size=len(batch_parameters)
        )
        # This is only accurat, because page size is as large as the metadata count.
        return database_cursor.rowcount
