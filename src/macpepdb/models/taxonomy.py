# internal imports
from __future__ import annotations
from enum import IntEnum, unique
from typing import List

# external imports
from psycopg2.extras import execute_values
from psycopg2.extensions import cursor as DatabaseCursor

@unique
class TaxonomyRank(IntEnum):
    """
    Defines all types of taxonomy ranks.
    Enum names created by rank.upper().replace(" ", "_")
    """
    BIOTYPE             = 0
    CLADE               = 1
    CLASS               = 2
    COHORT              = 3
    FAMILY              = 4
    FORMA               = 5
    FORMA_SPECIALIS     = 6
    GENOTYPE            = 7
    GENUS               = 8
    INFRACLASS          = 9
    INFRAORDER          = 10
    ISOLATE             = 11
    KINGDOM             = 12
    MORPH               = 13
    NO_RANK             = 14
    ORDER               = 15
    PARVORDER           = 16
    PATHOGROUP          = 17
    PHYLUM              = 18
    SECTION             = 19
    SERIES              = 20
    SEROGROUP           = 21
    SEROTYPE            = 22
    SPECIES             = 23
    SPECIES_GROUP       = 24
    SPECIES_SUBGROUP    = 25
    STRAIN              = 26
    SUBCLASS            = 27
    SUBCOHORT           = 28
    SUBFAMILY           = 29
    SUBGENUS            = 30
    SUBKINGDOM          = 31
    SUBORDER            = 32
    SUBPHYLUM           = 33
    SUBSECTION          = 34
    SUBSPECIES          = 35
    SUBTRIBE            = 36
    SUBVARIETY          = 37
    SUPERCLASS          = 38
    SUPERFAMILY         = 39
    SUPERKINGDOM        = 40
    SUPERORDER          = 41
    SUPERPHYLUM         = 42
    TRIBE               = 43
    VARIETAS            = 44

    def __str__(self):
        return self.name.lower().replace("_", " ")

    @classmethod
    def from_string(cls, rank: str):
        """
        Returns TaxonomyRank by string

        Parameters
        ----------
        rank : str
            Name of the rank

        Returns
        -------
        TaxonomyRank
        """
        rank = rank.upper().replace(" ", "_")
        if rank in cls.__members__:
            return cls.__members__[rank]

class Taxonomy:
    """
    Defines a taxonomy with id, parant id (to build tree), name and rank.
    
    Parameters
    ----------
    id : id
        ID
    parent_id : id
        Parent ID
    name : str
        Name
    rank : TaxonomyRank
        Rank
    """

    TABLE_NAME = "taxonomies"
    """Database table name
    """

    def __init__(self, id: int, parent_id: int, name: str, rank: TaxonomyRank):
        self.id = id
        self.parent_id = parent_id
        self.name = name
        self.rank = rank

    def parent(self, database_cursor):
        """
        Selects the parent.

        Parameters
        ----------
        database_cursor
            Active database cursor

        Returns
        -------
        Taxonomy
        """
        PARENT_QUERY = "SELECT id, parent_id, name, rank WHERE id = %s;"
        database_cursor.execute(
            PARENT_QUERY,
            (self.parent_id,)
        )
        row = database_cursor.fetchone()
        return Taxonomy(
            row[0],
            row[1],
            row[2],
            TaxonomyRank(row[3])
        )

    def __hash__(self):
        """
        Implements the ability to use as key in dictionaries and sets.
        """
        return hash(self.id)
    
    def __eq__(self, other):
        """
        Implements the equals operator.
        According to the Python documentation this should be implemented if __hash__() is implemented.
        """
        return self.id == other.id


    @staticmethod
    def insert(database_cursor, taxonomy: Taxonomy):
        """
        Inserts a taxonomy into the database.

        Parameters
        ----------
        database_cursor
            Database cursor
        taxonomy : Taxonomy
            Taxonomy to insert
        """
        INSERT_QUERY = f"INSERT INTO {Taxonomy.TABLE_NAME} (id, parent_id, name, rank) VALUES (%s, %s, %s, %s);"
        database_cursor.execute(
            INSERT_QUERY,
            (
                taxonomy.id,
                taxonomy.parent_id,
                taxonomy.name,
                taxonomy.rank.value
            )
        )

    @classmethod
    def bulk_insert(cls, database_cursor, taxonomies: list) -> int:
        """
        Efficiently inserts multiple taxonomies.

        database_cursor
            Database cursor with open transaction.
        taxonomies : List[Taxonomy]
            Taxonomies for bulk insert.
        """
        BULK_INSERT_QUERY = (
            f"INSERT INTO {cls.TABLE_NAME} (id, parent_id, name, rank) "
            "VALUES %s ON CONFLICT DO NOTHING;"
        )
        # Bulk insert the new peptides
        execute_values(
            database_cursor,
            BULK_INSERT_QUERY,
            [
                (
                    taxonomy.id,
                    taxonomy.parent_id,
                    taxonomy.name,
                    taxonomy.rank.value
                ) for taxonomy in taxonomies
            ],
            page_size=len(taxonomies)
        )
        # rowcount is only accurate, because the page size is as high as the number of inserted data. If the page size would be smaller rowcount would only return the rowcount of the last processed page.
        return database_cursor.rowcount
        
    @classmethod
    def select(cls, database_cursor, select_conditions: tuple = ("", []), fetchall: bool = False):
        """
        Selects one or many taxonomies.

        Parameters
        ----------
        database_cursor : 
            Active database cursor
        select_conditions : Tuple[str, List[Any]]
             A tupel with the where statement (without WHERE) and a list of parameters, e.g. ("id = %s", [1])
        fetchall : bool
             Indicates if multiple rows should be fetched

        Returns
        -------
        Taxonomy or list of taxonomies
        """
        select_query = f"SELECT id, parent_id, name, rank FROM {cls.TABLE_NAME}"
        if len(select_conditions) == 2 and len(select_conditions[0]):
            select_query += f" WHERE {select_conditions[0]}"
        select_query += ";"
        database_cursor.execute(select_query, select_conditions[1])
        
        if fetchall:
            return [cls(row[0], row[1], row[2], row[3]) for row in database_cursor.fetchall()]
        else:
            row = database_cursor.fetchone()
            if row:
                return cls(row[0], row[1], row[2], row[3])
            else:
                return None

    def sub_species(self, database_cursor: DatabaseCursor) -> List[Taxonomy]:
        """
        Returns all sub taxonomies with rank TaxonomyRank.SPECIES including itself if itself has rank TaxonomyRank.SPECIES.

        Parameters
        ----------
        database_cursor : DatabaseCursor
            Database cursor

        Returns
        -------
        List[Taxonomy]
            List of taxonomies including the self
        """
        recursive_subspecies_id_query = (
            "WITH RECURSIVE subtaxonomies AS ("
                "SELECT id, parent_id, name, rank "
                f"FROM {self.__class__.TABLE_NAME} "
                "WHERE id = %s "
                "UNION " 
                    "SELECT t.id, t.parent_id, t.name, t.rank "
                    f"FROM {self.__class__.TABLE_NAME} t "
                    "INNER JOIN subtaxonomies s ON s.id = t.parent_id "
            f") SELECT id, parent_id, name, rank FROM subtaxonomies WHERE rank = %s;"
        )
        database_cursor.execute(recursive_subspecies_id_query, (self.id, TaxonomyRank.SPECIES.value))
        return [self.__class__(row[0], row[1], row[2], row[3]) for row in database_cursor.fetchall()]
