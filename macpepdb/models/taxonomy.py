from __future__ import annotations
from enum import IntEnum, unique

@unique
class TaxonomyRank(IntEnum):
    # Enum names created by rank.upper().replace(" ", "_")
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
        rank = rank.upper().replace("_", " ")
        if rank in cls.__members__:
            return cls.__members__[rank]

class Taxonomy:
    TABLE_NAME = "taxonomies"

    def __init__(self, id: int, parent_id: int, name: str, rank: TaxonomyRank):
        self.id = id
        self.parent_id = parent_id
        self.name = name
        self.rank = rank

    def parent(self, database_cursor):
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

    # This method is implemented to use bot accessions as hash when a protein merge is stored in a hashable collection (Set, Dictionary, ...)
    def __hash__(self):
        return hash(self.id)
    
    # According to the Python documentation this should be implemented if __hash__() is implemented.
    def __eq__(self, other):
        return self.id == other.id


    @staticmethod
    def insert(database_cursor, taxonomy: Taxonomy):
        """
        @param database_cursor Database cursor
        @param taxonomy Taxonomy to insert
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