# std imports
from dataclasses import dataclass
from typing import Any, Iterable, List, Optional

# internal imports
from macpepdb.models.peptide_metadata import PeptideMetadata

@dataclass
class MetadataCondition:
    """
    Class to define metadata conditions for filtering peptides.
    """

    __slots__ = [
        "__is_swiss_prot",
        "__is_trembl",
        "__taxonomy_ids",
        "__unique_taxonomy_ids",
        "__proteome_id"
    ]

    __is_swiss_prot: Optional[bool]
    __is_trembl: Optional[bool]
    __taxonomy_ids: Optional[List[int]]
    __unique_taxonomy_ids: Optional[List[int]]
    __proteome_id: Optional[str]

    def __init__(self):
        self.__is_swiss_prot = None
        self.__is_trembl = None
        self.__taxonomy_ids = None
        self.__unique_taxonomy_ids = None
        self.__proteome_id = None

    @property
    def is_swiss_prot(self) -> Optional[bool]:
        return self.__is_swiss_prot

    @property
    def is_trembl(self) -> Optional[bool]:
        return self.__is_trembl

    @property
    def taxonomy_ids(self) -> Optional[List[int]]:
        return self.__taxonomy_ids

    @property
    def unique_taxonomy_ids(self) -> Optional[List[int]]:
        return self.__unique_taxonomy_ids

    @property
    def proteome_id(self) -> Optional[str]:
        return self.__proteome_id

    @is_swiss_prot.setter
    def is_swiss_prot(self, value: Optional[bool]):
        self.__is_swiss_prot = value

    @is_trembl.setter
    def is_trembl(self, value: Optional[bool]):
        self.__is_trembl = value

    @taxonomy_ids.setter
    def taxonomy_ids(self, value: Optional[List[int]]):
        self.__taxonomy_ids = value

    @unique_taxonomy_ids.setter
    def unique_taxonomy_ids(self, value: Optional[List[int]]):
        self.__unique_taxonomy_ids = value

    @proteome_id.setter
    def proteome_id(self, value: Optional[str]):
        self.__proteome_id = value


    def validate(self, metadata: PeptideMetadata) -> bool:
        if self.__is_swiss_prot is not None and metadata.is_swiss_prot != self.__is_swiss_prot:
            return False
        if self.__is_trembl is not None and metadata.is_trembl != self.__is_trembl:
            return False
        if self.__taxonomy_ids is not None and not self.__class__.is_intersecting(self.__taxonomy_ids, metadata.taxonomy_ids):
            return False
        if self.__unique_taxonomy_ids is not None and not self.__class__.is_intersecting(self.__unique_taxonomy_ids, metadata.unique_taxonomy_ids):
            return False
        if self.__proteome_id is not None and self.__proteome_id not in metadata.proteome_ids:
            return False
        return True

    def has_conditions(self) -> bool:
        """
        Checks if metadata conditions exists

        Returns
        -------
        bool
            False if no check is needed (not metadata condition)
        """
        return self.__is_swiss_prot is not None \
            or self.__is_trembl is not None \
            or self.__taxonomy_ids is not None \
            or self.__unique_taxonomy_ids is not None \
            or self.__proteome_id is not None \
        
    @classmethod
    def is_intersecting(cls, iterable_x: Iterable[Any], iterable_y: Iterable[Any]) -> bool:
        """
        Checks if two iterables are intersecting by checking if one element of iterable x is contained by iterable y.

        Arguments
        ---------
        iterable_x: Iterable[Any]
            List with elements
        iterable_y: Iterable[Any]
            List with elements
        
        Returns
        -------
        True if intersect
        """
        for x_element in iterable_x:
            if x_element in iterable_y:
                return True
        return False
