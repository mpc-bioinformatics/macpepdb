import re

from ...models import peptide
from ...models import protein

class DigestEnzyme:
    def __init__(self, name: str = "Abstract Digest Enzym", shortcut: str = "", regex: re = r".", max_number_of_missed_cleavages: int = 0, minimum_peptide_length: int = 0, maximum_peptide_length: int = 1):
        self.__name = name
        self.__shortcut = shortcut
        self.__regex = regex
        self.__max_number_of_missed_cleavages = max_number_of_missed_cleavages
        self.__minimum_peptide_length = minimum_peptide_length
        self.__maximum_peptide_length = maximum_peptide_length
        self.__peptide_range = range(minimum_peptide_length, maximum_peptide_length + 1) # upper border of range is excluded, so plus

    @property
    def max_number_of_missed_cleavages(self):
        return self.__max_number_of_missed_cleavages

    @property
    def minimum_peptide_length(self):
        return self.__minimum_peptide_length

    @property
    def maximum_peptide_length(self):
        return self.__maximum_peptide_length


    def digest(self, protein: 'protein.Protein') -> list:
        peptides = set()
        # Split protein sequence on every cleavage position
        protein_parts = re.split(self.__regex, protein.sequence)
        # Start with every part
        for part_index in range(0, len(protein_parts)):
            # Check if end of protein_parts is reached before the last missed cleavage (prevent overflow)
            last_part_to_add = min(
                part_index + self.__max_number_of_missed_cleavages + 1,
                len(protein_parts)
            )
            peptide_sequence = ""
            for missed_cleavage in range(part_index, last_part_to_add):
                peptide_sequence += protein_parts[missed_cleavage]
                if len(peptide_sequence) in self.__peptide_range and peptide_sequence.count('X') == 0:
                    peptides.add(peptide.Peptide(peptide_sequence, missed_cleavage - part_index))

        return list(peptides)

    @classmethod
    def get_enzyme_by_name(cls, name: str):
        # to prevent cyclic imports, import enzyms here not at the top
        from .trypsin import Trypsin
        if name.lower() == Trypsin.NAME.lower():
            return Trypsin
        # elif name.lower() == OtherEnzym.NAME.lower()
        raise NameError("Unknown enzym '{}'".format(name))

    @classmethod
    def get_known_enzymes(cls):
        from .trypsin import Trypsin
        return [
            Trypsin.NAME.lower()
        ]