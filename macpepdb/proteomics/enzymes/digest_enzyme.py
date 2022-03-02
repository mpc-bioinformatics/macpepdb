# std imports
import re
from typing import ClassVar

# internal imports
from macpepdb.models import peptide
from macpepdb.models import protein
from macpepdb.proteomics.amino_acid import X as UnknwonAminoAcid, REPLACEABLE_AMBIGIOUS_AMINO_ACID_LOOKUP

class DigestEnzyme:
    """
    Defines a enzyme for protein digestion

    Parameters
    ----------
    name: str
        Name
    shortcut: str
        Shortcut
    cleavage_regex: re
        cleavage_regex to find cleavage sites
    max_number_of_missed_cleavages : int
        Maximum number of missed cleavages
    minimum_peptide_length : int
        Minimum peptide length
    maximum_peptide_length : int
        Maxiumum peptide length
    """

    MISSED_CLEAVAGE_REGEX: ClassVar[re.Pattern] = re.compile(r"[A-Z](?!$)")
    """Regex to count missed cleavages
    """
    
    def __init__(self, name: str = "Abstract Digest Enzym", shortcut: str = "", cleavage_regex: re = r".", max_number_of_missed_cleavages: int = 0, minimum_peptide_length: int = 0, maximum_peptide_length: int = 1):
        self.__name = name
        self.__shortcut = shortcut
        self.__cleavage_regex = cleavage_regex
        self.__max_number_of_missed_cleavages = max_number_of_missed_cleavages
        self.__minimum_peptide_length = minimum_peptide_length
        self.__maximum_peptide_length = maximum_peptide_length
        self.__peptide_range = range(minimum_peptide_length, maximum_peptide_length + 1) # upper border of range is excluded, so plus

    @property
    def max_number_of_missed_cleavages(self):
        """
        Returns
        -------
        Maximal number of missed cleavages.
        """
        return self.__max_number_of_missed_cleavages

    @property
    def minimum_peptide_length(self):
        """
        Returns
        -------
        Minimum peptide length.
        """
        return self.__minimum_peptide_length

    @property
    def maximum_peptide_length(self):
        """
        Returns
        -------
        Maximal peptide length.
        """
        return self.__maximum_peptide_length


    def digest(self, protein: 'protein.Protein') -> list:
        """
        Digests a protein.

        Parameters
        ----------
        protein : Protein
            Protein to digest

        Returns
        -------
        List of peptides.
        """
        peptides = set()
        # Split protein sequence on every cleavage position
        protein_parts = re.split(self.__cleavage_regex, protein.sequence)
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
                if len(peptide_sequence) in self.__peptide_range and not UnknwonAminoAcid.one_letter_code in peptide_sequence:
                    peptides.add(peptide.Peptide(peptide_sequence, missed_cleavage - part_index))
                    if self.__class__.is_sequence_containing_replaceable_ambigous_amino_acids(peptide_sequence):
                        # If there is a replaceable ambigous amino acid within the sequence, calculate each sequence combination of the actual amino acids
                        # Note: Some protein sequences in SwissProt and TrEMBL contain ambigous amino acids (B, Z). In most cases B and Z are denoted with the average mass of their encoded amino acids (D, N and E, Q).
                        # The average mass makes it difficult to create precise queries for these sequences in MaCPepDB. Therefor we calculates each differentiated version of the ambigous sequence and store it with the differentiated masses.
                        # This works only, when the actual amino acids have distinct masses like for B and Z, therefore we have to tolerate Js.
                        differentiated_sequences = self.__class__.differentiate_ambigous_sequences(peptide_sequence)
                        for sequence in differentiated_sequences:
                            peptides.add(peptide.Peptide(sequence, missed_cleavage - part_index))
        return list(peptides)

    @classmethod
    def get_enzyme_by_name(cls, name: str):
        """
        Returns a enzyme by name.

        Argument
        ========

        Returns
        -------
        DigestEnzym
        """
        # to prevent cyclic imports, import enzyms here not at the top
        from .trypsin import Trypsin
        if name.lower() == Trypsin.NAME.lower():
            return Trypsin
        # elif name.lower() == OtherEnzym.NAME.lower()
        raise NameError("Unknown enzym '{}'".format(name))

    @classmethod
    def get_known_enzymes(cls):
        """
        Returns a list of DigestEnzym names.

        Returns
        -------
        List of names.
        """
        from .trypsin import Trypsin
        return [
            Trypsin.NAME.lower()
        ]

    @classmethod
    def is_sequence_containing_replaceable_ambigous_amino_acids(cls, sequence: str) -> bool:
        """
        Determines if a sequence contains ambigous amino acids.

        Parameters
        ----------
        sequence : str
            Amino acid sequence.
        
        Returns
        
        -------
        True if the sequence contains ambigous amino acids, otherwise False
        """
        for one_letter_code in REPLACEABLE_AMBIGIOUS_AMINO_ACID_LOOKUP.keys():
            if one_letter_code in sequence:
                return True
        return False
     

    @classmethod
    def differentiate_ambigous_sequences(cls, ambigous_sequence: str, ) -> set:
        """
        Calculates all possible combinations of an ambigous sequence.

        Parameters
        ----------
        ambigous_sequence : str
            Amino acid sequence with ambigous amino acids.

        Returns
        Set of sequences
        """
        differentiated_sequences = set()
        cls.__differentiate_ambigous_sequences(ambigous_sequence, differentiated_sequences)
        return differentiated_sequences

    @classmethod
    def __differentiate_ambigous_sequences(cls, sequence: str, differentiated_sequences: set, position: int = 0):
        """
        Recursivly calculates all possible differentiate combinations of ambigous sequence.

        Parameters
        ----------
        sequence : str
            Current sequence
        differentiated_sequences : Set[str]
            A set() where the differentiated sequences where stored.
        position : int
            current position in the sequence
        """
        # If recursion reached the end of the sequence add it to the set of sequences
        if position == len(sequence):
            differentiated_sequences.add(sequence)
            return
        if not sequence[position] in REPLACEABLE_AMBIGIOUS_AMINO_ACID_LOOKUP:
            # If the amino acid on the current position is not a replaceable ambigous amino acid, pass the unchanged sequence to to next recursion level
            cls.__differentiate_ambigous_sequences(sequence, differentiated_sequences, position + 1)
        else:
            # If the amino acid on the current position is a replaceable ambigous amino acid create a new sequence where the current ambigous amino acid is sequentially replaced by its actual amino acids.
            # Than pass the new sequence to the next recursion level.
            for replacement_amino_acid in REPLACEABLE_AMBIGIOUS_AMINO_ACID_LOOKUP[sequence[position]]:
                new_sequence = sequence[:position] + replacement_amino_acid.one_letter_code + sequence[position + 1:]
                cls.__differentiate_ambigous_sequences(new_sequence, differentiated_sequences, position + 1)

    @classmethod
    def count_missed_cleavages(cls, sequence: str) -> int:
        """
        Counts missed cleavage in sequence.

        Parameters
        ----------
        sequence : str
            Amino acid sequence

        Returns
        -------
        int
            Number of missed cleavages
        """
        return len(cls.MISSED_CLEAVAGE_REGEX.findall(sequence))
