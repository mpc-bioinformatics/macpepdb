# std imports
from collections import defaultdict
from dataclasses import dataclass
from typing import List, Dict, Optional, Set

# internal imports
from macpepdb.models.peptide_base import PeptideBase
from macpepdb.proteomics.modification import Modification
from macpepdb.proteomics.modification_collection import ModificationCollection
from macpepdb.proteomics.mass.precursor_range import PrecursorRange

@dataclass
class VariableModificationCombination:
    """
    Helper class for PeptideMassValidator. Represents a possible combination of variable modifications.

    Parameters
    -------
    modifications
        List of variable modifications
    """

    __slots__ = [
        "__modifications",
        "__anywhere_modifications",
        "__n_terminus_modification",
        "__c_terminus_modification",
        "__hash"
    ]

    __modifications: List[Modification]
    __anywhere_modifications: List[Modification]
    __n_terminus_modification: Optional[Modification]
    __c_terminus_modification: Optional[Modification]
    __hash: int

    def __init__(self, modifications: List[Modification]):
        self.__modifications = sorted(modifications, key=lambda mod: mod.accession)
        self.__anywhere_modifications = []
        self.__n_terminus_modification = None
        self.__c_terminus_modification = None
        for mod in modifications:
            if mod.is_position_anywhere:
                self.__anywhere_modifications.append(mod)
            elif mod.is_position_n_terminus:
                self.__n_terminus_modification = mod
            elif mod.is_position_c_terminus:
                self.__c_terminus_modification = mod

        self.__hash = hash("".join(f"{mod.accession}{str(mod.delta)}{mod.amino_acid.name}" for mod in self.__modifications))

    @property
    def modifications(self) -> List[Modification]:
        """
        Returns
        -------
        List[Modification]
            Modifications
        """
        return self.__modifications

    @property
    def hash(self) -> int:
        """
        Returns
        -------
        int
            Hash of sorted modification accesions
        """
        return self.__hash

    def __hash__(self) -> int:
        return self.__hash

    def __eq__(self, other):
        return self.__hash == other.hash

    def check_peptide_fits(self, peptide: PeptideBase, modification_slots: List[Optional[Modification]]) -> bool:
        """
        Try to apply this combination of variable modifications to the given peptide.

        Parameters
        ----------
        peptide : PeptideBase
            A peptide

        Returns
        -------
        bool
            True if modification combination is applicable
        """
        # Check if n_terminus position is present and fit
        if self.__n_terminus_modification is not None \
            and peptide.sequence[0] == self.__n_terminus_modification.amino_acid.one_letter_code:
            modification_slots[0] = self.__n_terminus_modification
        # Check if c_terminus position is present and fit
        if self.__c_terminus_modification is not None \
            and peptide.sequence[-1] == self.__c_terminus_modification.amino_acid.one_letter_code:
            modification_slots[-1] = self.__c_terminus_modification
        # Check for every anywhere position if it has a palce
        for mod in self.__anywhere_modifications:
            found_position = False
            for amino_acid_idx, amino_acid in enumerate(peptide.sequence):
                # If modification fits amino acid and modification slot is free set it
                if amino_acid == mod.amino_acid.one_letter_code \
                    and modification_slots[amino_acid_idx] is None:
                    modification_slots[amino_acid_idx] = mod
                    found_position = True
                    break
            if not found_position:
                return False
        
        return True

class PeptideMassValidator:
    """
    Checks if the mass of a peptide fits into the mass/modification-combination of the given ModificationCollection and PrecursorRange.
    Mainly used or testing.

    Parameters
    ----------
    modification_collection : ModificationCollection
        Collection of modifications
    maximum_number_of_variable_modifications : int
        Maximum number of variable modifications per peptide
    precursor_range : PrecursorRange
        Precursor / mass range
    """

    def __init__(self, modification_collection: ModificationCollection, maximum_number_of_variable_modifications: int, precursor_range: PrecursorRange):
        self.__modification_collection = modification_collection
        self.__maximum_number_of_variable_modifications = maximum_number_of_variable_modifications
        self.__precursor_range = precursor_range

        self.__static_anywhere_modifications = {
            modification.amino_acid.one_letter_code: modification 
            for modification in self.__modification_collection.static_anywhere_modifications
        }
        self.__static_n_terminus_modification = self.__modification_collection.static_n_terminus_modification
        self.__static_c_terminus_modification = self.__modification_collection.static_c_terminus_modification


        self.__applied_variable_modifications_delta: Dict[int, Set[VariableModificationCombination]] =  defaultdict(set)

        self.__create_applied_modifications_matrix()


    @property
    def precursor_range(self) -> PrecursorRange:
        """
        Returns
        -------
        PrecursorRange
            Precursor range
        """
        return self.__precursor_range

    @precursor_range.setter
    def precursor_range(self, precursor_range: PrecursorRange):
        """
        Set precursor_nrange

        Parameters
        ----------
        precursor_range : PrecursorRange
            New precursor range
        """
        self.__precursor_range = precursor_range

    def __create_applied_modifications_matrix(self):
        applied_modifications_comb = [-1] * self.__maximum_number_of_variable_modifications
        self.__create_applied_modifications_matrix_recursive(applied_modifications_comb, 0)


    def __create_applied_modifications_matrix_recursive(self, applied_modifications_comb: List[int], comb_pos: int):
        for modification_index in range(-1, len(self.__modification_collection.variable_modifications)):
            # Check if modification is valid to add (no double n-/c-terminuns modifications)
            if not self.__check_if_modification_is_valid_to_add(applied_modifications_comb, modification_index):
                continue
            # Add modification
            applied_modifications_comb[comb_pos] = modification_index
            # If position in combination is not at the end of combination add the next modification
            if comb_pos < len(applied_modifications_comb) - 1:
                self.__create_applied_modifications_matrix_recursive(applied_modifications_comb, comb_pos + 1)
            else:
                self.__applied_variable_modifications_delta[
                    self.__get_modification_combination_delta(applied_modifications_comb)
                ].add(VariableModificationCombination(
                    [self.__modification_collection.variable_modifications[idx] for idx in applied_modifications_comb if idx >= 0]
                ))

    def __get_modification_combination_delta(self, applied_modifications_comb: List[int]) -> int:
        """
        Calculates the total delta of given modifications

        Parameters
        ----------
        applied_modifications_comb : List[int]
            List of variable modifications indexes.

        Returns
        -------
        int
            Mass delta
        """
        mass = 0
        for idx in applied_modifications_comb:
            if idx >= 0:
                mass += self.__modification_collection.variable_modifications[idx].delta
        return mass

    def __check_if_modification_is_valid_to_add(self, applied_modifications_comb: List[int], modification_index: int) -> bool:
        """
        Checks if a modification is valid to add it to the modification combination. E.g. no redundant terminus modifications

        Parameters
        ----------
        applied_modifications_comb : List[int]
            List of variable modification indexes
        modification_index : int
            Index of modification to be added

        Returns
        -------
        bool
            Ture if modification can be added
        """
        if modification_index < 0:
            return True
        modification_to_add = self.__modification_collection.variable_modifications[modification_index]
        if modification_to_add.is_position_anywhere:
            return True
        elif modification_to_add.is_position_n_terminus:
            for idx in applied_modifications_comb:
                if idx >= 0 and self.__modification_collection.variable_modifications[idx].is_position_n_terminus:
                    return False
            return True
        elif modification_to_add.is_position_c_terminus:
            for idx in applied_modifications_comb:
                if idx >= 0 and self.__modification_collection.variable_modifications[idx].is_position_c_terminus:
                    return False
            return True

    def __create_sequence_with_modification_markers(self, sequence: str, static_n_terminus_modificaton: Optional[Modification], 
        static_c_terminus_modificaton: Optional[Modification], modification_slots: List[Optional[Modification]]) -> str:
        """
        Returns the amino acid sequence with modification indicators, e.g. Peptide CIYLMVVMIYLTHA
        [static n-terminus modification delta].C[s:delta of static modification of C]IYLMVVM[v: delta of variable modification of M]IYLTHA.[static c-terminus modification delta]

        Parameters
        ----------
        sequence : str
            Amino acid sequence
        static_n_terminus_modificaton : Optional[Modification]
            Static n terminus modification
        static_c_terminus_modificaton : Optional[Modification]
            Static c terminuns modification
        modification_slots : List[Option[Modification]]
            List od applied modifications

        Returns
        -------
        str
            Sequence with modification markers
        """
        sequence_with_modification_markers = ""
        if static_n_terminus_modificaton is not None:
            sequence_with_modification_markers += f"[{static_n_terminus_modificaton.delta}]."
        
        for amino_acid_idx, amino_acid in enumerate(sequence):
            sequence_with_modification_markers += amino_acid
            if modification_slots[amino_acid_idx] is not None:
                mod = modification_slots[amino_acid_idx]
                sequence_with_modification_markers += f"[{'v' if mod.is_variable else 's'}:{mod.delta}]"

        if static_c_terminus_modificaton is not None:
            sequence_with_modification_markers += f".[{static_c_terminus_modificaton.delta}]"

        return sequence_with_modification_markers


    def validate(self, peptide: PeptideBase, add_sequence_with_modification_markers: bool = False) -> bool:
        """
        Validates if the given peptides matches the precursor.

        Parameters
        ----------
        peptide : PeptideBase
            _description_

        Returns
        -------
        bool
            _description_
        """
        mass = peptide.mass
        static_n_terminus_modificaton = None
        static_c_terminus_modificaton = None
        modification_slots = [None] * len(peptide.sequence)

        # Calculate mass inclusive static modifications
        if self.__static_n_terminus_modification is not None:
            static_n_terminus_modificaton = self.__static_n_terminus_modification
            mass += self.__static_n_terminus_modification.delta
        if self.__static_c_terminus_modification is not None:
            static_c_terminus_modificaton = self.__static_c_terminus_modification
            mass += self.__static_c_terminus_modification.delta
        for idx, one_letter_code in enumerate(peptide.sequence):
            if one_letter_code in self.__static_anywhere_modifications:
                modification_slots[idx] = self.__static_anywhere_modifications[one_letter_code]
                mass += self.__static_anywhere_modifications[one_letter_code].delta

        for delta, mod_combinations in self.__applied_variable_modifications_delta.items():
            # Check if mass + delta fits into the precursor range
            if mass + delta in self.__precursor_range:
                for combination in mod_combinations:
                    # Check if combination can be applied to peptide
                    if combination.check_peptide_fits(peptide, modification_slots):
                        # Add sequence with modifications if necessary
                        if add_sequence_with_modification_markers:
                            peptide.sequence_with_modification_markers = self.__create_sequence_with_modification_markers(
                                peptide.sequence,
                                static_n_terminus_modificaton,
                                static_c_terminus_modificaton,
                                modification_slots
                            )
                        return True
                    else:
                        modification_slots = [slot for slot in modification_slots if slot is None or slot.is_static]
        if add_sequence_with_modification_markers:
            peptide.sequence_with_modification_markers = self.__create_sequence_with_modification_markers(
                peptide.sequence,
                static_n_terminus_modificaton,
                static_c_terminus_modificaton,
                modification_slots
            )
        return False
                    