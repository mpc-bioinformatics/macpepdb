from .proteomics.modification_collection import ModificationCollection
from .proteomics.modification import Modification
from .proteomics.mass.precursor_range import PrecursorRange

class PeptideMassValidator:
    def __init__(self, modification_collection: ModificationCollection, maximum_number_of_variable_modifications: int, precursor_range: PrecursorRange):
        self.__modification_collection = modification_collection
        self.__maximum_number_of_variable_modifications = maximum_number_of_variable_modifications
        self.__precursor_range = precursor_range
        self.__static_modifications = {}
        self.__variable_modifications = {}
        self.__variable_n_terminus_modifications = {}
        self.__variable_c_terminus_modifications = {}

        # Create lookup table for static non-terminus modifications
        for modification in self.__modification_collection.static:
            if not modification.is_terminus_modification:
                self.__static_modifications[modification.amino_acid.one_letter_code] = modification

        # Create lookup tables for variable modifications
        for modification in self.__modification_collection.variable:
            if not modification.is_terminus_modification:
                self.__variable_modifications.setdefault(modification.amino_acid.one_letter_code, []).append(modification)
            elif modification.is_position_n_terminus:
                self.__variable_n_terminus_modifications.setdefault(modification.amino_acid.one_letter_code, []).append(modification)
            elif modification.is_position_c_terminus:
                self.__variable_c_terminus_modifications.setdefault(modification.amino_acid.one_letter_code, []).append(modification)

        # Set to actual values in self.validate()
        self.__peptide = None
        self.__applied_modifications = None

    def set_precursor_range(self, precursor_range: PrecursorRange):
        """
        Set new precursor range
        @param precursor_range Instance of PrecursorRange
        """
        self.__precursor_range = precursor_range

    def set_maximum_number_of_variable_modifications(self, maximum_number_of_variable_modifications: int):
        """
        Sets a new limit for variable modifications
        @param maximum_number_of_variable_modifications New limit
        """
        self.__maximum_number_of_variable_modifications = maximum_number_of_variable_modifications

    def __current_peptide_mass(self) -> int:
        """
        Returns the current peptide mass inclusively the applied modifications
        @return int Returns the mass of the modified peptide
        """
        modifications_delta_sum = sum(modification.delta for modification in self.__applied_modifications if modification)
        return self.__peptide.mass + modifications_delta_sum 

    def __validate(self, current_mod_idx: int, number_of_variable_modifications: int) -> bool:
        """
        Recursivly try all variable modifications until 
        @param current_mod_idx Current modification index
        @param number_of_variable_modifications Current number of variable modifications
        @return bool True if the peptide + a combindation of modifications matches the precursor range 
        """
        # Return True if mass matches
        if self.__current_peptide_mass() in self.__precursor_range:
            return True
        # Return False if index is out of bounds
        if current_mod_idx == len(self.__applied_modifications):
            return False
        # Skip modification adjustments if current modification is static and go to next 
        if self.__applied_modifications[current_mod_idx] and self.__applied_modifications[current_mod_idx].is_static:
            return self.__validate(current_mod_idx + 1, number_of_variable_modifications)
        # Return false if no more variable modification is allowed
        if number_of_variable_modifications == self.__maximum_number_of_variable_modifications:
            return False

        # Try all variable non-terminus modification
        if 0 < current_mod_idx < len(self.__applied_modifications) - 1:
            amino_acid = self.__peptide.sequence[current_mod_idx - 1]
            for modification in self.__variable_modifications.get(amino_acid, []):
                self.__applied_modifications[current_mod_idx] = modification
                if self.__validate(current_mod_idx + 1, number_of_variable_modifications + 1):
                    return True
        # Try all variable n-terminus modification
        elif current_mod_idx == 0:
            amino_acid = self.__peptide.sequence[0]
            for modification in self.__variable_n_terminus_modifications.get(amino_acid, []):
                self.__applied_modifications[current_mod_idx] = modification
                if self.__validate(current_mod_idx + 1, number_of_variable_modifications + 1):
                    return True
        # Try all variable c-terminus modifications
        elif current_mod_idx == len(self.__applied_modifications) - 1:
            amino_acid = self.__peptide.sequence[-1]
            for modification in self.__variable_c_terminus_modifications.get(amino_acid, []):
                self.__applied_modifications[current_mod_idx] = modification
                if self.__validate(current_mod_idx + 1, number_of_variable_modifications + 1):
                    return True

        # Remove this variable modification
        self.__applied_modifications[current_mod_idx] = None
        # Try variable modification on next amino acid
        return self.__validate(current_mod_idx + 1, number_of_variable_modifications)



    def validate(self, peptide) -> bool:
        """
        Checks if the given peptide can be modified to match the precursor range.
        @param peptide Object based on PeptideBase
        @return boolean True if the peptide + a combination of modification matches the precursor, if not False
        """

        # Set peptide
        self.__peptide = peptide
        self.__applied_modifications = [None] * (peptide.length + 2) # +2 for n- and c-terminus

        # Set static modifications
        for idx, amino_acid in enumerate(peptide.sequence):
            self.__applied_modifications[idx + 1] = self.__static_modifications.get(amino_acid, None)

        if self.__modification_collection.static_n_terminus_modifications:
            self.__applied_modifications[0] = self.__modification_collection.static_n_terminus_modifications

        if self.__modification_collection.static_c_terminus_modifications:
            self.__applied_modifications[-1] = self.__modification_collection.static_c_terminus_modifications

        is_valid = self.__validate(0, 0)

        # Reset variables
        self.__peptide = None
        self.__applied_modifications = None
        self.__number_of_variable_modifications = 0

        return is_valid


