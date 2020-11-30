import pathlib
import csv
import operator
from math import floor

from sqlalchemy import or_, and_, between

from ..proteomics.modification_collection import ModificationCollection
from ..proteomics.mass.precursor_range import PrecursorRange

# Only a simple helper class to build the combination matrix in ModifiedPeptideQueryFilter
class ModificationCounter:
    def __init__(self, modification):
        self.modification = modification
        self.count = 0

class ModificatioCombination:
    def __init__(self, modification_counters: list):
        # Dict with amino acid one letter code as key and [count, is_static] as value, e.g.: "c": [4, True]
        self.__amino_acid_occurences = {}
        # Used terminus modification of the form [modification, is_applied], e.g.: [Modification, False]
        self.__n_terminus_modification = None
        self.__c_terminus_modification = None
        # Sum of the weight delta by the applied modifications
        self.__delta_sum = 0

        for counter in modification_counters:
            # If counter is for static non terminus modification or a variable modification with a count higher than zero. Add it to the occurances
            if counter.modification.is_static and not counter.modification.is_terminus_modification or (counter.modification.is_variable and counter.count):
                if counter.modification.amino_acid.one_letter_code in self.__amino_acid_occurences:
                    self.__amino_acid_occurences[counter.modification.amino_acid.one_letter_code][0] += counter.count
                else:
                    self.__amino_acid_occurences[counter.modification.amino_acid.one_letter_code] = [counter.count, counter.modification.is_static]

                if counter.modification.is_position_n_terminus:
                    self.__n_terminus_modification = [counter.modification, counter.count > 0]
                if counter.modification.is_position_c_terminus:
                    self.__c_terminus_modification = [counter.modification, counter.count > 0]
                
                self.__delta_sum += counter.count * counter.modification.delta
            elif counter.modification.is_static and counter.modification.is_terminus_modification:
                self.__delta_sum += counter.count * counter.modification.delta


    def to_where_clause(self, peptide_class, precursor: int, lower_precursor_tolerance_ppm: int, upper_precursor_tolerance_ppm: int):
        conditions = []
        for one_letter_code, count_and_type in self.__amino_acid_occurences.items():
            # Build the column name, e.g. a_count
            column_name = "{}_count".format(one_letter_code.lower())
            # Get the related column attribute by string
            column = getattr(peptide_class, column_name)
            # Add condition. In case the modification is fix we want equals the amino acid, if the modification is variable it could be more.
            if count_and_type[1]:
                conditions.append(column == count_and_type[0])
            else:
                conditions.append(column >= count_and_type[0])
        
        # For n- and c-terminus modification check if there is one:
        # * If the modification is variable and applied check for presence (absence is uninteresting because it is a variable modification)
        if self.__n_terminus_modification and (self.__n_terminus_modification[0].is_variable and self.__n_terminus_modification[1]):
            conditions.append(peptide_class.n_terminus == self.__n_terminus_modification[0].amino_acid.one_letter_code)

        if self.__c_terminus_modification and (self.__c_terminus_modification[0].is_variable and self.__c_terminus_modification[1]):
            conditions.append(peptide_class.c_terminus == self.__c_terminus_modification[0].amino_acid.one_letter_code)

        # Add the weight between condition
        precursor_range = PrecursorRange(precursor - self.__delta_sum, lower_precursor_tolerance_ppm, upper_precursor_tolerance_ppm)
        conditions.append(between(peptide_class.weight, precursor_range.lower_limit, precursor_range.upper_limit))
        # Concat all conditions with "and"
        return and_(*conditions)


class ModifiedPeptideWhereClauseBuilder:
    def __init__(self, modification_collection: ModificationCollection, precursor: int, lower_precursor_tolerance_ppm: int, upper_precursor_tolerance_ppm: int, variable_modification_maximum: int):
        self.__precursor = precursor
        self.__lower_precursor_tolerance_ppm = lower_precursor_tolerance_ppm
        self.__upper_precursor_tolerance_ppm = upper_precursor_tolerance_ppm
        self.__modification_counter = [ModificationCounter(modification) for modification in modification_collection.all]
        self.__modification_combinations = []
        self.__build_combinations(0, self.__precursor, variable_modification_maximum, False, False, False, False)


    def __build_combinations(self, counter_idx: int, remaining_precursor_mass: int, free_variable_modifications: int, is_n_terminus_used: bool, is_n_terminal_residue_used: bool, is_c_terminus_used: bool, is_c_terminal_residue_used: bool):
        # Exit method if index is greater than number of counters
        if counter_idx >= len(self.__modification_counter):
            return None

        # Save acid and mod for cleaner code
        modification = self.__modification_counter[counter_idx].modification

        # How many of the current acids with modification will apply to the current precursor mass
        if modification.is_static and not modification.is_terminus_modification:
            # Static vars apply any time the acid occues (except termini)
            mod_max_count = floor(remaining_precursor_mass / modification.mono_mass)
        elif modification.is_variable and not modification.is_terminus_modification:
            # Variable modifications apply only if there is "space" for more variable modifications
            # and the number of used modification, regardless of variable or not, has to fit
            mod_max_count = min(
                floor(remaining_precursor_mass / modification.mono_mass),
                free_variable_modifications
            )
        elif modification.is_static and modification.is_position_n_terminus:
            # 0 if terminus is alredy in use, 1 if the terminus is free
            # and the modification has to fit
            mod_max_count = min(
                0 if is_n_terminus_used else 1,
                floor(remaining_precursor_mass / modification.mono_mass)
            )
        elif modification.is_static and modification.is_position_c_terminus:
            # 0 if terminus is alredy in use, 1 if the terminus is free
            # and the modification has to fit
            mod_max_count = min(
                0 if is_c_terminus_used else 1,
                floor(remaining_precursor_mass / modification.mono_mass)
            )
        elif modification.is_variable and modification.is_position_n_terminus:
            # 0 if terminus is already in use or there is no more space for a variable modification left
            # and the modification has to fit
            mod_max_count = min(
                0 if is_n_terminal_residue_used or not free_variable_modifications else 1,
                floor(remaining_precursor_mass / modification.mono_mass)
            )
        elif modification.is_variable and modification.is_position_c_terminus:
            # 0 if terminus is already in use or there is no more space for a variable modification left
            # and the modification has to fit
            mod_max_count = min(
                0 if is_c_terminal_residue_used or not free_variable_modifications else 1,
                floor(remaining_precursor_mass / modification.mono_mass)
            )


        is_precursor_reached = False

        # Increase the counter for the current modification until maximum is reached
        # +1 to make range inclusive
        for count in range(mod_max_count + 1):
            # Reset all following modification counts to zero
            for i in range(counter_idx + 1, len(self.__modification_counter)):
                self.__modification_counter[i].count = 0

            # Calculate the remaining precursor mass for the following modifications
            next_remaining_precursor = remaining_precursor_mass - (modification.mono_mass * count)

            # Check if remaining precursor has space for more modifications
            if next_remaining_precursor > 0:
                # Assign the count for the current mod to the counter
                self.__modification_counter[counter_idx].count = count

                # calculate free_variable_modifications for next iteration
                if modification.is_static:
                    next_free_variable_modifications = free_variable_modifications
                else:
                    next_free_variable_modifications = free_variable_modifications - count


                # set next terminus/temrinal residue used to last value
                next_is_n_terminal_residue_used = is_n_terminal_residue_used
                next_is_c_terminal_residue_used = is_c_terminal_residue_used
                next_is_n_terminus_used = is_n_terminus_used
                next_is_c_terminus_used = is_c_terminus_used
                # check if n-terminal-residue is used
                if modification.is_variable and modification.is_position_n_terminus and count:
                    next_is_n_terminal_residue_used = True
                # check if c-terminal-residue is used
                elif modification.is_variable and modification.is_position_c_terminus and count:
                    next_is_c_terminal_residue_used = True
                # check if n-terminus is used
                elif modification.is_static and modification.is_position_n_terminus and count:
                    next_is_n_terminus_used = True
                # check if c-terminus is used
                elif modification.is_static and modification.is_position_c_terminus and count:
                    next_is_c_terminus_used = True


                # start the next iteration
                self.__build_combinations(counter_idx + 1, next_remaining_precursor, next_free_variable_modifications, next_is_n_terminus_used, next_is_n_terminal_residue_used, next_is_c_terminus_used, next_is_c_terminal_residue_used)
            else:
                is_precursor_reached = True

            # Add current counter state to matrix, if getcurrent counter is last counter or precursor is reached
            if counter_idx == len(self.__modification_counter) - 1 or is_precursor_reached:
                self.__modification_combinations.append(ModificatioCombination(self.__modification_counter))
                # line = []
                # for counter in self.__modification_counter:
                #     line.append("{} => {}".format(counter.modification.amino_acid.one_letter_code, counter.count))
                # print("; ".join(line))


            # Stop iteration if precursor is reached
            if is_precursor_reached:
                break


    def build(self, peptide_class):
        if len(self.__modification_combinations):
            conditions = []
            for combination in self.__modification_combinations:
                conditions.append(combination.to_where_clause(peptide_class, self.__precursor, self.__lower_precursor_tolerance_ppm, self.__upper_precursor_tolerance_ppm))
            # Concatenate all conditions with "or"
            return or_(*conditions)
        else:
            precursor_range = PrecursorRange(self.__precursor, self.__lower_precursor_tolerance_ppm, self.__upper_precursor_tolerance_ppm)
            return between(peptide_class.weight, precursor_range.lower_limit, precursor_range.upper_limit)