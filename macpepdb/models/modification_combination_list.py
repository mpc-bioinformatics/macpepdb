from __future__ import annotations
import pathlib
import csv
from math import floor

from ..proteomics.modification_collection import ModificationCollection
from ..proteomics.mass.precursor_range import PrecursorRange

# Only a simple helper class to build the combination matrix in ModifiedPeptideQueryFilter
class ModificationCounter:
    def __init__(self, modification):
        self.modification = modification
        self.count = 0

class ModificatioCombination:
    """
    This class keeps all information for querying peptides with a specific modification combination.
    The property `precursor_range` holds the PrecursorRange (mass tolerance) for the query and `column_conditions` contains a dictionary,
    where the key is the column name and the value is a tupel which keeps the operator and the value, e.g. `(operator, value)`
    """
    def __init__(self, modification_counters: list, precursor: int, lower_precursor_tolerance_ppm: int, upper_precursor_tolerance_ppm: int):
        self.__column_conditions = {}
        self.__precursor_range = (0, 0)
        self.__calculate_columns_and_precursor_range(modification_counters, precursor, lower_precursor_tolerance_ppm, upper_precursor_tolerance_ppm)

    @property
    def column_conditions(self):
        return self.__column_conditions

    @property
    def precursor_range(self):
        return self.__precursor_range

    def __calculate_columns_and_precursor_range(self, modification_counters: list, precursor: int, lower_precursor_tolerance_ppm: int, upper_precursor_tolerance_ppm: int):
        # Dict with amino acid one letter code as key and [count, is_static] as value, e.g.: "c": [4, True]
        amino_acid_occurences = {}
        # Used terminus modification of the form [modification, is_applied], e.g.: [Modification, False]
        n_terminus_modification = None
        c_terminus_modification = None
        # Sum of the weight delta by the applied modifications
        delta_sum = 0

        for counter in modification_counters:
            # If counter is for static non terminus modification or a variable modification with a count higher than zero. Add it to the occurances
            if counter.modification.is_static and not counter.modification.is_terminus_modification or (counter.modification.is_variable and counter.count):
                if counter.modification.amino_acid.one_letter_code in amino_acid_occurences:
                    amino_acid_occurences[counter.modification.amino_acid.one_letter_code][0] += counter.count
                else:
                    amino_acid_occurences[counter.modification.amino_acid.one_letter_code] = [counter.count, counter.modification.is_static]

                if counter.modification.is_position_n_terminus:
                    n_terminus_modification = [counter.modification, counter.count > 0]
                if counter.modification.is_position_c_terminus:
                    c_terminus_modification = [counter.modification, counter.count > 0]
                
                delta_sum += counter.count * counter.modification.delta
            elif counter.modification.is_static and counter.modification.is_terminus_modification:
                delta_sum += counter.count * counter.modification.delta

        for one_letter_code, count_and_type in amino_acid_occurences.items():
            # Build the column name, e.g. a_count
            column_name = "{}_count".format(one_letter_code.lower())
            # Add condition. In case the modification is fix we want equals the amino acid, if the modification is variable it could be more.
            self.__column_conditions[column_name] = (
                '=' if count_and_type[1] else '>=',
                count_and_type[0]
            )

        # For n- and c-terminus modification check if there is one:
        # * If the modification is variable and applied check for presence (absence is uninteresting because it is a variable modification)
        if n_terminus_modification and (n_terminus_modification[0].is_variable and n_terminus_modification[1]):
            self.__column_conditions['n_terminus'] = (
                '=',
                n_terminus_modification[0].amino_acid.one_letter_code
            )

        if c_terminus_modification and (c_terminus_modification[0].is_variable and c_terminus_modification[1]):
            self.__column_conditions['c_terminus'] = (
                '=',
                c_terminus_modification[0].amino_acid.one_letter_code
            )

        # Add the weight between condition
        self.__precursor_range = PrecursorRange(precursor - delta_sum, lower_precursor_tolerance_ppm, upper_precursor_tolerance_ppm)

    def to_sql(self) -> tuple:
        # Begin with weight
        query_string = "weight BETWEEN %s AND %s"
        values = [self.precursor_range.lower_limit, self.precursor_range.upper_limit]
        # Add each additional column condition
        for column_name, (operator, value) in self.column_conditions.items():
            query_string += " AND "
            query_string += f"{column_name} {operator} %s"
            values.append(value)
        return (query_string, values)

class ModificatioCombinationList:
    """
    This class creates all possible combination (ModificationCombination) from a given ModificationCollection for querying peptides.
    The class itself is iterable.
    """
    def __init__(self, modification_collection: ModificationCollection, precursor: int, lower_precursor_tolerance_ppm: int, upper_precursor_tolerance_ppm: int, variable_modification_maximum: int):
        self.__precursor = precursor
        self.__lower_precursor_tolerance_ppm = lower_precursor_tolerance_ppm
        self.__upper_precursor_tolerance_ppm = upper_precursor_tolerance_ppm
        self.__modification_counter = [ModificationCounter(modification) for modification in modification_collection.all]
        self.__modification_combinations = []
        self.__modification_combinations_iter = None
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
                self.__modification_combinations.append(ModificatioCombination(self.__modification_counter, self.__precursor, self.__lower_precursor_tolerance_ppm, self.__upper_precursor_tolerance_ppm))


            # Stop iteration if precursor is reached
            if is_precursor_reached:
                break

    def __iter__(self) -> ModificatioCombinationList:
        self.__modification_combinations_iter = iter(self.__modification_combinations)
        return self

    def __next__(self) -> ModificatioCombination:
        return next(self.__modification_combinations_iter)

    def __len__(self):
        return len(self.__modification_combinations)

    def to_sql(self) -> tuple:
        if len(self):
            query_string = []
            values = []
            for combination in self:
                combination_query_string, combination_values = combination.to_sql()
                query_string.append(combination_query_string)
                values += combination_values
            return (" OR ".join(query_string), values)
        else:
            precursor_range = PrecursorRange(self.__precursor, self.__lower_precursor_tolerance_ppm, self.__upper_precursor_tolerance_ppm)
            return ("weight BETWEEN %s AND %s", [precursor_range.lower_limit, precursor_range.upper_limit])