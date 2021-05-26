import pathlib
import csv

from .modification import Modification

class ModificationLimitError(BaseException):
    pass

class InvalidModificationCombinationError(BaseException):
    pass


class ModificationCollection:
    # Due to Comet this limitations applied
    MAX_VARIABLE_MODIFICATIONS = 9
    MAX_STATIC_N_TERMINUS_MODIFICATIONS = 1
    MAX_STATIC_C_TERMINUS_MODIFICATIONS = 1

    def __init__(self, modifications: list):
        variable_modification_counter = 0
        static_n_terminus_counter = 0
        static_c_terminus_counter = 0
        self.__modifications = modifications
        self.__variable_modifications = []
        self.__static_modifications = []
        self.__static_n_terminus_modifications = None
        self.__static_c_terminus_modifications = None
        self.__variable_n_terminus_modifications = []
        self.__variable_c_terminus_modifications = []
        # Sort the given modifications by their position into different array for direct access
        for modification in self.__modifications:
            if modification.is_static and not modification.is_terminus_modification:
                self.__static_modifications.append(modification)
            elif modification.is_static and modification.is_position_n_terminus:
                self.__static_modifications.append(modification)
                static_n_terminus_counter += 1
                self.__static_n_terminus_modifications = modification
            elif modification.is_static and modification.is_position_c_terminus:
                self.__static_modifications.append(modification)
                static_c_terminus_counter += 1
                self.__static_c_terminus_modifications = modification
            elif modification.is_variable:
                self.__variable_modifications.append(modification)
                variable_modification_counter += 1
                if modification.is_position_n_terminus:
                    self.__variable_n_terminus_modifications.append(modification)
                elif modification.is_position_c_terminus:
                    self.__variable_c_terminus_modifications.append(modification)

        # !!! TODO check if there are variable and static modifications for the same amino acids
        for static_mod in self.__static_modifications:
            for variable_mod in self.__variable_modifications:
                if static_mod.amino_acid == variable_mod.amino_acid:
                    raise InvalidModificationCombinationError("Static and variable modification for the same amino acid found:\n{}\n\n{}".format(static_mod, variable_mod))

        if variable_modification_counter > self.__class__.MAX_VARIABLE_MODIFICATIONS:
            raise ModificationLimitError("Only {} variable modifications are allowed.".format(self.__class__.MAX_VARIABLE_MODIFICATIONS))

        if static_n_terminus_counter > self.__class__.MAX_STATIC_N_TERMINUS_MODIFICATIONS:
            raise ModificationLimitError("Only {} static n-terminus modification are allowed.".format(self.__class__.MAX_STATIC_N_TERMINUS_MODIFICATIONS))

        if static_c_terminus_counter > self.__class__.MAX_STATIC_C_TERMINUS_MODIFICATIONS:
            raise ModificationLimitError("Only {} static c-terminus modification are allowed.".format(self.__class__.MAX_STATIC_C_TERMINUS_MODIFICATIONS))

    def __len__(self):
        return len(self.__modifications)

    @property
    def all(self):
        return self.__modifications

    @property
    def variable(self):
        return self.__variable_modifications

    @property
    def static(self):
        return self.__static_modifications

    @property
    def static_n_terminus_modifications(self):
        return self.__static_n_terminus_modifications

    @property
    def static_c_terminus_modifications(self):
        return self.__static_c_terminus_modifications

    @classmethod
    def read_from_csv_file(cls, csv_file_path: pathlib.Path):
        return cls(Modification.read_from_csv_file(csv_file_path))

        