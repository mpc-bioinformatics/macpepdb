import pathlib
import csv
from enum import IntEnum, unique

from .mass.convert import to_int as mass_to_int, to_float as mass_to_float
from .amino_acid import AminoAcid

@unique
class ModificationPosition(IntEnum):
    ANYWHERE = 1
    N_TERMINUS = 2
    C_TERMINUS = 3

    def __str__(self):
        return self.name.lower() # pylint: disable=no-member

    @classmethod
    def from_string(cls, position: str):
        return cls[position.upper()]

class Modification:
    # String to is static map
    IS_STATIC_STRING_LOOKUP = {
        "static": True,
        "variable": False
    }

    # Is static to string map
    IS_STATIC_BOOL_LOOKUP = {
        True: "static",
        False: "variable"
    }

    # Peff notation keys
    PEFF_KEY_UNIMOD = "ModResUnimod"
    PEFF_KEY_PSI = "ModResPsi"
    PEFF_KEY_OTHER = "ModRes"

    def __init__(self, accession: str, name: str, amino_acid: AminoAcid, delta: int, is_static: bool, position: ModificationPosition):
        self.__accession = accession
        self.__name = name
        self.__amino_acid = amino_acid
        self.__delta = delta
        self.__is_static = is_static
        self.__position = position

    @property
    def accession(self):
        return self.__accession

    @property
    def name(self):
        return self.__name

    @property
    def amino_acid(self):
        return self.__amino_acid

    @property
    def delta(self):
        return self.__delta

    @property
    def is_static(self):
        return self.__is_static

    @property
    def is_variable(self):
        return not self.__is_static

    @property
    def position(self):
        return self.__position

    # Returns the modifications delta mass + the amino acid mono mass
    @property
    def mono_mass(self):
        return self.__amino_acid.mono_mass + self.__delta

    @property
    def peff_key(self):
        if self.__accession.upper().startswith("UNIMOD:"):
            return self.__class__.PEFF_KEY_UNIMOD
        elif self.__accession.upper().startswith("MOD:"):
            return self.__class__.PEFF_KEY_PSI
        else:
            return self.__class__.PEFF_KEY_OTHER

    @classmethod
    def read_from_csv_file(cls, csv_file_path: pathlib.Path):
        modifications = []
        with csv_file_path.open("r") as csv_file:
            csv_reader = csv.reader(csv_file)
            # Omit header
            next(csv_reader)
            for row in csv_reader:
                modifications.append(
                    Modification(row[0], row[1], AminoAcid.get_by_one_letter_code(row[2]), mass_to_int(float(row[3])), cls.string_to_is_static(row[4].lower()), ModificationPosition.from_string(row[5]))
                )
        return modifications

    # Returns the modifications as comet parameter. If the the modification is static, the result is ready to use.
    # If the modification is variable or terminus modification, the resulting string contains a placeholder ('|i') for the variable modification index and another one ('|v') for the variable modification maximum.
    # Remember Comet allows only 9 variable modifications in total.
    def to_comet_parameter(self):
        if self.is_static and not self.is_terminus_modification:
            return "add_{}_{} = {}".format(self.amino_acid.one_letter_code, self.amino_acid.name.lower().replace(" ", "_"), mass_to_float(self.delta))
        elif self.is_static and self.is_position_n_terminus:
            return "add_Nterm_peptide = {}".format(mass_to_float(self.delta))
        elif self.is_static and self.is_position_c_terminus:
            return "add_Cterm_peptide = {}".format(mass_to_float(self.delta))
        elif self.is_variable and not self.is_terminus_modification:
            return "variable_mod0|i = {} {} 0 |v -1 2 0 0.0".format(
                mass_to_float(self.delta),
                self.amino_acid.one_letter_code
            )
        elif self.is_variable and self.is_position_n_terminus:
            return "variable_mod0|i = {} {} 0 |v 0 2 0 0.0".format(
                mass_to_float(self.delta),
                self.amino_acid.one_letter_code
            )
        elif self.is_variable and self.is_position_c_terminus:
            return "variable_mod0|i = {} {} 0 |v 0 3 0 0.0".format(
                mass_to_float(self.delta),
                self.amino_acid.one_letter_code
            )
                
    @property
    def is_position_anywhere(self):
        return self.__position == ModificationPosition.ANYWHERE
    
    @property
    def is_position_n_terminus(self):
        return self.__position == ModificationPosition.N_TERMINUS

    @property
    def is_position_c_terminus(self):
        return self.__position == ModificationPosition.C_TERMINUS

    @property
    def is_terminus_modification(self):
        return self.is_position_n_terminus or self.is_position_c_terminus

    def __str__(self):
        return "accession:  {}\nname:       {}\namino_acid: {}\ndelta:      {}\nstatic?:    {}\nposition:   {}".format(
            self.__accession,
            self.__name,
            self.__amino_acid.one_letter_code,
            self.__delta,
            "True" if self.__is_static else "False",
            str(self.__position)
        )

    @property
    def is_none_modification(self):
        return False

    @classmethod
    def string_to_is_static(cls, is_static: str):
        try:
            return cls.IS_STATIC_STRING_LOOKUP[is_static]
        except KeyError:
            raise AttributeError("Is static {} is unknown.".format(is_static))

    @classmethod
    def is_static_to_string(cls, is_static: bool):
        try:
            return cls.IS_STATIC_BOOL_LOOKUP[is_static]
        except KeyError:
            raise AttributeError("Is static is unknown.")

    def is_static_as_string(self):
        return self.__class__.is_static_to_string(self.__is_static)

    @classmethod
    def from_dict(cls, attributes: dict):
        return Modification(
            attributes["accession"],
            attributes["name"],
            AminoAcid.get_by_one_letter_code(attributes["amino_acid"]),
            attributes["delta"],
            attributes["is_static"],
            ModificationPosition.from_string(attributes["position"])
        )

