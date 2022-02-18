# std imports
from __future__ import annotations
from enum import IntEnum, unique
import csv
import pathlib
from typing import List, ClassVar

# internal imports
from macpepdb.proteomics.mass.convert import to_int as mass_to_int, to_float as mass_to_float
from macpepdb.proteomics.amino_acid import AminoAcid

@unique
class ModificationPosition(IntEnum):
    """
    Enum of possible modification positions.
    """
    ANYWHERE = 1
    N_TERMINUS = 2
    C_TERMINUS = 3

    def __str__(self) -> str:
        return self.name.lower() # pylint: disable=no-member

    @classmethod
    def from_string(cls, position: str) -> ModificationPosition:
        """
        Returns the position

        Parameters
        ----------
        position : str
            Name of the position (anywhere, n_terminus, c_terminus)

        Returns
        -------
        Modification position
        """
        return cls[position.upper()]

class Modification:
    """
    Defines a posttranslational amino acid modification (PTM)

    Parameters
    ----------
    accession : str
        Accession, e.g. UNIMOD:1
    name : str
        Human readable name
    amino_acid : AminoAcid
        Amino acid to modified
    delta : int
        Mass change
    is_static : bool
        Indicates if this is static or variable
    position : ModificationPosition
        Possible position of the modification
    """

    # String to is static map
    IS_STATIC_STRING_LOOKUP: ClassVar[dict] = {
        "static": True,
        "variable": False
    }
    """Lookup is static status to bool
    """

    # Is static to string map
    IS_STATIC_BOOL_LOOKUP: ClassVar[dict] = {
        True: "static",
        False: "variable"
    }
    """ Lookup is static status to string
    """

    PEFF_KEY_UNIMOD: ClassVar[str] = "ModResUnimod"
    """PEFF key for Unimod
    """
    PEFF_KEY_PSI: ClassVar[str] = "ModResPsi"
    """Peff key fopr PSI
    """

    PEFF_KEY_OTHER: ClassVar[str] = "ModRes"
    """PEFF key for others
    """

    def __init__(self, accession: str, name: str, amino_acid: AminoAcid, delta: int, is_static: bool, position: ModificationPosition):
        self.__accession = accession
        self.__name = name
        self.__amino_acid = amino_acid
        self.__delta = delta
        self.__is_static = is_static
        self.__position = position

    @property
    def accession(self) -> str:
        """
        Returns
        -------
        Modification accession
        """
        return self.__accession

    @property
    def name(self) -> str:
        """
        Returns
        -------
        Human readable modification name
        """
        return self.__name

    @property
    def amino_acid(self) -> str:
        """
        Returns
        -------
        Targeted amino acid
        """
        return self.__amino_acid

    @property
    def delta(self) -> int:
        """
        Returns
        -------
        Mass change applied by this modification
        """
        return self.__delta

    @property
    def is_static(self) -> bool:
        """
        Returns
        -------
        True if modification is static
        """
        return self.__is_static

    @property
    def is_variable(self)  -> bool:
        """
        Returns
        -------
        True if modification is variable
        """
        return not self.__is_static

    @property
    def position(self) -> ModificationPosition:
        """
        Returns
        -------
        Position of the modification
        """
        return self.__position

    # Returns the modifications delta mass + the amino acid mono mass
    @property
    def mono_mass(self) -> int:
        """
        Returns
        -------
        Mono mass of the amino acid including the mass change abblied by the modification.
        """
        return self.__amino_acid.mono_mass + self.__delta

    @property
    def peff_key(self) -> str:
        """
        Returns
        -------
        The modifications PEFF-key
        """
        if self.__accession.upper().startswith("UNIMOD:"):
            return self.__class__.PEFF_KEY_UNIMOD
        elif self.__accession.upper().startswith("MOD:"):
            return self.__class__.PEFF_KEY_PSI
        else:
            return self.__class__.PEFF_KEY_OTHER

    @classmethod
    def read_from_csv_file(cls, csv_file_path: pathlib.Path) -> List[Modification]:
        """
        Reads modifications from CSV file.

        Parameters
        ----------
        csv_file_path : pathlib.Path
            Path of the CSV-file

        Returns
        -------
        List of modifications
        """
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

    def to_comet_parameter(self) -> str:
        """
        Returns the modifications as comet parameter. If the the modification is static, the result is ready to use.
        If the modification is variable or terminus modification, the resulting string contains a placeholder ('|i') for the variable modification index and another one ('|v') for the variable modification maximum.
        Remember Comet allows only 9 variable modifications in total.

        Returns
        -------
        Mostly ready to use Comet configuration string of the modification
        """
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
    def is_position_anywhere(self) -> bool:
        """
        Returns
        -------
        True if the modification position is "anywhere"
        """
        return self.__position == ModificationPosition.ANYWHERE
    
    @property
    def is_position_n_terminus(self) -> bool:
        """
        Returns
        -------
        True if the modification position is "n_terminus"
        """
        return self.__position == ModificationPosition.N_TERMINUS

    @property
    def is_position_c_terminus(self) -> bool:
        """
        Returns
        -------
        True if the modification position is "c_terminus"
        """
        return self.__position == ModificationPosition.C_TERMINUS

    @property
    def is_terminus_modification(self) -> bool:
        """
        Returns
        -------
        True if the modification position is at a terminus
        """
        return self.is_position_n_terminus or self.is_position_c_terminus

    def __str__(self) -> str:
        return "accession:  {}\nname:       {}\namino_acid: {}\ndelta:      {}\nstatic?:    {}\nposition:   {}".format(
            self.__accession,
            self.__name,
            self.__amino_acid.one_letter_code,
            self.__delta,
            "True" if self.__is_static else "False",
            str(self.__position)
        )

    @property
    def is_none_modification(self) -> bool:
        """
        Depracted: Will be removed in a feature release.

        Returns
        -------
        False
        """
        return False

    @classmethod
    def string_to_is_static(cls, is_static: str) -> bool:
        """
        'Converts' string to boolean

        Paramters
        ---------
        is_static : str
            String which indicates if the modification is static or variable (possible values: 'static', 'variable')

        Returns
        -------
        True if is_static is 'static', False if is_static is 'variable'

        Raises
        ------
        AttributeError
            If is_static does not match 'static' or 'variable'
        """
        try:
            return cls.IS_STATIC_STRING_LOOKUP[is_static]
        except KeyError:
            raise AttributeError("Is static {} is unknown.".format(is_static))

    @classmethod
    def is_static_to_string(cls, is_static: bool) -> str:
        """
        'Converts' is bool to 'static' or 'variable'

        Parameters
        ----------
        is_static : bool

        Returns
        -------
        'static' if is_static was True or 'variable' if is_static was False

        Raises
        ------
        AttributeError
            If is_static does not match Ture or False
        """
        try:
            return cls.IS_STATIC_BOOL_LOOKUP[is_static]
        except KeyError:
            raise AttributeError("Is static is unknown.")

    def is_static_as_string(self) -> str:
        """
        See class method `Modification.string_to_is_static()`
        """
        return self.__class__.is_static_to_string(self.__is_static)

    @classmethod
    def from_dict(cls, attributes: dict) -> Modification:
        """
        Creates modification from dictionary.

        Parameters
        ----------
        attributes : dict
            Dictionary of attributes with key: 'accession', 'name', 'amino_acid', 'delta', 'is_static' and 'position'

        Returns
        -------
        Modification
        """
        return Modification(
            attributes["accession"],
            attributes["name"],
            AminoAcid.get_by_one_letter_code(attributes["amino_acid"]),
            attributes["delta"],
            attributes["is_static"],
            ModificationPosition.from_string(attributes["position"])
        )

