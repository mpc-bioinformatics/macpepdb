# std import
import re
from typing import ClassVar

# internal imports
from macpepdb.proteomics.enzymes.digest_enzyme import DigestEnzyme

class Trypsin(DigestEnzyme):
    """
    Defines Trypsin which cuts amino acid sequences on each K or R if not followed by P.

    Parameters
    ----------
    max_number_of_missed_cleavages : int
        Maximum number of missed cleavages
    minimum_peptide_length : int
        Minimum peptide length
    maximum_peptide_length : int
        Maxiumum peptide length
    """

    NAME = "Trypsin"
    """Enzym name
    """

    SHORTCUT = "try"
    """Enzym shortcut
    """

    CLEAVAGE_REGEX: ClassVar[re.Pattern] = re.compile(r"(?<=[KR])(?!P)")
    """Regex for finding cleavage positions
    """

    MISSED_CLEAVAGE_REGEX: ClassVar[re.Pattern] = re.compile(r"(R|K)(?!($|P))")
    """Regex to count missed cleavages: R or L not followed by P or end of string
    """

    def __init__(self, max_number_of_missed_cleavages = 0, minimum_peptide_length = 0, maximum_peptide_length = 1):
        super().__init__(
            self.NAME,
            self.SHORTCUT,
            self.CLEAVAGE_REGEX,
            max_number_of_missed_cleavages,
            minimum_peptide_length,
            maximum_peptide_length
        )