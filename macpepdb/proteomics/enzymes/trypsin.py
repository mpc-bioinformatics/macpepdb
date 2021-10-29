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

    REGEX = r"(?<=[KR])(?!P)"
    """Regex for finding cleavage positions
    """

    def __init__(self, max_number_of_missed_cleavages = 0, minimum_peptide_length = 0, maximum_peptide_length = 1):
        super().__init__(
            self.NAME,
            self.SHORTCUT,
            self.REGEX,
            max_number_of_missed_cleavages,
            minimum_peptide_length,
            maximum_peptide_length
        )