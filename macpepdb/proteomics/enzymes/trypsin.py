# internal imports
from macpepdb.proteomics.enzymes.digest_enzyme import DigestEnzyme

class Trypsin(DigestEnzyme):
    NAME = "Trypsin"
    SHORTCUT = "try"
    REGEX = r"(?<=[KR])(?!P)"


    def __init__(self, max_number_of_missed_cleavages = 0, minimum_peptide_length = 0, maximum_peptide_length = 1):
        super().__init__(
            self.NAME,
            self.SHORTCUT,
            self.REGEX,
            max_number_of_missed_cleavages,
            minimum_peptide_length,
            maximum_peptide_length
        )