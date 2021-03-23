from .mass.convert import to_int as mass_to_int

AMINO_ACIDS_FOR_COUNTING = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "Y"]

class AminoAcid:
    def __init__(self, name: str, one_letter_code: str, three_letter_code: str, chemical_formula: str, mono_mass: float, average_mass: float):
        self.name = name
        self.one_letter_code = one_letter_code
        self.three_letter_code = three_letter_code
        self.chemical_formula = chemical_formula
        self.mono_mass = mass_to_int(mono_mass)
        self.average_mass = mass_to_int(average_mass)

    # Returns an amino acid by one letter code.
    # If the given one letter code is unknown Unknown Amino Acid (X) is returned.
    @classmethod
    def get_by_one_letter_code(cls, one_letter_code: str):
        try:
            return eval(one_letter_code.upper())
        except NameError:
            return X

    # Returns Tryptophan (W)
    @classmethod
    def get_haviest(cls):
        return W

    # Returns Glycine (G)
    @classmethod
    def get_lightest(cls):
        return G

    @staticmethod
    def get_unknown():
        return X

    @staticmethod
    def all() -> tuple:
        """
        Returns a tuple with all known amino acids
        """
        return KNOWN_AMINO_ACIDS


# Standard amino acids
# https://proteomicsresource.washington.edu/protocols06/masses.php
A = AminoAcid("Alanine", 'A', "Ala", "C3H5ON", 71.037113805, 71.0788)
C = AminoAcid("Cysteine", 'C', "Cys", "C3H5ONS", 103.009184505, 103.1388)
D = AminoAcid("Aspartic acid", 'D', "Asp", "C4H5O3N", 115.026943065, 115.0886)
E = AminoAcid("Glutamic acid", 'E', "Glu", "C5H7O3N", 129.042593135, 129.1155)
F = AminoAcid("Phenylalanine", 'F', "Phe", "C9H9ON", 147.068413945, 147.1766)
G = AminoAcid("Glycine", 'G', "Gly", "C2H3ON", 57.021463735, 57.0519)
H = AminoAcid("Histidine", 'H', "His", "C6H7ON3", 137.058911875, 137.1411)
I = AminoAcid("Isoleucine", 'I', "Ile", "C6H11ON", 113.084064015, 113.1594)
K = AminoAcid("Lysine", 'K', "Lys", "C6H12ON2", 128.094963050, 128.1741)
L = AminoAcid("Leucine", 'L', "Leu", "C6H11ON", 113.084064015, 113.1594)
M = AminoAcid("Methionine", 'M', "Met", "C5H9ONS", 131.040484645, 131.1926)
N = AminoAcid("Asparagine", 'N', "Asn", "C4H6O2N2", 114.042927470, 114.1038)
O = AminoAcid("Pyrrolysine", 'O', "Pyl", "C12H19N3O2", 237.147726925, 237.29816)
P = AminoAcid("Proline", 'P', "Pro", "C5H7ON", 97.052763875, 97.1167)
Q = AminoAcid("Glutamine", 'Q', "Gln", "C5H8O2N2", 128.05857754, 128.1307)
R = AminoAcid("Arginine", 'R', "Arg", "C6H12ON4", 156.101111050, 156.1875)
S = AminoAcid("Serine", 'S', "Ser", "C3H5O2N", 87.032028435, 87.0782)
T = AminoAcid("Threonine", 'T', "Thr", "C4H7O2N", 101.047678505, 101.1051)
U = AminoAcid("Selenocysteine", 'U', "SeC", "C3H5NOSe", 150.953633405, 150.0379)
V = AminoAcid("Valine", 'V', "Val", "C5H9ON", 99.068413945, 99.1326)
W = AminoAcid("Tryptophan", 'W', "Trp", "C11H10ON2", 186.079312980, 186.2132)
Y = AminoAcid("Tyrosine", 'Y', "Tyr", "C9H9O2N", 163.063328575, 163.1760)
# Special amino acids
## Some Search Engines and Databases used the X Amino Acid for unknown amino acids
X = AminoAcid("Unknown Amino Acid", 'X', "Xaa", "Unknown", 0.0, 0.0)

# Tuple with containing all standard amino acids and X (unknown)
KNOWN_AMINO_ACIDS = (
    A,
    C,
    D,
    E,
    F,
    G,
    H,
    I,
    K,
    L,
    M,
    N,
    O,
    P,
    Q,
    R,
    S,
    T,
    U,
    V,
    W,
    Y,
    X
)