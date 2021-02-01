from sqlalchemy import Column, BigInteger, String, Integer, SmallInteger, CHAR
from sqlalchemy import orm

from ..proteomics.neutral_loss import NeutralLoss
from ..proteomics.amino_acid import AminoAcid
from ..proteomics.amino_acid import AMINO_ACIDS_FOR_COUNTING

# This class is only a super class acutal peptide classes e.g. Peptide
class PeptideBase(object):
    FASTA_HEADER_PREFIX = ">PEPTIDE_"


    id = Column(BigInteger, primary_key = True)
    # Attentions: Do not change the sequence if the peptide is stored in a hashabale collection, because sequence is used as hash-key,
    # therefore a change of the sequence will result in undetectability of the peptide within the collection.
    sequence = Column(String)
    length = Column(SmallInteger)
    number_of_missed_cleavages = Column(SmallInteger)
    weight = Column(Integer)
    a_count = Column(SmallInteger)
    c_count = Column(SmallInteger)
    d_count = Column(SmallInteger)
    e_count = Column(SmallInteger)
    f_count = Column(SmallInteger)
    g_count = Column(SmallInteger)
    h_count = Column(SmallInteger)
    i_count = Column(SmallInteger)
    k_count = Column(SmallInteger)
    l_count = Column(SmallInteger)
    m_count = Column(SmallInteger)
    n_count = Column(SmallInteger)
    o_count = Column(SmallInteger)
    p_count = Column(SmallInteger)
    q_count = Column(SmallInteger)
    r_count = Column(SmallInteger)
    s_count = Column(SmallInteger)
    t_count = Column(SmallInteger)
    u_count = Column(SmallInteger)
    v_count = Column(SmallInteger)
    w_count = Column(SmallInteger)
    y_count = Column(SmallInteger)
    n_terminus = Column(CHAR(1))
    c_terminus = Column(CHAR(1))


    def __init__(self, sequence: str, number_of_missed_cleavages: int):
        self.peff_notation_of_modifications = ""
        if len(sequence):
            self.set_sequence(sequence, number_of_missed_cleavages)

    @orm.reconstructor
    def init_on_load(self):
        self.peff_notation_of_modifications = ""

    # Sets new sequence and missed cleavages and sets new weight, length, amino acid counts and termini according to the new sequence
    def set_sequence(self, sequence: str, number_of_missed_cleavages: int):
        self.sequence = sequence.upper()
        self.length = len(self.sequence)
        self.number_of_missed_cleavages = number_of_missed_cleavages
        self.weight = self.__class__.calculate_weight(self.sequence)
        self.__count_amino_acids()
        self.n_terminus = self.sequence[0]
        self.c_terminus = self.sequence[len(self.sequence) - 1]

    # Calculats the weight of a sequence
    @classmethod
    def calculate_weight(cls, sequence: str) -> int:
        weight = NeutralLoss.get_by_name("H2O").mono_mass
        for amino_acid_one_letter_code in sequence:
            weight += AminoAcid.get_by_one_letter_code(amino_acid_one_letter_code).mono_mass
        return weight

    def __count_amino_acids(self):
        for one_letter_code in AMINO_ACIDS_FOR_COUNTING:
            # instead of 21 ifs we could use `exec("self.{}_count = {}".format(one_letter_code.lower(), self.sequence.count(one_letter_code)))`
            # but creating and executing code is slow
            if one_letter_code == 'A':
                self.a_count = self.sequence.count(one_letter_code)
            elif one_letter_code == 'C':
                self.c_count = self.sequence.count(one_letter_code)
            elif one_letter_code == 'D':
                self.d_count = self.sequence.count(one_letter_code)
            elif one_letter_code == 'E':
                self.e_count = self.sequence.count(one_letter_code)
            elif one_letter_code == 'F':
                self.f_count = self.sequence.count(one_letter_code)
            elif one_letter_code == 'G':
                self.g_count = self.sequence.count(one_letter_code)
            elif one_letter_code == 'H':
                self.h_count = self.sequence.count(one_letter_code)
            elif one_letter_code == 'I':
                self.i_count = self.sequence.count(one_letter_code)
            elif one_letter_code == 'K':
                self.k_count = self.sequence.count(one_letter_code)
            elif one_letter_code == 'L':
                self.l_count = self.sequence.count(one_letter_code)
            elif one_letter_code == 'M':
                self.m_count = self.sequence.count(one_letter_code)
            elif one_letter_code == 'N':
                self.n_count = self.sequence.count(one_letter_code)
            elif one_letter_code == 'O':
                self.o_count = self.sequence.count(one_letter_code)
            elif one_letter_code == 'P':
                self.p_count = self.sequence.count(one_letter_code)
            elif one_letter_code == 'Q':
                self.q_count = self.sequence.count(one_letter_code)
            elif one_letter_code == 'R':
                self.r_count = self.sequence.count(one_letter_code)
            elif one_letter_code == 'S':
                self.s_count = self.sequence.count(one_letter_code)
            elif one_letter_code == 'T':
                self.t_count = self.sequence.count(one_letter_code)
            elif one_letter_code == 'U':
                self.u_count = self.sequence.count(one_letter_code)
            elif one_letter_code == 'V':
                self.v_count = self.sequence.count(one_letter_code)
            elif one_letter_code == 'W':
                self.w_count = self.sequence.count(one_letter_code)
            elif one_letter_code == 'Y':
                self.y_count = self.sequence.count(one_letter_code)

    # This method is implemented to make sure only the sequence is used as hash when a protein is stored in a hashable collection (Set, Dictionary, ...)
    def __hash__(self):
        return hash(self.sequence)

    # According to the Python documentation this should be implemented if __hash__() is implemented.
    def __eq__(self, other):
        return self.sequence == other.sequence

    def __get_fasta_header(self):
        if len(self.peff_notation_of_modifications) > 0:
            return "{} {}".format(self.__class__.FASTA_HEADER_PREFIX, self.peff_notation_of_modifications)
        else:
            return self.__class__.FASTA_HEADER_PREFIX

    def to_fasta_entry(self):
        return "{}\n{}\n".format(self.__get_fasta_header(), self.sequence)

    def insert_values(self) -> dict:
        return {
            "sequence": self.sequence,
            "length": self.length,
            "number_of_missed_cleavages": self.number_of_missed_cleavages,
            "weight": self.weight,
            "a_count": self.a_count,
            "c_count": self.c_count,
            "d_count": self.d_count,
            "e_count": self.e_count,
            "f_count": self.f_count,
            "g_count": self.g_count,
            "h_count": self.h_count,
            "i_count": self.i_count,
            "k_count": self.k_count,
            "l_count": self.l_count,
            "m_count": self.m_count,
            "n_count": self.n_count,
            "o_count": self.o_count,
            "p_count": self.p_count,
            "q_count": self.q_count,
            "r_count": self.r_count,
            "s_count": self.s_count,
            "t_count": self.t_count,
            "u_count": self.u_count,
            "v_count": self.v_count,
            "w_count": self.w_count,
            "y_count": self.y_count,
            "n_terminus": self.n_terminus,
            "c_terminus": self.c_terminus 
        }

    def to_dict(self):
        dictionary = self.insert_values()
        dictionary["id"] = self.id
        dictionary["peff_notation_of_modifications"] = self.peff_notation_of_modifications
        return dictionary

    @classmethod
    def from_dict(cls, attributes: dict):
        peptide = cls("", 0)
        peptide.id = attributes["id"]
        peptide.sequence = attributes["sequence"]
        peptide.length = attributes["length"]
        peptide.number_of_missed_cleavages = attributes["number_of_missed_cleavages"]
        peptide.weight = attributes["weight"]
        peptide.a_count = attributes["a_count"]
        peptide.c_count = attributes["c_count"]
        peptide.d_count = attributes["d_count"]
        peptide.e_count = attributes["e_count"]
        peptide.f_count = attributes["f_count"]
        peptide.g_count = attributes["g_count"]
        peptide.h_count = attributes["h_count"]
        peptide.i_count = attributes["i_count"]
        peptide.k_count = attributes["k_count"]
        peptide.l_count = attributes["l_count"]
        peptide.m_count = attributes["m_count"]
        peptide.n_count = attributes["n_count"]
        peptide.o_count = attributes["o_count"]
        peptide.p_count = attributes["p_count"]
        peptide.q_count = attributes["q_count"]
        peptide.r_count = attributes["r_count"]
        peptide.s_count = attributes["s_count"]
        peptide.t_count = attributes["t_count"]
        peptide.u_count = attributes["u_count"]
        peptide.v_count = attributes["v_count"]
        peptide.w_count = attributes["w_count"]
        peptide.y_count = attributes["y_count"]
        peptide.n_terminus = attributes["n_terminus"]
        peptide.c_terminus = attributes["c_terminus"]
        peptide.peff_notation_of_modifications = attributes["peff_notation_of_modifications"]
        return peptide
