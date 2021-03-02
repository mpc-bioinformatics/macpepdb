from psycopg2.extras import execute_values

from ..proteomics.neutral_loss import NeutralLoss
from ..proteomics.amino_acid import AminoAcid
from ..proteomics.amino_acid import AMINO_ACIDS_FOR_COUNTING

# This class is only a super class acutal peptide classes e.g. Peptide
class PeptideBase:
    TABLE_NAME = 'peptide_base'
    FASTA_HEADER_PREFIX = ">PEPTIDE_"

    def __init__(self, sequence: str, number_of_missed_cleavages: int, id = None):
        self.__id = None
        self.peff_notation_of_modifications = ""
        if len(sequence):
            self.set_sequence(sequence, number_of_missed_cleavages)

    @property
    def id(self):
        return self.__id;

    # Sets new sequence and missed cleavages and sets new weight, length, amino acid counts and termini according to the new sequence
    def set_sequence(self, sequence: str, number_of_missed_cleavages: int):
        self.sequence = sequence.upper()
        self.length = len(self.sequence)
        self.number_of_missed_cleavages = number_of_missed_cleavages
        self.weight = self.__class__.calculate_weight(self.sequence)
        self.a_count = 0
        self.c_count = 0
        self.d_count = 0
        self.e_count = 0
        self.f_count = 0
        self.g_count = 0
        self.h_count = 0
        self.i_count = 0
        self.k_count = 0
        self.l_count = 0
        self.m_count = 0
        self.n_count = 0
        self.o_count = 0
        self.p_count = 0
        self.q_count = 0
        self.r_count = 0
        self.s_count = 0
        self.t_count = 0
        self.u_count = 0
        self.v_count = 0
        self.w_count = 0
        self.y_count = 0
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

    @classmethod
    def select(cls, database_cursor, select_conditions: tuple = ("", []), fetchall: bool = False):
        """
        @param database_cursor
        @param select_conditions A tupel with the where statement (without WHERE) and a list of parameters, e.g. ("accession = %s AND taxonomy_id = %s",["Q257X2", 6909])
        @param fetchall Indicates if multiple rows should be fetched
        @return Ppetide or list of proteins
        """
        select_query = f"SELECT id, sequence, number_of_missed_cleavages FROM {cls.TABLE_NAME}"
        if len(select_conditions) == 2 and len(select_conditions[0]):
            select_query += f" WHERE {select_conditions[0]}"
        select_query += ";"
        database_cursor.execute(select_query, select_conditions[1])
        
        if fetchall:
            return [cls(row[1], row[2], row[0]) for row in database_cursor.fetchall()]
        else:
            row = database_cursor.fetchall()
            if row:
                return cls(row[1], row[2], row[0])
            else:
                return None

    @classmethod
    def insert(cls, database_cursor, peptide) -> int:
        """
        @param database_cursor Database cursor with open transaction.
        @param peptide Peptide to insert
        @return int ID of inserted peptide
        """
        INSERT_QUERY = (
            f"INSERT INTO {cls.TABLE_NAME} (sequence, length, number_of_missed_cleavages, weight, a_count, c_count, d_count, e_count, f_count, g_count, h_count, i_count, k_count, l_count, m_count, n_count, o_count, p_count, q_count, r_count, s_count, t_count, u_count, v_count, w_count, y_count, n_terminus, c_terminus) "
            "VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s) "
            "RETURNING id;"
        )
        database_cursor.execute(
            INSERT_QUERY,
            (
                peptide.sequence,
                peptide.length,
                peptide.number_of_missed_cleavages,
                peptide.weight,
                peptide.a_count,
                peptide.c_count,
                peptide.d_count,
                peptide.e_count,
                peptide.f_count,
                peptide.g_count,
                peptide.h_count,
                peptide.i_count,
                peptide.k_count,
                peptide.l_count,
                peptide.m_count,
                peptide.n_count,
                peptide.o_count,
                peptide.p_count,
                peptide.q_count,
                peptide.r_count,
                peptide.s_count,
                peptide.t_count,
                peptide.u_count,
                peptide.v_count,
                peptide.w_count,
                peptide.y_count,
                peptide.n_terminus,
                peptide.c_terminus
            )
        )

    @classmethod
    def bulk_insert(cls, database_cursor, peptides: list) -> list:
        """
        @param database_cursor Database cursor with open transaction.
        @param peptides Peptides for bulk insert. Make sure they do not exitst before.
        @return list List of the new peptide IDs in the same order as the peptides where inserted.
        """
        BULK_INSERT_QUERY = (
            f"INSERT INTO {cls.TABLE_NAME} (sequence, length, number_of_missed_cleavages, weight, a_count, c_count, d_count, e_count, f_count, g_count, h_count, i_count, k_count, l_count, m_count, n_count, o_count, p_count, q_count, r_count, s_count, t_count, u_count, v_count, w_count, y_count, n_terminus, c_terminus) "
            "VALUES %s "
            "RETURNING id;"
        )
        # Bulk insert the new peptides
        id_rows = execute_values(
            database_cursor,
            BULK_INSERT_QUERY,
            [
                (
                    peptide.sequence,
                    peptide.length,
                    peptide.number_of_missed_cleavages,
                    peptide.weight,
                    peptide.a_count,
                    peptide.c_count,
                    peptide.d_count,
                    peptide.e_count,
                    peptide.f_count,
                    peptide.g_count,
                    peptide.h_count,
                    peptide.i_count,
                    peptide.k_count,
                    peptide.l_count,
                    peptide.m_count,
                    peptide.n_count,
                    peptide.o_count,
                    peptide.p_count,
                    peptide.q_count,
                    peptide.r_count,
                    peptide.s_count,
                    peptide.t_count,
                    peptide.u_count,
                    peptide.v_count,
                    peptide.w_count,
                    peptide.y_count,
                    peptide.n_terminus,
                    peptide.c_terminus
                ) for peptide in peptides
            ],
            fetch=True
        )
        return [row[0] for row in id_rows]
