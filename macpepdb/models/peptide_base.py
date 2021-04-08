from collections import Counter
from psycopg2.extras import execute_values

from ..proteomics.neutral_loss import H2O
from ..proteomics.amino_acid import AminoAcid
from ..proteomics.amino_acid import AMINO_ACIDS_FOR_COUNTING

# This class is only a super class acutal peptide classes e.g. Peptide
class PeptideBase:
    TABLE_NAME = 'peptide_base'
    FASTA_HEADER_PREFIX = ">PEPTIDE_"

    def __init__(self, sequence: str, number_of_missed_cleavages: int):
        self.__sequence = sequence.upper()
        self.__number_of_missed_cleavages = number_of_missed_cleavages
        self.__weight = self.__class__.calculate_weight(self.__sequence)
        # On demand values
        self.__amino_acid_counter = None

    @property
    def sequence(self):
        return self.__sequence

    @property
    def weight(self):
        return self.__weight

    @property
    def number_of_missed_cleavages(self):
        return self.__number_of_missed_cleavages

    @property
    def a_count(self):
        if not self.__amino_acid_counter:
            self.__count_amino_acids()
        return self.__amino_acid_counter['A']

    @property
    def b_count(self):
        if not self.__amino_acid_counter:
            self.__count_amino_acids()
        return self.__amino_acid_counter['B']

    @property
    def c_count(self):
        if not self.__amino_acid_counter:
            self.__count_amino_acids()
        return self.__amino_acid_counter['C']

    @property
    def d_count(self):
        if not self.__amino_acid_counter:
            self.__count_amino_acids()
        return self.__amino_acid_counter['D']

    @property
    def e_count(self):
        if not self.__amino_acid_counter:
            self.__count_amino_acids()
        return self.__amino_acid_counter['E']

    @property
    def f_count(self):
        if not self.__amino_acid_counter:
            self.__count_amino_acids()
        return self.__amino_acid_counter['F']

    @property
    def g_count(self):
        if not self.__amino_acid_counter:
            self.__count_amino_acids()
        return self.__amino_acid_counter['G']

    @property
    def h_count(self):
        if not self.__amino_acid_counter:
            self.__count_amino_acids()
        return self.__amino_acid_counter['H']

    @property
    def i_count(self):
        if not self.__amino_acid_counter:
            self.__count_amino_acids()
        return self.__amino_acid_counter['I']

    @property
    def j_count(self):
        if not self.__amino_acid_counter:
            self.__count_amino_acids()
        return self.__amino_acid_counter['J']

    @property
    def k_count(self):
        if not self.__amino_acid_counter:
            self.__count_amino_acids()
        return self.__amino_acid_counter['K']

    @property
    def l_count(self):
        if not self.__amino_acid_counter:
            self.__count_amino_acids()
        return self.__amino_acid_counter['L']

    @property
    def m_count(self):
        if not self.__amino_acid_counter:
            self.__count_amino_acids()
        return self.__amino_acid_counter['M']

    @property
    def n_count(self):
        if not self.__amino_acid_counter:
            self.__count_amino_acids()
        return self.__amino_acid_counter['N']

    @property
    def o_count(self):
        if not self.__amino_acid_counter:
            self.__count_amino_acids()
        return self.__amino_acid_counter['O']

    @property
    def p_count(self):
        if not self.__amino_acid_counter:
            self.__count_amino_acids()
        return self.__amino_acid_counter['P']

    @property
    def q_count(self):
        if not self.__amino_acid_counter:
            self.__count_amino_acids()
        return self.__amino_acid_counter['Q']

    @property
    def r_count(self):
        if not self.__amino_acid_counter:
            self.__count_amino_acids()
        return self.__amino_acid_counter['R']

    @property
    def s_count(self):
        if not self.__amino_acid_counter:
            self.__count_amino_acids()
        return self.__amino_acid_counter['S']

    @property
    def t_count(self):
        if not self.__amino_acid_counter:
            self.__count_amino_acids()
        return self.__amino_acid_counter['T']

    @property
    def u_count(self):
        if not self.__amino_acid_counter:
            self.__count_amino_acids()
        return self.__amino_acid_counter['U']

    @property
    def v_count(self):
        if not self.__amino_acid_counter:
            self.__count_amino_acids()
        return self.__amino_acid_counter['V']

    @property
    def w_count(self):
        if not self.__amino_acid_counter:
            self.__count_amino_acids()
        return self.__amino_acid_counter['W']

    @property
    def y_count(self):
        if not self.__amino_acid_counter:
            self.__count_amino_acids()
        return self.__amino_acid_counter['Y']

    @property
    def z_count(self):
        if not self.__amino_acid_counter:
            self.__count_amino_acids()
        return self.__amino_acid_counter['Z']

    @property
    def n_terminus(self):
        return self.sequence[0]

    @property
    def c_terminus(self):
        return self.sequence[-1]

    @property
    def length(self):
        return len(self.sequence)

    # Calculats the weight of a sequence
    @classmethod
    def calculate_weight(cls, sequence: str) -> int:
        weight = H2O.mono_mass
        for amino_acid_one_letter_code in sequence:
            weight += AminoAcid.get_by_one_letter_code(amino_acid_one_letter_code).mono_mass
        return weight

    def __count_amino_acids(self):
        self.__amino_acid_counter = Counter(self.__sequence)

    # This method is implemented to make sure only the sequence is used as hash when a protein is stored in a hashable collection (Set, Dictionary, ...)
    def __hash__(self):
        return hash((self.weight, self.sequence))

    # According to the Python documentation this should be implemented if __hash__() is implemented.
    def __eq__(self, other):
        return self.sequence == other.sequence

    def __get_fasta_header(self):
        return self.__class__.FASTA_HEADER_PREFIX

    def to_fasta_entry(self):
        return "{}\n{}\n".format(self.__get_fasta_header(), self.sequence)

    def to_dict(self):
        return {
            'sequence': self.sequence,
            'number_of_missed_cleavages': self.number_of_missed_cleavages,
            'weight': self.weight
        }

    @classmethod
    def from_dict(cls, attributes: dict):
        return cls(attributes["sequence"], attributes["number_of_missed_cleavages"])

    @classmethod
    def select(cls, database_cursor, select_conditions: tuple = ("", []), fetchall: bool = False):
        """
        @param database_cursor
        @param select_conditions A tupel with the where statement (without WHERE) and a list of parameters, e.g. ("accession = %s AND taxonomy_id = %s",["Q257X2", 6909])
        @param fetchall Indicates if multiple rows should be fetched
        @return Ppetide or list of proteins
        """
        select_query = f"SELECT sequence, number_of_missed_cleavages FROM {cls.TABLE_NAME}"
        if len(select_conditions) == 2 and len(select_conditions[0]):
            select_query += f" WHERE {select_conditions[0]}"
        select_query += ";"
        database_cursor.execute(select_query, select_conditions[1])
        
        if fetchall:
            return [cls(row[0], row[1]) for row in database_cursor.fetchall()]
        else:
            row = database_cursor.fetchall()
            if row:
                return cls(row[0], row[1])
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
            f"INSERT INTO {cls.TABLE_NAME} (weight, sequence, length, number_of_missed_cleavages, a_count, b_count, c_count, d_count, e_count, f_count, g_count, h_count, i_count, j_count, k_count, l_count, m_count, n_count, o_count, p_count, q_count, r_count, s_count, t_count, u_count, v_count, w_count, y_count, z_count, n_terminus, c_terminus) "
            "VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s);"
        )
        database_cursor.execute(
            INSERT_QUERY,
            (
                peptide.weight,
                peptide.sequence,
                peptide.length,
                peptide.number_of_missed_cleavages,
                peptide.a_count,
                peptide.b_count,
                peptide.c_count,
                peptide.d_count,
                peptide.e_count,
                peptide.f_count,
                peptide.g_count,
                peptide.h_count,
                peptide.i_count,
                peptide.j_count,
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
                peptide.z_count,
                peptide.n_terminus,
                peptide.c_terminus
            )
        )

    @classmethod
    def bulk_insert(cls, database_cursor, peptides: list):
        """
        @param database_cursor Database cursor with open transaction.
        @param peptides Peptides for bulk insert. Make sure they do not exitst before.
        @return list List of the new peptide IDs in the same order as the peptides where inserted.
        """
        BULK_INSERT_QUERY = (
            f"INSERT INTO {cls.TABLE_NAME} (weight, sequence, length, number_of_missed_cleavages, a_count, b_count, c_count, d_count, e_count, f_count, g_count, h_count, i_count, j_count, k_count, l_count, m_count, n_count, o_count, p_count, q_count, r_count, s_count, t_count, u_count, v_count, w_count, y_count, z_count, n_terminus, c_terminus) "
            "VALUES %s;"
        )
        # Bulk insert the new peptides
        execute_values(
            database_cursor,
            BULK_INSERT_QUERY,
            [
                (
                    peptide.weight,
                    peptide.sequence,
                    peptide.length,
                    peptide.number_of_missed_cleavages,
                    peptide.a_count,
                    peptide.b_count,
                    peptide.c_count,
                    peptide.d_count,
                    peptide.e_count,
                    peptide.f_count,
                    peptide.g_count,
                    peptide.h_count,
                    peptide.i_count,
                    peptide.j_count,
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
                    peptide.z_count,
                    peptide.n_terminus,
                    peptide.c_terminus
                ) for peptide in peptides
            ]
        )
