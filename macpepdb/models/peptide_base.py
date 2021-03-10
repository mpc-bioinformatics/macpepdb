from psycopg2.extras import execute_values

from ..proteomics.neutral_loss import NeutralLoss
from ..proteomics.amino_acid import AminoAcid
from ..proteomics.amino_acid import AMINO_ACIDS_FOR_COUNTING

# This class is only a super class acutal peptide classes e.g. Peptide
class PeptideBase:
    TABLE_NAME = 'peptide_base'
    FASTA_HEADER_PREFIX = ">PEPTIDE_"

    def __init__(self, sequence: str, number_of_missed_cleavages: int, id = None):
        self.__id = id
        self.__sequence = sequence.upper()
        self.__number_of_missed_cleavages = number_of_missed_cleavages
        self.__weight = None
        
    @property
    def id(self):
        return self.__id;

    @property
    def sequence(self):
        return self.__sequence

    @property
    def weight(self):
        if not self.__weight:
            self.__weight = self.__class__.calculate_weight(self.__sequence)
        return self.__weight

    @property
    def number_of_missed_cleavages(self):
        return self.__number_of_missed_cleavages

    @property
    def a_count(self):
        return self.sequence.count('A')

    @property
    def c_count(self):
        return self.sequence.count('C')

    @property
    def d_count(self):
        return self.sequence.count('D')

    @property
    def e_count(self):
        return self.sequence.count('E')

    @property
    def f_count(self):
        return self.sequence.count('F')

    @property
    def g_count(self):
        return self.sequence.count('G')

    @property
    def h_count(self):
        return self.sequence.count('H')

    @property
    def i_count(self):
        return self.sequence.count('I')

    @property
    def k_count(self):
        return self.sequence.count('K')

    @property
    def l_count(self):
        return self.sequence.count('L')

    @property
    def m_count(self):
        return self.sequence.count('M')

    @property
    def n_count(self):
        return self.sequence.count('N')

    @property
    def o_count(self):
        return self.sequence.count('O')

    @property
    def p_count(self):
        return self.sequence.count('P')

    @property
    def q_count(self):
        return self.sequence.count('Q')

    @property
    def r_count(self):
        return self.sequence.count('R')

    @property
    def s_count(self):
        return self.sequence.count('S')

    @property
    def t_count(self):
        return self.sequence.count('T')

    @property
    def u_count(self):
        return self.sequence.count('U')

    @property
    def v_count(self):
        return self.sequence.count('V')

    @property
    def w_count(self):
        return self.sequence.count('W')

    @property
    def y_count(self):
        return self.sequence.count('Y')

    @property
    def n_terminus(self):
        return self.sequence[0]

    @property
    def c_terminus(self):
        return self.sequence[-1]

    @property
    def length(self):
        return len(self.sequence)

    @property
    def partition_index(self):
        return self.__class__.get_parition_index(self.weight)

    # Calculats the weight of a sequence
    @classmethod
    def calculate_weight(cls, sequence: str) -> int:
        weight = NeutralLoss.get_by_name("H2O").mono_mass
        for amino_acid_one_letter_code in sequence:
            weight += AminoAcid.get_by_one_letter_code(amino_acid_one_letter_code).mono_mass
        return weight

    # This method is implemented to make sure only the sequence is used as hash when a protein is stored in a hashable collection (Set, Dictionary, ...)
    def __hash__(self):
        return hash(self.sequence)

    # According to the Python documentation this should be implemented if __hash__() is implemented.
    def __eq__(self, other):
        return self.sequence == other.sequence

    def __get_fasta_header(self):
        return self.__class__.FASTA_HEADER_PREFIX

    def to_fasta_entry(self):
        return "{}\n{}\n".format(self.__get_fasta_header(), self.sequence)

    def to_dict(self):
        return {
            'id': self.id,
            'sequence': self.sequence,
            'number_of_missed_cleavages': attributes["number_of_missed_cleavages"],
            'weight': self.weight
        }

    @classmethod
    def from_dict(cls, attributes: dict):
        return cls(attributes["sequence"], attributes["number_of_missed_cleavages"], attributes["id"])

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
