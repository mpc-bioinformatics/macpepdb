# std imports
from __future__ import annotations
from collections import Counter
from typing import Iterator, Optional, Union, List

# external imports
from psycopg2.extras import execute_values

# internal imports
from macpepdb.database.query_helpers.where_condition import WhereCondition
from macpepdb.proteomics.neutral_loss import H2O
from macpepdb.proteomics.amino_acid import AminoAcid

# This class is only a super class acutal peptide classes e.g. Peptide
class PeptideBase:
    """
    Base class for peptides. Keeps track of the sequence, mass, partition, amino acid counts of the peptide.

    Parameters
    ----------
    sequence : str
        Amino acid sequence
    number_of_missed_cleavages : int
        Number of missed cleavages
    """

    TABLE_NAME = 'peptide_base'
    """Name of the database table
    """
    
    FASTA_HEADER_PREFIX = ">PEPTIDE_"
    """Prefix for FASTA entries
    """

    PARTITONS = [
        853375237186,
        913477000136,
        955581465254,
        990535912010,
        1023509756866,
        1053679782375,
        1083451590056,
        1111525800894,
        1138631937691,
        1165675382443,
        1192670339929,
        1219667320195,
        1246666985860,
        1273670022365,
        1300685629579,
        1327695660223,
        1354681590068,
        1381654766255,
        1408643015600,
        1435684426855,
        1462744978267,
        1489818752295,
        1517750791935,
        1544867476661,
        1572829376619,
        1600811024379,
        1628809309848,
        1656831879069,
        1684878435320,
        1712960725178,
        1741885980286,
        1770825370650,
        1799070927569,
        1828830164540,
        1857929484664,
        1887032644581,
        1916995590399,
        1947000281810,
        1977053105308,
        2007947259995,
        2038113914664,
        2069105157658,
        2100148008735,
        2132041948331,
        2163249620300,
        2195290439278,
        2228096092405,
        2260246719816,
        2293254264779,
        2327102968698,
        2360273997212,
        2394317199394,
        2429197537790,
        2464165251412,
        2499250636001,
        2535173373290,
        2571256509291,
        2607459566240,
        2645235990340,
        2682438102505,
        2721166179430,
        2759481758269,
        2799346386970,
        2839371161917,
        2879572973290,
        2921417049495,
        2963464399148,
        3006439953053,
        3049644331975,
        3094535326835,
        3139670804689,
        3186521785254,
        3233657690608,
        3282521316303,
        3331771080889,
        3382807132061,
        3435662986988,
        3489629955549,
        3544787184542,
        3601800785747,
        3660712973701,
        3720915378753,
        3783794471564,
        3848092172594,
        3915841222060,
        3985795414949,
        4058165839496,
        4134194750959,
        4214018777167,
        4297086872145,
        4384179603378,
        4476133915516,
        4573190668592,
        4676247289029,
        4786251270318,
        4904412656668,
        5033600802562,
        5181449296042,
        5370697848017,
        11182769343501  # 60 times Tryptophan
    ]
    """Upper limit of peptide partition (partition are balanced)
    """

    def __init__(self, sequence: str, number_of_missed_cleavages: int):
        self.__sequence = sequence.upper()
        self.__number_of_missed_cleavages = number_of_missed_cleavages
        self.__mass = self.__class__.calculate_mass(self.__sequence)
        self.__partition = self.__class__.get_partition(self.__mass)
        # On demand values
        self.__amino_acid_counter = None

    @property
    def sequence(self):
        return self.__sequence

    @property
    def mass(self):
        """
        Returns
        -------
        Peptide mass
        """
        return self.__mass

    @property
    def partition(self):
        """
        Returns
        -------
        Peptide partition
        """
        return self.__partition

    @property
    def number_of_missed_cleavages(self):
        """
        Returns
        -------
        Number of missed cleavages
        """
        return self.__number_of_missed_cleavages

    @property
    def a_count(self):
        """
        Returns
        -------
        Number of A in the sequence
        """
        if not self.__amino_acid_counter:
            self.__count_amino_acids()
        return self.__amino_acid_counter['A']

    @property
    def b_count(self):
        """
        Returns
        -------
        Number of B in the sequence
        """
        if not self.__amino_acid_counter:
            self.__count_amino_acids()
        return self.__amino_acid_counter['B']

    @property
    def c_count(self):
        """
        Returns
        -------
        Number of C in the sequence
        """
        if not self.__amino_acid_counter:
            self.__count_amino_acids()
        return self.__amino_acid_counter['C']

    @property
    def d_count(self):
        """
        Returns
        -------
        Number of D in the sequence
        """
        if not self.__amino_acid_counter:
            self.__count_amino_acids()
        return self.__amino_acid_counter['D']

    @property
    def e_count(self):
        """
        Returns
        -------
        Number of E in the sequence
        """
        if not self.__amino_acid_counter:
            self.__count_amino_acids()
        return self.__amino_acid_counter['E']

    @property
    def f_count(self):
        """
        Returns
        -------
        Number of F in the sequence
        """
        if not self.__amino_acid_counter:
            self.__count_amino_acids()
        return self.__amino_acid_counter['F']

    @property
    def g_count(self):
        """
        Returns
        -------
        Number of G in the sequence
        """
        if not self.__amino_acid_counter:
            self.__count_amino_acids()
        return self.__amino_acid_counter['G']

    @property
    def h_count(self):
        """
        Returns
        -------
        Number of H in the sequence
        """
        if not self.__amino_acid_counter:
            self.__count_amino_acids()
        return self.__amino_acid_counter['H']

    @property
    def i_count(self):
        """
        Returns
        -------
        Number of I in the sequence
        """
        if not self.__amino_acid_counter:
            self.__count_amino_acids()
        return self.__amino_acid_counter['I']

    @property
    def j_count(self):
        """
        Returns
        -------
        Number of J in the sequence
        """
        if not self.__amino_acid_counter:
            self.__count_amino_acids()
        return self.__amino_acid_counter['J']

    @property
    def k_count(self):
        """
        Returns
        -------
        Number of K in the sequence
        """
        if not self.__amino_acid_counter:
            self.__count_amino_acids()
        return self.__amino_acid_counter['K']

    @property
    def l_count(self):
        """
        Returns
        -------
        Number of L in the sequence
        """
        if not self.__amino_acid_counter:
            self.__count_amino_acids()
        return self.__amino_acid_counter['L']

    @property
    def m_count(self):
        """
        Returns
        -------
        Number of M in the sequence
        """
        if not self.__amino_acid_counter:
            self.__count_amino_acids()
        return self.__amino_acid_counter['M']

    @property
    def n_count(self):
        """
        Returns
        -------
        Number of N in the sequence
        """
        if not self.__amino_acid_counter:
            self.__count_amino_acids()
        return self.__amino_acid_counter['N']

    @property
    def o_count(self):
        """
        Returns
        -------
        Number of O in the sequence
        """
        if not self.__amino_acid_counter:
            self.__count_amino_acids()
        return self.__amino_acid_counter['O']

    @property
    def p_count(self):
        """
        Returns
        -------
        Number of P in the sequence
        """
        if not self.__amino_acid_counter:
            self.__count_amino_acids()
        return self.__amino_acid_counter['P']

    @property
    def q_count(self):
        """
        Returns
        -------
        Number of Q in the sequence
        """
        if not self.__amino_acid_counter:
            self.__count_amino_acids()
        return self.__amino_acid_counter['Q']

    @property
    def r_count(self):
        """
        Returns
        -------
        Number of R in the sequence
        """
        if not self.__amino_acid_counter:
            self.__count_amino_acids()
        return self.__amino_acid_counter['R']

    @property
    def s_count(self):
        """
        Returns
        -------
        Number of S in the sequence
        """
        if not self.__amino_acid_counter:
            self.__count_amino_acids()
        return self.__amino_acid_counter['S']

    @property
    def t_count(self):
        """
        Returns
        -------
        Number of T in the sequence
        """
        if not self.__amino_acid_counter:
            self.__count_amino_acids()
        return self.__amino_acid_counter['T']

    @property
    def u_count(self):
        """
        Returns
        -------
        Number of U in the sequence
        """
        if not self.__amino_acid_counter:
            self.__count_amino_acids()
        return self.__amino_acid_counter['U']

    @property
    def v_count(self):
        """
        Returns
        -------
        Number of V in the sequence
        """
        if not self.__amino_acid_counter:
            self.__count_amino_acids()
        return self.__amino_acid_counter['V']

    @property
    def w_count(self):
        """
        Returns
        -------
        Number of W in the sequence
        """
        if not self.__amino_acid_counter:
            self.__count_amino_acids()
        return self.__amino_acid_counter['W']

    @property
    def y_count(self):
        """
        Returns
        -------
        Number of Y in the sequence
        """
        if not self.__amino_acid_counter:
            self.__count_amino_acids()
        return self.__amino_acid_counter['Y']

    @property
    def z_count(self):
        """
        Returns
        -------
        Number of Z in the sequence
        """
        if not self.__amino_acid_counter:
            self.__count_amino_acids()
        return self.__amino_acid_counter['Z']

    @property
    def n_terminus(self):
        """
        Returns
        -------
        First amino acid in the sequence.
        """
        return self.sequence[0]

    @property
    def c_terminus(self):
        """
        Returns
        -------
        Last amino acid in the sequence.
        """
        return self.sequence[-1]

    def get_n_terminus_ascii_dec(self) -> int:
        """
        Returns
        -------
        ASCII decimal represenation of the first amino acid in the sequence.
        """
        return ord(self.sequence[0])

    def get_c_terminus_ascii_dec(self) -> int:
        """
        Returns
        -------
        ASCII decimal represenation of the last amino acid in the sequence.
        """
        return ord(self.sequence[-1])

    @property
    def length(self):
        """
        Returns
        -------
        Length of the sequence.
        """
        return len(self.sequence)

    @classmethod
    def calculate_mass(cls, sequence: str) -> int:
        """
        Calculates the mass of the given amino acid sequence.

        Parameters
        ----------
        sequence: str
            Amino acid sequence
        """
        mass = H2O.mono_mass
        for amino_acid_one_letter_code in sequence:
            mass += AminoAcid.get_by_one_letter_code(amino_acid_one_letter_code).mono_mass
        return mass

    def __count_amino_acids(self):
        """
        Counts the amino acid in the peptides sequence.
        """
        self.__amino_acid_counter = Counter(self.__sequence)

    def __hash__(self):
        """
        Implements the ability to use peptides as key in sets and dictionaries.
        """
        return hash((self.mass, self.sequence))

    def __eq__(self, other):
        """
        Implements the equals operator.
        According to the Python documentation this should be implemented if __hash__() is implemented.
        """
        return self.sequence == other.sequence

    def __get_fasta_header(self):
        return self.__class__.FASTA_HEADER_PREFIX

    def to_fasta_entry(self):
        """
        Creates a fasta entry for the peptides.
        """
        return "{}\n{}\n".format(self.__get_fasta_header(), self.sequence)


    @classmethod
    def select(cls, database_cursor, where_condition: Optional[WhereCondition] = None, 
        order_by: Optional[str] = None, fetchall: bool = False, stream: bool = False) -> Optional[Union[PeptideBase, List[PeptideBase], Iterator[PeptideBase]]]:
        """
        Selects peptides.

        Parameters
        ----------
        database_cursor
            Database cursor
        where_condition : WhereCondition
            Where condition (optional)
        order_by : str
            Order by instruction, e.g `mass DESC` (optional)
        fetchall : bool
            Indicates if multiple rows should be fetched
        stream : bool
            If true, a generator is returned which yields all matching PeptideBase records

        Returns
        -------
        None, Petide, list of peptides or generator which yield peptides
        """
        select_query = f"SELECT sequence, number_of_missed_cleavages FROM {cls.TABLE_NAME}"
        select_values = []
        if where_condition is not None:
            select_query += f" WHERE {where_condition.get_condition_str()}"
            select_values += where_condition.values
        if order_by is not None:
            select_query += f" ORDER BY {order_by}"
        select_query += ";"
        database_cursor.execute(select_query, select_values)

        if not stream:
            if fetchall:
                return [cls(row[0], row[1]) for row in database_cursor.fetchall()]
            else:
                row = database_cursor.fetchone()
                if row:
                    return cls(row[0], row[1])
                else:
                    return None
        else:
            def gen():
                for row in database_cursor:
                    yield cls(row[0], row[1])
            return gen()

    @classmethod
    def insert(cls, database_cursor, peptide):
        """
        Inserts the given peptide.

        Arguments
        database_cursor
            Database cursor with open transaction.
        peptide : Peptide
            Peptide to insert
        """
        INSERT_QUERY = (
            f"INSERT INTO {cls.TABLE_NAME} (partition, mass, sequence, length, number_of_missed_cleavages, a_count, b_count, c_count, d_count, e_count, f_count, g_count, h_count, i_count, j_count, k_count, l_count, m_count, n_count, o_count, p_count, q_count, r_count, s_count, t_count, u_count, v_count, w_count, y_count, z_count, n_terminus, c_terminus) "
            "VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s) ON CONFLICT DO NOTHING;"
        )
        database_cursor.execute(
            INSERT_QUERY,
            (
                peptide.partition,
                peptide.mass,
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
                peptide.get_n_terminus_ascii_dec(),
                peptide.get_c_terminus_ascii_dec()
            )
        )
        return database_cursor.rowcount

    @classmethod
    def bulk_insert(cls, database_cursor, peptides: list) -> int:
        """
        Bulk inserts many peptides efficiently.

        Arguemnts
        =========
        database_cursor
            Database cursor with open transaction.
        peptides : List[Peptide]
            Peptides for bulk insert.
        """
        BULK_INSERT_QUERY = (
            f"INSERT INTO {cls.TABLE_NAME} (partition, mass, sequence, length, number_of_missed_cleavages, a_count, b_count, c_count, d_count, e_count, f_count, g_count, h_count, i_count, j_count, k_count, l_count, m_count, n_count, o_count, p_count, q_count, r_count, s_count, t_count, u_count, v_count, w_count, y_count, z_count, n_terminus, c_terminus) "
            "VALUES %s ON CONFLICT DO NOTHING;"
        )
        # Bulk insert the new peptides
        execute_values(
            database_cursor,
            BULK_INSERT_QUERY,
            [
                (
                    peptide.partition,
                    peptide.mass,
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
                    peptide.get_n_terminus_ascii_dec(),
                    peptide.get_c_terminus_ascii_dec()
                ) for peptide in peptides
            ],
            page_size=len(peptides)
        )
        # rowcount is only accurate, because the page size is as high as the number of inserted data. If the page size would be smaller rowcount would only return the rowcount of the last processed page.
        return database_cursor.rowcount

    @classmethod
    def get_partition(cls, mass: int) -> int:
        """
        Returns the partition id for the given mass.

        Parameters
        ----------
        mass : int
            Peptide mass

        Returns
        -------
        -1 if mass is not in any partition.
        """
        for idx in range(len(cls.PARTITONS)):
            if mass < cls.PARTITONS[idx]:
                return idx
        return -1
