from macpepdb.database.indexes.abstract_index import AbstractIndex
from macpepdb.database.indexes.column_definition import ColumnDefinition

class PostTranslationalModificationSearchIndex(AbstractIndex):
    AMINO_ACID_COUNT_OPERATOR = ">= %s"
    AMINO_ACID_COUNT_VALUE = (0,)
    TERMINI_OPERATOR = "!= %s"
    TERMINI_VALUE = ("0",)
    COLUMN_CONDITIONS_TEMPLATE = [
        ColumnDefinition("partition",   ">= %s",                    (0,)),
        ColumnDefinition("mass",        ">= %s",                    (0,)),
        ColumnDefinition("m_count",     AMINO_ACID_COUNT_OPERATOR,  AMINO_ACID_COUNT_VALUE),
        ColumnDefinition("c_count",     AMINO_ACID_COUNT_OPERATOR,  AMINO_ACID_COUNT_VALUE),
        ColumnDefinition("s_count",     AMINO_ACID_COUNT_OPERATOR,  AMINO_ACID_COUNT_VALUE),
        ColumnDefinition("t_count",     AMINO_ACID_COUNT_OPERATOR,  AMINO_ACID_COUNT_VALUE),
        ColumnDefinition("y_count",     AMINO_ACID_COUNT_OPERATOR,  AMINO_ACID_COUNT_VALUE),
        ColumnDefinition("n_terminus",  TERMINI_OPERATOR,           TERMINI_VALUE),
        ColumnDefinition("c_terminus",  TERMINI_OPERATOR,           TERMINI_VALUE),
        ColumnDefinition("r_count",     AMINO_ACID_COUNT_OPERATOR,  AMINO_ACID_COUNT_VALUE),
        ColumnDefinition("k_count",     AMINO_ACID_COUNT_OPERATOR,  AMINO_ACID_COUNT_VALUE),
        ColumnDefinition("a_count",     AMINO_ACID_COUNT_OPERATOR,  AMINO_ACID_COUNT_VALUE),
        ColumnDefinition("b_count",     AMINO_ACID_COUNT_OPERATOR,  AMINO_ACID_COUNT_VALUE),
        ColumnDefinition("d_count",     AMINO_ACID_COUNT_OPERATOR,  AMINO_ACID_COUNT_VALUE),
        ColumnDefinition("e_count",     AMINO_ACID_COUNT_OPERATOR,  AMINO_ACID_COUNT_VALUE),
        ColumnDefinition("f_count",     AMINO_ACID_COUNT_OPERATOR,  AMINO_ACID_COUNT_VALUE),
        ColumnDefinition("g_count",     AMINO_ACID_COUNT_OPERATOR,  AMINO_ACID_COUNT_VALUE),
        ColumnDefinition("h_count",     AMINO_ACID_COUNT_OPERATOR,  AMINO_ACID_COUNT_VALUE),
        ColumnDefinition("i_count",     AMINO_ACID_COUNT_OPERATOR,  AMINO_ACID_COUNT_VALUE),
        ColumnDefinition("j_count",     AMINO_ACID_COUNT_OPERATOR,  AMINO_ACID_COUNT_VALUE),
        ColumnDefinition("l_count",     AMINO_ACID_COUNT_OPERATOR,  AMINO_ACID_COUNT_VALUE),
        ColumnDefinition("n_count",     AMINO_ACID_COUNT_OPERATOR,  AMINO_ACID_COUNT_VALUE),
        ColumnDefinition("o_count",     AMINO_ACID_COUNT_OPERATOR,  AMINO_ACID_COUNT_VALUE),
        ColumnDefinition("p_count",     AMINO_ACID_COUNT_OPERATOR,  AMINO_ACID_COUNT_VALUE),
        ColumnDefinition("q_count",     AMINO_ACID_COUNT_OPERATOR,  AMINO_ACID_COUNT_VALUE),
        ColumnDefinition("u_count",     AMINO_ACID_COUNT_OPERATOR,  AMINO_ACID_COUNT_VALUE),
        ColumnDefinition("v_count",     AMINO_ACID_COUNT_OPERATOR,  AMINO_ACID_COUNT_VALUE),
        ColumnDefinition("w_count",     AMINO_ACID_COUNT_OPERATOR,  AMINO_ACID_COUNT_VALUE),
        ColumnDefinition("z_count",     AMINO_ACID_COUNT_OPERATOR,  AMINO_ACID_COUNT_VALUE)
    ]