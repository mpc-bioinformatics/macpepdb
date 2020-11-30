import unittest
import pathlib
import re

from trypperdb.proteomics.modification_collection import ModificationCollection
from trypperdb.models.modified_peptide_where_clause_builder import ModifiedPeptideWhereClauseBuilder
from trypperdb.proteomics.mass.convert import to_int as mass_to_int, thomson_to_dalton
from trypperdb.models.peptide import Peptide

from .abstract_database_test_case import AbstractDatabaseTestCase

# Values from scan 3224 in test file
MASS_TO_CHARGE_RATIO = 442.970306396484
CHARGE = 3
PRECURSOR_TOLERANCE = 5
VARIABLE_MODIFICATION_MAXIMUM = 3

PEPTIDES_FOR_UNMODIFIED_SEARCH = {
    'too_low': [
        # 10 ppm
        'GJJVJGJJFJJR',	    # 1325.879879215
        'JFJAVJJAVVJR',	    # 1325.879879215
        'JFJJJVGXJJGJR',	# 1325.879879215
        'JFJJJVGGJJXJR' 	# 1325.879879215
    ],
    'matching': [
        # first matching
        'KVKJJJJQMJK',	    # 1325.88325005
        # some matching between
        'JMJKJKNJJJK',	    # 1325.883250051
        'KKAJJJJCJJAK',	    # 1325.883250052
        'JTTJJJJJSJVK',	    # 1325.889775244
        'FJJJGJJRJJR',	    # 1325.891112585
        # last matching
        'AJJJJVRFJJR'	    # 1325.891112585
    ],
    'too_high': [
        # 10 ppm
        'KAJJTSKJJJVK',	    # 1325.901008612
        'STJJKJJKGJJK',	    # 1325.901008612
        'JJKJTGJJKSJK',	    # 1325.901008612
        'JVJGKKTJJJTK'	    # 1325.901008612
    ]
}

# precursor
# 1325887444084 == 1325.887444084
PEPTIDES_FOR_MODIFIED_SEARCH = {
    'too_low': [
        # 25 ppm
        'JJVCAJKRJJK',	    # 1268.836634178
        'JJJJGCJKJKR',	    # 1268.836634178
        'RKJGJJJCJJK',	    # 1268.836634178
        'KJJVAJJJRCK'	    # 1268.836634178
    ],
    'matching': [
        # first matching
        'VVJKKJJJCJK',	    # 1268.861786317
        # some matching between
        'MAKVVJJJAJKK',     # 1325.883250051
        'JJMJJJJKTJR',	    # 1325.883250053
        'JKRPJKFJAJK',	    # 1325.891112582
        'JRJAFVJJJJR',	    # 1325.891112585
        # last matching
        'AJJJJVRFJJR'	    # 1325.891112585
    ],
    'too_high': [
        # 15 ppm
        'KAJJTSKJJJVK',     # 1325.901008612
        'STJJKJJKGJJK',     # 1325.901008612
        'JJKJTGJJKSJK',     # 1325.901008612
        'JVJGKKTJJJTK'      # 1325.901008612
    ]
}

class ModifiedPeptideWhereClauseBuilderTestCase(AbstractDatabaseTestCase):
    WEIGHT_TOLERANCE_REGEX = re.compile(r"BETWEEN\W(?P<lower>\d+)\WAND\W(?P<upper>\d+)(\W|$)")

    def test_without_modifications(self):
        # Add peptides to database
        session = self.session_factory()
        for key in PEPTIDES_FOR_UNMODIFIED_SEARCH.keys():
            for sequence in PEPTIDES_FOR_UNMODIFIED_SEARCH[key]:
                peptide = Peptide(sequence, 0)
                session.add(peptide)
        session.commit()
        session.close()

        modification_collection = ModificationCollection([])
        precursor = mass_to_int(thomson_to_dalton(MASS_TO_CHARGE_RATIO, CHARGE))

        # Create fresh session
        session = self.session_factory()
        builder =  ModifiedPeptideWhereClauseBuilder(modification_collection, precursor, PRECURSOR_TOLERANCE, PRECURSOR_TOLERANCE, VARIABLE_MODIFICATION_MAXIMUM)
        where_clause = builder.build(Peptide)

        where_clause_string = str(where_clause.compile(compile_kwargs={"literal_binds": True}))
        matches = re.findall(self.__class__.WEIGHT_TOLERANCE_REGEX, where_clause_string)
        # Without modifications there is only one between-condition.
        self.assertEqual(len(matches), 1)

        peptides = session.query(Peptide).filter(where_clause).all()

        # Check if only matching peptides were found
        self.assertEqual(len(peptides), len(PEPTIDES_FOR_UNMODIFIED_SEARCH['matching']))
        for peptide in peptides:
            self.assertIn(peptide.sequence, PEPTIDES_FOR_UNMODIFIED_SEARCH['matching'])

    def test_with_modifications(self):# lower hit
        # Add peptides to database
        session = self.session_factory()
        for key in PEPTIDES_FOR_MODIFIED_SEARCH.keys():
            for sequence in PEPTIDES_FOR_MODIFIED_SEARCH[key]:
                peptide = Peptide(sequence, 0)
                session.add(peptide)
        session.commit()
        session.close()

        csv_file_path = pathlib.Path("./test_files/modifications.csv")
        modification_collection = ModificationCollection.read_from_csv_file(csv_file_path)
        
        precursor = mass_to_int(thomson_to_dalton(MASS_TO_CHARGE_RATIO, CHARGE))

        builder =  ModifiedPeptideWhereClauseBuilder(modification_collection, precursor, PRECURSOR_TOLERANCE, PRECURSOR_TOLERANCE, VARIABLE_MODIFICATION_MAXIMUM)

        where_clause = builder.build(Peptide)
        # Create fresh session
        session = self.session_factory()
        builder =  ModifiedPeptideWhereClauseBuilder(modification_collection, precursor, PRECURSOR_TOLERANCE, PRECURSOR_TOLERANCE, VARIABLE_MODIFICATION_MAXIMUM)
        where_clause = builder.build(Peptide)
        peptides = session.query(Peptide).filter(where_clause).all()

        # Check if only matching peptides were found
        self.assertEqual(len(peptides), len(PEPTIDES_FOR_MODIFIED_SEARCH['matching']))
        for peptide in peptides:
            self.assertIn(peptide.sequence, PEPTIDES_FOR_MODIFIED_SEARCH['matching'])