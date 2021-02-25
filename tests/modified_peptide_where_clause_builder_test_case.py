import unittest
import pathlib
import re
import os

from macpepdb.proteomics.modification_collection import ModificationCollection
from macpepdb.models.modified_peptide_where_clause_builder import ModifiedPeptideWhereClauseBuilder
from macpepdb.proteomics.mass.convert import to_int as mass_to_int, thomson_to_dalton
from macpepdb.models.peptide import Peptide
from macpepdb.tasks.digestion import Digestion
from macpepdb.proteomics.enzymes.digest_enzyme import DigestEnzyme
from macpepdb.peptide_mass_validator import PeptideMassValidator
from macpepdb.proteomics.mass.precursor_range import PrecursorRange

from .abstract_database_test_case import AbstractDatabaseTestCase

# Values from scan 3224 in test file
MASS_TO_CHARGE_RATIO = 442.970306396484
CHARGE = 3
PRECURSOR_TOLERANCE = 5
VARIABLE_MODIFICATION_MAXIMUM = 3

PEPTIDES_FOR_UNMODIFIED_SEARCH = {
    'too_low': [
        # 10 ppm
        'GLLVIGILFIIR',	    # 1325.879879215
        'LFIAVIIAVVLR',	    # 1325.879879215
        'IFLILVGXLLGIR',	# 1325.879879215
        'LFILLVGGILXIR' 	# 1325.879879215
    ],
    'matching': [
        # first matching
        'KVKLLLLQMLK',	    # 1325.88325005
        # some matching between
        'IMLKIKNLLLK',	    # 1325.883250051
        'KKAIIILCLIAK',	    # 1325.883250052
        'LTTIIILISIVK',	    # 1325.889775244
        'FLLLGIIRLIR',	    # 1325.891112585
        # last matching
        'AILLLVRFLIR'	    # 1325.891112585
    ],
    'too_high': [
        # 10 ppm
        'KAIITSKILLVK',	    # 1325.901008612
        'STLLKLLKGILK',	    # 1325.901008612
        'ILKLTGLIKSIK',	    # 1325.901008612
        'LVLGKKTILITK'	    # 1325.901008612
    ]
}

# precursor
# 1325887444084 == 1325.887444084
PEPTIDES_FOR_MODIFIED_SEARCH = {
    'too_low': [
        # 25 ppm
        'IIVCALKRILK',	    # 1268.836634178
        'ILIIGCIKIKR',	    # 1268.836634178
        'RKLGIIICLLK',	    # 1268.836634178
        'KIIVALLLRCK'	    # 1268.836634178
    ],
    'matching': [
        # first matching
        'VVLKKIIICIK',	    # 1268.861786317
        # some matching between
        'MAKVVLLIAIKK',     # 1325.883250051
        'ILMIIILKTIR',	    # 1325.883250053
        'LKRPLKFLAIK',	    # 1325.891112582
        'IRLAFVIIIIR',	    # 1325.891112585
        # last matching
        'ALLIIVRFLLR'	    # 1325.891112585
    ],
    'too_high': [
        # 15 ppm
        'KALITSKLLIVK',     # 1325.901008612
        'STILKILKGLLK',     # 1325.901008612
        'ILKLTGLLKSLK',     # 1325.901008612
        'IVIGKKTLLLTK'      # 1325.901008612
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

    def test_with_real_data(self):
        VARIABLE_MODIFICATION_LIMIT = 2
        # Mass of MFPVTJEDTEGNVJTVSPPCYGFJQJR
        PRECURSOR = 3025492916648
        PRECURSOR_TOLERANCE = 20
        # With the given mass, tolerance and modifications 3 peptides shoud be found

        protein_file_path = pathlib.Path('./test_files/UP000002006.txt')
        modifications_file_path = pathlib.Path('./test_files/modifications.csv')

        digestion = Digestion(
            [protein_file_path],
            pathlib.Path("./digestion_test.log"),
            pathlib.Path("./digestion_test.unprocessible_proteins.txt"),
            pathlib.Path("./digestion_test.statistics.csv"),
            5,
            4,
            "Trypsin",
            2,
            6   ,
            50
        )

        digestion.digest_to_database(os.getenv("TEST_MACPEPDB_URL"))

        modification_collection = ModificationCollection.read_from_csv_file(modifications_file_path)
        peptide_mass_validator = PeptideMassValidator(modification_collection, VARIABLE_MODIFICATION_LIMIT, PrecursorRange(PRECURSOR, PRECURSOR_TOLERANCE, PRECURSOR_TOLERANCE))

        validated_matching_peptide_sequences = set()

        session = self.session_factory()
        # Run through all peptidesv (in batches of 100 peptides) and check which match the precursor and modification requirements
        for peptide in session.query(Peptide).yield_per(1000):
            if peptide_mass_validator.validate(peptide):
                validated_matching_peptide_sequences.add(peptide.sequence)

        # Close and delete old session and create a new one so there is no cache
        session.close()
        del session
        session = self.session_factory()

        builder =  ModifiedPeptideWhereClauseBuilder(modification_collection, PRECURSOR, PRECURSOR_TOLERANCE, PRECURSOR_TOLERANCE, VARIABLE_MODIFICATION_LIMIT)
        where_clause = builder.build(Peptide)
        peptides = session.query(Peptide).filter(where_clause).all()

        queried_matching_peptide_sequences = set()
        for peptide in peptides:
            queried_matching_peptide_sequences.add(peptide.sequence)

        # Check length of both manually validated set and the queried set
        self.assertEqual(len(validated_matching_peptide_sequences), len(queried_matching_peptide_sequences))

        # Cross check if peptide from one set is in the other set
        for sequence in queried_matching_peptide_sequences:
            self.assertIn(sequence, validated_matching_peptide_sequences)

        for sequence in validated_matching_peptide_sequences:
            self.assertIn(sequence, queried_matching_peptide_sequences)



