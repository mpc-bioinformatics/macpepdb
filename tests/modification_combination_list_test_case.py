# std imports
import pathlib
import re
import os

# internal imports
from macpepdb.models.modification_combination_list import ModificationCombinationList
from macpepdb.models.peptide import Peptide
from macpepdb.peptide_mass_validator import PeptideMassValidator
from macpepdb.proteomics.mass.convert import to_int as mass_to_int, thomson_to_dalton
from macpepdb.proteomics.mass.precursor_range import PrecursorRange
from macpepdb.proteomics.modification_collection import ModificationCollection
from macpepdb.tasks.database_maintenance.database_maintenance import DatabaseMaintenance


# test imports
from tests.abstract_database_test_case import AbstractDatabaseTestCase
from tests.database_maintenance_workdir import DatabaseMaintenanceWorkdir

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

class ModifiedPeptideWhereClauseBuilderTestCase(AbstractDatabaseTestCase, DatabaseMaintenanceWorkdir):
    MASS_TOLERANCE_REGEX = re.compile(r"mass\WBETWEEN\W(?P<lower>\d+)\WAND\W(?P<upper>\d+)(\W|$)")

    def test_without_modifications(self):
        # Add peptides to database
        with self.database_connection:
            with self.database_connection.cursor() as database_cursor:
                for key in PEPTIDES_FOR_UNMODIFIED_SEARCH.keys():
                    Peptide.bulk_insert(database_cursor, [Peptide(sequence, 0) for sequence in PEPTIDES_FOR_UNMODIFIED_SEARCH[key]])

        modification_collection = ModificationCollection([])
        precursor = mass_to_int(thomson_to_dalton(MASS_TO_CHARGE_RATIO, CHARGE))

        with self.database_connection:
            with self.database_connection.cursor() as database_cursor:
                modification_combination_list =  ModificationCombinationList(modification_collection, precursor, PRECURSOR_TOLERANCE, PRECURSOR_TOLERANCE, VARIABLE_MODIFICATION_MAXIMUM)

                select_conditions = modification_combination_list.to_sql()
                select_conditions_string = database_cursor.mogrify(select_conditions[0], select_conditions[1]).decode('utf-8')
                matches = re.findall(self.__class__.MASS_TOLERANCE_REGEX, select_conditions_string)
                # Without modifications there is only one between-condition.
                self.assertEqual(len(matches), 1)

                peptides = Peptide.select(database_cursor, select_conditions, True)

                # Check if only matching peptides were found
                self.assertEqual(len(peptides), len(PEPTIDES_FOR_UNMODIFIED_SEARCH['matching']))
                for peptide in peptides:
                    self.assertIn(peptide.sequence, PEPTIDES_FOR_UNMODIFIED_SEARCH['matching'])

    def test_with_modifications(self):# lower hit
        # Add peptides to database
        with self.database_connection:
            with self.database_connection.cursor() as database_cursor:
                for key in PEPTIDES_FOR_MODIFIED_SEARCH.keys():
                    Peptide.bulk_insert(database_cursor, [Peptide(sequence, 0) for sequence in PEPTIDES_FOR_MODIFIED_SEARCH[key]])

        csv_file_path = pathlib.Path("./test_files/modifications.csv")
        modification_collection = ModificationCollection.read_from_csv_file(csv_file_path)
        
        precursor = mass_to_int(thomson_to_dalton(MASS_TO_CHARGE_RATIO, CHARGE))

        with self.database_connection:
            with self.database_connection.cursor() as database_cursor:
                modification_combination_list =  ModificationCombinationList(modification_collection, precursor, PRECURSOR_TOLERANCE, PRECURSOR_TOLERANCE, VARIABLE_MODIFICATION_MAXIMUM)
                select_conditions = modification_combination_list.to_sql()
                peptides = Peptide.select(database_cursor, select_conditions, True)

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

        modifications_file_path = pathlib.Path('./test_files/modifications.csv')

        work_dir = pathlib.Path(f"./tmp/{self.id()}")
        test_files_path = pathlib.Path('./test_files')
        protein_data_test_file_path = test_files_path.joinpath('proteins.txt')
        

        self.prepare_workdir(work_dir, test_files_path, protein_data_test_file_path)

        maintenance = DatabaseMaintenance(
            os.getenv("TEST_MACPEPDB_URL"),
            work_dir,
            4,
            5,
            'Trypsin',
            2,
            5,
            40
        )

        maintenance.start()

        modification_collection = ModificationCollection.read_from_csv_file(modifications_file_path)
        peptide_mass_validator = PeptideMassValidator(modification_collection, VARIABLE_MODIFICATION_LIMIT, PrecursorRange(PRECURSOR, PRECURSOR_TOLERANCE, PRECURSOR_TOLERANCE))

        validated_matching_peptide_sequences = set()

        with self.database_connection:
            with self.database_connection.cursor() as database_cursor:
                # Run through all peptides (in batches of 1000 peptides) and check which matchs the precursor and modification requirements
                database_cursor.execute(f"SELECT sequence, number_of_missed_cleavages FROM {Peptide.TABLE_NAME};")
                while True:
                    rows = database_cursor.fetchmany(1000)
                    if not len(rows):
                        break
                    for row in rows:
                        peptide = Peptide(row[0], row[1])
                        if peptide_mass_validator.validate(peptide):
                            validated_matching_peptide_sequences.add(peptide.sequence)

        with self.database_connection:
            with self.database_connection.cursor() as database_cursor:
                modification_combination_list =  ModificationCombinationList(modification_collection, PRECURSOR, PRECURSOR_TOLERANCE, PRECURSOR_TOLERANCE, VARIABLE_MODIFICATION_LIMIT)
                select_conditions = modification_combination_list.to_sql()
                peptides = Peptide.select(database_cursor, select_conditions, True)

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



