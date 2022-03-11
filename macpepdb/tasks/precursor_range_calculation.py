# std imports
from dataclasses import replace
import pathlib

# internal imports
from macpepdb.proteomics.mass.precursor_range import PrecursorRange
from macpepdb.proteomics.mass.convert import to_float as mass_to_float, to_int as mass_to_int
from macpepdb.proteomics.modification_collection import ModificationCollection
from macpepdb.models.modification_combination_list import ModificationCombinationList
from macpepdb.database.indexes.post_translational_modification_search_index import PostTranslationalModificationSearchIndex as PtmIndex


class PrecursorRangeCalculation:
    """
    Calculates a precursor / mass range from a mass and upper / lower tolerance.

    Parameters
    ----------
    precursor : int
        Precursor / mass
    lower_tolerance_ppm : int
        Lower tolerance (ppm)
    upper_tolerance_ppm : int
        Upper tolerance (ppm)
    modification_collection: ModificationCollection
        Collection of modifications
    max_number_of_variable_modifications: int
        Number of maximum variable modifications per peptide.
    """

    def __init__(self, precursor: int, lower_precursor_tolerance_ppm: int, upper_precursor_tolerance_ppm: int, modification_collection: ModificationCollection, max_number_of_variable_modifications: int):
        self.__precursor_range = PrecursorRange(precursor, lower_precursor_tolerance_ppm, upper_precursor_tolerance_ppm)
        self.__max_number_of_variable_modifications = max_number_of_variable_modifications
        self.__modification_collection_list = ModificationCombinationList(
            modification_collection,
            self.__precursor_range.precursor,
            upper_precursor_tolerance_ppm,
            lower_precursor_tolerance_ppm,
            self.__max_number_of_variable_modifications
        )

    def __str__(self):
        result = ""
        for combination in self.__modification_collection_list:
            mass = ""
            amino_acid_counts = []

            for condition in combination.where_conditions:
                if condition.column_name == "partition":
                    continue
                output = f"{condition.column_name} {condition.operator}"
                if condition.column_name == "mass":
                    mass = output % tuple([mass_to_float(mass) for mass in condition.values])
                else:
                    amino_acid_counts.append(output % condition.values)
            amino_acid_counts.sort()
            result += f"{mass} ({' & '.join(amino_acid_counts)})"
            result += "\n"
        return result .replace(" BETWEEN", ":").replace("AND", "-")

    @classmethod
    def start_from_comand_line(cls, args):
        """
        Starts a precursor range calculation with the arguments from the CLI.

        Parameters
        ----------
        args
            Arguments from the CLI parser
        """
        calculation = cls(
            mass_to_int(args.precursor),
            args.lower_precursor_tolerance,
            args.upper_precursor_tolerance,
            ModificationCollection.read_from_csv_file(pathlib.Path(args.modifications)) if args.modifications is not None else ModificationCollection([]),
            args.max_variable_modifications
        )
        print(calculation)

    @classmethod
    def comand_line_arguments(cls, subparsers):
        """
        Defines the CLI parameters for the precursor range calculation in the given subparser.

        Parameters
        ----------
        subparser : argparse._SubParsersAction
            Subparser of main CLI parser
        """
        parser = subparsers.add_parser('precursor-range', help="Calculates the precursor range for the given precursor, tolerance and PTMs.")
        parser.add_argument("--precursor", "-p", type=float, required=True, help="Precursor")
        parser.add_argument("--lower-precursor-tolerance", "-l", type=int, required=True, help="Lower precursor tolerance (ppm)")
        parser.add_argument("--upper-precursor-tolerance", "-u", type=int, required=True, help="Upper precursor tolerance (ppm)")
        parser.add_argument("--modifications", "-m", type=str, required=False, help="CSV with post translational modifications.")
        parser.add_argument("--max-variable-modifications", "-v", type=int, default=3, required=False, help="Maximum number of variable modifications per peptide.")
        parser.set_defaults(func=cls.start_from_comand_line)
