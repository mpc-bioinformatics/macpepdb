# internal imports
from macpepdb.proteomics.mass.precursor_range import PrecursorRange
from macpepdb.proteomics.mass.convert import to_float as mass_to_float, to_int as mass_to_int


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
    """

    def __init__(self, precursor: int, lower_precursor_tolerance_ppm: int, upper_precursor_tolerance_ppm: int):
        self.__precursor_range = PrecursorRange(precursor, lower_precursor_tolerance_ppm, upper_precursor_tolerance_ppm)

    @property
    def precursor_range(self):
        """
        Returns
        -------
        Precursor range
        """
        return self.__precursor_range

    def __str__(self):
        return "{} - {}".format(mass_to_float(self.precursor_range.lower_limit), mass_to_float(self.precursor_range.upper_limit))

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
            args.upper_precursor_tolerance
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
        parser = subparsers.add_parser('precursor-range', help="Calculates the precursor range for the given precursor and tolerance.")
        parser.add_argument("--precursor", "-p", type=float, required=True, help="Precursor")
        parser.add_argument("--lower-precursor-tolerance", "-l", type=int, required=True, help="Lower precursor tolerance (ppm)")
        parser.add_argument("--upper-precursor-tolerance", "-u", type=int, required=True, help="Upper precursor tolerance (ppm)")
        parser.set_defaults(func=cls.start_from_comand_line)
