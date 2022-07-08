# std imports
import argparse

# internal imports
from macpepdb.tasks.database_maintenance.database_maintenance import DatabaseMaintenance
from macpepdb.tasks.precursor_range_calculation import PrecursorRangeCalculation
from macpepdb.tasks.statistics import Statistics
from macpepdb.web.cli import ComandLineInterface as WebComandLineInterface

class ComandLineInterface():
    """
    Defines the command line interface by adding the `command_line_arguments` of each class in the `task`-module
    """

    def __init__(self):
        self.__parser = argparse.ArgumentParser(description='Create a large peptide database')
        subparsers = self.__parser.add_subparsers()
        DatabaseMaintenance.comand_line_arguments(subparsers)
        PrecursorRangeCalculation.comand_line_arguments(subparsers)
        Statistics.command_line_arguments(subparsers)
        WebComandLineInterface.add_cli_arguments(subparsers)

    def start(self):
        """
        Starts the command line parsing
        """
        args = self.__parser.parse_args()
        args.func(args)