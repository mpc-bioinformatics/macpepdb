import argparse

from macpep_db.tasks.digestion import Digestion
from .tasks.precursor_range_calculation import PrecursorRangeCalculation
from .tasks.taxonomy_maintenance import TaxonomyMaintenance
from .tasks.statistics import Statistics

class ComandLineInterface():
    def __init__(self):
        self.__parser = argparse.ArgumentParser(description='Create a large peptide database')
        subparsers = self.__parser.add_subparsers()
        Digestion.comand_line_arguments(subparsers)
        PrecursorRangeCalculation.comand_line_arguments(subparsers)
        TaxonomyMaintenance.comand_line_arguments(subparsers)
        Statistics.command_line_arguments(subparsers)

    def start(self):
        args = self.__parser.parse_args()
        args.func(args)