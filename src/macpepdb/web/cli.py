# std imports
import argparse

# internal imports
from macpepdb.web import server
from macpepdb.web.utility.configuration import Configuration


class ComandLineInterface:
    def add_cli_arguments(subparsers: argparse._SubParsersAction):
        parser = subparsers.add_parser("web", help="WebAPI for MaCPepDB")
        parser.set_defaults(func=lambda args: parser.print_help())
        subparsers = parser.add_subparsers()
        server.add_cli_arguments(subparsers)
        Configuration.add_cli_arguments(subparsers)

