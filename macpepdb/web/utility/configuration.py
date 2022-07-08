# std imports
from __future__ import annotations
import argparse
from enum import Enum
import os
from pathlib import Path
from typing import ClassVar, Dict, Any, Optional, Type

# 3rd party imports
from yaml import load as yaml_load, Loader as YamlLoader

class Environment(Enum):
    production = "production"
    development = "development"

    def __str__(self) -> str:
        return f"{self.value}"

    @classmethod
    def from_str(cls, value: str) -> Environment:
        return cls[value.lower()]


class Configuration():
    DEFAULT_CONFIG: ClassVar[str] = """
# Interface to use (127.0.0.1 = all interfaces)
interface: 127.0.0.1
# Port to use
port: 3000
# An aritrary ascii string to sign sessions etc. Make sure to back it up!
secret: "development"
# Debug outputs
debug: true
# Database which contains the proteins and peptides
macpepdb:
  # Database url: 
  url: postgresql://postgres:developer@127.0.0.1:5433/macpepdb_dev
  # Number of connections
  pool_size: 1 # Bjoern (the application server) is single threaded, so it is useless to support more then one connection
matomo:
  enabled: false
  url: ""
  site_id: 1
  auth_token: ""
"""
    __values: ClassVar[Dict[str, Any]] = {}
    __environment: ClassVar[Environment] = Environment.production

    @staticmethod
    def values() -> Dict[str, Any]:
        return Configuration.__values
    
    @staticmethod
    def environment() -> Environment:
        return Configuration.__environment

    @classmethod
    def initialize(cls, config_path: Optional[Path] = None, environment: Optional[Environment] = None):
        if environment is None:
            environment = Environment.from_str(os.getenv('MACPEPDB_WEB_ENV', str(Environment.development)))

        config = yaml_load(cls.DEFAULT_CONFIG, Loader=YamlLoader)
        if config_path:
            with config_path.open("r") as config_file:
                new_config = yaml_load(config_file.read(), Loader=YamlLoader)
                config = Configuration._merge_dicts_recursively(new_config, config)
        Configuration._validate_config(config)

        Configuration.__values = config
        Configuration.__environment = environment

    @classmethod
    def _validate_config(cls, config: Dict[str, Any]):
        """
        Checks config values

        Parameters
        ----------
        config : Dict[str, Any]
            Config

        Raises
        ------
        KeyError
            If key path not found
        """
        try:
            cls._validate_type(config['debug'], bool, 'boolean', 'debug')
            cls._validate_type(config['interface'], str, 'ip string', 'interface')
            cls._validate_type(config['port'], int, 'integer', 'port')
            cls._validate_ascii_string(config['secret'], 'secret')
            cls._validate_psql_url(config['macpepdb']['url'], 'macpepdb.url')
            cls._validate_type(config['macpepdb']['pool_size'], int, 'integer', 'macpepdb.pool_size')
        except KeyError as key_error:
            raise KeyError(f"The configuration key {key_error} is missing.")

    @staticmethod
    def _validate_psql_url(url: Any, key_path: str) -> bool:
        """
        Checks if url is Postgresql URL

        Parameters
        ----------
        url : Any
            Candidate to check
        key_path : str
            Path if keys to value

        Returns
        -------
        bool

        Raises
        ------
        TypeError
            If not a postgres url
        """
        Configuration._validate_type(url, str, 'string', key_path)
        if not url.startswith('postgresql://'):
            raise TypeError(f"{key_path} must start with 'postgresql://'.")
        return True

    @staticmethod
    def _validate_type(value: Any, expected_type: Type, expected_type_as_str: str, key_path: str) -> bool:
        """
        Checks if the given string is an ASCII string

        Parameters
        ----------
        value : Any
            Candidate to check
        expected_type : Type
            Expected type
        expected_type_as_str : str
            Excpected type as string
        key_path : str
            Path if keys to value

        Returns
        -------
        bool

        Raises
        ------
        TypeError
           If value is of expected type
        """
        if not isinstance(value, expected_type):
            raise TypeError(f"Configuration key '{key_path}' is not of type {expected_type_as_str}.")
        return True

    @staticmethod
    def _validate_ascii_string(value: Any, key_path: str) -> bool:
        """
        Checks if the given value is an ASCII string

        Parameters
        ----------
        value : Any
            Candidate to check
        key_path : str
            Path if keys to value

        Returns
        -------
        bool

        Raises
        ------
        TypeError
            If value is not an ASCII string
        """
        Configuration._validate_type(value, str, 'string', key_path)
        if not all(ord(char) < 128 for char in value):
            raise TypeError(f"Configuration key '{key_path}' contains non ascii character.")
        return True

    @staticmethod
    def _merge_dicts_recursively(source: Dict[str, Any], destination: Dict[str, Any]) -> Dict[str, Any]:
        """
        Merges source dictionary into destintion dictionary.

        Parameters
        ----------
        source : Dict[str, Any]
            Source dictionary
        destination : Dict[str, Any]
            Destination dictionary

        Returns
        -------
        Dict[str, Any]
            Merged dictionary
        """
        for key, value in source.items():
            if isinstance(value, dict):
                node = destination.setdefault(key, {})
                Configuration._merge_dicts_recursively(value, node)
            else:
                destination[key] = value

        return destination

    @classmethod
    def write_default_config(cls, config_path: Path):
        """
        Writes default config to file

        Parameters
        ----------
        config_path : Path
            Path to config file
        """
        with config_path.open("w") as config_file:
            config_file.write(cls.DEFAULT_CONFIG.strip())

    @staticmethod
    def write_config_file_by_cli(cli_args):
        """
        Write

        Parameters
        ----------
        cli_args : Any
            CLI arguments
        """
        Configuration.write_default_config(
            Path(cli_args.path)
        )


    @staticmethod
    def add_cli_arguments(subparsers: argparse._SubParsersAction):
        """
        Defines the CLI parameters for the web server.

        Parameters
        ----------
        subparser : argparse._SubParsersAction
            Subparser of main CLI parser
        """
        parser = subparsers.add_parser("write-config-file", help="Configuratio")
        parser.add_argument("path", help="Output file")
        parser.set_defaults(func=Configuration.write_config_file_by_cli)