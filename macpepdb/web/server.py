# std imports
import argparse
import json
from pathlib import Path
import re
import traceback
from threading import Thread
from typing import Optional, Union, Dict, Any, List

# 3rd part imports
from flask import Flask, g as request_store, request
from werkzeug.exceptions import HTTPException

# internal imports
from macpepdb.web.database.connection_pool import ConnectionPool
from macpepdb.web.utility.configuration import Configuration, Environment
from macpepdb.web.utility.headers.cross_origin_resource_sharing import add_allow_cors_headers
from macpepdb.web.utility.matomo import track_request as matomo_track_request

app: Optional[Flask] = None
"""Flask application used for serving request
"""

db_pool: Optional[ConnectionPool] = None 
"""Connection pool for web applications database.
Use `get_database_connection` to get a database connection which is valid during a request request.
It will be put back automatically.
However, if you want to use a database connection in a generator for streaming,
you have to manually get and put back the connection.
The best way to deal with it in a generator, is to use a try/catch-block
and put the connection back when GeneratorExit is thrown or in the finally-block.
"""

macpepdb_pool: Optional[ConnectionPool] = None
"""Connection pool for MaCPepDB database.
Use `get_database_connection` to get a database connection which is valid during a request request.
It will be put back automatically.
However, if you want to use a database connection in a generator for streaming,
you have to manually get and put back the connection.
The best way to deal with it in a generator, is to use a try/catch-block
and put the connection back when GeneratorExit is thrown or in the finally-block.
"""

def get_database_connection():
    """
    Takes a database connection from the pool, stores it in the request_store and returns it.
    It is automatically returned to the database pool after the requests is finished.
    """
    if "database_connection" not in request_store:
        request_store.database_connection = macpepdb_pool.getconn() # pylint: disable=assigning-non-slot
    return request_store.database_connection


def initialize_config(config_file: Optional[Union[Path, str]] = None, environment: Optional[Union[Environment, str]] = None):
    """
    Parses the config_file and environment if necessary and initializes the configuration.

    Parameters
    ----------
    config_file : Optional[Union[Path, str]], optional
        Path to config file, by default None
    environment : Optional[Union[Environment, str]], optional
        Environment as str, by default None
    """
    if config_file is not None and isinstance(config_file, str):
        config_file = Path(config_file)
    if environment is not None and isinstance(environment, str):
        environment = Environment.from_str(environment)

    Configuration.initialize(
        config_file,
        environment
    )

def get_app(config_file: Optional[Union[Path, str]] = None, environment: Optional[Union[Environment, str]] = None) -> Flask:
    """
    Build and set Flask app as global 

    Parameters
    ----------
    config_file : Optional[Union[Path, str]], optional
        Path to config file, by default None
    environment : Optional[Union[Environment, str]], optional
        Environment as str, by default None

    Returns
    -------
    Flask
        Flask app.
    """
    global app
    global macpepdb_pool

    initialize_config(config_file, environment)

    app = Flask('app')
    # Default Flask parameter
    app.config.update(
        ENV =  Configuration.environment().name,
        DEBUG =  Configuration.values()['debug'],
        SECRET_KEY = bytes( Configuration.values()['secret'], "ascii")
    )

    macpepdb_pool = ConnectionPool(1, Configuration.values()['macpepdb']['pool_size'], Configuration.values()['macpepdb']['url'])

    @app.before_request
    def track_request():
        if  Configuration.values()["matomo"]["enabled"]:
            track_thread = Thread(target=matomo_track_request, args=(
                request.headers.get("User-Agent", ""),
                request.remote_addr,
                request.headers.get("Referer", ""),
                request.headers.get("Accept-Language", ""),
                request.headers.get("Host", ""),
                request.full_path,
                request.query_string,
                request.url.startswith("https"),
                Configuration.values()["matomo"]["url"],
                Configuration.values()["matomo"]["site_id"],
                Configuration.values()["matomo"]["auth_token"], 
                app,
                Configuration.values()["debug"]
            ))
            track_thread.start()
            request_store.track_thread = track_thread 

    @app.teardown_appcontext
    def return_database_connection_to_pool(exception=None):
        """
        Returns the database connection to the pool
        """
        database_connection = request_store.pop("database_connection", None)
        if database_connection:
            macpepdb_pool.putconn(database_connection)

    @app.teardown_appcontext
    def wait_for_track_request(exception=None):
        track_thread = request_store.pop("track_thread", None)
        if track_thread:
            track_thread.join()

    @app.errorhandler(Exception)
    def handle_exception(e):
        response = None
        # pass through HTTP errors
        if isinstance(e, HTTPException):
            # Return JSON instead of HTML for HTTP errors
            # start with the correct headers and status code from the error
            response = e.get_response()
            # replace the body with JSON
            response.data = json.dumps({
                "errors": {
                    "general": e.description
                }
            })
            response.content_type = "application/json"
        else:
            response = app.response_class(
                response=json.dumps({
                    "errors": {
                        "general": str(e)
                    }
                }),
                status=500,
                mimetype='application/json'
            )
        if  Configuration.values()['debug']:
            app.logger.error(traceback.format_exc())
            response = add_allow_cors_headers(response)
        return response

    if  Configuration.values()['debug']:
        @app.after_request
        def add_cors_header_in_development_mode(response):
            return add_allow_cors_headers(response)


    # Import controllers.
    # Do not move this import to the top of the files. Each controller uses 'app' to build the routes.
    # Some controllers also import the connection pools.
    from macpepdb.web.routes import register_routes
    register_routes(app)

    bind: str = f"{Configuration.values()['interface']}:{ Configuration.values()['port']}"
    print(f"Start MaCPepDB webinterface in {Configuration.environment().name} mode on {bind}")
    return app

def start(config_file: Optional[Union[Path, str]] = None, environment: Optional[Union[Environment, str]] = None):
    """
    Starts Flask app with integrated server. Good for development.

    Parameters
    ----------
    config_file : Optional[Union[Path, str]], optional
        Path to config file, by default None
    environment : Optional[Union[Environment, str]], optional
        Environment as str, by default None
    """
    app: Flask = get_app(config_file, environment)
    app.run(
        Configuration.values()['interface'],
        Configuration.values()['port']
    )

def start_by_cli(cli_args):
    """
    Starts server with CLI arguments

    Parameters
    ----------
    cli_args : Any
        CLI arguments from argparse
    """
    config_path: Optional[Path] = Path(cli_args.config).absolute() if cli_args.config is not None else None
    environment: Optional[Environment] = Environment.from_str(cli_args.environment) if cli_args.environment is not None else None
    
    # If gunicorn is not set, just start the app
    # Otherwise build and sprint gunicorn parameters.
    if cli_args.gunicorn is None:
        start(
            config_path,
            environment
        )
    else:
        initialize_config(config_path, environment)

        config_path_str: str = f"\"{config_path}\"" if config_path is not None else "None"
        environment_str: str = f"\"{environment}\"" if environment is not None else "None"
        gunicorn_args: str = cli_args.gunicorn
        
        # Add bind from config if not defined in guncicorn
        if not "-b" in gunicorn_args and not "--bind" in gunicorn_args:
            gunicorn_args += f" -b {Configuration.values()['interface']}:{Configuration.values()['port']} "
        print(
            f"{gunicorn_args} 'macpepdb.web.server:get_app(config_file={config_path_str}, environment={environment_str})'"
        )


def add_cli_arguments(subparsers: argparse._SubParsersAction):
    """
    Defines the CLI parameters for the web server.

    Parameters
    ----------
    subparser : argparse._SubParsersAction
        Subparser of main CLI parser
    """
    parser = subparsers.add_parser("serve", help="Starts webserver")
    parser.add_argument("--environment", "-e", default=None, choices=[env.value for env in Environment], help="Environment to run")
    parser.add_argument("--config", "-c", required=False, default=None, help="Optional config file")
    parser.add_argument(
        "--gunicorn",
        nargs='?',
        type=str,
        default=None,                       # When not present 
        const="",                           # When present without value
        help=(
            "Can be used as flag or to pass arguments for Gunicorn webserver. "
            "If this is used (even without value), it returns a Gunicorn arguments string "
            "to start MaCPepDB using Gunicorn with the given config file and environment. "
            "Just pass the string to Gunicorn binary. "
            "If no Gunicorn bind option is added (-b|--bind) the interface and port of the config will be used."
        )
    )
    parser.set_defaults(func=start_by_cli)

