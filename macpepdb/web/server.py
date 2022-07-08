# std imports
import argparse
import json
import os
from pathlib import Path
import traceback
from threading import Thread
from typing import Optional

# 3rd part imports
import bjoern
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


def start(config_file: Optional[Path] = None, environment: Optional[Environment] = None):
    global app
    global macpepdb_pool

    Configuration.initialize(
        config_file,
        environment
    )

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

    print(f"Start MaCPepDB webinterface in {Configuration.environment().name} mode on { Configuration.values()['interface']}:{ Configuration.values()['port']}")

    if Configuration.environment() == Environment.production or os.getenv("USE_BJOERN", "false") == "true":
        bjoern.run(app,  Configuration.values()['interface'],  Configuration.values()['port'])
    else:
        app.run( Configuration.values()['interface'],  Configuration.values()['port'])

def start_by_cli(cli_args):
    start(
        Path(cli_args.config) if cli_args.config is not None else None,
        Environment.from_str(cli_args.environment) if cli_args.environment is not None else None
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
    parser.set_defaults(func=start_by_cli)

