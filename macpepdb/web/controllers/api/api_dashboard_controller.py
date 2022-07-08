# 3rd party imports
from flask import jsonify
from macpepdb.models.maintenance_information import MaintenanceInformation

# internal imports
from macpepdb.web.server import app, get_database_connection
from macpepdb.web.controllers.application_controller import ApplicationController

class ApiDashboardController(ApplicationController):
    @staticmethod
    def citus():
        database_connection = get_database_connection()
        with database_connection.cursor() as database_cursor:
            database_cursor.execute("SELECT count(*) from pg_dist_node");
            number_of_nodes = database_cursor.fetchone()[0]

            database_cursor.execute("SELECT * FROM get_rebalance_progress()");
            rebalance_job_rows = database_cursor.fetchall()
            finished_rebalance_job_rows = list(filter(lambda row: row[8] == 2, rebalance_job_rows))
            running_rebalance_job_rows = list(filter(lambda row: row[8] == 1, rebalance_job_rows))

            return jsonify({
                "number_of_nodes": number_of_nodes,
                "number_of_rebalance_jobs": len(rebalance_job_rows),
                "number_of_finished_rebalance_jobs": len(finished_rebalance_job_rows),
                "number_of_running_rebalance_jobs": len(running_rebalance_job_rows)
            })

    @staticmethod
    def maintenance():
        database_connection = get_database_connection()
        with database_connection.cursor() as database_cursor:
            comment = MaintenanceInformation.select(database_cursor, MaintenanceInformation.COMMENT_KEY)
            database_status = MaintenanceInformation.select(database_cursor, MaintenanceInformation.DATABASE_STATUS_KEY)
            digestion_parameter = MaintenanceInformation.select(database_cursor, MaintenanceInformation.DIGESTION_PARAMTERS_KEY)

            return jsonify({
                "comment": comment.values[comment] if comment is not None else None,
                "database_status": database_status.values,
                "digestion_parameters": digestion_parameter.values
            })
