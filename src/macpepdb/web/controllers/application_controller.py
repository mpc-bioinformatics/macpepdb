import datetime
import io

from flask import render_template
import matplotlib.pyplot as plt

from macpepdb.models.maintenance_information import MaintenanceInformation
from macpepdb.tasks.statistics import Statistics

from macpepdb.web.server import app, get_database_connection

class ApplicationController:
    pass