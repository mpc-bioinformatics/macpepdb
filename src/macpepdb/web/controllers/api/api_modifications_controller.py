from flask import jsonify

from macpepdb.proteomics.modification import ModificationPosition

from macpepdb.web.server import app
from macpepdb.web.controllers.application_controller import ApplicationController

class ApiModificationsController(ApplicationController):
    @staticmethod
    def modification_positions():
        return jsonify({
            "modification_positions": sorted([str(position) for position in ModificationPosition])
        })