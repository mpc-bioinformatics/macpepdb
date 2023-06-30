from collections import defaultdict
from sqlalchemy.orm import sessionmaker, selectinload, exc as sqlalchemy_exceptions
from flask import json, request, jsonify, url_for

from macpepdb.models.taxonomy import Taxonomy
from macpepdb.models.taxonomy_merge import TaxonomyMerge

from macpepdb.web.server import app, get_database_connection
from macpepdb.web.controllers.application_controller import ApplicationController

class ApiTaxonomiesController(ApplicationController):
    @staticmethod
    def search():
        data = request.get_json()
        errors = defaultdict(list)

        if not "query" in data:
            errors["query"].append("must be present")
            if not isinstance(data["query"], str) or not isinstance(data["query"], int):
                errors["query"].append("must be a string or an integer")
        
        if not len(errors):
            condition = None
            if isinstance(data["query"], str):
                query = data["query"].replace("*", "%")
                condition = ("name LIKE %s", [query])
            else: 
                condition = ("id = %s", [data["query"]])

            database_connection = get_database_connection()
            with database_connection.cursor() as database_cursor:
                taxonomies = Taxonomy.select(database_cursor, condition, fetchall=True)

                if isinstance(data["query"], int) and not len(taxonomies):
                    taxonomy_merge = TaxonomyMerge.select(database_cursor, ("source_id = %s", [data["query"]]))
                    if taxonomy_merge:
                        taxonomies = Taxonomy.select(database_cursor, ("id = %s", [taxonomy_merge.target_id]), fetchall=True)

                response = [{
                    "id": taxonomy.id,
                    "name": taxonomy.name
                } for taxonomy in taxonomies]

            return jsonify(response)
        else:
            return jsonify({
                "errors": errors
            }), 422

    @staticmethod
    def show(id):
        database_connection = get_database_connection()
        with database_connection.cursor() as database_cursor:
            taxonomy = Taxonomy.select(database_cursor, ("id = %s", [id]))

            if not taxonomy:
                taxonomy_merge = TaxonomyMerge.select(database_cursor, ("source_id = %s", [id]))
                if taxonomy_merge:
                    taxonomy = Taxonomy.select(database_cursor, ("id = %s", [taxonomy_merge.target_id]))


            response = None
            if taxonomy:
                response = {
                    "id": taxonomy.id,
                    "name": taxonomy.name,
                    "parent": taxonomy.parent_id,
                    "rank": taxonomy.rank,
                }

        if response:
            return jsonify(response)
        else:
            return jsonify(["not found"]), 422

    @staticmethod
    def sub_species(id):
        database_connection = get_database_connection()
        with database_connection.cursor() as database_cursor:
            taxonomy = Taxonomy.select(database_cursor, ("id = %s", [id]))

            if taxonomy is None:
                taxonomy_merge = TaxonomyMerge.select(database_cursor, ("source_id = %s", [id]))
                if taxonomy_merge:
                    taxonomy = Taxonomy.select(database_cursor, ("id = %s", [taxonomy_merge.target_id]))

            if taxonomy is not None:
                return jsonify({
                    "sub_species": [
                        {
                            "id": sub_taxonomy.id,
                            "name": sub_taxonomy.name,
                            "parent": sub_taxonomy.parent_id,
                            "rank": sub_taxonomy.rank,
                        } for sub_taxonomy in taxonomy.sub_species(database_cursor)
                    ]
                })
        return jsonify(["not found"]), 422

    @staticmethod
    def by_ids():
        data = request.get_json()
        errors = defaultdict(list)

        if not "ids" in data:
            errors["ids"].append("cannot be missing")
        
        if "ids" in data and isinstance(data["ids"], list):
            for idx, id in enumerate(data["ids"]):
                if not isinstance(id, int):
                    errors[f"ids[{idx}]"].append(f"must be an integer")
        else:
            errors["ids"].append("must be a array")

        if len(errors):
            return jsonify({
                "errors": errors
            })

        database_connection = get_database_connection()
        with database_connection.cursor() as database_cursor:
            return jsonify({
                "taxonomies": [{
                    "id": taxonomy.id,
                    "name": taxonomy.name
                } for taxonomy in Taxonomy.select(database_cursor, ("id = ANY(%s)", [data["ids"]]), True)]
            })