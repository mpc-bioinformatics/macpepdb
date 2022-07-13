from typing import Iterator
from flask import request, jsonify, Response

from macpepdb.database.query_helpers.where_condition import WhereCondition
from macpepdb.proteomics.amino_acid import AminoAcid
from macpepdb.models.protein import Protein
from macpepdb.models.peptide import Peptide
from macpepdb.models.protein_peptide_association import ProteinPeptideAssociation
from macpepdb.web.server import app, get_database_connection
from macpepdb.web.controllers.application_controller import ApplicationController
from macpepdb.web.controllers.api.api_digestion_controller import ApiDigestionController

class ApiProteinsController(ApplicationController):
    @staticmethod
    def show(accession: str):
        accession = accession.upper()

        database_connection = get_database_connection()

        with database_connection.cursor() as database_cursor:
            protein = Protein.select(
                database_cursor, 
                WhereCondition(
                    ["accession = %s"], 
                    [accession]
                ),
                False
            )

            if protein:
                return Response(
                    protein.to_json(),
                    content_type="application/json"
                )
            else:               
                return jsonify({
                    "errors": {
                        "accession": ["not found"]
                    }
                }), 404

    @staticmethod
    def peptides(accession: str):

        database_connection = get_database_connection()
        with database_connection.cursor() as database_cursor:
                peptides = Peptide.select(
                    database_cursor,
                    WhereCondition(
                        [f"(peps.partition, peps.mass, peps.sequence) IN (SELECT partition, peptide_mass, peptide_sequence FROM {ProteinPeptideAssociation.TABLE_NAME} as ppa WHERE ppa.protein_accession = %s)"],
                        [accession]
                    ),
                    fetchall=True,
                    include_metadata= True
                )
                peptides.sort(key = lambda peptide: peptide.mass)

                def json_stream() -> Iterator[bytes]:
                    yield b"{\"peptides\": ["
                    for peptide_idx, peptide in enumerate(peptides):
                        if peptide_idx > 0:
                            yield b","
                        yield from peptide.to_json()
                    yield b"]}"
                    
                return Response(
                    json_stream(),
                    content_type="application/json"
                )



    @staticmethod
    def digest():
        """
        Digests the seqeunce of the given protein.
        """
        data = request.get_json()
        errors = ApiDigestionController.check_digestion_parameters(data)

        if not "accession" in data:
            errors["accession"].append("cannot be empty")

        peptides = []
        if len(errors) == 0:
            database_connection = get_database_connection()
            with database_connection.cursor() as database_cursor:
                protein = Protein.select(
                    database_cursor, 
                    WhereCondition(
                        ["accession = %s"], 
                        [data["accession"]]
                    ),
                    False
                )
                if protein:
                    peptides = list(filter(
                        lambda peptide: peptide.number_of_missed_cleavages <= data["maximum_number_of_missed_cleavages"] \
                            and data["minimum_peptide_length"] <= peptide.length <= data["maximum_peptide_length"],
                        protein.peptides(database_cursor)
                    ))
                    peptides.sort(key = lambda peptide: peptide.mass)
                else:
                    errors["accession"].append("not found")

        if len(errors) == 0:
            def json_stream():
                yield b"{\"peptides\": ["
                for peptide_idx, peptide in enumerate(peptides):
                    if peptide_idx > 0:
                        yield b","
                    yield from peptide.to_json()
                yield f"], \"count\": {len(peptides)}}}".encode("utf-8 ")
            return Response(
                json_stream(),
                content_type="application/json"
            )
        else:
            return jsonify({
                "errors": errors
            }), 422

    @staticmethod
    def amino_acids():
        return jsonify({
            "amino_acids": [{
                "one_letter_code": amino_acid.one_letter_code,
                "name": amino_acid.name
            } for amino_acid in AminoAcid.all()]
        })