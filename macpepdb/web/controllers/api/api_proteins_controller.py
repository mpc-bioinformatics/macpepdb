from flask import request, jsonify

from macpepdb.database.query_helpers.where_condition import WhereCondition
from macpepdb.proteomics.amino_acid import AminoAcid
from macpepdb.proteomics.mass.convert import to_float as mass_to_float
from macpepdb.models.protein import Protein
from macpepdb.models.peptide import Peptide
from macpepdb.models.protein_peptide_association import ProteinPeptideAssociation
from macpepdb.models.taxonomy import Taxonomy

from macpepdb.web.server import app, get_database_connection
from macpepdb.web.models.convert import peptide_to_dict, protein_to_dict
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
                response_body = protein_to_dict(protein)
                database_cursor.execute(f"SELECT name FROM {Taxonomy.TABLE_NAME} WHERE id = %s;", (protein.taxonomy_id,))
                taxonomy_name_row = database_cursor.fetchone()
                response_body["taxonomy_name"] = taxonomy_name_row[0] if taxonomy_name_row else None

                return jsonify(response_body)
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

                return jsonify({
                    "peptides": [{
                        "mass": mass_to_float(peptide.mass),
                        "sequence": peptide.sequence,
                        "number_of_missed_cleavages": peptide.number_of_missed_cleavages,
                        "is_swiss_prot": peptide.metadata.is_swiss_prot,
                        "is_trembl": peptide.metadata.is_trembl,
                        "taxonomy_ids": peptide.metadata.taxonomy_ids,
                        "unique_taxonomy_ids": peptide.metadata.unique_taxonomy_ids,
                        "proteome_ids": peptide.metadata.proteome_ids
                    } for peptide in peptides]
                })


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
                        lambda peptide: peptide.number_of_missed_cleavages <= data["maximum_number_of_missed_cleavages"] and data["minimum_peptide_length"] <= peptide.length <= data["maximum_peptide_length"],
                        protein.peptides(database_cursor)
                    ))
                    peptides.sort(key = lambda peptide: peptide.mass)
                else:
                    errors["accession"].append("not found")

        if len(errors) == 0:
            return jsonify({
                "peptides": [peptide_to_dict(peptide) for peptide in peptides],
                "count": len(peptides)
            })
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