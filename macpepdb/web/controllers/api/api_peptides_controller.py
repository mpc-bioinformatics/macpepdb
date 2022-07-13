# std imports
from collections import defaultdict
import itertools
from typing import ClassVar, Iterator

# 3rd party imports
from flask import request, jsonify, Response
from macpepdb.database.query_helpers.where_condition import WhereCondition
from macpepdb.models.protein import Protein
from macpepdb.models.protein_peptide_association import ProteinPeptideAssociation
from macpepdb.proteomics.mass.convert import to_float as mass_to_float
from macpepdb.proteomics.enzymes import get_digestion_enzyme_by_name

# internal imports
from macpepdb.web.server import app, get_database_connection, macpepdb_pool
from macpepdb.web.models.convert import peptide_to_dict, protein_to_dict
from macpepdb.web.controllers.api.api_abstract_peptide_controller import ApiAbstractPeptideController
from macpepdb.web.controllers.api.api_digestion_controller import ApiDigestionController
from macpepdb.web.models.peptide import Peptide


class ApiPeptidesController(ApiAbstractPeptideController):
    PEPTIDE_LOOKUP_CHUNKS: ClassVar[int] = 500

    @staticmethod
    def search(file_extension: str = None):
        return ApiAbstractPeptideController._search(request, file_extension)

    @staticmethod
    def show(sequence: str):
        is_reviewed = request.args.get("is_reviewed", None)
        if is_reviewed is not None:
            is_reviewed = bool(is_reviewed)
        sequence = sequence.upper()
        database_connection = get_database_connection()
        with database_connection.cursor() as database_cursor:
            peptide = Peptide(sequence, 0, None)
            peptide = Peptide.select(
                database_cursor,
                WhereCondition(
                    ["partition = %s", "AND", "mass = %s", "AND", "sequence = %s"],
                    [peptide.partition, peptide.mass, peptide.sequence]
                ),
                include_metadata=True
            )
            if peptide is None:
                return jsonify({
                    "errors": {
                        "sequence": ["not found"]
                    }
                }), 404
            # Return peptide if is_reviewed is not requested (None),
            # or is_reviewed is requested and True and metadata is_swiss_prot is also True
            # or is_reviewed is requested and False and metadata is_trembl is True
            if is_reviewed is None \
                or is_reviewed and peptide.metadata.is_swiss_prot \
                or not is_reviewed and peptide.metadata.peptide.metadata.is_trembl:
                return Response(
                    peptide.to_json(),
                    content_type="application/json"
                )

    @staticmethod
    def proteins(sequence: str):
        peptide = Peptide(sequence.upper(), 0)
        database_connection = get_database_connection()
        with database_connection.cursor() as database_cursor:
            proteins = Protein.select(
                database_cursor,
                WhereCondition(
                    [
                        f"accession = ANY(SELECT protein_accession FROM {ProteinPeptideAssociation.TABLE_NAME} as ppa WHERE ppa.partition = %s AND ppa.peptide_mass = %s AND ppa.peptide_sequence = %s)"
                    ],
                    [
                        peptide.partition,
                        peptide.mass,
                        peptide.sequence
                    ]
                ),
                True
            )

            reviewed_proteins_rows = []
            unreviewed_proteins_rows = []

            for protein in proteins:
                if protein.is_reviewed:
                    reviewed_proteins_rows.append(
                        protein_to_dict(protein)
                    )
                else:
                    unreviewed_proteins_rows.append(
                        protein_to_dict(protein)
                    )

            return jsonify({
                "reviewed_proteins": reviewed_proteins_rows,
                "unreviewed_proteins": unreviewed_proteins_rows
            })

    @staticmethod
    def sequence_mass(sequence):
        peptide = Peptide(sequence, 0)

        return jsonify({
            'mass': mass_to_float(peptide.mass)
        })


    @staticmethod
    def digest():
        """
        Digest a given peptide/sequence, search the resulting peptides in the database and return matching and not matching peptides in separate array.
        """
        data = request.get_json()
        errors = ApiDigestionController.check_digestion_parameters(data)

        if not "sequence" in data:
            errors["sequence"].append("cannot be empty")

        digestion_peptides = []
        database_peptides = []
        if len(errors) == 0:
            EnzymeClass = get_digestion_enzyme_by_name("trypsin")
            enzyme = EnzymeClass(data["maximum_number_of_missed_cleavages"], data["minimum_peptide_length"], data["maximum_peptide_length"])
            digestion_peptides = enzyme.digest(Protein("TMP", [], "TMP", "TMP", data["sequence"], [], [], False, 0))

            where_clause = " OR ".join(["mass = %s AND sequence = %s"] * len(digestion_peptides))
            where_values = list(itertools.chain.from_iterable([[peptide.mass, peptide.sequence] for peptide in digestion_peptides]))

            if "do_database_search" in data and isinstance(data["do_database_search"], bool) and data["do_database_search"]:
                database_connection = get_database_connection()           
                with database_connection.cursor() as database_cursor:
                    database_peptides = Peptide.select(database_cursor, (where_clause, where_values), fetchall=True)
                database_peptides.sort(key = lambda peptide: peptide.mass)
                digestion_peptides = [peptide for peptide in digestion_peptides if peptide not in database_peptides]


            digestion_peptides.sort(key = lambda peptide: peptide.mass)

        if len(errors) == 0:
            return jsonify({
                "database": [peptide_to_dict(peptide) for peptide in database_peptides],
                "digestion": [peptide_to_dict(peptide) for peptide in digestion_peptides],
                "count": len(database_peptides) +  len(digestion_peptides)
            })
        else:
            return jsonify({
                "errors": errors
            }), 422


    @staticmethod
    def sequence_lookup():
        """
        Check if the incoming peptide sequences exists in MaCPepDB

        Returns
        -------
        Response
            Flask response
        """
        data = request.json
        errors = defaultdict(list)

        if not "sequences" in data:
            errors["sequences"].append("cannot be empty")
        elif not isinstance(data["sequences"], list):
            errors["sequences"].append("must be a list")

        peptides = [Peptide(sequence, 0) for sequence in data["sequences"]]
        # Sort peptides by partition
        partitions = defaultdict(list)
        for peptide in peptides:
            partitions[peptide.partition].append(peptide)

        if len(errors) > 0:
            return jsonify({
                "errors": errors
            }), 422

        def generate_text_stream() -> Iterator[str]:
            database_connection = macpepdb_pool.getconn()
            try:
                with database_connection.cursor() as database_cursor:
                    database_cursor.itersize = ApiPeptidesController.PEPTIDE_LOOKUP_CHUNKS
                    # Request each partition in a separate query
                    for partition, part_peptides in partitions.items():
                        # Chunk the peptides in 500 
                        chunk_start = 0
                        while chunk_start < len(part_peptides):
                            chunk = part_peptides[chunk_start:(chunk_start + ApiPeptidesController.PEPTIDE_LOOKUP_CHUNKS)]
                            for peptide_idx, peptide in enumerate(Peptide.select(
                                database_cursor,
                                WhereCondition(
                                    [
                                        "partition = %s",
                                        "AND"
                                        "(mass, sequence) IN (" + ",".join(["(%s, %s)"] * len(chunk)) + ")"
                                    ],
                                    [partition] + list(itertools.chain(*[(peptide.mass, peptide.sequence) for peptide in chunk])),
                                ),
                                stream=True
                            )):
                                if peptide_idx > 0:
                                    yield "\n"
                                yield peptide.sequence
                            chunk_start += ApiPeptidesController.PEPTIDE_LOOKUP_CHUNKS

            finally:
                macpepdb_pool.putconn(database_connection)

        return Response(generate_text_stream(), content_type="text/plain")
