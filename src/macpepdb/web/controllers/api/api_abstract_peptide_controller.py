from __future__ import annotations
from collections import defaultdict
from dataclasses import dataclass
from enum import Enum, unique
import json
import math
import traceback
from typing import ByteString, Callable, Iterator, List, Iterable, Optional, Any

from flask import jsonify, Response
from macpepdb.database.query_helpers.where_condition import WhereCondition
from macpepdb.proteomics.mass.convert import to_int as mass_to_int
from macpepdb.proteomics.modification import Modification
from macpepdb.proteomics.modification_collection import ModificationCollection
from macpepdb.models.modification_combination_list import ModificationCombinationList
from macpepdb.models.taxonomy import Taxonomy
from macpepdb.models.peptide import Peptide
from macpepdb.helpers.metadata_condition import MetadataCondition
from macpepdb.web.server import get_database_connection, macpepdb_pool, app
from macpepdb.web.controllers.application_controller import ApplicationController

@unique
class OutputFormat(Enum):
    json = "application/json"
    stream = "application/octet-stream"
    text = "text/plain"
    csv = "text/csv"
    fasta = "text/fasta"

    def __str__(self) -> str:
        return f"{self.value}"

    @classmethod
    def from_name(cls, name: str) -> OutputFormat:
        return cls[name.lower()]

    @classmethod
    def from_value(cls, value: str) -> OutputFormat:
        lower_value = value.lower()
        for format in cls:
            if format.value == lower_value:
                return format
        raise KeyError(f"f{value} not found")

class ApiAbstractPeptideController(ApplicationController):
    SUPPORTED_ORDER_COLUMNS = ['mass', 'length', 'sequence', 'number_of_missed_cleavages']
    SUPPORTED_ORDER_DIRECTIONS = ['asc', 'desc']
    FASTA_SEQUENCE_LEN = 60

    PEPTIDE_QUERY_DEFAULT_COLUMNS = ["mass", "sequence", "number_of_missed_cleavages", "length", "is_swiss_prot", "is_trembl", "taxonomy_ids", "unique_taxonomy_ids", "proteome_ids"]

    @staticmethod
    def _search(request, file_extension: str):
        errors = defaultdict(list)
        data = None
        if request.headers.get("Content-Type", "") == "application/json":
            data = request.get_json()
        elif request.headers.get("Content-Type", "") == "application/x-www-form-urlencoded":
            # For use with classical form-tag. The JSON-formatted search parameters should be provided in the form parameter "search_params"
            data = json.loads(request.form.get("search_params", "{}"))

        include_count = False
        if 'include_count' in data and isinstance(data['include_count'], bool):
            include_count = data['include_count']

        order_by = None
        if 'order_by' in data:
            if isinstance(data['order_by'], str) and data['order_by'] in ApiAbstractPeptideController.SUPPORTED_ORDER_COLUMNS:
                order_by = data['order_by']
            else:
                errors["order_by"].append(f"must be a string with one of following values: {', '.join(ApiAbstractPeptideController.SUPPORTED_ORDER_COLUMNS)}")

        if 'order_direction' in data:
            if not isinstance(data['order_direction'], str) or not data['order_direction'] in ApiAbstractPeptideController.SUPPORTED_ORDER_DIRECTIONS:
                errors["order_direction"].append(f"'order_direction' must be a string with one of following values: {', '.join(ApiAbstractPeptideController.SUPPORTED_ORDER_DIRECTIONS)}")
            
        include_metadata = False
        if "include_metadata" in data:
            if isinstance(data["include_metadata"], bool):
                include_metadata = data["include_metadata"]
            else:
                errors["include_metadata"].append("must be a boolean")

        output_style = None
        if file_extension is not None:
            try:
                output_style = OutputFormat.from_name(file_extension)
            except KeyError:
                pass
        else:
            try:
                output_style = OutputFormat.from_value(request.headers.get("accept", default=""))
            except KeyError:
                output_style = OutputFormat.json

        # validate int attributes
        for attribute in  ["lower_precursor_tolerance_ppm", "upper_precursor_tolerance_ppm", "variable_modification_maximum"]:
            if attribute in data:
                if isinstance(data[attribute], int):
                    if data[attribute] < 0:
                        errors[attribute].append("not greater or equals 0")
                else:
                    errors[attribute].append("not an integer")
            else:    
                errors[attribute].append("cannot be empty")

        modifications = []
        if "modifications" in data:
            if isinstance(data["modifications"], list):
                for idx, modification_attributes in enumerate(data["modifications"]):
                    if isinstance(modification_attributes, dict):
                        accession_and_name = "onlinemod:{}".format(idx)
                        try:
                            modification_attributes['accession'] = accession_and_name
                            modification_attributes['name'] = accession_and_name
                            modification_attributes['delta'] = mass_to_int(modification_attributes['delta'])
                            modifications.append(Modification.from_dict(modification_attributes))
                        except Exception as e:
                            errors[f"modifications[{idx}]"].append("is invalid")
                    else:
                        errors[f"modifications[{idx}]"].append("not a dictionary")
            else:
                errors["modifications"].append("modifications has to be of type list")
        
        try:
            modification_collection = ModificationCollection(modifications)
        except Exception as e:
            errors["modifications"].append(f"{e}")

        database_connection = get_database_connection()
        if not len(errors):
            if "precursor" in data:
                if isinstance(data["precursor"], float) or isinstance(data["precursor"], int):

                    modification_combination_list = ModificationCombinationList(
                        modification_collection, 
                        mass_to_int(data["precursor"]),
                        data["lower_precursor_tolerance_ppm"],
                        data["upper_precursor_tolerance_ppm"],
                        data["variable_modification_maximum"]
                    )

                    metadata_condition = MetadataCondition()

                    # List of metadata conditions
                    if "taxonomy_id" in data:
                        if isinstance(data["taxonomy_id"], int):
                            with database_connection.cursor() as database_cursor:
                                taxonomy = Taxonomy.select(
                                    database_cursor,
                                    (
                                        "id = %s",
                                        (data["taxonomy_id"], )
                                    )
                                )
                                if taxonomy is not None:
                                    metadata_condition.taxonomy_ids = [sub.id for sub in taxonomy.sub_species(database_cursor)]
                                else:
                                    errors["taxonomy_id"].append("not found")

                        else:
                            errors["taxonomy_id"].append("must be an integer")

                    if "proteome_id" in data:
                        if isinstance(data["proteome_id"], str):
                            metadata_condition.proteome_id = data["proteome_id"]
                        else:
                            errors["proteome_id"].append("must be a string")

                    if "is_reviewed" in data:
                        if isinstance(data["is_reviewed"], bool):
                            if data["is_reviewed"]:
                                metadata_condition.is_swiss_prot = True
                            else:
                                metadata_condition.is_trembl = True
                        else:
                            errors["is_reviewed"].append("must be a boolean")

                    # Sort by `order_by`
                    order_by_instruction = "sequence ASC"
                    if order_by and not output_style == OutputFormat.text:
                        order_by_instruction = f"{order_by} {data['order_direction']}, sequence ASC"

                    # Note about offset and limit: It is much faster to fetch data from server and discard rows below the offset and stop the fetching when the limit is reached, instead of applying LIMIT and OFFSET directly to the query.
                    # Even on high offsets, which discards a lot of rows, this approach is faster.
                    # Curl shows the diffences: curl -o foo.json --header "Content-Type: application/json" --request POST --data '{"include_count":true,"offset":0,"limit":50,"modifications":[{"amino_acid":"C","position":"anywhere","is_static":true,"delta":57.021464}],"lower_precursor_tolerance_ppm":5,"upper_precursor_tolerance_ppm":5,"variable_modification_maximum":0,"order":true,"precursor":859.49506802369}' http://localhost:3000/api/peptides/search
                    # Applying OFFSET and LIMIT to query: 49 - 52 seconds
                    # Discarding rows which are below the offset and stop the fetching early: a few hundred miliseconds (not printed by curl).
                    offset = 0
                    limit = math.inf
                    if "limit" in data:
                        if isinstance(data["limit"], int):
                            limit = data["limit"]
                        else:
                            errors["limit"].append("must be an integer")
                    if "offset" in data:
                        if isinstance(data["offset"], int):
                            offset = data["offset"]
                        else:
                            errors["offset"].append("must be an integer")

                else:
                    errors["precursor"] = ["must be an integer or float"]
            else:
                errors["precursor"] = ["cannot be missing"]

        if len(errors):
            return jsonify({
                "errors": errors
            }), 422

        include_metadata = include_metadata or metadata_condition.has_conditions()

        peptide_conversion = lambda _, __: (b"",)       # lambda to convert peptide to output type
        delimiter = b""                                 # delimiter between each converted peptide
        pre_peptide_content = b""                       # content before peptide
        post_peptide_content = lambda _, __, ___: b""   # content after peptides

        if output_style == OutputFormat.json:
            peptide_conversion = lambda _, peptide: peptide.to_json()
            delimiter = b","
            pre_peptide_content = b"{\"peptides\":["
            post_peptide_content = lambda _, __, ___: b"]}"
            if include_count:
                post_peptide_content = lambda database_cursor, where_condition, metadata_condition: f"],\"count\":{Peptide.count_on_stream(database_cursor, where_condition, metadata_condition)}}}".encode("utf-8")
        elif output_style == OutputFormat.stream:
            peptide_conversion = lambda _, peptide: peptide.to_json()
            delimiter = b"\n"
        elif output_style == OutputFormat.fasta:
            peptide_conversion = lambda peptide_idx, peptide: peptide.to_fasta_entry(f"P{peptide_idx}".encode())
            delimiter = b"\n"
        elif output_style == OutputFormat.csv:
            peptide_conversion = lambda _, peptide: peptide.to_csv_row()
            delimiter = b"\n"
            pre_peptide_content = (
                ",".join(Peptide.CSV_HEADER).encode("utf-8") if not include_metadata else \
                ",".join(Peptide.CSV_HEADER + Peptide.METADATA_CSV_HEADER).encode("utf-8")
            ) + b"\n"
        elif output_style == OutputFormat.text:
            peptide_conversion = lambda _, peptide: peptide.to_plain_text()
            delimiter = b"\n"

        return Response(
            ApiAbstractPeptideController.stream(
                peptide_conversion,
                delimiter,
                pre_peptide_content,
                post_peptide_content,
                modification_combination_list.to_where_condition(),
                order_by_instruction,
                offset,
                limit,
                include_metadata,
                metadata_condition
            ),
            content_type=f"{output_style}; charset=utf-8"
        )

    @staticmethod
    def stream(peptide_conversion: Callable[[int, Peptide], Iterator[ByteString]], delimiter: ByteString, pre_peptide_content: ByteString, post_peptide_content: Callable[[Any, WhereCondition, MetadataCondition], ByteString],
        where_condition: WhereCondition, order_by_instruction: str, offset: int, limit: int, include_metadata: bool,
        metadata_condition: MetadataCondition) -> Iterable[ByteString]:
        """
        Queries peptides and yields content for a stream response.

        Parameters
        ----------
        peptide_conversion : Callable[[int, Peptide], Iterator[ByteString]]
            Function with peptide index and peptide as input and yields the given peptide as bytes string for the response.
        delimiter : ByteString
            Delimiter between peptides
        pre_peptide_content : ByteString
            _description_
        post_peptide_content : Callable[[Any, WhereCondition], ByteString]
            Function which yields content after peptides, e.g. to add the count.
            Arguments are is the database cursor and the where condition
        where_condition : WhereCondition
            WhereCondition for SQL query
        order_by_instruction : str
            OrderBy instruction for SQL query
        offset : int
            Search offset
        limit : int
            Return limit
        include_metadata : bool
            If ture metadata will be included.
        metadata_condition : MetadataCondition
            Conditions for metadata. If not empty it sets metadata to true.

        Yields
        ------
        Iterator[Iterable[ByteString]]
            Response body
        """
        # In a generator the response is already returned and the app context is teared down. So we can not use a database connection from the actual request handling.
        # Get a new one from the pool and return it when the generator ist stopped (GeneratorExit is thrown).
        database_connection = macpepdb_pool.getconn()
        # Determine if metadata checks are necessary
        do_metadata_checks = metadata_condition.has_conditions()
        try:
            yield pre_peptide_content
            with database_connection.cursor(name="peptide_search") as database_cursor:
                # Counter for written peptides necessary of manual limit offset handling
                written_peptides = 0
                matching_peptides = 0
                for peptide in Peptide.select(database_cursor, where_condition, order_by=order_by_instruction, include_metadata=include_metadata, stream=True):
                    if do_metadata_checks and not metadata_condition.validate(peptide.metadata):
                        continue
                    matching_peptides += 1
                    if matching_peptides <= offset:
                        continue
                    if written_peptides > 0:
                        yield delimiter
                    yield from peptide_conversion(matching_peptides, peptide)
                    written_peptides += 1
                    # Break peptide cursor loop if limit is hit
                    if written_peptides == limit:
                        break
            with database_connection.cursor() as database_cursor:
                yield post_peptide_content(database_cursor, where_condition, metadata_condition if do_metadata_checks else None)
        except BaseException as e:
            app.logger.error(f"steam throws err => {e}\n\ntraceback => {traceback.format_exc()}")
            raise e
        finally:
            macpepdb_pool.putconn(database_connection)
