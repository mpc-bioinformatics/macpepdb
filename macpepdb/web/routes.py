# 3rd party imports
from flask import Flask

# internal imports
from macpepdb.web.controllers.api.api_dashboard_controller import ApiDashboardController
from macpepdb.web.controllers.api.api_digestion_controller import ApiDigestionController
from macpepdb.web.controllers.api.api_documents_controller import ApiDocumentsController
from macpepdb.web.controllers.api.api_modifications_controller import ApiModificationsController
from macpepdb.web.controllers.api.api_peptides_controller import ApiPeptidesController
from macpepdb.web.controllers.api.api_proteins_controller import ApiProteinsController
from macpepdb.web.controllers.api.api_taxonomies_controller import ApiTaxonomiesController



def register_routes(app: Flask):
    """
    Register routes to controller

    Parameters
    ----------
    app : Flask
        Flask app
    """
    # If multiple function assigned to view_func of app.add_url_rule has the same name, regardless of the controller name,
    # you have to give them a distinct endpioint name.

    # Dashboard controller
    app.add_url_rule("/api/dashboard/status", view_func=ApiDashboardController.status)
    app.add_url_rule("/api/dashboard/maintenance", view_func=ApiDashboardController.maintenance)
    
    # Modifications controller
    app.add_url_rule("/api/modifications/positions", view_func=ApiModificationsController.modification_positions)

    # Proteins controller
    app.add_url_rule("/api/proteins/<string:accession>", view_func=ApiProteinsController.show, endpoint="protein_show")
    app.add_url_rule("/api/proteins/<string:accession>/peptides", view_func=ApiProteinsController.peptides)
    app.add_url_rule("/api/proteins/digest", view_func=ApiProteinsController.digest, methods=["POST"], endpoint="protein_digest")
    app.add_url_rule("/api/proteins/amino-acids", view_func=ApiProteinsController.amino_acids)

    # Peptides controller
    app.add_url_rule("/api/peptides/search", view_func=ApiPeptidesController.search, methods=["POST"], endpoint="peptide_search")
    app.add_url_rule("/api/peptides/search.<string:file_extension>", view_func=ApiPeptidesController.search, methods=["POST"], endpoint="peptide_search_file_ext")
    app.add_url_rule("/api/peptides/<string:sequence>", view_func=ApiPeptidesController.show, methods=["GET"], endpoint="peptide_show")
    app.add_url_rule("/api/peptides/<string:sequence>/proteins", view_func=ApiPeptidesController.proteins, methods=["GET"])
    app.add_url_rule("/api/peptides/mass/<string:sequence>", view_func=ApiPeptidesController.sequence_mass, methods=["GET"])
    app.add_url_rule("/api/peptides/digest", view_func=ApiPeptidesController.digest, methods=["POST"], endpoint="peptide_digest")
    app.add_url_rule("/api/peptides/lookup", view_func=ApiPeptidesController.sequence_lookup, methods=["POST"])

    # Taxonomy controller
    app.add_url_rule("/api/taxonomies/search", view_func=ApiTaxonomiesController.search, methods=["POST"], endpoint="taxonomy_search")
    app.add_url_rule("/api/taxonomies/<int:id>", view_func=ApiTaxonomiesController.show, endpoint="taxonomy_show")
    app.add_url_rule("/api/taxonomies/by/ids", view_func=ApiTaxonomiesController.by_ids, methods=["POST"])

    # Documents controller
    app.add_url_rule("/api/documents/20220314-macpepdb__increasing-performance.pdf", view_func=ApiDocumentsController.increasing_performance)
    app.add_url_rule("/api/documents/20220331-enhancement_of_macpepdb.pdf", view_func=ApiDocumentsController.enhancement)
    
