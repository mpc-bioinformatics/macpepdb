# std imports
from pathlib import Path
from typing import ClassVar

# 3rd party imports
from flask import send_file

class ApiDocumentsController:
    DOCUMENTS_FOLDER: ClassVar[Path] = Path(__file__).parent.parent.parent.parent
    INCREASING_PERFORMANCE_POSTER_PATH: ClassVar[Path] = DOCUMENTS_FOLDER.joinpath("documents/20220314-macpepdb__increasing-performance.pdf")
    INCREASING_PERFORMANCE_POSTER_PATH: ClassVar[Path] = DOCUMENTS_FOLDER.joinpath("documents/20220331-enhancement_of_macpepdb.pdf")

    @staticmethod
    def increasing_performance():
        return send_file(
            ApiDocumentsController.INCREASING_PERFORMANCE_POSTER_PATH,
            as_attachment=True,
            attachment_filename=ApiDocumentsController.INCREASING_PERFORMANCE_POSTER_PATH.name
        )

    @staticmethod
    def enhancement():
        return send_file(
            ApiDocumentsController.INCREASING_PERFORMANCE_POSTER_PATH,
            as_attachment=True,
            attachment_filename=ApiDocumentsController.INCREASING_PERFORMANCE_POSTER_PATH.name
        )