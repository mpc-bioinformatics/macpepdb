from collections import defaultdict
from typing import DefaultDict, List


class ApiDigestionController:
    """
    Controller for amino acid sequence digestion.
    """
    @staticmethod
    def check_digestion_parameters(data: dict) -> DefaultDict[str, List[str]]:
        """
        Checks if the given data dictionary contains the necessary parameters for a digestion.

        Arguments
        =========
        data : dict
            Request data as dictionary

        Return
        ======
        Dictionary with parameter name as key and a list of errors as value.
        """
        errors = defaultdict(list)
        for attribute in  ["maximum_number_of_missed_cleavages", "minimum_peptide_length", "maximum_peptide_length"]:
            if attribute in data:
                if not isinstance(data[attribute], int):
                    errors[attribute].append("must be an integer")
            else:
                errors[attribute].append("cannot be missing")

        if len(errors) == 0:
            if data["maximum_number_of_missed_cleavages"] < 0:
                errors["maximum_number_of_missed_cleavages"].append("must be greater or equals 0")

            for attribute in  ["minimum_peptide_length", "maximum_peptide_length"]:
                if data[attribute] < 1:
                    errors["minimum_peptide_length"].append("must be greater or equals 1")

            if data["maximum_peptide_length"] < data["minimum_peptide_length"]:
                minimum_peptide_length = data["minimum_peptide_length"]
                errors["minimum_peptide_length"].append(f"must be greater or equals {minimum_peptide_length} (minimum peptide length)")

        return errors
