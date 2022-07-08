from macpepdb.proteomics import mass
from macpepdb.proteomics.mass.convert import to_float as mass_to_float
from macpepdb.models.peptide import Peptide
from macpepdb.models.protein import Protein

def peptide_to_dict(peptide: Peptide) -> dict:
    result =  {
        "mass": mass_to_float(peptide.mass),
        "sequence": peptide.sequence,
        "length": peptide.length,
        "number_of_missed_cleavages": peptide.number_of_missed_cleavages
    }
    if peptide.metadata:
        result["is_swiss_prot"] =  peptide.metadata.is_swiss_prot
        result["is_trembl"] =  peptide.metadata.is_trembl
        result["taxonomy_ids"] =  peptide.metadata.taxonomy_ids
        result["unique_taxonomy_ids"] =  peptide.metadata.unique_taxonomy_ids
        result["proteome_ids"] =  peptide.metadata.proteome_ids
    return result


def protein_to_dict(protein: Protein) -> dict:
    return {
        "accession": protein.accession,
        "entry_name": protein.entry_name,
        "name": protein.name,
        "sequence": protein.sequence,
        "taxonomy_id": protein.taxonomy_id,
        "proteome_id": protein.proteome_id,
        "is_reviewed": protein.is_reviewed
    }