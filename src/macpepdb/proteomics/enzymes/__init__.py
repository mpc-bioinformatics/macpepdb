# std imports
from typing import List, Type

# internal imports
from macpepdb.proteomics.enzymes.digest_enzyme import DigestEnzyme
from macpepdb.proteomics.enzymes.trypsin import Trypsin


def get_digestion_enzyme_by_name(name: str) -> Type[DigestEnzyme]:
    """
    Returns a enzyme by name.

    Argument
    ========

    Returns
    -------
    DigestEnzym
    """
    # to prevent cyclic imports, import enzyms here not at the top
    if name.lower() == Trypsin.NAME.lower():
        return Trypsin
    # elif name.lower() == OtherEnzym.NAME.lower()
    raise NameError("Unknown enzym '{}'".format(name))

def get_known_digestion_enzymes() -> List[str]:
    """
    Returns a list of DigestEnzym names.

    Returns
    -------
    List of names.
    """
    return [
        Trypsin.NAME.lower()
    ]
