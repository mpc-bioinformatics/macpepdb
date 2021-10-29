MASS_CONVERT_FACTOR = 1000000000.0

HYDROGEN_MONO_MASS = 1.007825035


def to_int(mass: float) -> int:
    """
    Converts the given mass into the internal used integer representation.

    Parameters
    ----------
    mass : float
        Human readable mass

    Returns
    -------
    Integer mass
    """
    return int(mass * MASS_CONVERT_FACTOR)

def to_float(mass: int) -> float:
    """
    Converts the given mass into the a human readable float representation.

    Parameters
    ----------
    mass : int
        Internal integer representation

    Returns
    -------
    Human readable float
    """
    return mass / MASS_CONVERT_FACTOR

def thomson_to_dalton(thomson: float, charge: int) -> float:
    """
    Converts thomson (mass to charge) into Dalton.

    Parameters
    ----------
    thomson : float
        Mass to charge ratio.
    charge : int
        Charge

    Returns
    -------
    dalton
    """
    return thomson * charge - HYDROGEN_MONO_MASS * charge