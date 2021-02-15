MASS_CONVERT_FACTOR = 1000000000.0

HYDROGEN_MONO_MASS = 1.007825035


def to_int(mass: float) -> int:
    return int(mass * MASS_CONVERT_FACTOR)

def to_float(mass: int) -> float:
    return mass / MASS_CONVERT_FACTOR

def thomson_to_dalton(thomson: float, charge: int) -> float:
    return thomson * charge - HYDROGEN_MONO_MASS * charge