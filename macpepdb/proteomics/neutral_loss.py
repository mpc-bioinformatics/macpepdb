# internal imports
from macpepdb.proteomics.mass.convert import to_int as mass_to_int

class NeutralLoss:
    def __init__(self, name: str, mono_mass: float, average_mass: float):
        self.name = name
        self.mono_mass = mass_to_int(mono_mass)
        self.average_mass = mass_to_int(average_mass)

    # Returns Neutral loss by name
    @classmethod
    def get_by_name(cls, name: str):
        try:
            return eval(name.upper())
        except NameError:
            return NONE

H2O = NeutralLoss("H2O", 18.010564700, 18.015)
NONE = NeutralLoss("None", 0.0, 0.0)