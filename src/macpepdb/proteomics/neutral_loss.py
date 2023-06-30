# internal imports
from macpepdb.proteomics.mass.convert import to_int as mass_to_int

class NeutralLoss:
    """
    Defines a neutral loss

    Parameters
    ----------
    name : str
        Name of the neutral loss
    mono_mass : int
        Mono mass as human readable float
    average_mass : int
        Average as human readable float
    """

    def __init__(self, name: str, mono_mass: float, average_mass: float):
        self.name = name
        self.mono_mass = mass_to_int(mono_mass)
        self.average_mass = mass_to_int(average_mass)

    @classmethod
    def get_by_name(cls, name: str):
        """
        Returns a neutral loss.

        Parameters
        ----------
        name : str
            Name of the NeutralLoss

        Returns
        -------
        NeutralLoss
        """
        try:
            return eval(name.upper())
        except NameError:
            return NONE

H2O = NeutralLoss("H2O", 18.010564700, 18.015)
NONE = NeutralLoss("None", 0.0, 0.0)