class PrecursorRange:
    """
    Defines a precursor or mass range.

    Parameters
    ----------
    precursor : int
        Precursor / mass
    lower_tolerance_ppm : int
        Lower tolerance (ppm)
    upper_tolerance_ppm : int
        Upper tolerance (ppm)
    """

    def __init__(self, precursor: int, lower_tolerance_ppm: int, upper_tolerance_ppm: int):
        self.__precursor = precursor
        self.__lower_tolerance_ppm = lower_tolerance_ppm
        self.__upper_tolerance_ppm = upper_tolerance_ppm
        # Calculates the absolute precursor difference for the given tolerances (ppm) and add/subtract it from precursor
        self.__lower_limit = int(self.__precursor - int(self.__precursor / 1000000.0 * lower_tolerance_ppm))
        self.__upper_limit = int(self.__precursor + int(self.__precursor / 1000000.0 * upper_tolerance_ppm))

    @property
    def precursor(self):
        """
        Returns
        -------
        Precursor
        """
        return self.__precursor

    @property
    def lower_limit(self):
        """
        Returns
        -------
        Lower limit of the range
        """
        return self.__lower_limit

    @property
    def upper_limit(self):
        """
        Returns
        -------
        Upper limit of the range.
        """
        return self.__upper_limit

    def __contains__(self, value: int):
        """
        Implements 'in'-operator, e.g. value in self

        Parameters
        ----------
        value : int
            Mass in integer representation.

        Returns
        -------
        True if the value is between lower and upper limit (both including).
        """
        return self.__lower_limit <= value and value <= self.__upper_limit

    def __lt__(self, value: int):
        """
        Implements less-then operator, e.g. self < value

        Parameters
        ----------
        value : int
            Mass in integer representation.

        Returns
        -------
        True if the range is smaller then the given value
        """
        return self.__upper_limit < value

    def __le__(self, value: int):
        """
        Implements less-or-equals operator, e.g. self <= value
        
        Parameters
        ----------
        value : int
            Mass in integer representation.

        Returns
        -------
        True if the range is smaller or equals then the given value
        """
        return self.__lower_limit <= value

    def __gt__(self, value: int):
        """
        Implements greater-than operator, e.g. self > value

        Parameters
        ----------
        value : int
            Mass in integer representation.

        Returns
        -------
        True if the range is greater then the given value
        """
        return self.__lower_limit > value

    def __ge__(self, value: int):
        """
        Implements greater-or-equals operator, e.g. self >= value

        Parameters
        ----------
        value : int
            Mass in integer representation.

        Returns
        -------
        True if the range is greater or equals then the given value
        """
        return self.__upper_limit >= value

    def __eq__(self, value: int):
        """
        Implements equals, e.g. self == value

        Parameters
        ----------
        value : int
            Mass in integer representation.

        Returns
        -------
        True if the value is within the range (look: __contains__)
        """
        return self.__contains__(value)

    @property
    def highest_tolerance(self):
        """
        Returns
        -------
        Highest value of lower and upper tolerance
        """
        return max(self.__lower_tolerance_ppm, self.__upper_tolerance_ppm)
