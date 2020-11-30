class PrecursorRange:
    def __init__(self, precursor: int, lower_tolerance_ppm: int, upper_tolerance_ppm: int):
        self.__precursor = precursor
        self.__lower_tolerance_ppm = lower_tolerance_ppm
        self.__upper_tolerance_ppm = upper_tolerance_ppm
        # Calculates the absolute precursor difference for the given tolerances (ppm) and add/subtract it from precursor
        self.__lower_limit = int(self.__precursor - int(self.__precursor / 1000000.0 * lower_tolerance_ppm))
        self.__upper_limit = int(self.__precursor + int(self.__precursor / 1000000.0 * upper_tolerance_ppm))

    @property
    def precursor(self):
        return self.__precursor

    @property
    def lower_limit(self):
        return self.__lower_limit

    @property
    def upper_limit(self):
        return self.__upper_limit

    # Implements 'in'-operator, e.g. value in self
    def __contains__(self, value: int):
        return self.__lower_limit <= value and value <= self.__upper_limit

    # Implements less-then operator, e.g. self < value
    def __lt__(self, value: int):
        return self.__upper_limit < value

    # Implements less-or-equals operator, e.g. self <= value
    def __le__(self, value: int):
        return self.__lower_limit <= value

    # Implements greater-than operator, e.g. self > value
    def __gt__(self, value: int):
        return self.__lower_limit > value

    # Implements greater-or-equals operator, e.g. self >= value
    def __ge__(self, value: int):
        return self.__upper_limit >= value

    # Implements equals, e.g. self == value
    def __eq__(self, value: int):
        return self.__contains__(value)

    # Returns the highest value of lower and upper tolerance
    @property
    def highest_tolerance(self):
        return max(self.__lower_tolerance_ppm, self.__upper_tolerance_ppm)
