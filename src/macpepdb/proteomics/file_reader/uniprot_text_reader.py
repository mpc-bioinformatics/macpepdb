# std imports
import re
from datetime import datetime, timedelta

# internal imports
from macpepdb.models.protein import Protein

class UniprotTextReader():
    """
    Reads EMBL files from UniProt.
    Class instances are iterable. The iterator will return instances of Protein.

    Parameters
    ----------
    file : TextIO
        EMBL file as file-like text object
    """

    TAXONOMY_ID_REGEX = re.compile(r".*=(?P<taxonomy_id>\d+)")
    """Matches the taxonomy ID after '='
    """

    NAME_REGEX = re.compile(r"Full=(?P<name>.*?)(\{|;)")
    """Regex for fullname
    """

    WHITESPACE_REGEX = re.compile(r"\s")
    SERIAL_WHITESPACES_REGEX = re.compile(r"\s{2,}")
    """Just a precompiled regex for matching whitespaces
    """
    
    DT_MONTH_LOOKUP_TABLE = {
        "JAN": 1,
        "FEB": 2,
        "MAR": 3,
        "APR": 4,
        "MAY": 5,
        "JUN": 6,
        "JUL": 7,
        "AUG": 8,
        "SEP": 9,
        "OCT": 10,
        "NOV": 11,
        "DEC": 12
    }
    """Lookup for month number by name. So no locale change is necessary
    """
    

    SERIAL_WHITESPACES_REGEX = re.compile(r"\s{2,}")
    """Just a precompiled regex for matching 2 or more sequential whitespaces
    """

    def __init__(self, file):
        self.__file = file

    def __iter__(self):
        self.__file_iter = iter(self.__file)
        return self

    def __next__(self):
        entry_name = ""
        name = ""
        is_reviewed = False
        accessions = []
        taxonomy_id = None
        sequence = ""
        proteome_id = None
        last_update = "01-JAN-1970"

        
        while True:
            line = next(self.__file_iter)
            
            line = line.rstrip()

            if len(line) >= 2:
                if line.startswith("ID"):
                    entry_name, is_reviewed = self.__process_id(line[5:])
                elif line.startswith("AC"):
                    accessions += self.__process_ac(line[5:])
                elif line.startswith("OX"):
                    taxonomy_id = self.__process_ox(line[5:])
                elif line.startswith("DR"):
                    if line[5:].startswith("Proteomes;"):
                        proteome_id = self.__process_dr_proteoms(line[5:])
                # sequence starts with two whitespaces
                elif line.startswith("  "):
                    sequence += self.__process_sq_no_header(line)
                elif line.startswith("DE"):
                    if name == "" and line[5:].startswith("RecName") or line[5:].startswith("AltName") or line[5:].startswith("Sub"):
                        name = self.__process_de_name(line[5:])
                elif line.startswith("DT"):
                    last_update = line[5:16]
                elif line.startswith("//"):
                    primary_accession = accessions.pop(0)
                    return Protein(primary_accession, accessions, entry_name, name, sequence, taxonomy_id, proteome_id, is_reviewed, self.__dt_date_to_utc_timestamp(last_update))

    def __process_id(self, line):
        """
        Returns the the uniprot entry name and the review status (ture|false)

        Parameters
        ----------
        line : str
            Line from EMBL file

        Returns
        -------
        Tuple with entry name and review status.
        """
        # split line by multiple sequential whitespaces
        splitted_id_line = self.SERIAL_WHITESPACES_REGEX.split(line)
        return splitted_id_line[0], splitted_id_line[1] == "Reviewed;"

    def __process_ac(self, line):
        """
        Returns accessions
        
        Parameters
        ----------
        line : str
            Line from EMBL file

        Returns
        -------
        List of accession
        """
        # Split by whitespaces and return the first element without last char (';')
        return [accession[:-1] for accession in line.split()]

    def __process_ox(self, line):
        """
        Returns the taxonomy id

        Parameters
        ----------
        line : str
            Line from EMBL file

        Returns
        -------
        Taxonomy ID or None
        """
        matches = self.TAXONOMY_ID_REGEX.search(line)
        if matches and "taxonomy_id" in matches.groupdict():
            return int(matches.group("taxonomy_id"))
        else:
            return None

    def __process_dr_proteoms(self, line):
        """
        Returns the proteom id

        Parameters
        ----------
        line : str
            Line from EMBL file

        Returns
        -------
        Proteome ID 
        """
        # Split line by spaces and return the second element without last character (';')
        return line.split()[1][:-1]

    def __process_sq_no_header(self, line):
        """
        Returns the sequence without any whitespaces

        Parameters
        ----------
        line : str
            Line from EMBL file

        Returns
        -------
        Amino Acid Sequence
        """
        return self.WHITESPACE_REGEX.sub("", line)

    def __process_de_name(self, line):
        """
        Returns the value of the FullName attribute

        Parameters
        ----------
        line : str
            Line from EMBL file

        Returns
        -------
        Fullname
        """
        matches = self.NAME_REGEX.search(line)
        if matches:
            return matches["name"].strip()
        return ""

    def __dt_date_to_utc_timestamp(self, dt_date: str) -> int:
        """
        Calculate UTC timestamp, see: https://docs.python.org/3/library/datetime.html#datetime.datetime.timestamp 

        Arguments
        ---------
        dt_date : str
            Date in the form 01-JAN-1970

        Return
        ------
        UTC unix timestamp
        """
        dt_date = dt_date.upper()
        day, month, year = dt_date.split("-")
        date = datetime(int(year), self.__class__.DT_MONTH_LOOKUP_TABLE.get(month, 1), int(day))
        return (date - datetime(1970, 1, 1)) / timedelta(seconds=1)