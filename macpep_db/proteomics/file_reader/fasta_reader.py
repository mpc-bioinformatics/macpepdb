import re
from typing import IO
from ...models.protein import Protein
from ...models.protein_merge import ProteinMerge

from .file_reader import FileReader

# Reads fasta file entry (header, sequence) by entry. It implements an iterator so it can be used in for-loops
class FastaReader(FileReader):
    ##
    # Line starts with ">"
    # >
    # Followd by "sp" (Swissprot) or "tr" (TrEMBL) in our case the review status
    # (?P<review_status>sp|tr)
    # Followed by "|" and the accession
    # \|(?P<accession>.+?)
    # Followed by the "|" and the entry_name
    # \|(?P<entry_name>\w+)
    # Followed by one ore more whitespaces "\s+"
    # \s+
    # Followed by the name which is matched by a positive lookahead "(?=...)". The name ...
    # ... cannot start with uppercase letters followed by a equals sign "(?![A-Z]{2,}=)".
    # ... can contain anything: ".*".
    # ... stops by one or more whitespaces followed by uppercase letters and an equals sign (first attribute) OR by zero or more whitespaces and a line end: "(\s+[A-Z]{2,}=|\s*$))".
    # ... it can also occure only zero or one: "{0,1}"
    # The reason it is matched with a positive lookahead is the absence of a defined stop sign. 
    # On of it possible stops is the name of the following attributes, which is then matched and cannot be matched again.
    # To use a positive lookahead right after the entry name prevents this but needs 
    # (?=(?P<name>(?![A-Z]{2,}=).+?)(\s+[A-Z]{2,}=|\s*$)){0,1}
    # Next, the attributes are follow. Starting with upper case letters and a equals sign.
    # This are matched by positive lookaheads (?= ...), so they can occure in any abitrary order.
    # There are also two MaCPep DB specific attributes "TDBPI" the proteome id and "TDBPM" old protein accessions which are merged with this proteins.
    # They will only added for unprocessible proteins during digestion if the input file was in text format.
    # The attribute regexes start with anything (because they are anchored after the entry name and the name must be skipped, see the description for the name) followed by uppercase letters and an equals sign.
    # The content can vary, so can the stop condition.
    # Each attribute can occure only never or once: "{0,1}"
    # Taxonomy id: Contain only numbers "\d+"
    # (?=.*OX=(?P<taxonomy_id>\d+)){0,1}
    # Protein id (MaCPep DB-specific): Contain only numbers "\d+"
    # (?=.*TDBPI=(?P<proteome_id>UP\d+)){0,1}
    # Protein merges: Contain non-whitespace characters
    # (?=.*TDBPM=(?P<protein_merges>(\S+))){0,1}
    # Other attributes can be added in the same way. For example the Organism, which can contain anything, so it mus stop by the next attribute or on a line ending: "(\s[A-Z]{2,}=|\s*$)"
    # (?=.*OS=(?P<organism>.+?)(\s[A-Z]{2,}=|\s*$)){0,1}
    HEADER_REGEX = re.compile(r">(?P<review_status>sp|tr)\|(?P<accession>.+?)\|(?P<entry_name>\w+)\s+(?=(?P<name>(?![A-Z]{2,}=).+?)(\s+[A-Z]{2,}=|\s*$)){0,1}(?=.*OX=(?P<taxonomy_id>\d+)){0,1}(?=.*TDBPI=(?P<proteome_id>UP\d+)){0,1}(?=.*TDBPM=(?P<protein_merges>(\S+))){0,1}")

    def __init__(self, file: IO):
        super().__init__(file)
        self.__file_iter = None
        self.__last_read_protein_header = None
        self.__has_reached_eof = False


    def __iter__(self):
        self.__file_iter = iter(self.file)
        return self

    def __next__(self):
        # If the EOF is reached, close the file and stop the FastaReader
        if self.__has_reached_eof:
            raise StopIteration

        sequence = ""
        try:
            while True:
                # Read next line and strip leading and tailing whitespaces (newlines, etc.)
                line = next(self.__file_iter).strip()
                # Lines not starting with ">" are sequence lines
                if not line.startswith(">"):
                    sequence += line
                else:
                    # Lines starting with ">" are headers.
                    # If we read a header before ...
                    if self.__last_read_protein_header:
                        is_reviewed, accession, entry_name, name, taxonomy_id, proteome_id, accession_merges =  self.__get_information_from_header()
                        # ... initialize the protein with the sequence
                        protein = Protein(accession, entry_name, name, sequence, taxonomy_id, proteome_id, is_reviewed)
                        # store the current header for the next protein
                        self.__last_read_protein_header = line
                        # and return the new protein. Fasta files don't contain information about protein merges, so an empty array is returned
                        return protein, accession_merges
                    else:
                        # If no header is stored (at the beginnging of the file) we need to store the first header and read the sequence lines 
                        self.__last_read_protein_header = line
        # If next(self.__file_iter) raises a StopIteration, we've reached the EOF ...
        except StopIteration:
            # stores this for next call of this iterator
            self.__has_reached_eof = True
            is_reviewed, accession, entry_name, name, taxonomy_id, proteome_id, accession_merges =  self.__get_information_from_header()
            # and can build and return the last protein. Fasta files don't contain information about protein merges, so an empty array is returned
            return Protein(accession, entry_name, name, sequence, taxonomy_id, proteome_id, is_reviewed), accession_merges
            
    # Returns accession and taxonomy id from header
    def __get_information_from_header(self) -> tuple:
        is_reviewed = False
        accession = "None"
        name = ""
        entry_name = ""
        taxonomy_id = None
        proteome_id = None
        protein_merges = []

        matches = self.HEADER_REGEX.search(self.__last_read_protein_header)
        if matches:
            header_informations = matches.groupdict()
            if "review_status" in header_informations:
                is_reviewed = header_informations["review_status"] == "sp"
            if "accession" in header_informations:
                accession = header_informations["accession"]
            if "entry_name" in header_informations:
                entry_name = header_informations["entry_name"]
            if "name" in header_informations:
                name = header_informations["name"]
            if "taxonomy_id" in header_informations:
                taxonomy_id = int(header_informations["taxonomy_id"])
            # The following two attributes are not part of the actual fasta format, see regex description.
            if "proteome_id" in header_informations:
                proteome_id = header_informations["proteome_id"]
            if "protein_merges" in header_informations and header_informations["protein_merges"]:
                # Split commaseperated list
                accessions = header_informations["protein_merges"].split(",")
                protein_merges = [ProteinMerge(source_accession, accession) for source_accession in accessions]

        return is_reviewed, accession, entry_name, name, taxonomy_id, proteome_id, protein_merges

