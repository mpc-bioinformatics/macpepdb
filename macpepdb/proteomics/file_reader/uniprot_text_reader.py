import re

from macpepdb.models.protein import Protein
from macpepdb.models.protein_merge import ProteinMerge

class UniprotTextReader():
    TAXONOMY_ID_REGEX = re.compile(r".*=(?P<taxonomy_id>\d+)")
    NAME_REGEX = re.compile(r"Full=(?P<name>.*?)(\{|;)")
    WHITESPACE_REGEX = re.compile(r"\s")
    SERIAL_WHITESPACES_REGEX = re.compile(r"\s{2,}")


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

        
        while True:
            line = next(self.__file_iter)

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
            elif line.startswith("//"):
                accession = accessions.pop(0)
                return Protein(accession, entry_name, name, sequence, taxonomy_id, proteome_id, is_reviewed), self.__create_protein_merges(accessions, accession)


    # Returns the the uniprot entry name and the review status (ture|false)
    def __process_id(self, line):
        # split line by multiple sequential whitespaces
        splitted_id_line = self.SERIAL_WHITESPACES_REGEX.split(line)
        return splitted_id_line[0], splitted_id_line[1] == "Reviewed;"

    # Returns accessions
    def __process_ac(self, line):
        # Split by whitespaces and return the first element without last char (';')
        return [accession[:-1] for accession in line.split()]

    # Returns the taxonomy id
    def __process_ox(self, line):
        matches = self.TAXONOMY_ID_REGEX.search(line)
        if matches and "taxonomy_id" in matches.groupdict():
            return int(matches.group("taxonomy_id"))
        else:
            return None

    # Returns the proteom id
    def __process_dr_proteoms(self, line):
        # Split line by spaces and return the second element without last character (';')
        return line.split()[1][:-1]

    # Returns the sequence without any whitespaces
    def __process_sq_no_header(self, line):
        return self.WHITESPACE_REGEX.sub("", line)

    # returns the value of the FullName attribute
    def __process_de_name(self, line):
        matches = self.NAME_REGEX.search(line)
        if matches:
            return matches["name"].strip()
        return ""
            
    def __create_protein_merges(self, source_accessions, target_accession):
        return [ProteinMerge(source_accession, target_accession) for source_accession in source_accessions]