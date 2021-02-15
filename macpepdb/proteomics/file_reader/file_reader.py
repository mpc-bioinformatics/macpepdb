from typing import IO

class FileReader:
    def __init__(self, file: IO):
        self.__file = file

    @property
    def file(self):
        return self.__file

    @classmethod
    def get_reader_by_file_format(cls, format: str):
        from .fasta_reader import FastaReader
        if format.upper() == "FASTA":
            return FastaReader
        elif format.upper() == "TEXT":
            from .uniprot_text_reader import UniprotTextReader
            return UniprotTextReader
        # elif format.upper() == "SOME OTHER FORMAT"
        else:
            raise NameError("Unknown input format '{}'.".format(format))