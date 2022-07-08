# std imports
from typing import Iterator, ByteString

# 3rd party imports
from macpepdb.models.protein import Protein as BaseProtein

class Protein(BaseProtein):

    def to_json(self) -> Iterator[ByteString]:
        """
        Generator which yields the protein as a json formatted string

        Yields
        ------
        Iterator[ByteString]
            JSON formatted string.
        """
        yield b"{\"accession\":\""
        yield self.accession.encode("utf-8")
        yield b"\",\"entry_name\":"
        yield self.accession.encode("utf-8")
        yield b"\",\"name\":\""
        yield self.name.encode("utf-8")
        yield b"\",\"sequence\":\""
        yield self.sequence.encode("utf-8")
        yield b"\",\"taxonomy_id\":"
        yield str(self.taxonomy_id).encode("utf-8")
        yield b",\"proteome_id\":\""
        yield str(self.proteome_id).encode("utf-8")
        yield b"\",\"is_reviewed\":"
        yield b"true" if self.is_reviewed else b"false"
        yield b"}"
