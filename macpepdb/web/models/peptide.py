# std imports
import base64
from typing import ClassVar, Iterator, ByteString, List, Optional
import zlib

# 3rd party imports
from macpepdb.models.peptide import Peptide as BasePeptide
from macpepdb.proteomics.mass.convert import to_float as mass_to_float


class Peptide(BasePeptide):
    CSV_HEADER: ClassVar[List[str]] = ["mass", "sequence", "number_of_missed_cleavages"]
    """Header for CSV output
    """

    METADATA_CSV_HEADER: ClassVar[List[str]] = ["in_swiss_prot", "in_trembl", "taxonomy_ids", "unique_for_taxonomy_ids", "proteome_ids"]
    """Additional CSV header for metadata
    """

    @classmethod
    def __taxonomy_ids_to_json(cls, taxonomy_ids: List[int]) -> Iterator[ByteString]:
        """
        Generator which yields array of taxonomy IDs as json formatted string.

        Parameters
        ----------
        taxonomy_ids : List[int]
            Taxonomy IDs

        Yields
        ------
        Iterator[ByteString]
            JSON formatted string
        """
        yield b"["
        for idx, taxonomy_id in enumerate(taxonomy_ids):
            if idx > 0:
                yield b","
            yield str(taxonomy_id).encode("utf-8")
        yield b"]"

    def to_json(self) -> Iterator[ByteString]:
        """
        Generator which yields the peptides as a json formatted string

        Yields
        ------
        Iterator[ByteString]
            JSON formatted string.
        """
        yield b"{\"mass\":"
        yield str(mass_to_float(self.mass)).encode("utf-8")
        yield b",\"sequence\":\""
        yield self.sequence.encode("utf-8")
        yield b"\",\"length\":"
        yield str(self.length).encode("utf-8")
        yield b",\"number_of_missed_cleavages\":"
        yield str(self.number_of_missed_cleavages).encode("utf-8")
        yield b",\"metadata\":"
        if self.metadata is not None:
            yield b"{\"is_swiss_prot\":"
            yield b"true" if self.metadata.is_swiss_prot else b"false"
            yield b",\"is_trembl\":"
            yield b"true" if self.metadata.is_trembl else b"false"
            yield b",\"taxonomy_ids\":"
            for taxonomy_id_chunk in self.__class__.__taxonomy_ids_to_json(self.metadata.taxonomy_ids):
                yield taxonomy_id_chunk
            yield b",\"unique_taxonomy_ids\":"
            for taxonomy_id_chunk in self.__class__.__taxonomy_ids_to_json(self.metadata.unique_taxonomy_ids):
                yield taxonomy_id_chunk
            yield b",\"proteome_ids\":["
            for idx, proteome_id in enumerate(self.metadata.proteome_ids):
                if idx > 0:
                    yield b","
                yield b"\""
                yield proteome_id.encode("utf-8")
                yield b"\""
            yield b"]}"
        else:
            yield b"null"
        yield b"}"

    def to_fasta_entry(self, accession: Optional[ByteString] = None) -> Iterator[ByteString]:
        """
        Peptide as FASTA entry

        Parameters
        ----------
        accession : Optional[ByteString]
            Accession for fasta entry. If none, the base64 encoded, zlib compression of the encoded sequence is used.

        Yields
        ------
        Iterator[ByteString]
            FASTA entry
        """
        # Begin FASTA entry with '>macpepdb|' ...
        yield b">macpepdb|"
        encoded_seqeunce = self.sequence.encode()
        # Use 'P' + index or base64 encoded, zlib compression if sequence as accession
        yield accession if accession is not None else base64.b64encode(zlib.compress(encoded_seqeunce))
        yield b"|"
        # ... add the sequence in chunks of 60 amino acids ...
        for chunk_start in range(0, len(self.encoded_seqeunce), 60):
            yield b"\n"
            yield self.encoded_seqeunce[chunk_start : chunk_start+60]

    def to_csv_row(self) -> Iterator[ByteString]:
        """
        Yields the peptide as CSV row

        Yields
        ------
        Iterator[ByteString]
            CSV row with 2 or 7 columns, depending if the peptide was queried with metadata
        """
        yield str(mass_to_float(self.mass)).encode("utf-8")
        # At this point only string values are added to the cssv, so we quote them for better compatibility
        yield b",\""
        yield self.sequence.encode("utf-8")
        yield b"\","
        yield str(self.number_of_missed_cleavages).encode("utf-8")
        if self.metadata is not None:
            yield ",\""
            yield b"true" if self.metadata.is_swiss_prot else b"false"
            yield b"\",\""
            yield b"true" if self.metadata.is_trembl else b"false"
            yield b"\",\""
            yield ",".join([str(taxonomy_id) for taxonomy_id in self.metadata.taxonomy_ids]).encode()
            yield b"\",\""
            yield ",".join([str(taxonomy_id) for taxonomy_id in self.metadata.unique_taxonomy_ids]).encode()
            yield b"\",\""
            yield ",".join([f"{proteome_id}" for proteome_id in self.metadata.proteome_ids]).encode()
            yield b"\""

    def to_plain_text(self) -> Iterator[ByteString]:
        yield self.sequence.encode("utf-8")
