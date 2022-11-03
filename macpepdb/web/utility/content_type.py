# std imports
from __future__ import annotations
from enum import Enum, unique
from typing import Optional, Iterable, List

# 3rd party imports
from flask import Request

@unique
class ContentType(Enum):
    """
    Helper class to determine and represent accepted content types.
    """
    json = "application/json"
    stream = "application/octet-stream"
    text = "text/plain"
    csv = "text/csv"
    fasta = "text/fasta"

    def __str__(self) -> str:
        """
        Returns content type as mime type (string)

        Returns
        -------
        str
            Mimetype
        """
        return f"{self.value}"

    @classmethod
    def from_name(cls, name: str) -> ContentType:
        """
        Returns the content type by name.

        Parameters
        ----------
        name : str
            Short name (json, csv, ...) of content type

        Returns
        -------
        ContentType
            Content type

        Raises
        ------
        KeyError
            If not found
        """
        return cls[name.lower()]

    @classmethod
    def from_value(cls, value: str) -> ContentType:
        """
        Returns content type by value.

        Parameters
        ----------
        value : str
            Content type value, e.g. "application/json"

        Returns
        -------
        ContentType
            Content type

        Raises
        ------
        KeyError
            If not found
        """
        lower_value = value.lower()
        for format in cls:
            if format.value == lower_value:
                return format
        raise KeyError(f"f{value} not found")

    @classmethod
    def get_content_type_by_request(cls, request: Request, permitted_content_types: Optional[Iterable[ContentType]] = None, default: Optional[ContentType] = None) -> ContentType:
        """
        Tries to determine the contentn type from requests via accept header or file extension.

        Parameters
        ----------
        request : Request
            Incoming requests
        permitted_content_types : Optional[List[ContentType]], optional
            List of permitted content types, by default all.
        default : ContentType, optional
            Default content type in case no permitted one was found, by default ContentType.json

        Returns
        -------
        ContentType
            Found content type
        """
        if permitted_content_types is None:
            permitted_content_types = list(ContentType)
        if default is None:
            default = ContentType.json
        accept_header: Optional[str] = request.headers.get("Accept", None)
        if accept_header is not None:
            try:
                accepted_content_type: ContentType = ContentType.from_value(accept_header)
                if accepted_content_type in permitted_content_types:
                    return accepted_content_type
            except KeyError:
                pass
        last_path_segment: str = request.path.split("/")[-1]
        if "." in last_path_segment:
            file_extension: str = last_path_segment.split(".")[-1]
            try:
                requested_content_type: ContentType = ContentType.from_name(file_extension)
                if requested_content_type in permitted_content_types:
                    return requested_content_type
            except KeyError:
                pass
        return default