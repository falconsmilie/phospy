"""Kinase Library resource provenance."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

from phospy.provenance.immutability import (
    freeze_json_mapping,
    freeze_optional_json_mapping,
)
from phospy.provenance.models._shared import JsonValue
from phospy.provenance.models.tables import (
    TableFingerprint,
    _required_table_fingerprint_tuple,
)


@dataclass(frozen=True, slots=True)
class KinaseLibraryResourceProvenance:
    """Resolved provenance for local Kinase Library-style motif resources."""

    source_type: str
    source_name: str
    source_version: str
    license: str
    score_scale: str
    organisms: tuple[str, ...]
    sequence_window: Mapping[str, JsonValue]
    source_files: Mapping[str, JsonValue]
    table_fingerprints: tuple[TableFingerprint, ...]
    retrieved_at: str | None = None
    manifest: Mapping[str, JsonValue] | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "organisms", tuple(self.organisms))
        object.__setattr__(
            self,
            "sequence_window",
            freeze_json_mapping(
                self.sequence_window,
                field_name="kinase_library_resource_provenance.sequence_window",
            ),
        )
        object.__setattr__(
            self,
            "source_files",
            freeze_json_mapping(
                self.source_files,
                field_name="kinase_library_resource_provenance.source_files",
            ),
        )
        object.__setattr__(
            self,
            "table_fingerprints",
            _required_table_fingerprint_tuple(
                self.table_fingerprints,
                field_name="kinase_library_resource_provenance.table_fingerprints",
            ),
        )
        object.__setattr__(
            self,
            "manifest",
            freeze_optional_json_mapping(
                self.manifest,
                field_name="kinase_library_resource_provenance.manifest",
            ),
        )
