"""Reference manifest file and sequence-window value models."""

from __future__ import annotations

from dataclasses import dataclass

from phospy.provenance.models import JsonValue


@dataclass(frozen=True, slots=True)
class SequenceWindowDefinition:
    """Reference sequence-window definition for centralized site sequences."""

    upstream_residues: int
    downstream_residues: int
    central_residue_required: bool

    def to_payload(self) -> dict[str, JsonValue]:
        return {
            "upstream_residues": int(self.upstream_residues),
            "downstream_residues": int(self.downstream_residues),
            "central_residue_required": bool(self.central_residue_required),
        }


@dataclass(frozen=True, slots=True)
class ReferenceFileManifest:
    """Hash-verifiable file metadata for one file in a reference bundle."""

    relative_path: str
    role: str
    format: str
    sha256: str
    row_count: int | None = None
    column_names: tuple[str, ...] | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "relative_path", str(self.relative_path).strip())
        object.__setattr__(self, "role", str(self.role).strip())
        object.__setattr__(self, "format", str(self.format).strip())
        object.__setattr__(self, "sha256", str(self.sha256).strip())
        if self.row_count is not None:
            object.__setattr__(self, "row_count", int(self.row_count))
        if self.column_names is not None:
            object.__setattr__(
                self,
                "column_names",
                tuple(str(item) for item in self.column_names),
            )

    def to_payload(self) -> dict[str, JsonValue]:
        return {
            "relative_path": self.relative_path,
            "role": self.role,
            "format": self.format,
            "sha256": self.sha256,
            "row_count": self.row_count,
            "column_names": (
                None if self.column_names is None else list(self.column_names)
            ),
        }

    def to_source_file_payload(self) -> dict[str, JsonValue]:
        return {
            "path": self.relative_path,
            "relative_path": self.relative_path,
            "role": self.role,
            "format": self.format,
            "sha256": self.sha256,
            "row_count": self.row_count,
            "column_names": (
                None if self.column_names is None else list(self.column_names)
            ),
        }


def _source_file_key(file_manifest: ReferenceFileManifest) -> str:
    normalized = file_manifest.role.strip().lower().replace("-", "_").replace(" ", "_")
    return normalized or file_manifest.relative_path


__all__ = [
    "ReferenceFileManifest",
    "SequenceWindowDefinition",
]
