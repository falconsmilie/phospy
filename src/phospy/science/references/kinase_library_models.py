"""Kinase Library-style motif resource models."""

from __future__ import annotations

__phospy_contracts_facade_role__ = "science_owned_public_model"

import math
import re
from collections.abc import Sequence
from dataclasses import InitVar, dataclass
from datetime import date
from enum import Enum
from os import PathLike
from pathlib import Path
from typing import cast

import pandas as pd

from phospy.errors.validation import ReferenceValidationError
from phospy.frames.ownership import export_dataframe, own_dataframe
from phospy.provenance.models import KinaseLibraryResourceProvenance
from phospy.science.prediction.motif_scoring.models import AMINO_ACIDS
from phospy.science.references.models import (
    Organism,
    ReferenceBuildPath,
    SequenceWindowDefinition,
)

KinaseLibraryPath = str | Path | PathLike[str]

_COLUMN_TOKEN_PATTERN = re.compile(r"[^a-z0-9+-]+")

_ORGANISM_SPLIT_PATTERN = re.compile(r"[|,]+")


class KinaseLibraryResidueClass(str, Enum):
    """Supported Kinase Library phospho-acceptor lanes."""

    SER_THR = "ser_thr"
    TYR = "tyr"


@dataclass(frozen=True, slots=True)
class KinaseLibraryResourceLoadRequest:
    """Request for loading a local Kinase Library-style matrix resource.

    Metadata may be supplied here or repeated in the source table. Loader request
    values take precedence over source-table metadata.
    """

    path: ReferenceBuildPath
    source_name: str | None = None
    source_version: str | None = None
    retrieved_at: date | str | None = None
    license: str | None = None
    score_scale: str | None = None
    sequence_window: SequenceWindowDefinition | None = None
    organisms: tuple[Organism | str, ...] = ()


@dataclass(frozen=True, slots=True)
class KinaseLibraryMatrix:
    """Position-specific Kinase Library score table for one kinase lane."""

    kinase: str
    residue_class: KinaseLibraryResidueClass
    score_table: pd.DataFrame
    kinase_family: str | None = None
    kinase_group: str | None = None
    _assume_owned: InitVar[bool] = False

    def __post_init__(self, _assume_owned: bool) -> None:
        kinase = _require_text(self.kinase, field_name="kinase_library.kinase")
        residue_class = _coerce_residue_class(
            self.residue_class,
            field_name=f"kinase_library[{kinase}].residue_class",
        )
        score_table = own_dataframe(
            self.score_table,
            field_name=f"kinase_library[{kinase}].score_table",
            error_type=ReferenceValidationError,
            assume_owned=_assume_owned,
        )
        score_table = _normalise_score_table(
            score_table,
            context=f"kinase_library[{kinase}:{residue_class.value}]",
        )
        object.__setattr__(self, "kinase", kinase)
        object.__setattr__(self, "residue_class", residue_class)
        object.__setattr__(
            self,
            "kinase_family",
            _optional_text(
                self.kinase_family,
                field_name=f"kinase_library[{kinase}].kinase_family",
            ),
        )
        object.__setattr__(
            self,
            "kinase_group",
            _optional_text(
                self.kinase_group,
                field_name=f"kinase_library[{kinase}].kinase_group",
            ),
        )
        object.__setattr__(self, "score_table", score_table)

    def score_table_dataframe(self) -> pd.DataFrame:
        """Return a defensive snapshot of the raw position-specific score table."""

        return export_dataframe(self.score_table)


@dataclass(frozen=True, slots=True)
class KinaseLibraryResource:
    """Validated local Kinase Library-style matrix resource.

    Scores are preserved as provider-scale numeric values. This model does not
    interpret them as probabilities and does not connect them to workflow
    scoring.
    """

    matrices: tuple[KinaseLibraryMatrix, ...]
    source_name: str
    source_version: str
    score_scale: str
    sequence_window: SequenceWindowDefinition
    organisms: tuple[str, ...]
    license: str
    provenance: KinaseLibraryResourceProvenance
    retrieved_at: str | None = None

    def __post_init__(self) -> None:
        matrices = _coerce_matrix_tuple(self.matrices)
        source_name = _require_text(
            self.source_name,
            field_name="kinase_library.source_name",
        )
        source_version = _require_text(
            self.source_version,
            field_name="kinase_library.source_version",
        )
        score_scale = _require_text(
            self.score_scale,
            field_name="kinase_library.score_scale",
        )
        license_name = _require_text(
            self.license,
            field_name="kinase_library.license",
        )
        sequence_window = self.sequence_window
        if not isinstance(cast(object, sequence_window), SequenceWindowDefinition):
            raise ReferenceValidationError(
                "kinase_library.sequence_window must be SequenceWindowDefinition"
            )
        organisms = _coerce_organism_tokens(
            self.organisms,
            field_name="kinase_library.organisms",
        )
        retrieved_at = _optional_text(
            self.retrieved_at,
            field_name="kinase_library.retrieved_at",
        )
        provenance = self.provenance
        if not isinstance(cast(object, provenance), KinaseLibraryResourceProvenance):
            raise ReferenceValidationError(
                "kinase_library.provenance must be KinaseLibraryResourceProvenance"
            )

        object.__setattr__(self, "matrices", matrices)
        object.__setattr__(self, "source_name", source_name)
        object.__setattr__(self, "source_version", source_version)
        object.__setattr__(self, "score_scale", score_scale)
        object.__setattr__(self, "license", license_name)
        object.__setattr__(self, "organisms", organisms)
        object.__setattr__(self, "retrieved_at", retrieved_at)

        _validate_kinase_library_resource_contract(self)

    def matrix_for(
        self,
        kinase: str,
        residue_class: KinaseLibraryResidueClass | str,
    ) -> KinaseLibraryMatrix:
        """Return the matrix for one kinase and residue-class lane."""

        kinase_id = _require_text(kinase, field_name="kinase")
        residue = _coerce_residue_class(
            residue_class,
            field_name="residue_class",
        )
        for matrix in self.matrices:
            if matrix.kinase == kinase_id and matrix.residue_class is residue:
                return matrix
        raise KeyError(f"no Kinase Library matrix for {kinase_id}/{residue.value}")


def _validate_kinase_library_resource_contract(
    resource: KinaseLibraryResource,
) -> None:
    if not isinstance(resource, KinaseLibraryResource):
        raise ReferenceValidationError(
            "kinase_library resource must be KinaseLibraryResource"
        )
    _validate_sequence_window(resource.sequence_window)
    _validate_metadata(resource)
    _validate_matrices(
        matrices=resource.matrices,
        sequence_window=resource.sequence_window,
    )
    _validate_provenance(resource.provenance)


def _validate_metadata(resource: KinaseLibraryResource) -> None:
    for field_name, value in (
        ("source_name", resource.source_name),
        ("source_version", resource.source_version),
        ("score_scale", resource.score_scale),
        ("license", resource.license),
    ):
        if not isinstance(value, str) or not value.strip():
            raise ReferenceValidationError(
                f"kinase_library.{field_name} must be a non-empty string"
            )
    if not resource.organisms:
        raise ReferenceValidationError("kinase_library.organisms must not be empty")
    for organism in resource.organisms:
        if not isinstance(organism, str) or not organism.strip():
            raise ReferenceValidationError(
                "kinase_library.organisms must contain non-empty strings"
            )


def _validate_sequence_window(
    sequence_window: SequenceWindowDefinition,
) -> None:
    if not isinstance(cast(object, sequence_window), SequenceWindowDefinition):
        raise ReferenceValidationError(
            "kinase_library.sequence_window must be SequenceWindowDefinition"
        )
    if sequence_window.upstream_residues < 0:
        raise ReferenceValidationError(
            "kinase_library.sequence_window.upstream_residues must be >= 0"
        )
    if sequence_window.downstream_residues < 0:
        raise ReferenceValidationError(
            "kinase_library.sequence_window.downstream_residues must be >= 0"
        )
    if not isinstance(sequence_window.central_residue_required, bool):
        raise ReferenceValidationError(
            "kinase_library.sequence_window.central_residue_required must be bool"
        )


def _validate_matrices(
    *,
    matrices: tuple[KinaseLibraryMatrix, ...],
    sequence_window: SequenceWindowDefinition,
) -> None:
    if not matrices:
        raise ReferenceValidationError("kinase_library.matrices must not be empty")
    expected_positions = tuple(
        range(
            -int(sequence_window.upstream_residues),
            int(sequence_window.downstream_residues) + 1,
        )
    )
    seen_keys: set[tuple[str, str]] = set()
    for matrix in matrices:
        if not isinstance(matrix, KinaseLibraryMatrix):
            raise ReferenceValidationError(
                "kinase_library.matrices must contain KinaseLibraryMatrix values"
            )
        key = (matrix.kinase, matrix.residue_class.value)
        if key in seen_keys:
            raise ReferenceValidationError(
                "kinase_library.matrices contains duplicate kinase/residue_class "
                f"entry: {matrix.kinase}/{matrix.residue_class.value}"
            )
        seen_keys.add(key)
        _validate_score_table(
            matrix.score_table,
            expected_positions=expected_positions,
            context=(f"kinase_library[{matrix.kinase}:{matrix.residue_class.value}]"),
        )


def _validate_score_table(
    score_table: pd.DataFrame,
    *,
    expected_positions: tuple[int, ...],
    context: str,
) -> None:
    if score_table.empty:
        raise ReferenceValidationError(f"{context}.score_table must be non-empty")
    observed_positions = tuple(int(position) for position in score_table.columns)
    missing_positions = [
        position
        for position in expected_positions
        if position not in observed_positions
    ]
    unexpected_positions = [
        position
        for position in observed_positions
        if position not in expected_positions
    ]
    if missing_positions:
        raise ReferenceValidationError(
            f"{context}.score_table is missing required positions: "
            f"{_format_positions(missing_positions)}"
        )
    if unexpected_positions:
        raise ReferenceValidationError(
            f"{context}.score_table contains positions outside sequence_window: "
            f"{_format_positions(unexpected_positions)}"
        )
    if observed_positions != expected_positions:
        raise ReferenceValidationError(
            f"{context}.score_table positions must be ordered as "
            f"{_format_positions(list(expected_positions))}"
        )
    if score_table.isna().to_numpy().any():
        raise ReferenceValidationError(
            f"{context}.score_table contains missing score values"
        )
    values = score_table.to_numpy(dtype=float, copy=False)
    for value in values.ravel().tolist():
        if not math.isfinite(float(value)):
            raise ReferenceValidationError(
                f"{context}.score_table contains non-finite score values"
            )


def _validate_provenance(
    provenance: KinaseLibraryResourceProvenance,
) -> None:
    if not isinstance(cast(object, provenance), KinaseLibraryResourceProvenance):
        raise ReferenceValidationError(
            "kinase_library.provenance must be KinaseLibraryResourceProvenance"
        )
    if provenance.source_type != "local":
        raise ReferenceValidationError(
            "kinase_library.provenance.source_type must be 'local'"
        )
    if not provenance.source_files:
        raise ReferenceValidationError(
            "kinase_library.provenance.source_files must not be empty"
        )
    if not provenance.table_fingerprints:
        raise ReferenceValidationError(
            "kinase_library.provenance.table_fingerprints must not be empty"
        )


def _coerce_matrix_tuple(
    matrices: tuple[KinaseLibraryMatrix, ...],
) -> tuple[KinaseLibraryMatrix, ...]:
    if not isinstance(matrices, tuple):
        raise ReferenceValidationError("kinase_library.matrices must be a tuple")
    if not matrices:
        raise ReferenceValidationError("kinase_library.matrices must not be empty")
    for matrix in matrices:
        if not isinstance(matrix, KinaseLibraryMatrix):
            raise ReferenceValidationError(
                "kinase_library.matrices must contain KinaseLibraryMatrix values"
            )
    return matrices


def _normalise_score_table(score_table: pd.DataFrame, *, context: str) -> pd.DataFrame:
    if score_table.empty:
        raise ReferenceValidationError(f"{context}.score_table must be non-empty")
    normalized = score_table.copy(deep=True)
    amino_acids = [
        _coerce_amino_acid(value, field_name=f"{context}.score_table.index")
        for value in normalized.index.tolist()
    ]
    positions = [
        _coerce_position(value, field_name=f"{context}.score_table.columns")
        for value in normalized.columns.tolist()
    ]
    if len(set(amino_acids)) != len(amino_acids):
        raise ReferenceValidationError(
            f"{context}.score_table contains duplicate amino-acid rows"
        )
    if len(set(positions)) != len(positions):
        raise ReferenceValidationError(
            f"{context}.score_table contains duplicate position columns"
        )
    normalized.index = pd.Index(amino_acids, name="amino_acid")
    normalized.columns = pd.Index(positions, name="position")
    try:
        return normalized.astype(float)
    except ValueError as exc:
        raise ReferenceValidationError(
            f"{context}.score_table must contain numeric score values"
        ) from exc


def _coerce_residue_class(
    value: object,
    *,
    field_name: str,
) -> KinaseLibraryResidueClass:
    if isinstance(value, KinaseLibraryResidueClass):
        return value
    text = _require_text(value, field_name=field_name).lower()
    normalized = _COLUMN_TOKEN_PATTERN.sub("_", text).strip("_")
    if normalized in {"ser_thr", "s_t", "st", "serine_threonine"}:
        return KinaseLibraryResidueClass.SER_THR
    if normalized in {"tyr", "y", "tyrosine"}:
        return KinaseLibraryResidueClass.TYR
    raise ReferenceValidationError(
        f"{field_name} must be 'ser_thr' or 'tyr'; got {text!r}"
    )


def _coerce_amino_acid(value: object, *, field_name: str) -> str:
    amino_acid = _require_text(value, field_name=field_name).upper()
    if len(amino_acid) != 1 or amino_acid not in AMINO_ACIDS:
        raise ReferenceValidationError(
            f"{field_name} must be one supported one-letter amino-acid code"
        )
    return amino_acid


def _coerce_position(value: object, *, field_name: str) -> int:
    if isinstance(value, int) and not isinstance(value, bool):
        return int(value)
    text = _require_text(value, field_name=field_name)
    try:
        return int(text)
    except ValueError as exc:
        raise ReferenceValidationError(f"{field_name} must be an integer") from exc


def _coerce_organism_tokens(
    values: Sequence[object],
    *,
    field_name: str,
) -> tuple[str, ...]:
    tokens: list[str] = []
    for value in values:
        if isinstance(value, Organism):
            token = value.value
        else:
            token = _require_text(value, field_name=field_name).lower()
        for split_token in _split_organism_value(token):
            if split_token not in tokens:
                tokens.append(split_token)
    if not tokens:
        raise ReferenceValidationError(f"{field_name} must not be empty")
    return tuple(tokens)


def _split_organism_value(value: object) -> list[str]:
    text = _require_text(value, field_name="kinase_library.organisms").lower()
    return [
        token.strip() for token in _ORGANISM_SPLIT_PATTERN.split(text) if token.strip()
    ]


def _require_text(value: object, *, field_name: str) -> str:
    text = _optional_text(value, field_name=field_name)
    if text is None:
        raise ReferenceValidationError(f"{field_name} must be a non-empty string")
    return text


def _optional_text(value: object, *, field_name: str) -> str | None:
    if value is None:
        return None
    if isinstance(value, date):
        return value.isoformat()
    if value is pd.NA or value is pd.NaT:
        return None
    if isinstance(value, float) and pd.isna(value):
        return None
    text = str(value).strip()
    if text == "":
        return None
    return text


def _format_positions(positions: list[int]) -> str:
    if not positions:
        return "(none)"
    return ", ".join(str(position) for position in positions)


__all__ = [
    "KinaseLibraryMatrix",
    "KinaseLibraryPath",
    "KinaseLibraryResidueClass",
    "KinaseLibraryResource",
    "KinaseLibraryResourceLoadRequest",
]
