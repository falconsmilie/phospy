"""Kinase Library-style motif resource schema and local loader."""

from __future__ import annotations

import re
from collections.abc import Sequence
from dataclasses import InitVar, dataclass
from datetime import date
from enum import Enum
from os import PathLike
from pathlib import Path
from typing import cast

import pandas as pd

from phospy.errors.references import ReferenceResolutionError
from phospy.errors.validation import ReferenceValidationError
from phospy.frames.ownership import export_dataframe, own_dataframe
from phospy.io.bundles.reference_sources import ReferenceSourceTableReader
from phospy.provenance.hashing import fingerprint_table
from phospy.provenance.models import JsonValue, KinaseLibraryResourceProvenance
from phospy.provenance.references import fingerprint_local_reference_source_file
from phospy.science.prediction.motif_scoring.models import AMINO_ACIDS
from phospy.science.references.models import (
    Organism,
    ReferenceBuildPath,
    SequenceWindowDefinition,
)

KinaseLibraryPath = str | Path | PathLike[str]

_COLUMN_TOKEN_PATTERN = re.compile(r"[^a-z0-9+-]+")
_POSITION_COLUMN_PATTERN = re.compile(r"^(?:p|pos|position)?_?([+-]?\d+)$")
_ORGANISM_SPLIT_PATTERN = re.compile(r"[|,]+")
_TRUE_TOKENS = frozenset({"true", "t", "1", "yes", "y"})
_FALSE_TOKENS = frozenset({"false", "f", "0", "no", "n"})
_KINASE_ALIASES = (
    "kinase",
    "kinase_id",
    "kinase_identifier",
    "kinase_name",
    "gene",
    "gene_symbol",
)
_KINASE_FAMILY_ALIASES = ("kinase_family", "family")
_KINASE_GROUP_ALIASES = ("kinase_group", "group")
_RESIDUE_CLASS_ALIASES = (
    "residue_class",
    "residue_type",
    "substrate_residue_class",
    "phospho_acceptor_class",
)
_POSITION_ALIASES = ("position", "pos", "relative_position")
_AMINO_ACID_ALIASES = ("amino_acid", "aa", "residue")
_SCORE_ALIASES = ("score", "matrix_score", "value")
_SOURCE_NAME_ALIASES = ("source_name", "source")
_SOURCE_VERSION_ALIASES = ("source_version", "version")
_RETRIEVED_AT_ALIASES = ("retrieved_at", "downloaded_at")
_LICENSE_ALIASES = ("license", "licence")
_SCORE_SCALE_ALIASES = ("score_scale", "scale")
_ORGANISMS_ALIASES = ("organisms", "organism", "species", "applicable_organisms")
_UPSTREAM_ALIASES = ("upstream_residues", "upstream", "n_terminal_positions")
_DOWNSTREAM_ALIASES = ("downstream_residues", "downstream", "c_terminal_positions")
_CENTRAL_REQUIRED_ALIASES = (
    "central_residue_required",
    "requires_central_residue",
)
_LONG_REQUIRED_ALIASES = (
    _KINASE_ALIASES,
    _RESIDUE_CLASS_ALIASES,
    _POSITION_ALIASES,
    _AMINO_ACID_ALIASES,
    _SCORE_ALIASES,
)
_METADATA_ALIAS_GROUPS = (
    _KINASE_ALIASES,
    _KINASE_FAMILY_ALIASES,
    _KINASE_GROUP_ALIASES,
    _RESIDUE_CLASS_ALIASES,
    _POSITION_ALIASES,
    _AMINO_ACID_ALIASES,
    _SCORE_ALIASES,
    _SOURCE_NAME_ALIASES,
    _SOURCE_VERSION_ALIASES,
    _RETRIEVED_AT_ALIASES,
    _LICENSE_ALIASES,
    _SCORE_SCALE_ALIASES,
    _ORGANISMS_ALIASES,
    _UPSTREAM_ALIASES,
    _DOWNSTREAM_ALIASES,
    _CENTRAL_REQUIRED_ALIASES,
)


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

        from phospy.validation.references.kinase_library import (
            KinaseLibraryResourceValidator,
        )

        KinaseLibraryResourceValidator().run(self)

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


class KinaseLibraryResourceLoader:
    """Load local Kinase Library-style files into validated resource models."""

    def __init__(self, source_reader: ReferenceSourceTableReader | None = None) -> None:
        self._source_reader = source_reader or ReferenceSourceTableReader()

    def run(
        self,
        request: KinaseLibraryResourceLoadRequest | KinaseLibraryPath,
    ) -> KinaseLibraryResource:
        """Read a local matrix table and return a validated resource."""

        load_request = _coerce_load_request(request)
        path = Path(load_request.path)
        source = self._source_reader.run(
            path,
            field_name="kinase library resource path",
        )
        long_table = _normalise_source_table(source)
        source_name = _resolve_metadata_text(
            source,
            aliases=_SOURCE_NAME_ALIASES,
            override=load_request.source_name,
            field_name="kinase_library.source_name",
            required=True,
        )
        source_version = _resolve_metadata_text(
            source,
            aliases=_SOURCE_VERSION_ALIASES,
            override=load_request.source_version,
            field_name="kinase_library.source_version",
            required=True,
        )
        retrieved_at = _resolve_metadata_text(
            source,
            aliases=_RETRIEVED_AT_ALIASES,
            override=load_request.retrieved_at,
            field_name="kinase_library.retrieved_at",
            required=False,
        )
        license_name = _resolve_metadata_text(
            source,
            aliases=_LICENSE_ALIASES,
            override=load_request.license,
            field_name="kinase_library.license",
            required=True,
        )
        score_scale = _resolve_metadata_text(
            source,
            aliases=_SCORE_SCALE_ALIASES,
            override=load_request.score_scale,
            field_name="kinase_library.score_scale",
            required=True,
        )
        if source_name is None:
            raise ReferenceValidationError("kinase_library.source_name is required")
        if source_version is None:
            raise ReferenceValidationError("kinase_library.source_version is required")
        if license_name is None:
            raise ReferenceValidationError("kinase_library.license is required")
        if score_scale is None:
            raise ReferenceValidationError("kinase_library.score_scale is required")
        organisms = _resolve_organisms(
            source,
            override=load_request.organisms,
        )
        sequence_window = _resolve_sequence_window(
            source,
            positions=long_table.loc[:, "position"],
            override=load_request.sequence_window,
        )
        matrices = _build_matrices(long_table)
        provenance = _build_resource_provenance(
            path=path,
            matrices=matrices,
            source_name=source_name,
            source_version=source_version,
            retrieved_at=retrieved_at,
            license_name=license_name,
            score_scale=score_scale,
            sequence_window=sequence_window,
            organisms=organisms,
        )
        return KinaseLibraryResource(
            matrices=matrices,
            source_name=source_name,
            source_version=source_version,
            retrieved_at=retrieved_at,
            license=license_name,
            score_scale=score_scale,
            sequence_window=sequence_window,
            organisms=organisms,
            provenance=provenance,
        )


def load_kinase_library_resource(
    request: KinaseLibraryResourceLoadRequest | KinaseLibraryPath,
) -> KinaseLibraryResource:
    """Convenience wrapper for ``KinaseLibraryResourceLoader().run(...)``."""

    return KinaseLibraryResourceLoader().run(request)


def _coerce_load_request(
    request: KinaseLibraryResourceLoadRequest | KinaseLibraryPath,
) -> KinaseLibraryResourceLoadRequest:
    if isinstance(request, KinaseLibraryResourceLoadRequest):
        return request
    return KinaseLibraryResourceLoadRequest(path=request)


def _normalise_source_table(frame: pd.DataFrame) -> pd.DataFrame:
    if frame.empty:
        raise ReferenceValidationError("kinase library source table must be non-empty")
    if all(
        _find_column(frame, aliases) is not None for aliases in _LONG_REQUIRED_ALIASES
    ):
        return _normalise_long_source_table(frame)
    return _normalise_wide_source_table(frame)


def _normalise_long_source_table(frame: pd.DataFrame) -> pd.DataFrame:
    columns = {
        "kinase": _require_column(frame, _KINASE_ALIASES, semantic_name="kinase"),
        "residue_class": _require_column(
            frame,
            _RESIDUE_CLASS_ALIASES,
            semantic_name="residue_class",
        ),
        "position": _require_column(frame, _POSITION_ALIASES, semantic_name="position"),
        "amino_acid": _require_column(
            frame,
            _AMINO_ACID_ALIASES,
            semantic_name="amino_acid",
        ),
        "score": _require_column(frame, _SCORE_ALIASES, semantic_name="score"),
    }
    family_column = _find_column(frame, _KINASE_FAMILY_ALIASES)
    group_column = _find_column(frame, _KINASE_GROUP_ALIASES)
    normalized = pd.DataFrame(
        {
            "kinase": frame.loc[:, columns["kinase"]],
            "residue_class": frame.loc[:, columns["residue_class"]],
            "position": frame.loc[:, columns["position"]],
            "amino_acid": frame.loc[:, columns["amino_acid"]],
            "score": frame.loc[:, columns["score"]],
        }
    )
    normalized.loc[:, "kinase_family"] = (
        frame.loc[:, family_column] if family_column is not None else ""
    )
    normalized.loc[:, "kinase_group"] = (
        frame.loc[:, group_column] if group_column is not None else ""
    )
    return _clean_long_table(normalized)


def _normalise_wide_source_table(frame: pd.DataFrame) -> pd.DataFrame:
    kinase_column = _require_column(frame, _KINASE_ALIASES, semantic_name="kinase")
    residue_class_column = _require_column(
        frame,
        _RESIDUE_CLASS_ALIASES,
        semantic_name="residue_class",
    )
    amino_acid_column = _require_column(
        frame,
        _AMINO_ACID_ALIASES,
        semantic_name="amino_acid",
    )
    family_column = _find_column(frame, _KINASE_FAMILY_ALIASES)
    group_column = _find_column(frame, _KINASE_GROUP_ALIASES)
    position_columns = _find_wide_position_columns(frame)
    if not position_columns:
        raise ReferenceResolutionError(
            "kinase library source table must contain long columns "
            "(kinase, residue_class, position, amino_acid, score) or wide "
            "position columns such as p-1, p0, p1"
        )
    working = frame.copy(deep=True)
    if family_column is None:
        working.loc[:, "kinase_family"] = ""
        family_column = "kinase_family"
    if group_column is None:
        working.loc[:, "kinase_group"] = ""
        group_column = "kinase_group"
    melted = working.melt(
        id_vars=[
            kinase_column,
            residue_class_column,
            amino_acid_column,
            family_column,
            group_column,
        ],
        value_vars=[column for column, _ in position_columns],
        var_name="position",
        value_name="score",
    )
    position_lookup = {str(column): position for column, position in position_columns}
    normalized = pd.DataFrame(
        {
            "kinase": melted.loc[:, kinase_column],
            "residue_class": melted.loc[:, residue_class_column],
            "position": melted.loc[:, "position"].map(
                lambda value: position_lookup[str(value)]
            ),
            "amino_acid": melted.loc[:, amino_acid_column],
            "score": melted.loc[:, "score"],
            "kinase_family": melted.loc[:, family_column],
            "kinase_group": melted.loc[:, group_column],
        }
    )
    return _clean_long_table(normalized)


def _clean_long_table(frame: pd.DataFrame) -> pd.DataFrame:
    cleaned = pd.DataFrame(index=frame.index.copy())
    cleaned.loc[:, "kinase"] = [
        _require_text(value, field_name="kinase_library.kinase")
        for value in frame.loc[:, "kinase"].tolist()
    ]
    cleaned.loc[:, "residue_class"] = [
        _coerce_residue_class(
            value,
            field_name="kinase_library.residue_class",
        )
        for value in frame.loc[:, "residue_class"].tolist()
    ]
    cleaned.loc[:, "position"] = [
        _coerce_position(value, field_name="kinase_library.position")
        for value in frame.loc[:, "position"].tolist()
    ]
    cleaned.loc[:, "amino_acid"] = [
        _coerce_amino_acid(value, field_name="kinase_library.amino_acid")
        for value in frame.loc[:, "amino_acid"].tolist()
    ]
    cleaned.loc[:, "kinase_family"] = [
        _optional_text(value, field_name="kinase_library.kinase_family") or ""
        for value in frame.loc[:, "kinase_family"].tolist()
    ]
    cleaned.loc[:, "kinase_group"] = [
        _optional_text(value, field_name="kinase_library.kinase_group") or ""
        for value in frame.loc[:, "kinase_group"].tolist()
    ]
    try:
        cleaned.loc[:, "score"] = pd.to_numeric(frame.loc[:, "score"], errors="raise")
    except ValueError as exc:
        raise ReferenceValidationError(
            "kinase library score values must be numeric"
        ) from exc
    duplicated = cleaned.duplicated(
        subset=["kinase", "residue_class", "amino_acid", "position"],
        keep=False,
    )
    if duplicated.any():
        duplicate_rows = ", ".join(
            str(index) for index in cleaned.index[duplicated].tolist()[:5]
        )
        raise ReferenceValidationError(
            "kinase library source contains duplicate kinase/residue_class/"
            f"amino_acid/position rows at source row positions: {duplicate_rows}"
        )
    return cleaned


def _build_matrices(long_table: pd.DataFrame) -> tuple[KinaseLibraryMatrix, ...]:
    matrices: list[KinaseLibraryMatrix] = []
    grouped = long_table.groupby(["kinase", "residue_class"], sort=False)
    for (kinase, residue_class), group in grouped:
        family = _single_optional_metadata_value(
            group.loc[:, "kinase_family"],
            field_name=f"kinase_library[{kinase}].kinase_family",
        )
        kinase_group = _single_optional_metadata_value(
            group.loc[:, "kinase_group"],
            field_name=f"kinase_library[{kinase}].kinase_group",
        )
        score_table = group.pivot(
            index="amino_acid",
            columns="position",
            values="score",
        )
        score_table = score_table.loc[
            [
                amino_acid
                for amino_acid in AMINO_ACIDS
                if amino_acid in score_table.index
            ]
        ]
        score_table = score_table.loc[:, sorted(score_table.columns.tolist())]
        matrices.append(
            KinaseLibraryMatrix(
                kinase=str(kinase),
                kinase_family=family,
                kinase_group=kinase_group,
                residue_class=cast(KinaseLibraryResidueClass, residue_class),
                score_table=score_table,
                _assume_owned=True,
            )
        )
    return tuple(matrices)


def _build_resource_provenance(
    *,
    path: Path,
    matrices: tuple[KinaseLibraryMatrix, ...],
    source_name: str,
    source_version: str,
    retrieved_at: str | None,
    license_name: str,
    score_scale: str,
    sequence_window: SequenceWindowDefinition,
    organisms: tuple[str, ...],
) -> KinaseLibraryResourceProvenance:
    source_files: dict[str, JsonValue] = {
        "kinase_library": fingerprint_local_reference_source_file(
            path,
            role="kinase_library",
        )
    }
    manifest: dict[str, JsonValue] = {
        "resource_type": "kinase_library",
        "source_name": source_name,
        "source_version": source_version,
        "license": license_name,
        "score_scale": score_scale,
        "organisms": organisms,
        "sequence_window": sequence_window.to_payload(),
        "source_files": source_files,
    }
    if retrieved_at is not None:
        manifest["retrieved_at"] = retrieved_at
    table_fingerprints = (
        fingerprint_table(
            _matrix_inventory_frame(matrices),
            name="references.kinase_library.matrix_index",
        ),
        *(
            fingerprint_table(
                matrix.score_table,
                name=(
                    "references.kinase_library.score_table."
                    f"{_safe_table_name(matrix.kinase)}.{matrix.residue_class.value}"
                ),
            )
            for matrix in matrices
        ),
    )
    return KinaseLibraryResourceProvenance(
        source_type="local",
        source_name=source_name,
        source_version=source_version,
        retrieved_at=retrieved_at,
        license=license_name,
        score_scale=score_scale,
        organisms=organisms,
        sequence_window=sequence_window.to_payload(),
        source_files=source_files,
        manifest=manifest,
        table_fingerprints=table_fingerprints,
    )


def _matrix_inventory_frame(
    matrices: tuple[KinaseLibraryMatrix, ...],
) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "kinase": matrix.kinase,
                "residue_class": matrix.residue_class.value,
                "kinase_family": matrix.kinase_family,
                "kinase_group": matrix.kinase_group,
                "amino_acids": "|".join(str(item) for item in matrix.score_table.index),
                "positions": "|".join(
                    str(item) for item in matrix.score_table.columns.tolist()
                ),
                "rows": int(matrix.score_table.shape[0]),
                "columns": int(matrix.score_table.shape[1]),
            }
            for matrix in matrices
        ],
        columns=[
            "kinase",
            "residue_class",
            "kinase_family",
            "kinase_group",
            "amino_acids",
            "positions",
            "rows",
            "columns",
        ],
    )


def _resolve_sequence_window(
    frame: pd.DataFrame,
    *,
    positions: pd.Series,
    override: SequenceWindowDefinition | None,
) -> SequenceWindowDefinition:
    if override is not None:
        return override
    upstream_column = _find_column(frame, _UPSTREAM_ALIASES)
    downstream_column = _find_column(frame, _DOWNSTREAM_ALIASES)
    central_column = _find_column(frame, _CENTRAL_REQUIRED_ALIASES)
    if upstream_column is not None and downstream_column is not None:
        upstream = _coerce_metadata_int(
            _single_metadata_value(
                frame.loc[:, upstream_column],
                field_name="kinase_library.upstream_residues",
            ),
            field_name="kinase_library.upstream_residues",
        )
        downstream = _coerce_metadata_int(
            _single_metadata_value(
                frame.loc[:, downstream_column],
                field_name="kinase_library.downstream_residues",
            ),
            field_name="kinase_library.downstream_residues",
        )
        central_required = True
        if central_column is not None:
            central_required = _coerce_metadata_bool(
                _single_metadata_value(
                    frame.loc[:, central_column],
                    field_name="kinase_library.central_residue_required",
                ),
                field_name="kinase_library.central_residue_required",
            )
        return SequenceWindowDefinition(
            upstream_residues=upstream,
            downstream_residues=downstream,
            central_residue_required=central_required,
        )
    position_values = sorted({int(position) for position in positions.tolist()})
    if not position_values:
        raise ReferenceValidationError(
            "kinase library sequence_window cannot be inferred from empty positions"
        )
    min_position = min(position_values)
    max_position = max(position_values)
    if min_position > 0 or max_position < 0:
        raise ReferenceValidationError(
            "kinase library positions must include position 0 or explicit "
            "sequence_window metadata"
        )
    expected = list(range(min_position, max_position + 1))
    if position_values != expected:
        missing = sorted(set(expected).difference(position_values))
        raise ReferenceValidationError(
            "kinase library positions must be contiguous when sequence_window "
            f"is inferred; missing positions: {_format_positions(missing)}"
        )
    return SequenceWindowDefinition(
        upstream_residues=abs(min_position),
        downstream_residues=max_position,
        central_residue_required=True,
    )


def _resolve_metadata_text(
    frame: pd.DataFrame,
    *,
    aliases: tuple[str, ...],
    override: object,
    field_name: str,
    required: bool,
) -> str | None:
    if override is not None:
        if isinstance(override, date):
            return override.isoformat()
        if required:
            return _require_text(override, field_name=field_name)
        return _optional_text(override, field_name=field_name)
    column = _find_column(frame, aliases)
    if column is None:
        if required:
            raise ReferenceValidationError(f"{field_name} is required")
        return None
    value = _single_metadata_value(frame.loc[:, column], field_name=field_name)
    if required:
        return _require_text(value, field_name=field_name)
    return _optional_text(value, field_name=field_name)


def _resolve_organisms(
    frame: pd.DataFrame,
    *,
    override: tuple[Organism | str, ...],
) -> tuple[str, ...]:
    if override:
        return _coerce_organism_tokens(override, field_name="kinase_library.organisms")
    column = _find_column(frame, _ORGANISMS_ALIASES)
    if column is None:
        raise ReferenceValidationError("kinase_library.organisms is required")
    tokens: list[str] = []
    for value in _unique_nonblank_values(
        frame.loc[:, column],
        field_name="kinase_library.organisms",
    ):
        tokens.extend(_split_organism_value(value))
    return _coerce_organism_tokens(tuple(tokens), field_name="kinase_library.organisms")


def _single_metadata_value(series: pd.Series, *, field_name: str) -> str:
    values = _unique_nonblank_values(series, field_name=field_name)
    if not values:
        raise ReferenceValidationError(f"{field_name} is required")
    if len(values) > 1:
        preview = ", ".join(repr(value) for value in values[:5])
        raise ReferenceValidationError(
            f"{field_name} must have one value across the kinase library resource; "
            f"observed values: {preview}"
        )
    return values[0]


def _single_optional_metadata_value(
    series: pd.Series, *, field_name: str
) -> str | None:
    values = _unique_nonblank_values(series, field_name=field_name)
    if not values:
        return None
    if len(values) > 1:
        preview = ", ".join(repr(value) for value in values[:5])
        raise ReferenceValidationError(
            f"{field_name} must have one value per kinase matrix; observed values: "
            f"{preview}"
        )
    return values[0]


def _unique_nonblank_values(series: pd.Series, *, field_name: str) -> list[str]:
    values: list[str] = []
    for value in series.tolist():
        text = _optional_text(value, field_name=field_name)
        if text is None:
            continue
        if text not in values:
            values.append(text)
    return values


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
            f"{field_name} must be one canonical one-letter amino-acid code"
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


def _coerce_metadata_int(value: object, *, field_name: str) -> int:
    resolved = _coerce_position(value, field_name=field_name)
    if resolved < 0:
        raise ReferenceValidationError(f"{field_name} must be >= 0")
    return resolved


def _coerce_metadata_bool(value: object, *, field_name: str) -> bool:
    if isinstance(value, bool):
        return bool(value)
    text = _require_text(value, field_name=field_name).lower()
    if text in _TRUE_TOKENS:
        return True
    if text in _FALSE_TOKENS:
        return False
    raise ReferenceValidationError(f"{field_name} must be boolean")


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


def _require_column(
    frame: pd.DataFrame,
    aliases: tuple[str, ...],
    *,
    semantic_name: str,
) -> str:
    column = _find_column(frame, aliases)
    if column is not None:
        return column
    accepted = ", ".join(aliases)
    available = ", ".join(str(item) for item in frame.columns.tolist())
    if not available:
        available = "(none)"
    raise ReferenceResolutionError(
        f"kinase library source is missing required {semantic_name} column; "
        f"accepted aliases: {accepted}; available columns: {available}"
    )


def _find_column(frame: pd.DataFrame, aliases: tuple[str, ...]) -> str | None:
    keyed_columns = {
        _normalise_column_token(column): str(column)
        for column in frame.columns.tolist()
    }
    for alias in aliases:
        column = keyed_columns.get(alias)
        if column is not None:
            return column
    return None


def _find_wide_position_columns(frame: pd.DataFrame) -> list[tuple[str, int]]:
    metadata_tokens = {
        _normalise_column_token(alias)
        for aliases in _METADATA_ALIAS_GROUPS
        for alias in aliases
    }
    columns: list[tuple[str, int]] = []
    for column in frame.columns.tolist():
        token = _normalise_column_token(column)
        if token in metadata_tokens:
            continue
        match = _POSITION_COLUMN_PATTERN.match(token)
        if match is None:
            continue
        columns.append((str(column), int(match.group(1))))
    return columns


def _normalise_column_token(column: object) -> str:
    raw = str(column).strip().lower()
    return _COLUMN_TOKEN_PATTERN.sub("_", raw).strip("_")


def _safe_table_name(value: str) -> str:
    return _COLUMN_TOKEN_PATTERN.sub("_", value.lower()).strip("_") or "unknown"


def _format_positions(positions: list[int]) -> str:
    if not positions:
        return "(none)"
    return ", ".join(str(position) for position in positions)


__all__ = [
    "KinaseLibraryMatrix",
    "KinaseLibraryResidueClass",
    "KinaseLibraryResource",
    "KinaseLibraryResourceLoadRequest",
    "KinaseLibraryResourceLoader",
    "load_kinase_library_resource",
]
