"""Complete kinase-scoring state codec shared by reloadable bundles."""

from __future__ import annotations

from collections import Counter
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import TypeGuard, cast

import pandas as pd

from phospy.core.errors.input import PhosPyInputError
from phospy.core.errors.validation import PhosPyValidationError
from phospy.io.bundles._shared.primitives import (
    require_int,
    require_mapping,
    require_str,
    validate_json_safe_mapping,
)
from phospy.io.bundles._shared.tables import (
    read_optional_table,
    read_required_table,
    write_bundle_table,
    write_optional_bundle_table,
)
from phospy.kinase.scoring.motif import (
    MotifLibraryValidationResult,
    MotifLibraryValidationRow,
)
from phospy.kinase.scoring.results import KinaseScoringResult
from phospy.kinase.scoring.sequence_validation import (
    SEQUENCE_VALIDATION_STATUS_INVALID_SITE_ID,
    SEQUENCE_VALIDATION_STATUS_MISSING_SEQUENCE,
    SEQUENCE_VALIDATION_STATUS_NON_PHOSPHO_CENTRE_RESIDUE,
    SEQUENCE_VALIDATION_STATUS_OFF_CENTRE_SEQUENCE,
    SEQUENCE_VALIDATION_STATUS_SHORT_SEQUENCE,
    SEQUENCE_VALIDATION_STATUS_SITE_RESIDUE_MISMATCH,
    SEQUENCE_VALIDATION_STATUS_UNSUPPORTED_RESIDUE_CHARACTER,
    SEQUENCE_VALIDATION_STATUS_VALID,
    SequenceValidationResult,
    SequenceValidationRow,
    SequenceValidationStatus,
)
from phospy.provenance.immutability import thaw_json_mapping


@dataclass(frozen=True, slots=True)
class _ScoringTableSpec:
    key: str
    required: bool = False
    workflow_level: bool = False


_SCORING_TABLE_SPECS = (
    _ScoringTableSpec("profile_scores", required=True),
    _ScoringTableSpec("motif_scores"),
    _ScoringTableSpec("rank_weighted_fusion_scores"),
    _ScoringTableSpec("kinase_library_motif_scores"),
    _ScoringTableSpec("combined_profile_motif_scores"),
    _ScoringTableSpec("score_fusion_weights"),
    _ScoringTableSpec("score_source_matrix"),
    _ScoringTableSpec("score_source_summary"),
    _ScoringTableSpec("profile_score_diagnostics"),
    _ScoringTableSpec("kinase_library_site_diagnostics"),
    _ScoringTableSpec("kinase_library_kinase_diagnostics"),
    _ScoringTableSpec("substrate_contributions", workflow_level=True),
)
KINASE_SCORING_STATE_FIELDS = frozenset({"metadata", "table_metadata", "tables"})
KINASE_SCORING_TABLE_KEYS = frozenset(spec.key for spec in _SCORING_TABLE_SPECS)
KINASE_SCORING_OPTIONAL_TABLE_KEYS = frozenset(
    spec.key for spec in _SCORING_TABLE_SPECS if not spec.required
)

_SCORING_METADATA_FIELDS = frozenset(
    {
        "scoring_mode",
        "score_source",
        "score_scale",
        "score_scale_metadata",
        "profile_self_inclusion_policy",
        "motif_sequence_validation",
        "motif_library_validation",
    }
)
_SEQUENCE_VALIDATION_FIELDS = frozenset(
    {
        "total_sequences",
        "valid_sequences",
        "invalid_sequences",
        "short_sequences",
        "off_centre_sequences",
        "site_residue_mismatches",
        "invalid_site_ids",
        "non_phospho_centre_residues",
        "unsupported_residue_characters",
        "sequences_excluded_from_motif_scoring",
        "excluded_site_ids",
        "rows",
    }
)
_SEQUENCE_ROW_FIELDS = frozenset(
    {
        "site_id",
        "sequence",
        "status",
        "reason",
        "expected_centre_residue",
        "observed_centre_residue",
        "sequence_length",
    }
)
_MOTIF_LIBRARY_VALIDATION_FIELDS = frozenset(
    {
        "total_reference_sequences",
        "accepted_reference_sequences",
        "excluded_reference_sequences",
        "missing_sequences",
        "short_sequences",
        "off_centre_sequences",
        "site_residue_mismatches",
        "invalid_site_ids",
        "non_phospho_centre_residues",
        "unsupported_residue_characters",
        "sequences_excluded_from_motif_profile_construction",
        "expected_window_size",
        "supported_amino_acids",
        "accepted_window_length_policy",
        "unsupported_residue_policy",
        "excluded_reference_ids",
        "rows",
    }
)
_MOTIF_LIBRARY_ROW_FIELDS = frozenset(
    {
        "reference_id",
        "site_id",
        "kinase",
        "sequence",
        "status",
        "reason",
        "expected_centre_residue",
        "observed_centre_residue",
        "sequence_length",
    }
)
_SEQUENCE_STATUSES = frozenset(
    {
        SEQUENCE_VALIDATION_STATUS_VALID,
        SEQUENCE_VALIDATION_STATUS_MISSING_SEQUENCE,
        SEQUENCE_VALIDATION_STATUS_SHORT_SEQUENCE,
        SEQUENCE_VALIDATION_STATUS_OFF_CENTRE_SEQUENCE,
        SEQUENCE_VALIDATION_STATUS_SITE_RESIDUE_MISMATCH,
        SEQUENCE_VALIDATION_STATUS_INVALID_SITE_ID,
        SEQUENCE_VALIDATION_STATUS_NON_PHOSPHO_CENTRE_RESIDUE,
        SEQUENCE_VALIDATION_STATUS_UNSUPPORTED_RESIDUE_CHARACTER,
    }
)


@dataclass(frozen=True, slots=True)
class KinaseScoringStateMetadata:
    """Decoded non-tabular scoring semantics."""

    scoring_mode: str
    score_source: str
    score_scale: str
    score_scale_metadata: Mapping[str, object] | None
    profile_self_inclusion_policy: str
    motif_sequence_validation: SequenceValidationResult | None
    motif_library_validation: MotifLibraryValidationResult | None


@dataclass(frozen=True, slots=True)
class KinaseScoringStateManifest:
    """Parsed scoring subsection before bundle file reconstruction."""

    metadata: KinaseScoringStateMetadata
    table_metadata: Mapping[str, ScoringTableStorageMetadata | None]
    tables: Mapping[str, object]


@dataclass(frozen=True, slots=True)
class ScoringTableStorageMetadata:
    """Storage details CSV cannot preserve without manifest metadata."""

    index_name: str | None
    columns_name: str | None
    index_dtype: str
    column_dtypes: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class ReconstructedKinaseScoringState:
    """Reconstructed scoring result and workflow-level contribution table."""

    result: KinaseScoringResult
    substrate_contributions: pd.DataFrame | None


def write_kinase_scoring_tables(
    *,
    result: KinaseScoringResult,
    substrate_contributions: pd.DataFrame | None,
    bundle_root: Path,
    suffix: str,
    written: dict[str, Path],
) -> dict[str, object]:
    """Write every durable scoring table with explicit null entries."""

    table_values = _scoring_table_values(
        result=result,
        substrate_contributions=substrate_contributions,
    )
    entries: dict[str, object] = {}
    for spec, (_, table) in zip(_SCORING_TABLE_SPECS, table_values, strict=True):
        table_key = spec.key
        relative_path = Path("scoring") / f"{table_key}{suffix}"
        written_key = f"scoring.{table_key}"
        if spec.required:
            if table is None:
                raise PhosPyInputError("scoring_result.profile_scores is required")
            entries[table_key] = write_bundle_table(
                table=table,
                bundle_root=bundle_root,
                relative_path=relative_path,
                written=written,
                written_key=written_key,
            )
        else:
            entries[table_key] = write_optional_bundle_table(
                table=table,
                bundle_root=bundle_root,
                relative_path=relative_path,
                written=written,
                written_key=written_key,
            )
    return entries


def kinase_scoring_state_to_payload(
    *,
    result: KinaseScoringResult,
    substrate_contributions: pd.DataFrame | None,
    tables: Mapping[str, object],
) -> dict[str, object]:
    """Encode one complete scoring subsection for a bundle manifest."""

    return {
        "metadata": _scoring_metadata_to_payload(result),
        "table_metadata": _table_metadata_to_payload(
            _scoring_table_values(
                result=result,
                substrate_contributions=substrate_contributions,
            )
        ),
        "tables": dict(tables),
    }


def parse_kinase_scoring_state(
    value: object,
    *,
    field_name: str,
) -> KinaseScoringStateManifest:
    """Parse and validate a complete scoring manifest subsection."""

    payload = require_mapping(value, field_name=field_name)
    _require_exact_fields(
        payload, field_name=field_name, fields=KINASE_SCORING_STATE_FIELDS
    )
    metadata_payload = require_mapping(
        payload.get("metadata"), field_name=f"{field_name}.metadata"
    )
    tables = require_mapping(payload.get("tables"), field_name=f"{field_name}.tables")
    _require_exact_fields(
        tables,
        field_name=f"{field_name}.tables",
        fields=KINASE_SCORING_TABLE_KEYS,
    )
    return KinaseScoringStateManifest(
        metadata=_scoring_metadata_from_payload(
            metadata_payload,
            field_name=f"{field_name}.metadata",
        ),
        table_metadata=_table_metadata_from_payload(
            payload.get("table_metadata"),
            field_name=f"{field_name}.table_metadata",
        ),
        tables=tables,
    )


def reconstruct_kinase_scoring_state(
    *,
    bundle_root: Path,
    state: KinaseScoringStateManifest,
    field_name: str,
) -> ReconstructedKinaseScoringState:
    """Read and validate all tables before reconstructing the scoring model."""

    tables = state.tables
    loaded: dict[str, pd.DataFrame | None] = {}
    for table_key in sorted(KINASE_SCORING_TABLE_KEYS):
        table_field = f"{field_name}.tables.{table_key}"
        if table_key == "profile_scores":
            loaded[table_key] = read_required_table(
                bundle_root=bundle_root,
                tables=tables,
                table_key=table_key,
                field_name=table_field,
            )
        else:
            loaded[table_key] = read_optional_table(
                bundle_root=bundle_root,
                tables=tables,
                table_key=table_key,
                field_name=table_field,
            )
        loaded[table_key] = _restore_table_metadata(
            loaded[table_key],
            state.table_metadata[table_key],
            field_name=f"{field_name}.table_metadata.{table_key}",
        )
    metadata = state.metadata
    try:
        scoring_result = KinaseScoringResult(
            profile_scores=_required_loaded_table(loaded, "profile_scores"),
            motif_scores=loaded["motif_scores"],
            rank_weighted_fusion_scores=loaded["rank_weighted_fusion_scores"],
            kinase_library_motif_scores=loaded["kinase_library_motif_scores"],
            combined_profile_motif_scores=loaded["combined_profile_motif_scores"],
            score_fusion_weights=loaded["score_fusion_weights"],
            score_source_matrix=loaded["score_source_matrix"],
            score_source_summary=loaded["score_source_summary"],
            profile_score_diagnostics=loaded["profile_score_diagnostics"],
            kinase_library_site_diagnostics=loaded["kinase_library_site_diagnostics"],
            kinase_library_kinase_diagnostics=loaded[
                "kinase_library_kinase_diagnostics"
            ],
            scoring_mode=metadata.scoring_mode,
            score_source=metadata.score_source,
            score_scale=metadata.score_scale,
            score_scale_metadata=metadata.score_scale_metadata,
            profile_self_inclusion_policy=metadata.profile_self_inclusion_policy,
            motif_sequence_validation=metadata.motif_sequence_validation,
            motif_library_validation=metadata.motif_library_validation,
        )
    except (PhosPyValidationError, ValueError, TypeError) as exc:
        raise PhosPyInputError(f"{field_name} is invalid: {exc}") from exc
    return ReconstructedKinaseScoringState(
        result=scoring_result,
        substrate_contributions=loaded["substrate_contributions"],
    )


def _scoring_metadata_to_payload(result: KinaseScoringResult) -> dict[str, object]:
    return {
        "scoring_mode": result.scoring_mode,
        "score_source": result.score_source,
        "score_scale": result.score_scale,
        "score_scale_metadata": (
            None
            if result.score_scale_metadata is None
            else thaw_json_mapping(
                result.score_scale_metadata,
                field_name="scoring_result.score_scale_metadata",
            )
        ),
        "profile_self_inclusion_policy": str(result.profile_self_inclusion_policy),
        "motif_sequence_validation": _sequence_validation_to_payload(
            result.motif_sequence_validation
        ),
        "motif_library_validation": _motif_library_validation_to_payload(
            result.motif_library_validation
        ),
    }


def _scoring_table_values(
    *,
    result: KinaseScoringResult,
    substrate_contributions: pd.DataFrame | None,
) -> tuple[tuple[str, pd.DataFrame | None], ...]:
    return tuple(
        (
            spec.key,
            (
                substrate_contributions
                if spec.workflow_level
                else cast(pd.DataFrame | None, getattr(result, spec.key))
            ),
        )
        for spec in _SCORING_TABLE_SPECS
    )


def _table_metadata_to_payload(
    table_values: tuple[tuple[str, pd.DataFrame | None], ...],
) -> dict[str, object]:
    return {
        table_key: (
            None
            if table is None
            else {
                "index_name": _axis_name_to_payload(table.index.name),
                "columns_name": _axis_name_to_payload(table.columns.name),
                "index_dtype": str(table.index.dtype),
                "column_dtypes": [str(dtype) for dtype in table.dtypes],
            }
        )
        for table_key, table in table_values
    }


def _table_metadata_from_payload(
    value: object,
    *,
    field_name: str,
) -> dict[str, ScoringTableStorageMetadata | None]:
    payload = require_mapping(value, field_name=field_name)
    _require_exact_fields(
        payload,
        field_name=field_name,
        fields=KINASE_SCORING_TABLE_KEYS,
    )
    parsed: dict[str, ScoringTableStorageMetadata | None] = {}
    for table_key in sorted(KINASE_SCORING_TABLE_KEYS):
        raw_metadata = payload[table_key]
        if raw_metadata is None:
            parsed[table_key] = None
            continue
        metadata = require_mapping(raw_metadata, field_name=f"{field_name}.{table_key}")
        _require_exact_fields(
            metadata,
            field_name=f"{field_name}.{table_key}",
            fields=frozenset(
                {"index_name", "columns_name", "index_dtype", "column_dtypes"}
            ),
        )
        parsed[table_key] = ScoringTableStorageMetadata(
            index_name=_optional_string(
                metadata.get("index_name"),
                field_name=f"{field_name}.{table_key}.index_name",
            ),
            columns_name=_optional_string(
                metadata.get("columns_name"),
                field_name=f"{field_name}.{table_key}.columns_name",
            ),
            index_dtype=require_str(
                metadata.get("index_dtype"),
                field_name=f"{field_name}.{table_key}.index_dtype",
            ),
            column_dtypes=_string_tuple(
                metadata.get("column_dtypes"),
                field_name=f"{field_name}.{table_key}.column_dtypes",
            ),
        )
    return parsed


def _axis_name_to_payload(value: object) -> str | None:
    if value is None:
        return None
    if not isinstance(value, str):
        raise PhosPyInputError("scoring table axis names must be strings or null")
    return value


def _restore_table_metadata(
    table: pd.DataFrame | None,
    metadata: ScoringTableStorageMetadata | None,
    *,
    field_name: str,
) -> pd.DataFrame | None:
    if table is None:
        if metadata is not None:
            raise PhosPyInputError(
                f"{field_name} must be null when the table entry is null"
            )
        return None
    if metadata is None:
        raise PhosPyInputError(
            f"{field_name} must be populated when the table entry is populated"
        )
    if len(metadata.column_dtypes) != len(table.columns):
        raise PhosPyInputError(
            f"{field_name}.column_dtypes must contain one entry per table column"
        )
    try:
        restored = table.astype(
            dict(zip(table.columns, metadata.column_dtypes, strict=True))
        )
        restored.index = restored.index.astype(metadata.index_dtype)
    except (KeyError, TypeError, ValueError) as exc:
        raise PhosPyInputError(f"{field_name} contains invalid dtype metadata") from exc
    restored.index.name = metadata.index_name
    restored.columns.name = metadata.columns_name
    return restored


def _scoring_metadata_from_payload(
    payload: Mapping[str, object],
    *,
    field_name: str,
) -> KinaseScoringStateMetadata:
    _require_exact_fields(
        payload, field_name=field_name, fields=_SCORING_METADATA_FIELDS
    )
    raw_scale_metadata = payload.get("score_scale_metadata")
    return KinaseScoringStateMetadata(
        scoring_mode=require_str(
            payload.get("scoring_mode"), field_name=f"{field_name}.scoring_mode"
        ),
        score_source=require_str(
            payload.get("score_source"), field_name=f"{field_name}.score_source"
        ),
        score_scale=require_str(
            payload.get("score_scale"), field_name=f"{field_name}.score_scale"
        ),
        score_scale_metadata=(
            None
            if raw_scale_metadata is None
            else validate_json_safe_mapping(
                raw_scale_metadata,
                field_name=f"{field_name}.score_scale_metadata",
            )
        ),
        profile_self_inclusion_policy=require_str(
            payload.get("profile_self_inclusion_policy"),
            field_name=f"{field_name}.profile_self_inclusion_policy",
        ),
        motif_sequence_validation=_sequence_validation_from_payload(
            payload.get("motif_sequence_validation"),
            field_name=f"{field_name}.motif_sequence_validation",
        ),
        motif_library_validation=_motif_library_validation_from_payload(
            payload.get("motif_library_validation"),
            field_name=f"{field_name}.motif_library_validation",
        ),
    )


def _sequence_validation_to_payload(
    value: SequenceValidationResult | None,
) -> dict[str, object] | None:
    if value is None:
        return None
    return {
        "total_sequences": value.total_sequences,
        "valid_sequences": value.valid_sequences,
        "invalid_sequences": value.invalid_sequences,
        "short_sequences": value.short_sequences,
        "off_centre_sequences": value.off_centre_sequences,
        "site_residue_mismatches": value.site_residue_mismatches,
        "invalid_site_ids": value.invalid_site_ids,
        "non_phospho_centre_residues": value.non_phospho_centre_residues,
        "unsupported_residue_characters": value.unsupported_residue_characters,
        "sequences_excluded_from_motif_scoring": (
            value.sequences_excluded_from_motif_scoring
        ),
        "excluded_site_ids": list(value.excluded_site_ids),
        "rows": [_sequence_row_to_payload(row) for row in value.rows],
    }


def _sequence_row_to_payload(value: SequenceValidationRow) -> dict[str, object]:
    return {
        "site_id": value.site_id,
        "sequence": value.sequence,
        "status": value.status,
        "reason": value.reason,
        "expected_centre_residue": value.expected_centre_residue,
        "observed_centre_residue": value.observed_centre_residue,
        "sequence_length": value.sequence_length,
    }


def _sequence_validation_from_payload(
    value: object,
    *,
    field_name: str,
) -> SequenceValidationResult | None:
    if value is None:
        return None
    payload = require_mapping(value, field_name=field_name)
    _require_exact_fields(
        payload, field_name=field_name, fields=_SEQUENCE_VALIDATION_FIELDS
    )
    rows = tuple(
        _sequence_row_from_payload(row, field_name=f"{field_name}.rows[{index}]")
        for index, row in enumerate(
            _require_sequence(payload.get("rows"), field_name=f"{field_name}.rows")
        )
    )
    result = SequenceValidationResult(
        total_sequences=_require_nonnegative_int(
            payload, "total_sequences", field_name
        ),
        valid_sequences=_require_nonnegative_int(
            payload, "valid_sequences", field_name
        ),
        invalid_sequences=_require_nonnegative_int(
            payload, "invalid_sequences", field_name
        ),
        short_sequences=_require_nonnegative_int(
            payload, "short_sequences", field_name
        ),
        off_centre_sequences=_require_nonnegative_int(
            payload, "off_centre_sequences", field_name
        ),
        site_residue_mismatches=_require_nonnegative_int(
            payload, "site_residue_mismatches", field_name
        ),
        invalid_site_ids=_require_nonnegative_int(
            payload, "invalid_site_ids", field_name
        ),
        non_phospho_centre_residues=_require_nonnegative_int(
            payload, "non_phospho_centre_residues", field_name
        ),
        unsupported_residue_characters=_require_nonnegative_int(
            payload, "unsupported_residue_characters", field_name
        ),
        sequences_excluded_from_motif_scoring=_require_nonnegative_int(
            payload, "sequences_excluded_from_motif_scoring", field_name
        ),
        excluded_site_ids=_string_tuple(
            payload.get("excluded_site_ids"),
            field_name=f"{field_name}.excluded_site_ids",
        ),
        rows=rows,
    )
    _validate_sequence_summary(result, field_name=field_name)
    return result


def _sequence_row_from_payload(
    value: object, *, field_name: str
) -> SequenceValidationRow:
    payload = require_mapping(value, field_name=field_name)
    _require_exact_fields(payload, field_name=field_name, fields=_SEQUENCE_ROW_FIELDS)
    status = _require_status(payload.get("status"), field_name=f"{field_name}.status")
    return SequenceValidationRow(
        site_id=require_str(payload.get("site_id"), field_name=f"{field_name}.site_id"),
        sequence=_optional_string(
            payload.get("sequence"), field_name=f"{field_name}.sequence"
        ),
        status=status,
        reason=_optional_string(
            payload.get("reason"), field_name=f"{field_name}.reason"
        ),
        expected_centre_residue=_optional_string(
            payload.get("expected_centre_residue"),
            field_name=f"{field_name}.expected_centre_residue",
        ),
        observed_centre_residue=_optional_string(
            payload.get("observed_centre_residue"),
            field_name=f"{field_name}.observed_centre_residue",
        ),
        sequence_length=_optional_nonnegative_int(
            payload.get("sequence_length"), field_name=f"{field_name}.sequence_length"
        ),
    )


def _motif_library_validation_to_payload(
    value: MotifLibraryValidationResult | None,
) -> dict[str, object] | None:
    if value is None:
        return None
    return {
        "total_reference_sequences": value.total_reference_sequences,
        "accepted_reference_sequences": value.accepted_reference_sequences,
        "excluded_reference_sequences": value.excluded_reference_sequences,
        "missing_sequences": value.missing_sequences,
        "short_sequences": value.short_sequences,
        "off_centre_sequences": value.off_centre_sequences,
        "site_residue_mismatches": value.site_residue_mismatches,
        "invalid_site_ids": value.invalid_site_ids,
        "non_phospho_centre_residues": value.non_phospho_centre_residues,
        "unsupported_residue_characters": value.unsupported_residue_characters,
        "sequences_excluded_from_motif_profile_construction": value.sequences_excluded_from_motif_profile_construction,
        "expected_window_size": value.expected_window_size,
        "supported_amino_acids": list(value.supported_amino_acids),
        "accepted_window_length_policy": value.accepted_window_length_policy,
        "unsupported_residue_policy": value.unsupported_residue_policy,
        "excluded_reference_ids": list(value.excluded_reference_ids),
        "rows": [_motif_library_row_to_payload(row) for row in value.rows],
    }


def _motif_library_row_to_payload(
    value: MotifLibraryValidationRow,
) -> dict[str, object]:
    return {
        "reference_id": value.reference_id,
        "site_id": value.site_id,
        "kinase": value.kinase,
        "sequence": value.sequence,
        "status": value.status,
        "reason": value.reason,
        "expected_centre_residue": value.expected_centre_residue,
        "observed_centre_residue": value.observed_centre_residue,
        "sequence_length": value.sequence_length,
    }


def _motif_library_validation_from_payload(
    value: object,
    *,
    field_name: str,
) -> MotifLibraryValidationResult | None:
    if value is None:
        return None
    payload = require_mapping(value, field_name=field_name)
    _require_exact_fields(
        payload, field_name=field_name, fields=_MOTIF_LIBRARY_VALIDATION_FIELDS
    )
    rows = tuple(
        _motif_library_row_from_payload(row, field_name=f"{field_name}.rows[{index}]")
        for index, row in enumerate(
            _require_sequence(payload.get("rows"), field_name=f"{field_name}.rows")
        )
    )
    result = MotifLibraryValidationResult(
        total_reference_sequences=_require_nonnegative_int(
            payload, "total_reference_sequences", field_name
        ),
        accepted_reference_sequences=_require_nonnegative_int(
            payload, "accepted_reference_sequences", field_name
        ),
        excluded_reference_sequences=_require_nonnegative_int(
            payload, "excluded_reference_sequences", field_name
        ),
        missing_sequences=_require_nonnegative_int(
            payload, "missing_sequences", field_name
        ),
        short_sequences=_require_nonnegative_int(
            payload, "short_sequences", field_name
        ),
        off_centre_sequences=_require_nonnegative_int(
            payload, "off_centre_sequences", field_name
        ),
        site_residue_mismatches=_require_nonnegative_int(
            payload, "site_residue_mismatches", field_name
        ),
        invalid_site_ids=_require_nonnegative_int(
            payload, "invalid_site_ids", field_name
        ),
        non_phospho_centre_residues=_require_nonnegative_int(
            payload, "non_phospho_centre_residues", field_name
        ),
        unsupported_residue_characters=_require_nonnegative_int(
            payload, "unsupported_residue_characters", field_name
        ),
        sequences_excluded_from_motif_profile_construction=_require_nonnegative_int(
            payload, "sequences_excluded_from_motif_profile_construction", field_name
        ),
        expected_window_size=_require_positive_int(
            payload, "expected_window_size", field_name
        ),
        supported_amino_acids=_string_tuple(
            payload.get("supported_amino_acids"),
            field_name=f"{field_name}.supported_amino_acids",
        ),
        accepted_window_length_policy=require_str(
            payload.get("accepted_window_length_policy"),
            field_name=f"{field_name}.accepted_window_length_policy",
        ),
        unsupported_residue_policy=require_str(
            payload.get("unsupported_residue_policy"),
            field_name=f"{field_name}.unsupported_residue_policy",
        ),
        excluded_reference_ids=_string_tuple(
            payload.get("excluded_reference_ids"),
            field_name=f"{field_name}.excluded_reference_ids",
        ),
        rows=rows,
    )
    _validate_motif_library_summary(result, field_name=field_name)
    return result


def _motif_library_row_from_payload(
    value: object, *, field_name: str
) -> MotifLibraryValidationRow:
    payload = require_mapping(value, field_name=field_name)
    _require_exact_fields(
        payload, field_name=field_name, fields=_MOTIF_LIBRARY_ROW_FIELDS
    )
    return MotifLibraryValidationRow(
        reference_id=require_str(
            payload.get("reference_id"), field_name=f"{field_name}.reference_id"
        ),
        site_id=_optional_string(
            payload.get("site_id"), field_name=f"{field_name}.site_id"
        ),
        kinase=require_str(payload.get("kinase"), field_name=f"{field_name}.kinase"),
        sequence=_optional_string(
            payload.get("sequence"), field_name=f"{field_name}.sequence"
        ),
        status=_require_status(
            payload.get("status"), field_name=f"{field_name}.status"
        ),
        reason=_optional_string(
            payload.get("reason"), field_name=f"{field_name}.reason"
        ),
        expected_centre_residue=_optional_string(
            payload.get("expected_centre_residue"),
            field_name=f"{field_name}.expected_centre_residue",
        ),
        observed_centre_residue=_optional_string(
            payload.get("observed_centre_residue"),
            field_name=f"{field_name}.observed_centre_residue",
        ),
        sequence_length=_optional_nonnegative_int(
            payload.get("sequence_length"), field_name=f"{field_name}.sequence_length"
        ),
    )


def _validate_sequence_summary(
    value: SequenceValidationResult, *, field_name: str
) -> None:
    counts = Counter(row.status for row in value.rows)
    excluded = tuple(
        row.site_id
        for row in value.rows
        if row.status != SEQUENCE_VALIDATION_STATUS_VALID
    )
    expected = {
        "total_sequences": len(value.rows),
        "valid_sequences": counts[SEQUENCE_VALIDATION_STATUS_VALID],
        "invalid_sequences": len(value.rows) - counts[SEQUENCE_VALIDATION_STATUS_VALID],
        "short_sequences": counts[SEQUENCE_VALIDATION_STATUS_SHORT_SEQUENCE],
        "off_centre_sequences": counts[SEQUENCE_VALIDATION_STATUS_OFF_CENTRE_SEQUENCE],
        "site_residue_mismatches": counts[
            SEQUENCE_VALIDATION_STATUS_SITE_RESIDUE_MISMATCH
        ],
        "invalid_site_ids": counts[SEQUENCE_VALIDATION_STATUS_INVALID_SITE_ID],
        "non_phospho_centre_residues": counts[
            SEQUENCE_VALIDATION_STATUS_NON_PHOSPHO_CENTRE_RESIDUE
        ],
        "unsupported_residue_characters": counts[
            SEQUENCE_VALIDATION_STATUS_UNSUPPORTED_RESIDUE_CHARACTER
        ],
        "sequences_excluded_from_motif_scoring": len(excluded),
    }
    _require_summary_counts(value, expected, field_name=field_name)
    if value.excluded_site_ids != excluded:
        raise PhosPyInputError(
            f"{field_name}.excluded_site_ids must match excluded rows"
        )


def _validate_motif_library_summary(
    value: MotifLibraryValidationResult, *, field_name: str
) -> None:
    counts = Counter(row.status for row in value.rows)
    excluded = tuple(
        row.reference_id
        for row in value.rows
        if row.status != SEQUENCE_VALIDATION_STATUS_VALID
    )
    expected = {
        "total_reference_sequences": len(value.rows),
        "accepted_reference_sequences": counts[SEQUENCE_VALIDATION_STATUS_VALID],
        "excluded_reference_sequences": len(excluded),
        "missing_sequences": counts[SEQUENCE_VALIDATION_STATUS_MISSING_SEQUENCE],
        "short_sequences": counts[SEQUENCE_VALIDATION_STATUS_SHORT_SEQUENCE],
        "off_centre_sequences": counts[SEQUENCE_VALIDATION_STATUS_OFF_CENTRE_SEQUENCE],
        "site_residue_mismatches": counts[
            SEQUENCE_VALIDATION_STATUS_SITE_RESIDUE_MISMATCH
        ],
        "invalid_site_ids": counts[SEQUENCE_VALIDATION_STATUS_INVALID_SITE_ID],
        "non_phospho_centre_residues": counts[
            SEQUENCE_VALIDATION_STATUS_NON_PHOSPHO_CENTRE_RESIDUE
        ],
        "unsupported_residue_characters": counts[
            SEQUENCE_VALIDATION_STATUS_UNSUPPORTED_RESIDUE_CHARACTER
        ],
        "sequences_excluded_from_motif_profile_construction": len(excluded),
    }
    _require_summary_counts(value, expected, field_name=field_name)
    if value.excluded_reference_ids != excluded:
        raise PhosPyInputError(
            f"{field_name}.excluded_reference_ids must match excluded rows"
        )


def _require_summary_counts(
    value: object, expected: Mapping[str, int], *, field_name: str
) -> None:
    for name, expected_value in expected.items():
        if getattr(value, name) != expected_value:
            raise PhosPyInputError(f"{field_name}.{name} is inconsistent with rows")


def _required_loaded_table(
    tables: Mapping[str, pd.DataFrame | None], key: str
) -> pd.DataFrame:
    table = tables[key]
    if table is None:
        raise PhosPyInputError(f"scoring table {key} is required")
    return table


def _require_exact_fields(
    payload: Mapping[str, object], *, field_name: str, fields: frozenset[str]
) -> None:
    unknown = sorted(str(key) for key in payload if str(key) not in fields)
    if unknown:
        raise PhosPyInputError(
            f"{field_name} contains unsupported field(s): {', '.join(unknown)}"
        )
    missing = sorted(field for field in fields if field not in payload)
    if missing:
        raise PhosPyInputError(
            f"{field_name} is missing required field(s): {', '.join(missing)}"
        )


def _require_sequence(value: object, *, field_name: str) -> Sequence[object]:
    if not _is_sequence(value):
        raise PhosPyInputError(f"{field_name} must be an array")
    return value


def _is_sequence(value: object) -> TypeGuard[Sequence[object]]:
    return isinstance(value, Sequence) and not isinstance(
        value, (str, bytes, bytearray)
    )


def _string_tuple(value: object, *, field_name: str) -> tuple[str, ...]:
    return tuple(
        require_str(item, field_name=f"{field_name}[{index}]")
        for index, item in enumerate(_require_sequence(value, field_name=field_name))
    )


def _optional_string(value: object, *, field_name: str) -> str | None:
    if value is None:
        return None
    if not isinstance(value, str):
        raise PhosPyInputError(f"{field_name} must be a string or null")
    return value


def _require_status(value: object, *, field_name: str) -> SequenceValidationStatus:
    status = require_str(value, field_name=field_name)
    if status not in _SEQUENCE_STATUSES:
        raise PhosPyInputError(
            f"{field_name} must be one of: {', '.join(sorted(_SEQUENCE_STATUSES))}"
        )
    return cast(SequenceValidationStatus, status)


def _require_nonnegative_int(
    payload: Mapping[str, object], key: str, field_name: str
) -> int:
    value = require_int(payload.get(key), field_name=f"{field_name}.{key}")
    if value < 0:
        raise PhosPyInputError(f"{field_name}.{key} must be >= 0")
    return value


def _require_positive_int(
    payload: Mapping[str, object], key: str, field_name: str
) -> int:
    value = _require_nonnegative_int(payload, key, field_name)
    if value == 0:
        raise PhosPyInputError(f"{field_name}.{key} must be > 0")
    return value


def _optional_nonnegative_int(value: object, *, field_name: str) -> int | None:
    if value is None:
        return None
    parsed = require_int(value, field_name=field_name)
    if parsed < 0:
        raise PhosPyInputError(f"{field_name} must be >= 0")
    return parsed


__all__ = [
    "KINASE_SCORING_OPTIONAL_TABLE_KEYS",
    "KINASE_SCORING_STATE_FIELDS",
    "KINASE_SCORING_TABLE_KEYS",
    "KinaseScoringStateManifest",
    "ReconstructedKinaseScoringState",
    "kinase_scoring_state_to_payload",
    "parse_kinase_scoring_state",
    "reconstruct_kinase_scoring_state",
    "write_kinase_scoring_tables",
]
