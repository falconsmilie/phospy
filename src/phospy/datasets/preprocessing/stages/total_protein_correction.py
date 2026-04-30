"""Total-protein correction stage for dataset preprocessing."""

from __future__ import annotations

from dataclasses import dataclass, replace

import pandas as pd

from phospy.api.configs import (
    DATASET_INTENSITY_TRANSFORM_POLICY_LOG2,
    DATASET_TOTAL_PROTEIN_CORRECTION_IDENTITY_MODE_DIRECT,
    DATASET_TOTAL_PROTEIN_CORRECTION_IDENTITY_MODE_MAPPING_TABLE,
    DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_NONE,
    DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_SUBTRACT_LOG_TOTAL,
    DATASET_TOTAL_PROTEIN_CORRECTION_UNMATCHED_POLICY_ALLOW_UNCORRECTED,
    DATASET_TOTAL_PROTEIN_CORRECTION_UNMATCHED_POLICY_ERROR,
)
from phospy.datasets.preprocessing.models import (
    DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION,
    PreprocessingStageResult,
    PreprocessingState,
    TotalProteinCorrectionIdentityPolicy,
)
from phospy.datasets.processing_state import (
    TOTAL_PROTEIN_CORRECTION_DIAGNOSTICS_SCHEMA_VERSION_V1,
)
from phospy.errors.input import PhosPyInputError
from phospy.provenance.hashing import hash_table

_INDEX_IDENTITY_KEY = "__index__"
_GENE_SYMBOL_COLUMN = "gene_symbol"


@dataclass(frozen=True, slots=True)
class _ResolvedIdentityMapping:
    phosphosite_to_total_row: dict[str, str]
    corrected_phosphosite_rows: tuple[str, ...]
    uncorrected_phosphosite_rows: tuple[str, ...]
    used_total_rows: tuple[str, ...]
    unused_total_rows: tuple[str, ...]
    total_rows_used_by_multiple_phosphosites: tuple[str, ...]
    gene_symbol_matching_used: bool
    gene_symbol_warning: str | None


class TotalProteinCorrectionStage:
    """Apply log-scale phospho-to-total subtraction when requested."""

    stage_key = DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION

    def run(self, state: PreprocessingState) -> PreprocessingStageResult:
        requested_policy = state.plan.total_protein_correction_policy
        identity_policy = state.plan.total_protein_correction_identity_policy
        if requested_policy == DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_NONE:
            return PreprocessingStageResult(
                state=state,
                diagnostics={
                    "dropped_row_ids": (),
                    "dropped_row_count": 0,
                    "imputed_cell_count": 0,
                    "imputed_row_ids": (),
                    "notes": "stage executed",
                    "diagnostics": {
                        "diagnostics_schema_version": (
                            TOTAL_PROTEIN_CORRECTION_DIAGNOSTICS_SCHEMA_VERSION_V1
                        ),
                        "policy": str(requested_policy),
                        "requested_policy": str(requested_policy),
                        "resolved_policy": str(requested_policy),
                        "identity_mode": str(identity_policy.mode),
                        "phosphosite_key": str(identity_policy.phosphosite_key),
                        "total_protein_key": str(identity_policy.total_protein_key),
                        "duplicate_policy": str(identity_policy.duplicate_policy),
                        "unmatched_policy": str(identity_policy.unmatched_policy),
                    },
                },
            )
        if (
            requested_policy
            != DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_SUBTRACT_LOG_TOTAL
        ):
            raise PhosPyInputError(
                "dataset build request preprocessing_config contains an unsupported "
                "total_protein_correction.policy"
            )
        if (
            state.plan.intensity_transform_policy
            != DATASET_INTENSITY_TRANSFORM_POLICY_LOG2
        ):
            raise PhosPyInputError(
                "dataset build request "
                f"preprocessing_config.total_protein_correction.policy={requested_policy!r} "
                "requires log2-scale phospho and total values. Configure "
                "preprocessing_config.intensity_transform.policy='log2', or disable "
                "total-protein correction."
            )
        if state.total is None:
            raise PhosPyInputError(
                "dataset build request preprocessing_config.total_protein_correction."
                f"policy={requested_policy!r} requires total input data"
            )
        if not state.total.columns.equals(state.phospho.columns):
            raise PhosPyInputError(
                "dataset build request preprocessing total/protein correction requires "
                "total columns to exactly match phospho columns"
            )
        _require_numeric_matrix(state.phospho, field_name="phospho")
        _require_numeric_matrix(state.total, field_name="total")

        mapping = _resolve_identity_mapping(
            state=state,
            identity_policy=identity_policy,
        )
        corrected_rows = mapping.corrected_phosphosite_rows
        if not corrected_rows:
            raise PhosPyInputError(
                "dataset build request preprocessing total/protein correction found "
                "no phosphosite rows that resolved to a total-protein row under the "
                "configured identity policy"
            )

        corrected = state.phospho.copy(deep=True)
        corrected_index = corrected.index
        for phosphosite_row_id in corrected_rows:
            total_row_id = mapping.phosphosite_to_total_row[phosphosite_row_id]
            corrected.loc[phosphosite_row_id, :] = (
                state.phospho.loc[phosphosite_row_id, :]
                - state.total.loc[total_row_id, :]
            ).to_numpy(copy=False)
        corrected.index = corrected_index.copy()

        diagnostics = _build_diagnostics(
            requested_policy=requested_policy,
            identity_policy=identity_policy,
            mapping=mapping,
            phospho=state.phospho,
            total=state.total,
            corrected=corrected,
        )
        next_state = replace(state, phospho=corrected)
        return PreprocessingStageResult(
            state=next_state,
            diagnostics={
                "dropped_row_ids": (),
                "dropped_row_count": 0,
                "imputed_cell_count": 0,
                "imputed_row_ids": (),
                "notes": "stage executed",
                "diagnostics": diagnostics,
            },
        )


def _build_diagnostics(
    *,
    requested_policy: str,
    identity_policy: TotalProteinCorrectionIdentityPolicy,
    mapping: _ResolvedIdentityMapping,
    phospho: pd.DataFrame,
    total: pd.DataFrame,
    corrected: pd.DataFrame,
) -> dict[str, object]:
    corrected_rows = int(len(mapping.corrected_phosphosite_rows))
    uncorrected_rows = int(len(mapping.uncorrected_phosphosite_rows))
    unused_total_rows = int(len(mapping.unused_total_rows))
    shared_total_rows = int(len(mapping.total_rows_used_by_multiple_phosphosites))
    diagnostics: dict[str, object] = {
        "diagnostics_schema_version": TOTAL_PROTEIN_CORRECTION_DIAGNOSTICS_SCHEMA_VERSION_V1,
        "policy": str(requested_policy),
        "requested_policy": str(requested_policy),
        "resolved_policy": str(requested_policy),
        "formula": "log2_phospho - log2_total",
        "requires_log_scale": True,
        "input_scale": "log2",
        "output_scale": "log2_ratio",
        "quantitative_meaning": "phospho_total_log_ratio",
        "matched_rows": corrected_rows,
        "identity_mode": str(identity_policy.mode),
        "phosphosite_key": str(identity_policy.phosphosite_key),
        "total_protein_key": str(identity_policy.total_protein_key),
        "mapping_phosphosite_key": identity_policy.mapping_phosphosite_key,
        "mapping_total_protein_key": identity_policy.mapping_total_protein_key,
        "mapping_table_fingerprint": identity_policy.mapping_table_fingerprint,
        "duplicate_policy": str(identity_policy.duplicate_policy),
        "unmatched_policy": str(identity_policy.unmatched_policy),
        "phosphosite_row_count": int(phospho.shape[0]),
        "total_protein_row_count": int(total.shape[0]),
        "corrected_row_count": corrected_rows,
        "uncorrected_row_count": uncorrected_rows,
        "unused_total_protein_row_count": unused_total_rows,
        "total_rows_used_by_multiple_phosphosites": shared_total_rows,
        "unmatched_phosphosite_row_ids": list(mapping.uncorrected_phosphosite_rows),
        "unused_total_protein_row_ids": list(mapping.unused_total_rows),
        "gene_symbol_matching_used": bool(mapping.gene_symbol_matching_used),
        "gene_symbol_identity_warning": mapping.gene_symbol_warning,
        "total_table_hash": hash_table(
            total,
            name="total_protein_correction.total",
        ),
        "input_phospho_hash": hash_table(
            phospho,
            name="total_protein_correction.input.phospho",
        ),
        "output_phospho_hash": hash_table(
            corrected,
            name="total_protein_correction.output.phospho",
        ),
    }
    return diagnostics


def _resolve_identity_mapping(
    *,
    state: PreprocessingState,
    identity_policy: TotalProteinCorrectionIdentityPolicy,
) -> _ResolvedIdentityMapping:
    total_key = _resolve_total_key_series(
        total=state.total,
        key_name=identity_policy.total_protein_key,
    )
    phosphosite_key = _resolve_phosphosite_key_series(
        site_metadata=state.site_metadata,
        phospho=state.phospho,
        key_name=identity_policy.phosphosite_key,
    )

    gene_symbol_matching_used = _is_gene_symbol_identity_policy(identity_policy)
    normalize_gene = bool(gene_symbol_matching_used)
    normalized_total_key = _normalize_identifier_series(
        total_key, gene_symbol=normalize_gene
    )
    normalized_phosphosite_key = _normalize_identifier_series(
        phosphosite_key, gene_symbol=normalize_gene
    )
    _require_non_empty_keys(
        normalized_total_key,
        field_name=(
            "dataset build request preprocessing total/protein correction total "
            f"identity key '{identity_policy.total_protein_key}'"
        ),
    )
    _require_non_empty_keys(
        normalized_phosphosite_key,
        field_name=(
            "dataset build request preprocessing total/protein correction phosphosite "
            f"identity key '{identity_policy.phosphosite_key}'"
        ),
    )
    duplicate_total_key_mask = normalized_total_key.duplicated(keep=False)
    if bool(duplicate_total_key_mask.any()):
        preview = _unique_preview(normalized_total_key.loc[duplicate_total_key_mask])
        raise PhosPyInputError(
            "dataset build request preprocessing total/protein correction requires "
            "unique resolved total-protein identity keys; duplicate keys: "
            f"{preview}"
        )

    total_key_to_row = {
        str(key): str(row_id) for row_id, key in normalized_total_key.items()
    }

    mapping_mode = identity_policy.mode
    if mapping_mode == DATASET_TOTAL_PROTEIN_CORRECTION_IDENTITY_MODE_DIRECT:
        row_mapping = _resolve_direct_mapping(
            normalized_phosphosite_key=normalized_phosphosite_key,
            total_key_to_row=total_key_to_row,
        )
    elif mapping_mode == DATASET_TOTAL_PROTEIN_CORRECTION_IDENTITY_MODE_MAPPING_TABLE:
        row_mapping = _resolve_mapping_table_mode_mapping(
            phospho=state.phospho,
            normalized_phosphosite_key=normalized_phosphosite_key,
            total_key_to_row=total_key_to_row,
            identity_policy=identity_policy,
        )
    else:  # pragma: no cover - guarded by config validation
        raise PhosPyInputError(
            "dataset build request preprocessing total/protein correction identity "
            f"mode {mapping_mode!r} is unsupported"
        )

    corrected_rows = tuple(sorted(row_mapping))
    uncorrected_rows = tuple(
        row_id
        for row_id in state.phospho.index.astype(str).tolist()
        if row_id not in row_mapping
    )
    if uncorrected_rows and (
        identity_policy.unmatched_policy
        == DATASET_TOTAL_PROTEIN_CORRECTION_UNMATCHED_POLICY_ERROR
    ):
        preview = ", ".join(uncorrected_rows[:5])
        total_rows = int(len(state.phospho.index))
        unmatched_rows = int(len(uncorrected_rows))
        percentage = (unmatched_rows / total_rows) * 100.0 if total_rows else 0.0
        raise PhosPyInputError(
            "dataset build request preprocessing total/protein correction "
            "requires complete phosphosite-to-total mapping under the configured "
            f"identity policy but found {unmatched_rows} unmatched phosphosite rows "
            f"out of {total_rows} ({percentage:.1f}%): {preview}"
        )
    if uncorrected_rows and (
        identity_policy.unmatched_policy
        != DATASET_TOTAL_PROTEIN_CORRECTION_UNMATCHED_POLICY_ALLOW_UNCORRECTED
    ):
        raise PhosPyInputError(
            "dataset build request preprocessing total/protein correction has "
            "unmatched phosphosite rows but unmatched_policy does not support "
            "retaining uncorrected rows"
        )

    used_total_rows_counter: dict[str, int] = {}
    for total_row_id in row_mapping.values():
        used_total_rows_counter[total_row_id] = (
            used_total_rows_counter.get(total_row_id, 0) + 1
        )
    used_total_rows = tuple(sorted(used_total_rows_counter))
    if state.total is None:
        raise PhosPyInputError(
            "total matrix is required for total protein correction identity resolution"
        )
    all_total_rows = tuple(
        str(label) for label in state.total.index.astype(str).tolist()
    )
    unused_total_rows = tuple(
        row_id for row_id in all_total_rows if row_id not in used_total_rows_counter
    )
    shared_total_rows = tuple(
        row_id
        for row_id, count in sorted(used_total_rows_counter.items())
        if int(count) > 1
    )
    gene_symbol_warning = None
    if gene_symbol_matching_used:
        gene_symbol_warning = (
            "Gene-symbol identity matching is a convenience policy and may be "
            "biologically ambiguous for isoform- or protein-group-specific "
            "datasets."
        )
    return _ResolvedIdentityMapping(
        phosphosite_to_total_row=row_mapping,
        corrected_phosphosite_rows=corrected_rows,
        uncorrected_phosphosite_rows=uncorrected_rows,
        used_total_rows=used_total_rows,
        unused_total_rows=unused_total_rows,
        total_rows_used_by_multiple_phosphosites=shared_total_rows,
        gene_symbol_matching_used=gene_symbol_matching_used,
        gene_symbol_warning=gene_symbol_warning,
    )


def _resolve_direct_mapping(
    *,
    normalized_phosphosite_key: pd.Series,
    total_key_to_row: dict[str, str],
) -> dict[str, str]:
    mapping: dict[str, str] = {}
    for row_id, key in normalized_phosphosite_key.items():
        total_row_id = total_key_to_row.get(str(key))
        if total_row_id is None:
            continue
        mapping[str(row_id)] = total_row_id
    return mapping


def _resolve_mapping_table_mode_mapping(
    *,
    phospho: pd.DataFrame,
    normalized_phosphosite_key: pd.Series,
    total_key_to_row: dict[str, str],
    identity_policy: TotalProteinCorrectionIdentityPolicy,
) -> dict[str, str]:
    mapping_rows = identity_policy.mapping_table
    if mapping_rows is None:
        raise PhosPyInputError(
            "dataset build request preprocessing total/protein correction "
            "identity.mapping_table is required when identity.mode='mapping_table'"
        )
    phosphosite_key_to_rows: dict[str, list[str]] = {}
    for row_id, key in normalized_phosphosite_key.items():
        resolved_key = str(key)
        if resolved_key not in phosphosite_key_to_rows:
            phosphosite_key_to_rows[resolved_key] = []
        phosphosite_key_to_rows[resolved_key].append(str(row_id))
    mapping_table = pd.DataFrame(
        mapping_rows,
        columns=pd.Index(["mapping_phosphosite_key", "mapping_total_protein_key"]),
    )
    mapping_table.loc[:, "mapping_phosphosite_key"] = (
        mapping_table.loc[:, "mapping_phosphosite_key"].astype("string").str.strip()
    )
    mapping_table.loc[:, "mapping_total_protein_key"] = (
        mapping_table.loc[:, "mapping_total_protein_key"].astype("string").str.strip()
    )
    _require_non_empty_keys(
        mapping_table.loc[:, "mapping_phosphosite_key"],
        field_name=(
            "dataset build request preprocessing total/protein correction identity."
            "mapping_table phosphosite key"
        ),
    )
    _require_non_empty_keys(
        mapping_table.loc[:, "mapping_total_protein_key"],
        field_name=(
            "dataset build request preprocessing total/protein correction identity."
            "mapping_table total key"
        ),
    )
    duplicate_mapping_rows = mapping_table.duplicated(keep=False)
    if bool(duplicate_mapping_rows.any()):
        preview = ", ".join(
            [
                f"{left}->{right}"
                for left, right in mapping_table.loc[duplicate_mapping_rows, :]
                .head(5)
                .itertuples(index=False)
            ]
        )
        raise PhosPyInputError(
            "dataset build request preprocessing total/protein correction "
            "identity.mapping_table contains duplicate phosphosite-to-total mapping "
            f"rows: {preview}"
        )
    ambiguous_mapping = (
        mapping_table.groupby("mapping_phosphosite_key", sort=False)[
            "mapping_total_protein_key"
        ]
        .nunique(dropna=False)
        .loc[lambda series: series > 1]
    )
    if not ambiguous_mapping.empty:
        preview = ", ".join(str(item) for item in ambiguous_mapping.index.tolist()[:5])
        raise PhosPyInputError(
            "dataset build request preprocessing total/protein correction "
            "identity.mapping_table maps one phosphosite key to multiple total-"
            f"protein keys; ambiguous phosphosite keys: {preview}"
        )

    unknown_phosphosite = sorted(
        {
            str(key)
            for key in mapping_table.loc[:, "mapping_phosphosite_key"].tolist()
            if str(key) not in phosphosite_key_to_rows
        }
    )
    if unknown_phosphosite:
        preview = ", ".join(unknown_phosphosite[:5])
        raise PhosPyInputError(
            "dataset build request preprocessing total/protein correction "
            "identity.mapping_table references unknown phosphosite keys: "
            f"{preview}"
        )
    unknown_total = sorted(
        {
            str(key)
            for key in mapping_table.loc[:, "mapping_total_protein_key"].tolist()
            if str(key) not in total_key_to_row
        }
    )
    if unknown_total:
        preview = ", ".join(unknown_total[:5])
        raise PhosPyInputError(
            "dataset build request preprocessing total/protein correction "
            "identity.mapping_table references unknown total-protein keys: "
            f"{preview}"
        )

    mapping: dict[str, str] = {}
    for phosphosite_key, total_key in mapping_table.itertuples(index=False):
        resolved_total_row = total_key_to_row[str(total_key)]
        for phosphosite_row_id in phosphosite_key_to_rows[str(phosphosite_key)]:
            previous = mapping.get(str(phosphosite_row_id))
            if previous is not None and previous != resolved_total_row:
                raise PhosPyInputError(
                    "dataset build request preprocessing total/protein correction "
                    "identity.mapping_table resolves one phosphosite row to multiple "
                    "total-protein rows"
                )
            mapping[str(phosphosite_row_id)] = resolved_total_row
    if len(mapping) > int(len(phospho.index)):
        raise PhosPyInputError(
            "dataset build request preprocessing total/protein correction mapping "
            "resolved to more phosphosite rows than exist in the phospho matrix"
        )
    return mapping


def _resolve_total_key_series(
    *, total: pd.DataFrame | None, key_name: str
) -> pd.Series:
    if total is None:  # pragma: no cover - guarded earlier
        raise PhosPyInputError(
            "dataset build request preprocessing total/protein correction requires "
            "total input data"
        )
    resolved_key = str(key_name).strip()
    if resolved_key == _INDEX_IDENTITY_KEY:
        return pd.Series(total.index, index=total.index, dtype="string")
    index_name = None if total.index.name is None else str(total.index.name).strip()
    if index_name == resolved_key:
        return pd.Series(total.index, index=total.index, dtype="string")
    raise PhosPyInputError(
        "dataset build request preprocessing total/protein correction cannot resolve "
        f"total identity key {resolved_key!r}. This total matrix currently exposes "
        "identity via total.index. Set identity.total_protein_key='__index__' or "
        "set total.index.name to the configured key."
    )


def _resolve_phosphosite_key_series(
    *,
    site_metadata: pd.DataFrame,
    phospho: pd.DataFrame,
    key_name: str,
) -> pd.Series:
    resolved_key = str(key_name).strip()
    if resolved_key == _INDEX_IDENTITY_KEY:
        return pd.Series(phospho.index, index=phospho.index, dtype="string")
    if resolved_key not in site_metadata.columns:
        raise PhosPyInputError(
            "dataset build request preprocessing total/protein correction cannot "
            f"resolve phosphosite identity key {resolved_key!r} from site_metadata"
        )
    return pd.Series(
        site_metadata.loc[:, resolved_key], index=phospho.index, dtype="string"
    )


def _require_non_empty_keys(series: pd.Series, *, field_name: str) -> None:
    normalized = series.astype("string").str.strip()
    invalid_mask = normalized.isna() | (normalized == "")
    if bool(invalid_mask.any()):
        invalid_index = [
            str(index_value)
            for index_value, is_invalid in zip(
                series.index.tolist(),
                invalid_mask.tolist(),
                strict=True,
            )
            if bool(is_invalid)
        ]
        preview = ", ".join(invalid_index[:5])
        raise PhosPyInputError(
            f"{field_name} contains null/empty identifiers at rows: {preview}"
        )


def _normalize_identifier_series(series: pd.Series, *, gene_symbol: bool) -> pd.Series:
    normalized = series.astype("string").str.strip()
    if gene_symbol:
        return normalized.str.upper()
    return normalized


def _is_gene_symbol_identity_policy(
    identity_policy: TotalProteinCorrectionIdentityPolicy,
) -> bool:
    phosphosite_key = str(identity_policy.phosphosite_key).strip().lower()
    total_key = str(identity_policy.total_protein_key).strip().lower()
    if (
        identity_policy.mode
        == DATASET_TOTAL_PROTEIN_CORRECTION_IDENTITY_MODE_MAPPING_TABLE
    ):
        mapping_total_key = (
            ""
            if identity_policy.mapping_total_protein_key is None
            else str(identity_policy.mapping_total_protein_key).strip().lower()
        )
        mapping_site_key = (
            ""
            if identity_policy.mapping_phosphosite_key is None
            else str(identity_policy.mapping_phosphosite_key).strip().lower()
        )
        return "gene_symbol" in {
            phosphosite_key,
            total_key,
            mapping_site_key,
            mapping_total_key,
        }
    return "gene_symbol" in {phosphosite_key, total_key}


def _unique_preview(series: pd.Series, *, limit: int = 5) -> str:
    values = [
        str(item) for item in series.dropna().astype(str).unique().tolist()[:limit]
    ]
    return ", ".join(values)


def _require_numeric_matrix(frame: pd.DataFrame, *, field_name: str) -> None:
    non_numeric_columns = [
        str(column)
        for column in frame.columns
        if (
            not pd.api.types.is_numeric_dtype(frame[column])
            or pd.api.types.is_bool_dtype(frame[column])
        )
    ]
    if non_numeric_columns:
        raise PhosPyInputError(
            "dataset build request preprocessing total/protein correction requires "
            f"numeric non-boolean {field_name} columns; invalid columns: "
            + ", ".join(non_numeric_columns)
        )


__all__ = ["TotalProteinCorrectionStage"]
