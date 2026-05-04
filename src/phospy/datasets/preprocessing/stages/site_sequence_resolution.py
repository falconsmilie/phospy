"""Site-sequence resolution stage for dataset preprocessing."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import replace

import pandas as pd

from phospy.api.configs import (
    DATASET_SITE_SEQUENCE_CONFLICT_POLICY_ERROR,
    DATASET_SITE_SEQUENCE_CONFLICT_POLICY_PRESERVE_EXISTING,
    DATASET_SITE_SEQUENCE_CONFLICT_POLICY_REPLACE_EXISTING,
    DATASET_SITE_SEQUENCE_RESOLUTION_MODE_FILL_MISSING_ONLY,
    DATASET_SITE_SEQUENCE_RESOLUTION_MODE_REPLACE_EXISTING,
    DATASET_SITE_SEQUENCE_RESOLUTION_MODE_VALIDATE_EXISTING_ONLY,
)
from phospy.datasets.preprocessing.models import (
    DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION,
    PreprocessingStageResult,
    PreprocessingState,
)
from phospy.errors.input import PhosPyInputError, UnsupportedInputFormatError
from phospy.sequences import FastaProteinSequenceRepository
from phospy.sequences.resolver import (
    RESOLUTION_STATUS_RESOLVED,
    PhosphositeSequenceResolutionRequest,
    PhosphositeSequenceResolver,
)


class SiteSequenceResolutionStage:
    """Resolve site-sequence support from local FASTA for dataset preprocessing."""

    stage_key = DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION
    _RESOLVER_VERSION = "phospy.sequences.resolver.v1"

    def __init__(
        self,
        *,
        resolver: PhosphositeSequenceResolver | None = None,
    ) -> None:
        self._resolver = resolver or PhosphositeSequenceResolver()

    def run(self, state: PreprocessingState) -> PreprocessingStageResult:
        plan = state.plan
        if (
            not plan.site_sequence_resolution_enabled
            or plan.site_sequence_resolution_fasta_path is None
        ):
            return PreprocessingStageResult(
                state=state,
                diagnostics={
                    "dropped_row_ids": (),
                    "dropped_row_count": 0,
                    "imputed_cell_count": 0,
                    "imputed_row_ids": (),
                    "notes": "stage skipped",
                    "diagnostics": {"configured": False},
                },
            )

        accession_column = plan.site_sequence_resolution_accession_column
        site_column = plan.site_sequence_resolution_site_column
        self._require_columns(
            site_metadata=state.site_metadata,
            accession_column=accession_column,
            site_column=site_column,
        )
        repository = FastaProteinSequenceRepository.from_path(
            plan.site_sequence_resolution_fasta_path,
            source_label="dataset.site_sequence_resolution",
        )
        mode = plan.site_sequence_resolution_mode
        conflict_policy = plan.site_sequence_resolution_conflict_policy
        if (
            mode == DATASET_SITE_SEQUENCE_RESOLUTION_MODE_REPLACE_EXISTING
            and conflict_policy
            == DATASET_SITE_SEQUENCE_CONFLICT_POLICY_PRESERVE_EXISTING
        ):
            # Backward-compatible behavior for legacy plans that used mode-only
            # replacement semantics before explicit conflict-policy support.
            conflict_policy = DATASET_SITE_SEQUENCE_CONFLICT_POLICY_REPLACE_EXISTING
        flank_size = int(plan.site_sequence_resolution_flank_size)

        existing_site_sequence = (
            state.site_metadata.loc[:, "site_sequence"]
            if "site_sequence" in state.site_metadata.columns
            else pd.Series(pd.NA, index=state.site_metadata.index, dtype="string")
        )
        normalized_existing = _normalize_optional_site_sequence_series(
            existing_site_sequence
        )
        updated_site_sequence = existing_site_sequence.copy(deep=True)
        updated_site_metadata = state.site_metadata.copy(deep=True)

        resolved_site_count = 0
        unresolved_site_count = 0
        unresolved_counts_by_reason: dict[str, int] = defaultdict(int)
        existing_sequence_conflict_count = 0
        filled_missing_count = 0
        replaced_existing_count = 0
        preserved_existing_count = 0
        row_status: list[dict[str, object]] = []
        row_diagnostics: list[dict[str, object]] = []
        conflict_error_rows: list[dict[str, object]] = []
        fasta_source_path = repository.metadata.source_path
        fasta_sha256 = repository.metadata.sha256

        for row_index, row_id in enumerate(updated_site_metadata.index.tolist()):
            row_key = str(row_id)
            normalized_existing_value = normalized_existing.loc[row_id]
            has_existing = not bool(pd.isna(normalized_existing_value))
            site_id = row_key
            if (
                has_existing
                and mode == DATASET_SITE_SEQUENCE_RESOLUTION_MODE_FILL_MISSING_ONLY
            ):
                preserved_existing_count += 1
                status = "preserved_existing"
                reason = "existing site_sequence preserved (fill_missing_only)"
                row_status.append(
                    {"row_id": row_key, "status": status, "reason": reason}
                )
                row_diagnostics.append(
                    {
                        "row_index": int(row_index),
                        "row_id": row_key,
                        "site_id": site_id,
                        "status": status,
                        "existing_site_sequence": (
                            None
                            if pd.isna(normalized_existing_value)
                            else str(normalized_existing_value)
                        ),
                        "fasta_site_sequence": None,
                        "resolved_site_sequence": (
                            None
                            if pd.isna(normalized_existing_value)
                            else str(normalized_existing_value)
                        ),
                        "action": "preserve_existing",
                        "reason": reason,
                        "conflict_policy": conflict_policy,
                        "resolver_version": self._RESOLVER_VERSION,
                        "fasta_source_path": fasta_source_path,
                        "fasta_sha256": fasta_sha256,
                    }
                )
                continue
            if (
                not has_existing
            ) and mode == DATASET_SITE_SEQUENCE_RESOLUTION_MODE_VALIDATE_EXISTING_ONLY:
                status = "missing_existing_sequence"
                reason = "missing site_sequence preserved (validate_existing_only)"
                row_status.append(
                    {
                        "row_id": row_key,
                        "status": status,
                        "reason": reason,
                    }
                )
                row_diagnostics.append(
                    {
                        "row_index": int(row_index),
                        "row_id": row_key,
                        "site_id": site_id,
                        "status": status,
                        "existing_site_sequence": None,
                        "fasta_site_sequence": None,
                        "resolved_site_sequence": None,
                        "action": "skip_missing_existing",
                        "reason": reason,
                        "conflict_policy": conflict_policy,
                        "resolver_version": self._RESOLVER_VERSION,
                        "fasta_source_path": fasta_source_path,
                        "fasta_sha256": fasta_sha256,
                    }
                )
                continue

            resolution = self._resolver.run(
                PhosphositeSequenceResolutionRequest(
                    accession=updated_site_metadata.at[row_id, accession_column],
                    site_token=updated_site_metadata.at[row_id, site_column],
                    flank_size=flank_size,
                ),
                repository,
            )
            if resolution.status == RESOLUTION_STATUS_RESOLVED:
                resolved_site_count += 1
                resolved_sequence = (
                    None
                    if resolution.site_sequence is None
                    else resolution.site_sequence
                )
                if resolved_sequence is None:
                    unresolved_site_count += 1
                    unresolved_counts_by_reason["resolved_missing_sequence"] += 1
                    status = "resolved_missing_sequence"
                    reason = "resolver reported success without sequence output"
                    row_status.append(
                        {
                            "row_id": row_key,
                            "status": status,
                            "reason": reason,
                        }
                    )
                    row_diagnostics.append(
                        {
                            "row_index": int(row_index),
                            "row_id": row_key,
                            "site_id": site_id,
                            "status": status,
                            "existing_site_sequence": (
                                None
                                if pd.isna(normalized_existing_value)
                                else str(normalized_existing_value)
                            ),
                            "fasta_site_sequence": None,
                            "resolved_site_sequence": (
                                None
                                if pd.isna(normalized_existing_value)
                                else str(normalized_existing_value)
                            ),
                            "action": (
                                "preserve_existing" if has_existing else "unresolved"
                            ),
                            "reason": reason,
                            "conflict_policy": conflict_policy,
                            "resolver_version": self._RESOLVER_VERSION,
                            "fasta_source_path": fasta_source_path,
                            "fasta_sha256": fasta_sha256,
                        }
                    )
                    if has_existing:
                        preserved_existing_count += 1
                    continue

                if has_existing:
                    existing_text = str(normalized_existing_value)
                    if existing_text == resolved_sequence:
                        preserved_existing_count += 1
                        status = "resolved"
                        reason = "existing site_sequence validated against FASTA"
                        row_status.append(
                            {
                                "row_id": row_key,
                                "status": status,
                                "reason": reason,
                            }
                        )
                        row_diagnostics.append(
                            {
                                "row_index": int(row_index),
                                "row_id": row_key,
                                "site_id": site_id,
                                "status": status,
                                "existing_site_sequence": existing_text,
                                "fasta_site_sequence": resolved_sequence,
                                "resolved_site_sequence": existing_text,
                                "action": "validate_existing",
                                "reason": reason,
                                "conflict_policy": conflict_policy,
                                "resolver_version": self._RESOLVER_VERSION,
                                "fasta_source_path": fasta_source_path,
                                "fasta_sha256": fasta_sha256,
                            }
                        )
                        continue
                    existing_sequence_conflict_count += 1
                    if conflict_policy == DATASET_SITE_SEQUENCE_CONFLICT_POLICY_ERROR:
                        preserved_existing_count += 1
                        status = "existing_sequence_conflict"
                        reason = (
                            "existing site_sequence conflicts with FASTA-derived "
                            "sequence"
                        )
                        row_status.append(
                            {
                                "row_id": row_key,
                                "status": status,
                                "reason": reason,
                            }
                        )
                        conflict_row = {
                            "row_index": int(row_index),
                            "row_id": row_key,
                            "site_id": site_id,
                            "status": status,
                            "existing_site_sequence": existing_text,
                            "fasta_site_sequence": resolved_sequence,
                            "resolved_site_sequence": existing_text,
                            "action": "error",
                            "reason": reason,
                            "conflict_policy": conflict_policy,
                            "resolver_version": self._RESOLVER_VERSION,
                            "fasta_source_path": fasta_source_path,
                            "fasta_sha256": fasta_sha256,
                        }
                        row_diagnostics.append(conflict_row)
                        conflict_error_rows.append(conflict_row)
                        continue
                    if (
                        conflict_policy
                        == DATASET_SITE_SEQUENCE_CONFLICT_POLICY_REPLACE_EXISTING
                    ):
                        updated_site_sequence.at[row_id] = resolved_sequence
                        replaced_existing_count += 1
                        resolved_site_sequence = resolved_sequence
                        action = "replace_existing"
                    elif (
                        conflict_policy
                        == DATASET_SITE_SEQUENCE_CONFLICT_POLICY_PRESERVE_EXISTING
                    ):
                        preserved_existing_count += 1
                        resolved_site_sequence = existing_text
                        action = "preserve_existing"
                    else:  # pragma: no cover - defensive guard
                        raise PhosPyInputError(
                            "dataset preprocessing stage 'site_sequence_resolution' "
                            "received unsupported conflict policy: "
                            f"{conflict_policy!r}"
                        )
                    status = "existing_sequence_conflict"
                    reason = (
                        "existing site_sequence conflicts with FASTA-derived sequence"
                    )
                    row_status.append(
                        {
                            "row_id": row_key,
                            "status": status,
                            "reason": reason,
                        }
                    )
                    row_diagnostics.append(
                        {
                            "row_index": int(row_index),
                            "row_id": row_key,
                            "site_id": site_id,
                            "status": status,
                            "existing_site_sequence": existing_text,
                            "fasta_site_sequence": resolved_sequence,
                            "resolved_site_sequence": resolved_site_sequence,
                            "action": action,
                            "reason": reason,
                            "conflict_policy": conflict_policy,
                            "resolver_version": self._RESOLVER_VERSION,
                            "fasta_source_path": fasta_source_path,
                            "fasta_sha256": fasta_sha256,
                        }
                    )
                    continue

                updated_site_sequence.at[row_id] = resolved_sequence
                filled_missing_count += 1
                status = "resolved"
                reason = "missing site_sequence resolved from FASTA"
                row_status.append(
                    {
                        "row_id": row_key,
                        "status": status,
                        "reason": reason,
                    }
                )
                row_diagnostics.append(
                    {
                        "row_index": int(row_index),
                        "row_id": row_key,
                        "site_id": site_id,
                        "status": status,
                        "existing_site_sequence": None,
                        "fasta_site_sequence": resolved_sequence,
                        "resolved_site_sequence": resolved_sequence,
                        "action": "fill_missing",
                        "reason": reason,
                        "conflict_policy": conflict_policy,
                        "resolver_version": self._RESOLVER_VERSION,
                        "fasta_source_path": fasta_source_path,
                        "fasta_sha256": fasta_sha256,
                    }
                )
                continue

            unresolved_site_count += 1
            unresolved_reason = str(resolution.status)
            unresolved_counts_by_reason[unresolved_reason] += 1
            if has_existing:
                preserved_existing_count += 1
            status = unresolved_reason
            reason = resolution.reason
            row_status.append(
                {
                    "row_id": row_key,
                    "status": status,
                    "reason": reason,
                }
            )
            row_diagnostics.append(
                {
                    "row_index": int(row_index),
                    "row_id": row_key,
                    "site_id": site_id,
                    "status": status,
                    "existing_site_sequence": (
                        None
                        if pd.isna(normalized_existing_value)
                        else str(normalized_existing_value)
                    ),
                    "fasta_site_sequence": None,
                    "resolved_site_sequence": (
                        None
                        if pd.isna(normalized_existing_value)
                        else str(normalized_existing_value)
                    ),
                    "action": "preserve_existing" if has_existing else "unresolved",
                    "reason": reason,
                    "conflict_policy": conflict_policy,
                    "resolver_version": self._RESOLVER_VERSION,
                    "fasta_source_path": fasta_source_path,
                    "fasta_sha256": fasta_sha256,
                }
            )

        diagnostics = {
            "fasta_source_path": repository.metadata.source_path,
            "fasta_source_label": repository.metadata.source_label,
            "fasta_sha256": repository.metadata.sha256,
            "resolver_version": self._RESOLVER_VERSION,
            "flank_size": flank_size,
            "mode": mode,
            "conflict_policy": conflict_policy,
            "accession_column": accession_column,
            "site_column": site_column,
            "resolved_site_count": int(resolved_site_count),
            "unresolved_site_count": int(unresolved_site_count),
            "unresolved_counts_by_reason": {
                key: int(unresolved_counts_by_reason[key])
                for key in sorted(unresolved_counts_by_reason)
            },
            "existing_sequence_conflict_count": int(existing_sequence_conflict_count),
            "filled_missing_count": int(filled_missing_count),
            "replaced_existing_count": int(replaced_existing_count),
            "preserved_existing_count": int(preserved_existing_count),
            "row_status": row_status,
            "row_diagnostics": row_diagnostics,
        }
        if conflict_error_rows:
            raise PhosPyInputError(
                "dataset preprocessing stage 'site_sequence_resolution' detected "
                "conflicts between existing site_sequence and FASTA-derived values "
                "under conflict_policy='error'",
                diagnostics=diagnostics,
            )
        updated_site_metadata.loc[:, "site_sequence"] = updated_site_sequence.astype(
            "string"
        )
        next_state = replace(state, site_metadata=updated_site_metadata)
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

    @staticmethod
    def _require_columns(
        *,
        site_metadata: pd.DataFrame,
        accession_column: str,
        site_column: str,
    ) -> None:
        missing_columns = [
            column
            for column in (accession_column, site_column)
            if column not in site_metadata.columns
        ]
        if not missing_columns:
            return
        joined = ", ".join(missing_columns)
        raise UnsupportedInputFormatError(
            "dataset preprocessing stage 'site_sequence_resolution' requires "
            f"site_metadata column(s): {joined}"
        )


def _normalize_optional_site_sequence_series(column: pd.Series) -> pd.Series:
    as_string = column.astype("string").str.strip()
    missing = column.isna() | as_string.isna() | (as_string == "")
    return as_string.where(~missing, other=pd.NA)


__all__ = ["SiteSequenceResolutionStage"]
