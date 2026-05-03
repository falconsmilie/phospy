"""Site-sequence resolution stage for dataset preprocessing."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import replace

import pandas as pd

from phospy.api.configs import (
    DATASET_SITE_SEQUENCE_RESOLUTION_MODE_FILL_MISSING_ONLY,
    DATASET_SITE_SEQUENCE_RESOLUTION_MODE_REPLACE_EXISTING,
    DATASET_SITE_SEQUENCE_RESOLUTION_MODE_VALIDATE_EXISTING_ONLY,
)
from phospy.datasets.preprocessing.models import (
    DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION,
    PreprocessingStageResult,
    PreprocessingState,
)
from phospy.errors.input import UnsupportedInputFormatError
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
        flank_size = int(plan.site_sequence_resolution_flank_size)

        normalized_existing = _normalize_optional_site_sequence_series(
            state.site_metadata.loc[:, "site_sequence"]
            if "site_sequence" in state.site_metadata.columns
            else pd.Series(pd.NA, index=state.site_metadata.index, dtype="string")
        )
        updated_site_sequence = normalized_existing.copy()
        updated_site_metadata = state.site_metadata.copy(deep=True)

        resolved_site_count = 0
        unresolved_site_count = 0
        unresolved_counts_by_reason: dict[str, int] = defaultdict(int)
        existing_sequence_conflict_count = 0
        filled_missing_count = 0
        replaced_existing_count = 0
        preserved_existing_count = 0
        row_status: list[dict[str, object]] = []

        for row_id in updated_site_metadata.index.tolist():
            row_key = str(row_id)
            existing_value = updated_site_sequence.loc[row_id]
            has_existing = not bool(pd.isna(existing_value))
            if (
                has_existing
                and mode == DATASET_SITE_SEQUENCE_RESOLUTION_MODE_FILL_MISSING_ONLY
            ):
                preserved_existing_count += 1
                row_status.append(
                    {
                        "row_id": row_key,
                        "status": "preserved_existing",
                        "reason": "existing site_sequence preserved (fill_missing_only)",
                    }
                )
                continue
            if (
                not has_existing
            ) and mode == DATASET_SITE_SEQUENCE_RESOLUTION_MODE_VALIDATE_EXISTING_ONLY:
                row_status.append(
                    {
                        "row_id": row_key,
                        "status": "missing_existing_sequence",
                        "reason": "missing site_sequence preserved (validate_existing_only)",
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
                    row_status.append(
                        {
                            "row_id": row_key,
                            "status": "resolved_missing_sequence",
                            "reason": "resolver reported success without sequence output",
                        }
                    )
                    if has_existing:
                        preserved_existing_count += 1
                    continue

                if has_existing:
                    existing_text = str(existing_value)
                    if existing_text == resolved_sequence:
                        preserved_existing_count += 1
                        row_status.append(
                            {
                                "row_id": row_key,
                                "status": "resolved",
                                "reason": "existing site_sequence validated against FASTA",
                            }
                        )
                        continue
                    existing_sequence_conflict_count += 1
                    if mode == DATASET_SITE_SEQUENCE_RESOLUTION_MODE_REPLACE_EXISTING:
                        updated_site_sequence.at[row_id] = resolved_sequence
                        replaced_existing_count += 1
                    else:
                        preserved_existing_count += 1
                    row_status.append(
                        {
                            "row_id": row_key,
                            "status": "existing_sequence_conflict",
                            "reason": (
                                "existing site_sequence conflicts with FASTA-derived "
                                "sequence"
                            ),
                        }
                    )
                    continue

                updated_site_sequence.at[row_id] = resolved_sequence
                filled_missing_count += 1
                row_status.append(
                    {
                        "row_id": row_key,
                        "status": "resolved",
                        "reason": "missing site_sequence resolved from FASTA",
                    }
                )
                continue

            unresolved_site_count += 1
            unresolved_reason = str(resolution.status)
            unresolved_counts_by_reason[unresolved_reason] += 1
            if has_existing:
                preserved_existing_count += 1
            row_status.append(
                {
                    "row_id": row_key,
                    "status": unresolved_reason,
                    "reason": resolution.reason,
                }
            )

        updated_site_metadata.loc[:, "site_sequence"] = updated_site_sequence.astype(
            "string"
        )
        next_state = replace(state, site_metadata=updated_site_metadata)
        diagnostics = {
            "fasta_source_path": repository.metadata.source_path,
            "fasta_source_label": repository.metadata.source_label,
            "fasta_sha256": repository.metadata.sha256,
            "resolver_version": self._RESOLVER_VERSION,
            "flank_size": flank_size,
            "mode": mode,
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
        }
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
