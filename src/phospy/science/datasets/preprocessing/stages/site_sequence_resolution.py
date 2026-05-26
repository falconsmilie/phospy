"""Site-sequence resolution stage for dataset preprocessing."""

from __future__ import annotations

from dataclasses import replace

import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.science.datasets.preprocessing.models import (
    DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION,
    PreprocessingPlan,
    PreprocessingStageResult,
    PreprocessingState,
    PreprocessingStateTableKey,
)
from phospy.science.datasets.preprocessing.site_sequence import (
    SiteSequenceConflictResolver,
    SiteSequenceDiagnosticsBuilder,
    SiteSequenceMetadataUpdater,
    SiteSequenceReferenceLoader,
    SiteSequenceResolutionRequestBuilder,
)
from phospy.science.datasets.preprocessing.stage_contract import (
    PreprocessingStageContract,
)
from phospy.science.sequences.resolver import (
    RESOLUTION_STATUS_RESOLVED,
    PhosphositeSequenceResolver,
)


class SiteSequenceResolutionStage:
    """Resolve site-sequence support from local FASTA for dataset preprocessing."""

    stage_key = DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION
    _RESOLVER_VERSION = "phospy.science.sequences.resolver.v1"

    def __init__(
        self,
        *,
        resolver: PhosphositeSequenceResolver | None = None,
        reference_loader: SiteSequenceReferenceLoader | None = None,
        request_builder: SiteSequenceResolutionRequestBuilder | None = None,
        conflict_resolver: SiteSequenceConflictResolver | None = None,
    ) -> None:
        self._resolver = resolver or PhosphositeSequenceResolver()
        self._reference_loader = reference_loader or SiteSequenceReferenceLoader()
        self._request_builder = (
            request_builder or SiteSequenceResolutionRequestBuilder()
        )
        self._conflict_resolver = conflict_resolver or SiteSequenceConflictResolver()

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

        self._request_builder.require_columns(
            site_metadata=state.site_metadata,
            accession_column=plan.site_sequence_resolution_accession_column,
            site_column=plan.site_sequence_resolution_site_column,
        )
        context = self._reference_loader.load(plan)
        existing_site_sequence = (
            state.site_metadata.loc[:, "site_sequence"]
            if "site_sequence" in state.site_metadata.columns
            else pd.Series(pd.NA, index=state.site_metadata.index, dtype="string")
        )
        normalized_existing = _normalize_optional_site_sequence_series(
            existing_site_sequence
        )
        row_requests = self._request_builder.build(
            site_metadata=state.site_metadata,
            normalized_existing_site_sequence=normalized_existing,
            context=context,
        )
        diagnostics_builder = SiteSequenceDiagnosticsBuilder(
            context=context,
            resolver_version=self._RESOLVER_VERSION,
        )
        metadata_updater = SiteSequenceMetadataUpdater(
            site_metadata=state.site_metadata,
            existing_site_sequence=existing_site_sequence,
        )

        for row in row_requests:
            if row.skip_status is not None:
                diagnostics_builder.record_pre_resolution_skip(row)
                continue

            request = row.resolver_request
            if request is None:
                raise PhosPyInputError(
                    "dataset preprocessing stage 'site_sequence_resolution' produced "
                    "an invalid row request without resolver input"
                )
            resolution = self._resolver.run(request, context.repository)
            if resolution.status == RESOLUTION_STATUS_RESOLVED:
                if resolution.site_sequence is None:
                    diagnostics_builder.record_resolved_missing_sequence(row=row)
                    continue
                if row.has_existing:
                    existing_sequence = row.existing_site_sequence
                    if existing_sequence is None:
                        raise PhosPyInputError(
                            "dataset preprocessing stage 'site_sequence_resolution' "
                            "detected inconsistent existing site_sequence state"
                        )
                    conflict = self._conflict_resolver.resolve(
                        existing_sequence=existing_sequence,
                        fasta_sequence=resolution.site_sequence,
                        conflict_policy=context.conflict_policy,
                    )
                    if conflict.should_replace_existing:
                        metadata_updater.assign(
                            row_id=row.row_id,
                            site_sequence=resolution.site_sequence,
                        )
                    diagnostics_builder.record_existing_resolution(
                        row=row,
                        fasta_site_sequence=resolution.site_sequence,
                        conflict=conflict,
                    )
                    continue

                metadata_updater.assign(
                    row_id=row.row_id,
                    site_sequence=resolution.site_sequence,
                )
                diagnostics_builder.record_fill_missing(
                    row=row,
                    resolved_site_sequence=resolution.site_sequence,
                )
                continue

            diagnostics_builder.record_unresolved(
                row=row,
                status=str(resolution.status),
                reason=resolution.reason,
            )

        diagnostics = diagnostics_builder.build(
            accession_column=context.accession_column,
            site_column=context.site_column,
        )
        if diagnostics.has_conflict_errors:
            raise PhosPyInputError(
                "dataset preprocessing stage 'site_sequence_resolution' detected "
                "conflicts between existing site_sequence and FASTA-derived values "
                "under conflict_policy='error'",
                diagnostics=diagnostics.payload,
            )

        next_state = replace(state, site_metadata=metadata_updater.build())
        return PreprocessingStageResult(
            state=next_state,
            diagnostics={
                "dropped_row_ids": (),
                "dropped_row_count": 0,
                "imputed_cell_count": 0,
                "imputed_row_ids": (),
                "notes": "stage executed",
                "diagnostics": diagnostics.payload,
            },
        )


def _normalize_optional_site_sequence_series(column: pd.Series) -> pd.Series:
    as_string = column.astype("string").str.strip()
    missing = column.isna() | as_string.isna() | (as_string == "")
    return as_string.where(~missing, other=pd.NA)


def _resolve_operation(plan: PreprocessingPlan) -> str:
    return plan.site_sequence_resolution_mode.value


def _resolve_parameters(plan: PreprocessingPlan) -> dict[str, object]:
    return {
        "enabled": bool(plan.site_sequence_resolution_enabled),
        "fasta_path": plan.site_sequence_resolution_fasta_path,
        "mode": plan.site_sequence_resolution_mode.value,
        "conflict_policy": plan.site_sequence_resolution_conflict_policy.value,
        "flank_size": int(plan.site_sequence_resolution_flank_size),
        "accession_column": plan.site_sequence_resolution_accession_column,
        "site_column": plan.site_sequence_resolution_site_column,
    }


def _include_when_enabled(plan: PreprocessingPlan) -> bool:
    return bool(plan.site_sequence_resolution_enabled)


SITE_SEQUENCE_RESOLUTION_STAGE_CONTRACT = PreprocessingStageContract(
    stage_key=DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION,
    display_label=DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION,
    provenance_stage=DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION,
    operation_name=_resolve_operation,
    serialize_parameters=_resolve_parameters,
    consumed_input_tables=(PreprocessingStateTableKey.DATASET_SITE_METADATA,),
    produced_output_tables=(PreprocessingStateTableKey.DATASET_SITE_METADATA,),
    stage_factory=SiteSequenceResolutionStage,
    backend="phospy.science.sequences",
    include_when=_include_when_enabled,
    diagnostics_metadata={
        "diagnostics_schema_version": 1,
        "known_diagnostics_fields": (
            "configured",
            "mode",
            "flank_size",
            "fasta_source_path",
            "fasta_source_label",
            "fasta_sha256",
            "resolver_version",
            "resolved_site_count",
            "unresolved_site_count",
            "unresolved_counts_by_reason",
            "filled_missing_count",
            "replaced_existing_count",
            "preserved_existing_count",
            "existing_sequence_conflict_count",
            "conflict_policy",
            "accession_column",
            "site_column",
            "row_status",
            "row_diagnostics",
        ),
        "known_row_diagnostic_fields": (
            "row_index",
            "row_id",
            "site_id",
            "status",
            "existing_site_sequence",
            "fasta_site_sequence",
            "resolved_site_sequence",
            "action",
            "reason",
            "conflict_policy",
            "resolver_version",
            "fasta_source_path",
            "fasta_sha256",
        ),
    },
)


__all__ = ["SITE_SEQUENCE_RESOLUTION_STAGE_CONTRACT", "SiteSequenceResolutionStage"]
