"""Internal executor for the dataset builder path.

The public builder lane stays intentionally narrow: establish supported
intensity scale state after applying explicit builder preprocessing policy.
"""

from __future__ import annotations

from collections.abc import Mapping

import pandas as pd

from phospy.datasets.builders.contracts import (
    DatasetPreprocessorContract,
    InterpretedDatasetBuildRequest,
    PreprocessedDatasetBuildTables,
)
from phospy.datasets.builders.preprocessing import (
    DatasetPreprocessor,
    build_dataset_processing_state,
)
from phospy.datasets.builders.transformation_resolver import (
    DatasetIntensityScaleResolver,
)
from phospy.datasets.models import (
    AnalysisReadyPhosphoDataset,
    DatasetPreprocessingReport,
    SiteSequenceResolutionReport,
)
from phospy.datasets.preprocessing.models import (
    DATASET_PREPROCESSING_STAGE_COMPARISONS,
    DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM,
    DATASET_PREPROCESSING_STAGE_LOCALISATION,
    DATASET_PREPROCESSING_STAGE_MISSING_DATA,
    DATASET_PREPROCESSING_STAGE_NORMALISATION,
    DATASET_PREPROCESSING_STAGE_ORDER_DEFAULT,
    DATASET_PREPROCESSING_STAGE_SITE_MATRIX,
    DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION,
    DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION,
    PreprocessingPlan,
    PreprocessingStageExecution,
    PreprocessingStageOrderResolution,
    TotalProteinCorrectionIdentityPolicy,
)
from phospy.datasets.preprocessing.report_schema import (
    PreprocessingOperationRow,
    PreprocessingRowCountRow,
    comparison_group_stats_rows_from_dataframe,
    comparison_pair_stats_rows_from_dataframe,
    duplicate_site_resolution_rows_from_dataframe,
    metadata_conflict_rows_from_dataframe,
    operation_rows_from_dataframe,
    row_audit_rows_from_dataframe,
    row_count_rows_from_dataframe,
)
from phospy.errors.build import DatasetBuildError
from phospy.errors.input import PhosPyInputError
from phospy.errors.transformations import (
    TransformationStateEstablishmentError,
)
from phospy.policy_models import IntensityTransformPolicy
from phospy.provenance.environment import collect_environment_provenance
from phospy.provenance.hashing import fingerprint_optional_table
from phospy.provenance.models import (
    PREPROCESSING_STAGE_DETERMINISM_PURE,
    JsonValue,
    PreprocessingStageProvenance,
    RunProvenance,
    TableFingerprint,
)
from phospy.scientific_policies import (
    PreprocessingStageOrderPolicy,
    ScientificPolicyRecord,
    build_duplicate_site_resolution_policy,
)
from phospy.site_ids import canonicalize_site_components, canonicalize_site_identifier
from phospy.transformations.contracts import Transformer
from phospy.transformations.models import IntensityScaleKind
from phospy.transformations.transformers import IdentityTransformer

_FINAL_DATASET_STAGE = "final_dataset_construction"
_SUPPORTED_PREPROCESSING_STAGE_ORDER = (
    DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION,
    DATASET_PREPROCESSING_STAGE_LOCALISATION,
    DATASET_PREPROCESSING_STAGE_MISSING_DATA,
    DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM,
    DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION,
    DATASET_PREPROCESSING_STAGE_SITE_MATRIX,
    DATASET_PREPROCESSING_STAGE_NORMALISATION,
    DATASET_PREPROCESSING_STAGE_COMPARISONS,
)
_SITE_SEQUENCE_RESOLUTION_STAGE = DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION
_SITE_SEQUENCE_FAILURE_CATEGORY_MISSING_REFERENCE_SUPPORT = "missing_reference_support"
_SITE_SEQUENCE_FAILURE_CATEGORY_AMBIGUOUS_MAPPING = "ambiguous_mapping"
_SITE_SEQUENCE_FAILURE_CATEGORY_INVALID_METADATA = "invalid_metadata"
_SITE_SEQUENCE_CONFLICT_POLICY_NOT_APPLIED = "not_applied"


class DatasetBuildExecutor:
    """Construct `AnalysisReadyPhosphoDataset` from interpreted builder input.

    Default policy uses the identity transformer, which is a pass-through
    establisher for already-prepared quantitative matrices after internal
    preprocessing stages (including optional site-matrix construction).
    """

    def __init__(
        self,
        *,
        transformer: Transformer | None = None,
        intensity_scale_resolver: DatasetIntensityScaleResolver | None = None,
        preprocessor: DatasetPreprocessorContract | None = None,
    ) -> None:
        self._intensity_scale_resolver = (
            intensity_scale_resolver
            or DatasetIntensityScaleResolver(
                transformer=transformer or IdentityTransformer()
            )
        )
        self._preprocessor = preprocessor or DatasetPreprocessor()

    def run(
        self, request: InterpretedDatasetBuildRequest
    ) -> AnalysisReadyPhosphoDataset:
        preprocessed = self._preprocessor.run(
            phospho=request.phospho,
            site_metadata=request.site_metadata,
            sample_metadata=request.sample_metadata,
            total=request.total,
            plan=request.preprocessing_plan,
        )
        validated_site_metadata = _require_analysis_ready_site_sequence_support(
            site_metadata=preprocessed.site_metadata,
            preprocessing_trace=preprocessed.preprocessing_trace,
        )
        resolved = self._intensity_scale_resolver.run(
            phospho=preprocessed.phospho,
            total=preprocessed.total,
            expected_scale_kind=_resolve_expected_intensity_scale_kind(
                request.preprocessing_plan
            ),
        )
        if not resolved.intensity_scale_state.is_established:
            raise TransformationStateEstablishmentError(
                "intensity-scale resolver returned a non-established "
                "intensity scale state; this violates the dataset boundary "
                "contract"
            )
        processing_state = build_dataset_processing_state(
            plan=request.preprocessing_plan,
            intensity_scale_state=resolved.intensity_scale_state,
            explicit_quantitative_meaning=request.quantitative_meaning,
            preprocessing_trace=preprocessed.preprocessing_trace,
            final_phospho=resolved.phospho,
            final_site_metadata=validated_site_metadata,
            final_sample_metadata=preprocessed.sample_metadata,
        )
        intensity_scale_state = processing_state.intensity_scale
        quantitative_meaning = intensity_scale_state.quantity
        if quantitative_meaning is None:
            raise DatasetBuildError(
                "intensity-scale state is missing quantitative meaning"
            )
        report = _build_dataset_preprocessing_report(
            row_counts=preprocessed.preprocessing_row_counts,
            operations=preprocessed.preprocessing_operations,
            row_audit=preprocessed.row_audit,
            duplicate_site_resolution=preprocessed.duplicate_site_resolution,
            metadata_conflicts=preprocessed.metadata_conflicts,
            comparison_group_stats=preprocessed.comparison_group_stats,
            comparison_pair_stats=preprocessed.comparison_pair_stats,
            preprocessing_trace=preprocessed.preprocessing_trace,
            site_sequence_derivation=request.site_sequence_derivation,
            input_site_count=int(request.site_metadata.shape[0]),
            final_dataset_rows=int(len(resolved.phospho.index)),
            intensity_scale_label=intensity_scale_state.label,
            quantitative_meaning=quantitative_meaning.value,
        )
        provenance = _build_dataset_run_provenance(
            request=request,
            preprocessed=preprocessed,
            validated_site_metadata=validated_site_metadata,
            resolved_phospho=resolved.phospho,
            resolved_total=resolved.total,
            preprocessing_trace=preprocessed.preprocessing_trace,
            intensity_scale_label=intensity_scale_state.label,
            quantitative_meaning=quantitative_meaning.value,
        )
        return AnalysisReadyPhosphoDataset._from_owned(
            phospho=resolved.phospho,
            site_metadata=validated_site_metadata,
            sample_metadata=preprocessed.sample_metadata,
            total=resolved.total,
            comparisons=preprocessed.comparisons,
            organism=request.organism,
            intensity_scale_state=intensity_scale_state,
            processing_state=processing_state,
            preprocessing_report=report,
            provenance=provenance,
        )


def _require_analysis_ready_site_sequence_support(
    *,
    site_metadata: pd.DataFrame,
    preprocessing_trace: tuple[PreprocessingStageExecution, ...] | None,
) -> pd.DataFrame:
    if "site_sequence" not in site_metadata.columns:
        _raise_site_sequence_boundary_error(
            site_metadata=site_metadata,
            invalid_index=list(range(int(len(site_metadata.index)))),
            preprocessing_trace=preprocessing_trace,
            reason_override="site_sequence column is missing after preprocessing",
        )
    sequence_column = site_metadata.loc[:, "site_sequence"]
    missing_or_blank_mask = sequence_column.isna()
    missing_or_blank_mask = missing_or_blank_mask | (
        sequence_column.astype("string").str.strip().isna()
    )
    missing_or_blank_mask = missing_or_blank_mask | (
        sequence_column.astype("string").str.strip() == ""
    )
    non_string_mask = ~sequence_column.map(lambda value: isinstance(value, str))
    non_string_mask = non_string_mask & ~sequence_column.isna()
    invalid_mask = missing_or_blank_mask | non_string_mask
    if not bool(invalid_mask.any()):
        return site_metadata
    invalid_positions = [
        position
        for position, flagged in enumerate(invalid_mask.tolist())
        if bool(flagged)
    ]
    _raise_site_sequence_boundary_error(
        site_metadata=site_metadata,
        invalid_index=invalid_positions,
        preprocessing_trace=preprocessing_trace,
        reason_override=None,
    )
    return site_metadata


def _raise_site_sequence_boundary_error(
    *,
    site_metadata: pd.DataFrame,
    invalid_index: list[int],
    preprocessing_trace: tuple[PreprocessingStageExecution, ...] | None,
    reason_override: str | None,
) -> None:
    context_by_row_id = _resolve_site_sequence_failure_context(preprocessing_trace)
    row_messages: list[str] = []
    for position in invalid_index[:8]:
        row_id = str(site_metadata.index[position])
        row = site_metadata.iloc[position]
        gene_symbol = _optional_row_value(row, "gene_symbol")
        site_value = _optional_row_value(row, "site")
        protein_id = _optional_row_value(row, "protein_id")
        category, identity_reason = _classify_identity_failure(
            row_id=row_id,
            gene_symbol=gene_symbol,
            site=site_value,
        )
        context = context_by_row_id.get(row_id)
        if context is None:
            reason = (
                identity_reason
                if reason_override is None
                else f"{reason_override}; {identity_reason}"
            )
            reason_category = category
        else:
            context_reason, context_category = context
            if reason_override is None:
                reason = context_reason
            else:
                reason = f"{reason_override}; {context_reason}"
            reason_category = context_category
        row_messages.append(
            f"site_id={row_id!r}, gene_symbol={gene_symbol!r}, site={site_value!r}, "
            f"protein_id={protein_id!r}, reason={reason!r}, failure_category={reason_category!r}"
        )
    suffix = "" if len(invalid_index) <= 8 else f"; +{len(invalid_index) - 8} more rows"
    details = " | ".join(row_messages)
    raise PhosPyInputError(
        "dataset builder cannot construct AnalysisReadyPhosphoDataset because "
        "dataset.site_metadata.site_sequence is missing, blank, or invalid after "
        f"builder enrichment; unresolved_rows={len(invalid_index)}; row_context={details}{suffix}"
    )


def _resolve_site_sequence_failure_context(
    preprocessing_trace: tuple[PreprocessingStageExecution, ...] | None,
) -> dict[str, tuple[str, str]]:
    if preprocessing_trace is None:
        return {}
    contexts: dict[str, tuple[str, str]] = {}
    for stage in preprocessing_trace:
        if stage.stage != _SITE_SEQUENCE_RESOLUTION_STAGE:
            continue
        diagnostics = stage.diagnostics if stage.diagnostics is not None else {}
        row_diagnostics = diagnostics.get("row_diagnostics")
        if not isinstance(row_diagnostics, list):
            continue
        for row in row_diagnostics:
            if not isinstance(row, dict):
                continue
            row_id_value = row.get("row_id")
            if row_id_value is None:
                continue
            row_id = str(row_id_value)
            resolved_value = row.get("resolved_site_sequence")
            if (
                resolved_value is not None
                and str(resolved_value).strip() != ""
                and str(resolved_value).strip().lower() != "none"
            ):
                continue
            status = str(row.get("status", "")).strip().lower()
            reason = str(row.get("reason", "")).strip()
            category = _classify_stage_failure_category(status, reason)
            resolved_reason = (
                reason
                if reason
                else "site_sequence resolution failed during FASTA enrichment"
            )
            contexts[row_id] = (resolved_reason, category)
    return contexts


def _classify_stage_failure_category(status: str, reason: str) -> str:
    combined = f"{status} {reason}".strip().lower()
    if "ambiguous" in combined:
        return _SITE_SEQUENCE_FAILURE_CATEGORY_AMBIGUOUS_MAPPING
    invalid_tokens = (
        "missing_accession",
        "missing_existing_sequence",
        "invalid",
        "blank",
        "empty",
    )
    if any(token in combined for token in invalid_tokens):
        return _SITE_SEQUENCE_FAILURE_CATEGORY_INVALID_METADATA
    return _SITE_SEQUENCE_FAILURE_CATEGORY_MISSING_REFERENCE_SUPPORT


def _classify_identity_failure(
    *,
    row_id: str,
    gene_symbol: str | None,
    site: str | None,
) -> tuple[str, str]:
    metadata_site_id: str | None = None
    index_site_id: str | None = None
    metadata_error: str | None = None
    index_error: str | None = None
    if gene_symbol is not None and site is not None:
        try:
            metadata_site_id = canonicalize_site_components(
                gene_symbol,
                site,
                field_name=(
                    "dataset build request site_metadata.gene_symbol/site for "
                    "analysis-ready sequence enforcement"
                ),
                error_type=PhosPyInputError,
            )
        except PhosPyInputError as exc:
            metadata_error = str(exc)
    else:
        metadata_error = "missing gene_symbol/site metadata"
    try:
        if ";" in str(row_id):
            index_site_id = canonicalize_site_identifier(
                row_id,
                field_name=(
                    "dataset build request site_metadata.index for analysis-ready "
                    "sequence enforcement"
                ),
                error_type=PhosPyInputError,
            )
    except PhosPyInputError as exc:
        index_error = str(exc)
    if metadata_site_id is not None and index_site_id is not None:
        if metadata_site_id != index_site_id:
            return (
                _SITE_SEQUENCE_FAILURE_CATEGORY_AMBIGUOUS_MAPPING,
                "gene/site metadata and site index map to different canonical sites",
            )
        return (
            _SITE_SEQUENCE_FAILURE_CATEGORY_MISSING_REFERENCE_SUPPORT,
            "site_sequence is unresolved for canonical site despite coherent metadata",
        )
    if metadata_site_id is not None or index_site_id is not None:
        return (
            _SITE_SEQUENCE_FAILURE_CATEGORY_MISSING_REFERENCE_SUPPORT,
            "site_sequence is unresolved for canonical site despite partial identity metadata",
        )
    detail_chunks = [chunk for chunk in (metadata_error, index_error) if chunk]
    detail = (
        "; ".join(detail_chunks) if detail_chunks else "invalid site identity metadata"
    )
    return (_SITE_SEQUENCE_FAILURE_CATEGORY_INVALID_METADATA, detail)


def _optional_row_value(row: pd.Series, column_name: str) -> str | None:
    if column_name not in row.index:
        return None
    value = row[column_name]
    if pd.isna(value):
        return None
    if isinstance(value, str):
        stripped = value.strip()
        return stripped if stripped else None
    return str(value)


def _build_dataset_preprocessing_report(
    *,
    row_counts: pd.DataFrame | None,
    operations: pd.DataFrame | None,
    row_audit: pd.DataFrame | None,
    duplicate_site_resolution: pd.DataFrame | None,
    metadata_conflicts: pd.DataFrame | None,
    comparison_group_stats: pd.DataFrame | None,
    comparison_pair_stats: pd.DataFrame | None,
    preprocessing_trace: tuple[PreprocessingStageExecution, ...] | None,
    site_sequence_derivation: dict[str, object] | None,
    input_site_count: int,
    final_dataset_rows: int,
    intensity_scale_label: str,
    quantitative_meaning: str,
) -> DatasetPreprocessingReport:
    row_count_rows = list(row_count_rows_from_dataframe(row_counts))
    operation_rows = list(operation_rows_from_dataframe(operations))
    row_audit_rows = row_audit_rows_from_dataframe(row_audit)
    duplicate_site_resolution_rows = duplicate_site_resolution_rows_from_dataframe(
        duplicate_site_resolution
    )
    metadata_conflict_rows = metadata_conflict_rows_from_dataframe(metadata_conflicts)
    comparison_group_stats_rows = comparison_group_stats_rows_from_dataframe(
        comparison_group_stats
    )
    comparison_pair_stats_rows = comparison_pair_stats_rows_from_dataframe(
        comparison_pair_stats
    )

    row_count_rows.append(
        PreprocessingRowCountRow(
            stage=_FINAL_DATASET_STAGE,
            input_rows=final_dataset_rows,
            output_rows=final_dataset_rows,
            dropped_rows=0,
        )
    )
    if not operation_rows:
        final_step_order = 1
    else:
        final_step_order = int(max(row.step_order for row in operation_rows)) + 1
    operation_rows.append(
        PreprocessingOperationRow(
            step_order=final_step_order,
            stage=_FINAL_DATASET_STAGE,
            operation="construct_analysis_ready_dataset",
            parameters={
                "intensity_scale_label": intensity_scale_label,
                "quantitative_meaning": quantitative_meaning,
            },
            input_rows=final_dataset_rows,
            output_rows=final_dataset_rows,
            notes="analysis-ready dataset boundary construction",
        )
    )
    site_sequence_resolution = _build_site_sequence_resolution_report(
        preprocessing_trace=preprocessing_trace,
        site_sequence_derivation=site_sequence_derivation,
        total_sites=int(input_site_count),
        final_sequence_complete_sites=int(final_dataset_rows),
    )
    return DatasetPreprocessingReport.from_rows(
        row_count_rows=tuple(row_count_rows),
        operation_rows=tuple(operation_rows),
        row_audit_rows=row_audit_rows,
        duplicate_site_resolution_rows=duplicate_site_resolution_rows,
        metadata_conflict_rows=metadata_conflict_rows,
        comparison_group_stats_rows=comparison_group_stats_rows,
        comparison_pair_stats_rows=comparison_pair_stats_rows,
        site_sequence_resolution=site_sequence_resolution,
    )


def _build_site_sequence_resolution_report(
    *,
    preprocessing_trace: tuple[PreprocessingStageExecution, ...] | None,
    site_sequence_derivation: dict[str, object] | None,
    total_sites: int,
    final_sequence_complete_sites: int,
) -> SiteSequenceResolutionReport:
    stage_diagnostics = _resolve_site_sequence_stage_diagnostics(preprocessing_trace)
    if stage_diagnostics is not None:
        (
            provided_by_input,
            resolved_from_fasta,
            unresolved,
        ) = _summarize_stage_sequence_origins(stage_diagnostics)
        conflicts = _coerce_non_negative_int(
            stage_diagnostics.get("existing_sequence_conflict_count"),
            default=0,
        )
        conflict_policy = _resolve_conflict_policy(
            stage_diagnostics.get("conflict_policy")
        )
        return SiteSequenceResolutionReport(
            total_sites=int(max(total_sites, 0)),
            provided_by_input=int(max(provided_by_input, 0)),
            resolved_from_fasta=int(max(resolved_from_fasta, 0)),
            resolved_from_reference=0,
            unresolved=int(max(unresolved, 0)),
            conflicts=int(max(conflicts, 0)),
            conflict_policy=conflict_policy,
            final_sequence_complete_sites=int(max(final_sequence_complete_sites, 0)),
        )

    derivation = (
        {}
        if not isinstance(site_sequence_derivation, Mapping)
        else site_sequence_derivation
    )
    provided_by_input = _coerce_non_negative_int(
        derivation.get("provided_sequence_count"),
        default=0,
    )
    resolved_from_reference = _coerce_non_negative_int(
        derivation.get("derived_sequence_count"),
        default=0,
    )
    unresolved = _coerce_non_negative_int(
        derivation.get("unresolved_sequence_count"),
        default=max(total_sites - provided_by_input - resolved_from_reference, 0),
    )
    conflicts = _coerce_non_negative_int(
        derivation.get("existing_sequence_conflict_count"),
        default=0,
    )
    return SiteSequenceResolutionReport(
        total_sites=int(max(total_sites, 0)),
        provided_by_input=int(max(provided_by_input, 0)),
        resolved_from_fasta=0,
        resolved_from_reference=int(max(resolved_from_reference, 0)),
        unresolved=int(max(unresolved, 0)),
        conflicts=int(max(conflicts, 0)),
        conflict_policy=_SITE_SEQUENCE_CONFLICT_POLICY_NOT_APPLIED,
        final_sequence_complete_sites=int(max(final_sequence_complete_sites, 0)),
    )


def _resolve_site_sequence_stage_diagnostics(
    preprocessing_trace: tuple[PreprocessingStageExecution, ...] | None,
) -> dict[str, object] | None:
    if preprocessing_trace is None:
        return None
    for stage in preprocessing_trace:
        if stage.stage == _SITE_SEQUENCE_RESOLUTION_STAGE:
            return (
                {}
                if not isinstance(stage.diagnostics, Mapping)
                else dict(stage.diagnostics)
            )
    return None


def _summarize_stage_sequence_origins(
    stage_diagnostics: Mapping[str, object],
) -> tuple[int, int, int]:
    row_diagnostics = stage_diagnostics.get("row_diagnostics")
    if not isinstance(row_diagnostics, list):
        provided_by_input = _coerce_non_negative_int(
            stage_diagnostics.get("preserved_existing_count"),
            default=0,
        )
        resolved_from_fasta = _coerce_non_negative_int(
            stage_diagnostics.get("filled_missing_count"),
            default=0,
        ) + _coerce_non_negative_int(
            stage_diagnostics.get("replaced_existing_count"),
            default=0,
        )
        unresolved = _coerce_non_negative_int(
            stage_diagnostics.get("unresolved_site_count"),
            default=0,
        )
        return (provided_by_input, resolved_from_fasta, unresolved)

    provided_by_input = 0
    resolved_from_fasta = 0
    unresolved = 0
    for row in row_diagnostics:
        if not isinstance(row, Mapping):
            continue
        action = str(row.get("action", "")).strip().lower()
        existing_site_sequence = row.get("existing_site_sequence")
        resolved_site_sequence = row.get("resolved_site_sequence")
        has_existing = _has_resolved_site_sequence(existing_site_sequence)
        has_resolved = _has_resolved_site_sequence(resolved_site_sequence)

        if action in {"fill_missing", "replace_existing"} and has_resolved:
            resolved_from_fasta += 1
            continue
        if not has_resolved:
            unresolved += 1
            continue
        if has_existing:
            provided_by_input += 1
            continue
        if action in {"validate_existing", "preserve_existing"}:
            provided_by_input += 1
            continue
        resolved_from_fasta += 1
    return (provided_by_input, resolved_from_fasta, unresolved)


def _resolve_conflict_policy(value: object) -> str:
    if isinstance(value, str):
        normalized = value.strip()
        if normalized:
            return normalized
    return _SITE_SEQUENCE_CONFLICT_POLICY_NOT_APPLIED


def _has_resolved_site_sequence(value: object) -> bool:
    if value is None:
        return False
    text = str(value).strip()
    if not text:
        return False
    return text.lower() != "none"


def _coerce_non_negative_int(value: object, *, default: int) -> int:
    if isinstance(value, bool):
        return int(default)
    if isinstance(value, int):
        return max(int(value), 0)
    return int(max(default, 0))


def _resolve_expected_intensity_scale_kind(
    preprocessing_plan: PreprocessingPlan,
) -> IntensityScaleKind:
    if preprocessing_plan.intensity_transform_policy is IntensityTransformPolicy.LOG2:
        return IntensityScaleKind.LOG2
    return IntensityScaleKind.LINEAR


def _build_dataset_run_provenance(
    *,
    request: InterpretedDatasetBuildRequest,
    preprocessed: PreprocessedDatasetBuildTables,
    validated_site_metadata: pd.DataFrame,
    resolved_phospho: pd.DataFrame,
    resolved_total: pd.DataFrame | None,
    preprocessing_trace: tuple[PreprocessingStageExecution, ...] | None,
    intensity_scale_label: str,
    quantitative_meaning: str,
) -> RunProvenance:
    input_tables = _collect_fingerprints(
        (
            ("dataset.phospho", request.phospho),
            ("dataset.site_metadata", request.site_metadata),
            ("dataset.sample_metadata", request.sample_metadata),
            ("dataset.total", request.total),
        )
    )
    output_tables = _collect_fingerprints(
        (
            ("dataset.phospho", resolved_phospho),
            ("dataset.site_metadata", validated_site_metadata),
            ("dataset.sample_metadata", preprocessed.sample_metadata),
            ("dataset.total", resolved_total),
            ("dataset.comparisons", preprocessed.comparisons),
        )
    )
    return RunProvenance(
        environment=collect_environment_provenance(),
        input_tables=input_tables,
        preprocessing_stages=_stage_trace_to_provenance(preprocessing_trace),
        reference=None,
        workflow_name="dataset_builder",
        workflow_parameters={
            "preprocessing_plan": _preprocessing_plan_to_payload(
                request.preprocessing_plan
            ),
            "intensity_scale_label": intensity_scale_label,
            "quantitative_meaning": quantitative_meaning,
            "site_identifier_normalisation": (
                None
                if request.site_identifier_normalisation is None
                else request.site_identifier_normalisation.to_payload()
            ),
            "site_sequence_derivation": request.site_sequence_derivation,
        },
        random_state=None,
        random_seed_policy=None,
        output_tables=output_tables,
        scientific_policies=_dataset_scientific_policies(request.preprocessing_plan),
    )


def _collect_fingerprints(
    entries: tuple[tuple[str, pd.DataFrame | None], ...],
) -> tuple[TableFingerprint, ...]:
    fingerprints: list[TableFingerprint] = []
    for name, table in entries:
        fingerprint = fingerprint_optional_table(table, name=name)
        if fingerprint is None:
            continue
        fingerprints.append(fingerprint)
    return tuple(fingerprints)


def _stage_trace_to_provenance(
    trace: tuple[PreprocessingStageExecution, ...] | None,
) -> tuple[PreprocessingStageProvenance, ...]:
    if trace is None:
        return ()
    return tuple(
        PreprocessingStageProvenance(
            stage=item.stage,
            operation=item.operation,
            parameters=dict(item.parameters),
            input_shape=item.input_shape,
            output_shape=item.output_shape,
            input_hash=item.input_hash,
            output_hash=item.output_hash,
            phospho_input_hash=item.phospho_input_hash,
            phospho_output_hash=item.phospho_output_hash,
            dropped_row_ids=item.dropped_row_ids,
            dropped_row_count=int(item.dropped_row_count),
            schema_version=int(item.schema_version),
            consumed_input_tables=tuple(item.consumed_input_tables),
            produced_output_tables=tuple(item.produced_output_tables),
            backend=item.backend,
            random_seed=item.random_seed,
            determinism=(
                str(item.determinism).strip()
                if str(item.determinism).strip()
                else PREPROCESSING_STAGE_DETERMINISM_PURE
            ),
            is_deterministic=bool(item.is_deterministic),
            imputed_cell_count=int(item.imputed_cell_count),
            imputed_row_ids=item.imputed_row_ids,
            notes=item.notes,
            diagnostics=_to_json_mapping(item.diagnostics),
        )
        for item in trace
    )


def _to_json_mapping(values: Mapping[str, object]) -> dict[str, JsonValue]:
    return {str(key): _to_json_value(value) for key, value in values.items()}


def _to_json_value(value: object) -> JsonValue:
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, dict):
        return {str(key): _to_json_value(item) for key, item in value.items()}
    if isinstance(value, tuple):
        return [_to_json_value(item) for item in value]
    if isinstance(value, list):
        return [_to_json_value(item) for item in value]
    return str(value)


def _preprocessing_plan_to_payload(plan: PreprocessingPlan) -> dict[str, object]:
    payload: dict[str, object] = {
        "intensity_transform_policy": plan.intensity_transform_policy.value,
        "intensity_transform_pseudocount": float(plan.intensity_transform_pseudocount),
        "normalisation_policy": plan.normalisation_policy.value,
        "missing_data_policy": plan.missing_data_policy.value,
        "missing_data_min_observed_values": plan.missing_data_min_observed_values,
        "missing_data_q": plan.missing_data_q,
        "missing_data_width": plan.missing_data_width,
        "missing_data_seed": plan.missing_data_seed,
        "missing_data_k": plan.missing_data_k,
        "missing_data_distance": plan.missing_data_distance,
        "missing_data_max_missing_fraction_per_row": (
            plan.missing_data_max_missing_fraction_per_row
        ),
        "localisation_mode": plan.localisation_mode.value,
        "localisation_min_confidence": float(plan.localisation_min_confidence),
        "localisation_confidence_column": plan.localisation_confidence_column,
        "localisation_waiver_reason": plan.localisation_waiver_reason,
        "site_sequence_resolution_enabled": plan.site_sequence_resolution_enabled,
        "site_sequence_resolution_fasta_path": (
            plan.site_sequence_resolution_fasta_path
        ),
        "site_sequence_resolution_mode": plan.site_sequence_resolution_mode.value,
        "site_sequence_resolution_flank_size": int(
            plan.site_sequence_resolution_flank_size
        ),
        "site_sequence_resolution_accession_column": (
            plan.site_sequence_resolution_accession_column
        ),
        "site_sequence_resolution_site_column": (
            plan.site_sequence_resolution_site_column
        ),
        "total_protein_correction_policy": plan.total_protein_correction_policy.value,
        "total_protein_correction_identity_policy": (
            _total_correction_identity_policy_to_payload(
                plan.total_protein_correction_identity_policy
            )
        ),
        "site_matrix_policy": plan.site_matrix_policy.value,
        "comparison_building_policy": plan.comparison_building_policy.value,
        "site_matrix_duplicate_site_policy": plan.site_matrix_duplicate_site_policy.value,
        "site_matrix_missing_data_policy": plan.site_matrix_missing_data_policy.value,
        "site_matrix_minimum_observed_values": plan.site_matrix_minimum_observed_values,
        "comparison_sample_group_column": plan.comparison_sample_group_column,
        "comparison_pairs": (
            None if plan.comparison_pairs is None else list(plan.comparison_pairs)
        ),
        "ruv_readiness_enabled": bool(plan.ruv_readiness_enabled),
        "ruv_readiness_control_feature_column": (
            plan.ruv_readiness_control_feature_column
        ),
        "ruv_readiness_replicate_group_column": (
            plan.ruv_readiness_replicate_group_column
        ),
        "ruv_readiness_batch_column": plan.ruv_readiness_batch_column,
        "stage_order": list(plan.stage_order),
        "resolved_stage_order": _stage_order_resolution_to_payload(
            plan.stage_order_resolution
        ),
    }
    return payload


def _stage_order_resolution_to_payload(
    stage_order_resolution: tuple[PreprocessingStageOrderResolution, ...],
) -> list[dict[str, object]]:
    return [
        {
            "stage": item.stage,
            "order_index": int(item.order_index),
            "rationale": str(item.rationale),
        }
        for item in stage_order_resolution
    ]


def _total_correction_identity_policy_to_payload(
    policy: TotalProteinCorrectionIdentityPolicy,
) -> dict[str, object]:
    return {
        "mode": str(policy.mode),
        "matching_policy": str(policy.matching_policy),
        "phosphosite_key": policy.phosphosite_key,
        "total_protein_key": policy.total_protein_key,
        "mapping_phosphosite_key": policy.mapping_phosphosite_key,
        "mapping_total_protein_key": policy.mapping_total_protein_key,
        "mapping_table_fingerprint": policy.mapping_table_fingerprint,
        "mapping_table_row_count": (
            None if policy.mapping_table is None else int(len(policy.mapping_table))
        ),
        "duplicate_policy": str(policy.duplicate_policy),
        "unmatched_policy": str(policy.unmatched_policy),
    }


def _dataset_scientific_policies(
    preprocessing_plan: PreprocessingPlan,
) -> tuple[ScientificPolicyRecord, ...]:
    policies = [
        PreprocessingStageOrderPolicy(
            configured_stage_order=tuple(
                str(stage) for stage in preprocessing_plan.stage_order
            ),
            default_stage_order=tuple(
                str(stage) for stage in DATASET_PREPROCESSING_STAGE_ORDER_DEFAULT
            ),
            supported_stage_order=_SUPPORTED_PREPROCESSING_STAGE_ORDER,
        ).record,
    ]
    if DATASET_PREPROCESSING_STAGE_SITE_MATRIX in preprocessing_plan.stage_order:
        policies.append(
            build_duplicate_site_resolution_policy(
                duplicate_site_policy=(
                    preprocessing_plan.site_matrix_duplicate_site_policy.value
                )
            )
        )
    return tuple(policies)
