"""Internal executor for the dataset builder path.

The public builder lane stays intentionally narrow: establish supported
intensity scale state after applying explicit builder preprocessing policy.
"""

from __future__ import annotations

from dataclasses import asdict

import pandas as pd

from phospy.api.configs import DATASET_INTENSITY_TRANSFORM_POLICY_LOG2
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
)
from phospy.datasets.preprocessing.models import (
    PreprocessingPlan,
    PreprocessingStageExecution,
)
from phospy.errors.build import DatasetBuildError
from phospy.errors.input import PhosPyInputError
from phospy.errors.transformations import (
    PhosPyTransformationError,
    TransformationStateEstablishmentError,
)
from phospy.errors.validation import PhosPyValidationError
from phospy.provenance.environment import collect_environment_provenance
from phospy.provenance.hashing import fingerprint_optional_table
from phospy.provenance.models import (
    PreprocessingStageProvenance,
    RunProvenance,
    TableFingerprint,
)
from phospy.transformations.contracts import Transformer
from phospy.transformations.models import IntensityScaleKind
from phospy.transformations.transformers import IdentityTransformer

_ROW_COUNT_COLUMNS = ("stage", "input_rows", "output_rows", "dropped_rows")
_OPERATION_COLUMNS = (
    "step_order",
    "stage",
    "operation",
    "parameters",
    "input_rows",
    "output_rows",
    "notes",
)
_ROW_AUDIT_COLUMNS = (
    "stage",
    "action",
    "reason",
    "source_row_id",
    "site_id",
    "retained",
    "retained_row_id",
    "source_rows",
    "retained_row",
    "parameter_snapshot",
)
_DUPLICATE_SITE_RESOLUTION_COLUMNS = (
    "site_id",
    "source_row_id",
    "retained",
    "resolution_policy",
    "retained_reason",
    "dropped_reason",
    "observed_values",
    "mean_signal",
    "n_source_rows",
    "n_aggregated_rows",
    "source_protein_id",
    "source_gene_symbol",
    "source_site",
    "source_site_sequence",
    "metadata_conflict_detected",
)
_METADATA_CONFLICT_COLUMNS = (
    "site_id",
    "field",
    "values",
    "n_distinct_values",
    "source_row_ids",
)
_COMPARISON_GROUP_STATS_COLUMNS = (
    "site_id",
    "group",
    "n",
    "mean",
    "sd",
    "sem",
    "median",
    "min",
    "max",
    "sample_ids",
)
_COMPARISON_PAIR_STATS_COLUMNS = (
    "site_id",
    "comparison",
    "left_group",
    "right_group",
    "left_n",
    "right_n",
    "left_mean",
    "right_mean",
    "left_sd",
    "right_sd",
    "left_sem",
    "right_sem",
    "effect_size",
    "left_median",
    "right_median",
    "left_min",
    "right_min",
    "left_max",
    "right_max",
)
_FINAL_DATASET_STAGE = "final_dataset_construction"


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
        try:
            preprocessed = self._preprocessor.run(
                phospho=request.phospho,
                site_metadata=request.site_metadata,
                sample_metadata=request.sample_metadata,
                total=request.total,
                plan=request.preprocessing_plan,
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
            )
            report = _build_dataset_preprocessing_report(
                row_counts=preprocessed.preprocessing_row_counts,
                operations=preprocessed.preprocessing_operations,
                row_audit=preprocessed.row_audit,
                duplicate_site_resolution=preprocessed.duplicate_site_resolution,
                metadata_conflicts=preprocessed.metadata_conflicts,
                comparison_group_stats=preprocessed.comparison_group_stats,
                comparison_pair_stats=preprocessed.comparison_pair_stats,
                final_dataset_rows=int(len(resolved.phospho.index)),
                intensity_scale_label=resolved.intensity_scale_state.label,
            )
            provenance = _build_dataset_run_provenance(
                request=request,
                preprocessed=preprocessed,
                resolved_phospho=resolved.phospho,
                resolved_total=resolved.total,
                preprocessing_trace=preprocessed.preprocessing_trace,
                intensity_scale_label=resolved.intensity_scale_state.label,
            )
            return AnalysisReadyPhosphoDataset._from_owned(
                phospho=resolved.phospho,
                site_metadata=preprocessed.site_metadata,
                sample_metadata=preprocessed.sample_metadata,
                total=resolved.total,
                comparisons=preprocessed.comparisons,
                organism=request.organism,
                intensity_scale_state=resolved.intensity_scale_state,
                processing_state=processing_state,
                preprocessing_report=report,
                provenance=provenance,
            )
        except (
            PhosPyInputError,
            PhosPyTransformationError,
            PhosPyValidationError,
            DatasetBuildError,
        ):
            raise
        except Exception as exc:  # pragma: no cover - defensive boundary translation
            raise DatasetBuildError(
                "failed to construct AnalysisReadyPhosphoDataset from interpreted input"
            ) from exc


def _build_dataset_preprocessing_report(
    *,
    row_counts: pd.DataFrame | None,
    operations: pd.DataFrame | None,
    row_audit: pd.DataFrame | None,
    duplicate_site_resolution: pd.DataFrame | None,
    metadata_conflicts: pd.DataFrame | None,
    comparison_group_stats: pd.DataFrame | None,
    comparison_pair_stats: pd.DataFrame | None,
    final_dataset_rows: int,
    intensity_scale_label: str,
) -> DatasetPreprocessingReport:
    base_row_counts = (
        pd.DataFrame.from_records([], columns=_ROW_COUNT_COLUMNS)
        if row_counts is None
        else row_counts
    )
    base_operations = (
        pd.DataFrame.from_records([], columns=_OPERATION_COLUMNS)
        if operations is None
        else operations
    )
    base_row_audit = (
        pd.DataFrame.from_records([], columns=_ROW_AUDIT_COLUMNS)
        if row_audit is None
        else row_audit
    )
    base_duplicate_site_resolution = (
        pd.DataFrame.from_records([], columns=_DUPLICATE_SITE_RESOLUTION_COLUMNS)
        if duplicate_site_resolution is None
        else duplicate_site_resolution
    )
    base_metadata_conflicts = (
        pd.DataFrame.from_records([], columns=_METADATA_CONFLICT_COLUMNS)
        if metadata_conflicts is None
        else metadata_conflicts
    )
    base_comparison_group_stats = (
        pd.DataFrame.from_records([], columns=_COMPARISON_GROUP_STATS_COLUMNS)
        if comparison_group_stats is None
        else comparison_group_stats
    )
    base_comparison_pair_stats = (
        pd.DataFrame.from_records([], columns=_COMPARISON_PAIR_STATS_COLUMNS)
        if comparison_pair_stats is None
        else comparison_pair_stats
    )
    final_row_counts = pd.concat(
        [
            base_row_counts,
            pd.DataFrame.from_records(
                [
                    {
                        "stage": _FINAL_DATASET_STAGE,
                        "input_rows": final_dataset_rows,
                        "output_rows": final_dataset_rows,
                        "dropped_rows": 0,
                    }
                ],
                columns=_ROW_COUNT_COLUMNS,
            ),
        ],
        axis=0,
        ignore_index=True,
    )
    if base_operations.empty:
        final_step_order = 1
    else:
        final_step_order = int(base_operations.loc[:, "step_order"].max()) + 1
    final_operations = pd.concat(
        [
            base_operations,
            pd.DataFrame.from_records(
                [
                    {
                        "step_order": final_step_order,
                        "stage": _FINAL_DATASET_STAGE,
                        "operation": "construct_analysis_ready_dataset",
                        "parameters": {"intensity_scale_label": intensity_scale_label},
                        "input_rows": final_dataset_rows,
                        "output_rows": final_dataset_rows,
                        "notes": "analysis-ready dataset boundary construction",
                    }
                ],
                columns=_OPERATION_COLUMNS,
            ),
        ],
        axis=0,
        ignore_index=True,
    )
    return DatasetPreprocessingReport._from_owned(
        row_counts=final_row_counts,
        operations=final_operations,
        row_audit=base_row_audit,
        duplicate_site_resolution=base_duplicate_site_resolution,
        metadata_conflicts=base_metadata_conflicts,
        comparison_group_stats=base_comparison_group_stats,
        comparison_pair_stats=base_comparison_pair_stats,
    )


def _resolve_expected_intensity_scale_kind(
    preprocessing_plan: PreprocessingPlan,
) -> IntensityScaleKind:
    if (
        preprocessing_plan.intensity_transform_policy
        == DATASET_INTENSITY_TRANSFORM_POLICY_LOG2
    ):
        return IntensityScaleKind.LOG2
    return IntensityScaleKind.LINEAR


def _build_dataset_run_provenance(
    *,
    request: InterpretedDatasetBuildRequest,
    preprocessed: PreprocessedDatasetBuildTables,
    resolved_phospho: pd.DataFrame,
    resolved_total: pd.DataFrame | None,
    preprocessing_trace: tuple[PreprocessingStageExecution, ...] | None,
    intensity_scale_label: str,
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
            ("dataset.site_metadata", preprocessed.site_metadata),
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
            "preprocessing_plan": asdict(request.preprocessing_plan),
            "intensity_scale_label": intensity_scale_label,
        },
        random_state=None,
        random_seed_policy=None,
        output_tables=output_tables,
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
            dropped_row_ids=item.dropped_row_ids,
            dropped_row_count=int(item.dropped_row_count),
            imputed_cell_count=int(item.imputed_cell_count),
            imputed_row_ids=item.imputed_row_ids,
            notes=item.notes,
            diagnostics=dict(item.diagnostics),
        )
        for item in trace
    )
