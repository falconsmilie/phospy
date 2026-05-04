"""Central registry for preprocessing stage metadata and provenance shape."""

from __future__ import annotations

from collections.abc import Callable, Mapping
from dataclasses import dataclass, field

from phospy.datasets.preprocessing.models import (
    DATASET_PREPROCESSING_STAGE_COMPARISONS,
    DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM,
    DATASET_PREPROCESSING_STAGE_MISSING_DATA,
    DATASET_PREPROCESSING_STAGE_NORMALISATION,
    DATASET_PREPROCESSING_STAGE_SITE_MATRIX,
    DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION,
    DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION,
    PreprocessingPlan,
)
from phospy.errors.build import DatasetBuildError

_ParameterSerializer = Callable[[PreprocessingPlan], dict[str, object]]
_OperationResolver = Callable[[PreprocessingPlan], str]


def _always_include(_plan: PreprocessingPlan) -> bool:
    return True


@dataclass(frozen=True, slots=True)
class PreprocessingStageMetadata:
    """Single canonical metadata record for a preprocessing stage."""

    stage_key: str
    display_label: str
    provenance_stage: str
    operation_name: _OperationResolver
    serialize_parameters: _ParameterSerializer
    consumed_input_tables: tuple[str, ...]
    produced_output_tables: tuple[str, ...]
    backend: str | None = None
    include_in_builder_provenance: bool = True
    include_when: Callable[[PreprocessingPlan], bool] = field(default=_always_include)
    diagnostics_metadata: Mapping[str, object] = field(default_factory=dict)


def _resolve_total_protein_parameters(plan: PreprocessingPlan) -> dict[str, object]:
    identity = plan.total_protein_correction_identity_policy
    return {
        "total_protein_correction_policy": plan.total_protein_correction_policy.value,
        "identity_mode": identity.mode,
        "phosphosite_key": identity.phosphosite_key,
        "total_protein_key": identity.total_protein_key,
        "mapping_phosphosite_key": identity.mapping_phosphosite_key,
        "mapping_total_protein_key": identity.mapping_total_protein_key,
        "mapping_table_fingerprint": identity.mapping_table_fingerprint,
        "mapping_table_row_count": (
            None if identity.mapping_table is None else len(identity.mapping_table)
        ),
        "duplicate_policy": identity.duplicate_policy,
        "unmatched_policy": identity.unmatched_policy,
    }


def _resolve_intensity_transform_parameters(
    plan: PreprocessingPlan,
) -> dict[str, object]:
    return {"pseudocount": float(plan.intensity_transform_pseudocount)}


def _resolve_site_sequence_parameters(plan: PreprocessingPlan) -> dict[str, object]:
    return {
        "enabled": bool(plan.site_sequence_resolution_enabled),
        "fasta_path": plan.site_sequence_resolution_fasta_path,
        "mode": plan.site_sequence_resolution_mode,
        "conflict_policy": plan.site_sequence_resolution_conflict_policy,
        "flank_size": int(plan.site_sequence_resolution_flank_size),
        "accession_column": plan.site_sequence_resolution_accession_column,
        "site_column": plan.site_sequence_resolution_site_column,
    }


def _resolve_missing_data_parameters(plan: PreprocessingPlan) -> dict[str, object]:
    return {
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
    }


def _resolve_site_matrix_parameters(plan: PreprocessingPlan) -> dict[str, object]:
    return {
        "site_matrix_policy": plan.site_matrix_policy,
        "site_matrix_duplicate_site_policy": plan.site_matrix_duplicate_site_policy,
        "site_matrix_missing_data_policy": plan.site_matrix_missing_data_policy,
        "site_matrix_minimum_observed_values": plan.site_matrix_minimum_observed_values,
    }


def _resolve_comparisons_parameters(plan: PreprocessingPlan) -> dict[str, object]:
    return {
        "comparison_building_policy": plan.comparison_building_policy,
        "comparison_sample_group_column": plan.comparison_sample_group_column,
        "comparison_pairs": plan.comparison_pairs,
    }


def _resolve_normalisation_parameters(_plan: PreprocessingPlan) -> dict[str, object]:
    return {}


def _resolve_missing_data_operation(plan: PreprocessingPlan) -> str:
    return plan.missing_data_policy.value


def _resolve_intensity_transform_operation(plan: PreprocessingPlan) -> str:
    return plan.intensity_transform_policy


def _resolve_total_protein_operation(plan: PreprocessingPlan) -> str:
    return plan.total_protein_correction_policy.value


def _resolve_site_matrix_operation(plan: PreprocessingPlan) -> str:
    return plan.site_matrix_policy


def _resolve_normalisation_operation(plan: PreprocessingPlan) -> str:
    return plan.normalisation_policy


def _resolve_comparisons_operation(plan: PreprocessingPlan) -> str:
    return plan.comparison_building_policy


def _resolve_site_sequence_operation(plan: PreprocessingPlan) -> str:
    return plan.site_sequence_resolution_mode


def _include_when_site_sequence_enabled(plan: PreprocessingPlan) -> bool:
    return bool(plan.site_sequence_resolution_enabled)


PREPROCESSING_STAGE_REGISTRY: tuple[PreprocessingStageMetadata, ...] = (
    PreprocessingStageMetadata(
        stage_key=DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION,
        display_label=DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION,
        provenance_stage=DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION,
        operation_name=_resolve_site_sequence_operation,
        serialize_parameters=_resolve_site_sequence_parameters,
        consumed_input_tables=("dataset.site_metadata",),
        produced_output_tables=("dataset.site_metadata",),
        backend="phospy.sequences",
        include_when=_include_when_site_sequence_enabled,
        diagnostics_metadata={"diagnostics_version": 1},
    ),
    PreprocessingStageMetadata(
        stage_key=DATASET_PREPROCESSING_STAGE_MISSING_DATA,
        display_label=DATASET_PREPROCESSING_STAGE_MISSING_DATA,
        provenance_stage=DATASET_PREPROCESSING_STAGE_MISSING_DATA,
        operation_name=_resolve_missing_data_operation,
        serialize_parameters=_resolve_missing_data_parameters,
        consumed_input_tables=("dataset.phospho", "dataset.site_metadata"),
        produced_output_tables=(
            "dataset.phospho",
            "dataset.site_metadata",
            "report.row_audit",
        ),
        backend="pandas",
        diagnostics_metadata={"diagnostics_version": 1},
    ),
    PreprocessingStageMetadata(
        stage_key=DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM,
        display_label=DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM,
        provenance_stage=DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM,
        operation_name=_resolve_intensity_transform_operation,
        serialize_parameters=_resolve_intensity_transform_parameters,
        consumed_input_tables=("dataset.phospho", "dataset.total"),
        produced_output_tables=("dataset.phospho", "dataset.total"),
        backend="numpy",
    ),
    PreprocessingStageMetadata(
        stage_key=DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION,
        display_label=DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION,
        provenance_stage=DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION,
        operation_name=_resolve_total_protein_operation,
        serialize_parameters=_resolve_total_protein_parameters,
        consumed_input_tables=(
            "dataset.phospho",
            "dataset.total",
            "dataset.site_metadata",
        ),
        produced_output_tables=("dataset.phospho",),
        backend="pandas",
        diagnostics_metadata={"diagnostics_version": 1},
    ),
    PreprocessingStageMetadata(
        stage_key=DATASET_PREPROCESSING_STAGE_SITE_MATRIX,
        display_label=DATASET_PREPROCESSING_STAGE_SITE_MATRIX,
        provenance_stage=DATASET_PREPROCESSING_STAGE_SITE_MATRIX,
        operation_name=_resolve_site_matrix_operation,
        serialize_parameters=_resolve_site_matrix_parameters,
        consumed_input_tables=("dataset.phospho", "dataset.site_metadata"),
        produced_output_tables=(
            "dataset.phospho",
            "dataset.site_metadata",
            "report.duplicate_site_resolution",
            "report.metadata_conflicts",
            "report.row_audit",
        ),
        backend="pandas",
    ),
    PreprocessingStageMetadata(
        stage_key=DATASET_PREPROCESSING_STAGE_NORMALISATION,
        display_label=DATASET_PREPROCESSING_STAGE_NORMALISATION,
        provenance_stage=DATASET_PREPROCESSING_STAGE_NORMALISATION,
        operation_name=_resolve_normalisation_operation,
        serialize_parameters=_resolve_normalisation_parameters,
        consumed_input_tables=("dataset.phospho",),
        produced_output_tables=("dataset.phospho",),
        backend="numpy",
    ),
    PreprocessingStageMetadata(
        stage_key=DATASET_PREPROCESSING_STAGE_COMPARISONS,
        display_label=DATASET_PREPROCESSING_STAGE_COMPARISONS,
        provenance_stage=DATASET_PREPROCESSING_STAGE_COMPARISONS,
        operation_name=_resolve_comparisons_operation,
        serialize_parameters=_resolve_comparisons_parameters,
        consumed_input_tables=("dataset.phospho", "dataset.sample_metadata"),
        produced_output_tables=(
            "dataset.comparisons",
            "report.comparison_group_stats",
            "report.comparison_pair_stats",
        ),
        backend="pandas",
    ),
)

_STAGE_METADATA_BY_KEY = {
    metadata.stage_key: metadata for metadata in PREPROCESSING_STAGE_REGISTRY
}


def get_preprocessing_stage_metadata(stage_key: str) -> PreprocessingStageMetadata:
    """Return metadata for ``stage_key`` or raise with an actionable error."""

    metadata = _STAGE_METADATA_BY_KEY.get(stage_key)
    if metadata is not None:
        return metadata
    supported = ", ".join(sorted(_STAGE_METADATA_BY_KEY))
    raise DatasetBuildError(
        "dataset preprocessing stage metadata is not registered for stage "
        f"{stage_key!r}; supported stages: {supported}"
    )


def resolve_builder_provenance_stage_order(
    plan: PreprocessingPlan,
) -> tuple[PreprocessingStageMetadata, ...]:
    """Return canonical stage metadata sequence used by builder provenance tables."""

    return tuple(
        metadata
        for metadata in PREPROCESSING_STAGE_REGISTRY
        if metadata.include_in_builder_provenance and metadata.include_when(plan)
    )


def list_registered_preprocessing_stages() -> tuple[PreprocessingStageMetadata, ...]:
    """Return the canonical preprocessing stage metadata registry."""

    return PREPROCESSING_STAGE_REGISTRY


def merge_preprocessing_stage_metadata(
    overrides: tuple[PreprocessingStageMetadata, ...] | None,
) -> dict[str, PreprocessingStageMetadata]:
    """Return default metadata map merged with optional custom overrides."""

    merged = dict(_STAGE_METADATA_BY_KEY)
    if overrides is None:
        return merged
    for metadata in overrides:
        merged[metadata.stage_key] = metadata
    return merged


__all__ = [
    "PreprocessingStageMetadata",
    "get_preprocessing_stage_metadata",
    "list_registered_preprocessing_stages",
    "merge_preprocessing_stage_metadata",
    "resolve_builder_provenance_stage_order",
]
