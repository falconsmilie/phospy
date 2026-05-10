"""Central registry for preprocessing stage metadata and provenance shape."""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field

from phospy.datasets._processing_state.json_contracts import (
    MISSING_DATA_DIAGNOSTICS_SCHEMA_VERSION_V1,
    TOTAL_PROTEIN_CORRECTION_DIAGNOSTICS_SCHEMA_VERSION_V1,
    V1_KNOWN_MISSING_DATA_DIAGNOSTICS_FIELDS,
    V1_KNOWN_TOTAL_PROTEIN_DIAGNOSTICS_FIELDS,
)
from phospy.datasets.preprocessing.models import (
    DATASET_PREPROCESSING_STAGE_COMPARISONS,
    DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM,
    DATASET_PREPROCESSING_STAGE_MISSING_DATA,
    DATASET_PREPROCESSING_STAGE_NORMALISATION,
    DATASET_PREPROCESSING_STAGE_SITE_MATRIX,
    DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION,
    DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION,
    PREPROCESSING_STATE_TABLE_KEYS,
    PreprocessingPlan,
    PreprocessingStage,
    PreprocessingStateTableKey,
)
from phospy.datasets.preprocessing.stages.comparisons import ComparisonsStage
from phospy.datasets.preprocessing.stages.intensity_transform import (
    IntensityTransformStage,
)
from phospy.datasets.preprocessing.stages.missing_data import MissingDataStage
from phospy.datasets.preprocessing.stages.normalisation import NormalisationStage
from phospy.datasets.preprocessing.stages.site_matrix import SiteMatrixStage
from phospy.datasets.preprocessing.stages.site_sequence_resolution import (
    SiteSequenceResolutionStage,
)
from phospy.datasets.preprocessing.stages.total_protein_correction import (
    TotalProteinCorrectionStage,
)
from phospy.errors.build import DatasetBuildError

_ParameterSerializer = Callable[[PreprocessingPlan], dict[str, object]]
_OperationResolver = Callable[[PreprocessingPlan], str]
_StageFactory = Callable[[], PreprocessingStage]


def _always_include(_plan: PreprocessingPlan) -> bool:
    return True


@dataclass(frozen=True, slots=True)
class PreprocessingStageMetadata:
    """Single canonical metadata record for a preprocessing stage."""

    stage_key: str
    display_label: str
    operation_name: _OperationResolver
    serialize_parameters: _ParameterSerializer
    consumed_input_tables: tuple[PreprocessingStateTableKey, ...]
    produced_output_tables: tuple[PreprocessingStateTableKey, ...]
    stage_factory: _StageFactory | None = None
    provenance_stage: str | None = None
    backend: str | None = None
    include_in_builder_provenance: bool = True
    include_when: Callable[[PreprocessingPlan], bool] = field(default=_always_include)
    diagnostics_metadata: Mapping[str, object] = field(default_factory=dict)

    def __post_init__(self) -> None:
        normalized_stage_key = str(self.stage_key).strip()
        object.__setattr__(self, "stage_key", normalized_stage_key)
        object.__setattr__(
            self,
            "consumed_input_tables",
            _normalize_stage_table_keys(
                stage_key=normalized_stage_key,
                table_keys=self.consumed_input_tables,
                role="consumed_input_tables",
            ),
        )
        object.__setattr__(
            self,
            "produced_output_tables",
            _normalize_stage_table_keys(
                stage_key=normalized_stage_key,
                table_keys=self.produced_output_tables,
                role="produced_output_tables",
            ),
        )

    @property
    def provenance_stage_key(self) -> str:
        if self.provenance_stage is None:
            return self.stage_key
        normalized = self.provenance_stage.strip()
        if not normalized:
            return self.stage_key
        return normalized


_SITE_SEQUENCE_DIAGNOSTICS_FIELDS = (
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
)
_SITE_SEQUENCE_ROW_DIAGNOSTIC_FIELDS = (
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
)
_INTENSITY_TRANSFORM_DIAGNOSTICS_FIELDS = (
    "policy",
    "pseudocount",
    "affected_matrices",
    "input_phospho_hash",
    "output_phospho_hash",
    "input_total_hash",
    "output_total_hash",
)
_NORMALISATION_DIAGNOSTICS_FIELDS = (
    "policy",
    "affected_columns",
    "input_phospho_hash",
    "output_phospho_hash",
    "note",
)
_COMPARISONS_DIAGNOSTICS_FIELDS = (
    "policy",
    "sample_group_column",
    "resolved_comparison_pairs",
    "group_labels",
    "output_comparison_hash",
)
_SITE_MATRIX_DIAGNOSTICS_FIELDS = (
    "dropped_missing_sequence_row_ids",
    "dropped_incomplete_row_ids",
    "dropped_row_ids",
    "duplicate_site_policy",
    "missing_data_policy",
    "required_observed_count",
    "final_constructed_site_ids",
    "duplicate_aggregation",
    "duplicate_site_decisions",
)


def _resolve_total_protein_parameters(plan: PreprocessingPlan) -> dict[str, object]:
    identity = plan.total_protein_correction_identity_policy
    return {
        "total_protein_correction_policy": plan.total_protein_correction_policy.value,
        "identity_mode": str(identity.mode),
        "identity_matching_policy": str(identity.matching_policy),
        "phosphosite_key": identity.phosphosite_key,
        "total_protein_key": identity.total_protein_key,
        "mapping_phosphosite_key": identity.mapping_phosphosite_key,
        "mapping_total_protein_key": identity.mapping_total_protein_key,
        "mapping_table_fingerprint": identity.mapping_table_fingerprint,
        "mapping_table_row_count": (
            None if identity.mapping_table is None else len(identity.mapping_table)
        ),
        "duplicate_policy": str(identity.duplicate_policy),
        "unmatched_policy": str(identity.unmatched_policy),
    }


def _resolve_intensity_transform_parameters(
    plan: PreprocessingPlan,
) -> dict[str, object]:
    return {"pseudocount": float(plan.intensity_transform_pseudocount)}


def _resolve_site_sequence_parameters(plan: PreprocessingPlan) -> dict[str, object]:
    return {
        "enabled": bool(plan.site_sequence_resolution_enabled),
        "fasta_path": plan.site_sequence_resolution_fasta_path,
        "mode": plan.site_sequence_resolution_mode.value,
        "conflict_policy": plan.site_sequence_resolution_conflict_policy.value,
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
        "site_matrix_policy": plan.site_matrix_policy.value,
        "site_matrix_duplicate_site_policy": plan.site_matrix_duplicate_site_policy.value,
        "site_matrix_missing_data_policy": plan.site_matrix_missing_data_policy.value,
        "site_matrix_minimum_observed_values": plan.site_matrix_minimum_observed_values,
    }


def _resolve_comparisons_parameters(plan: PreprocessingPlan) -> dict[str, object]:
    return {
        "comparison_building_policy": plan.comparison_building_policy.value,
        "comparison_sample_group_column": plan.comparison_sample_group_column,
        "comparison_pairs": plan.comparison_pairs,
    }


def _resolve_normalisation_parameters(_plan: PreprocessingPlan) -> dict[str, object]:
    return {}


def _resolve_missing_data_operation(plan: PreprocessingPlan) -> str:
    return plan.missing_data_policy.value


def _resolve_intensity_transform_operation(plan: PreprocessingPlan) -> str:
    return plan.intensity_transform_policy.value


def _resolve_total_protein_operation(plan: PreprocessingPlan) -> str:
    return plan.total_protein_correction_policy.value


def _resolve_site_matrix_operation(plan: PreprocessingPlan) -> str:
    return plan.site_matrix_policy.value


def _resolve_normalisation_operation(plan: PreprocessingPlan) -> str:
    return plan.normalisation_policy.value


def _resolve_comparisons_operation(plan: PreprocessingPlan) -> str:
    return plan.comparison_building_policy.value


def _resolve_site_sequence_operation(plan: PreprocessingPlan) -> str:
    return plan.site_sequence_resolution_mode.value


def _include_when_site_sequence_enabled(plan: PreprocessingPlan) -> bool:
    return bool(plan.site_sequence_resolution_enabled)


def _normalize_stage_table_keys(
    *,
    stage_key: str,
    table_keys: tuple[PreprocessingStateTableKey, ...],
    role: str,
) -> tuple[PreprocessingStateTableKey, ...]:
    normalized: list[PreprocessingStateTableKey] = []
    for index, table_key in enumerate(table_keys):
        if isinstance(table_key, PreprocessingStateTableKey):
            normalized.append(table_key)
            continue
        if not isinstance(table_key, str):
            raise DatasetBuildError(
                "dataset preprocessing stage metadata contains non-string table key: "
                f"stage={stage_key or '<empty>'!r}, field={role!r}[{index}], "
                f"got {table_key!r} ({type(table_key).__name__})"
            )
        try:
            normalized.append(PreprocessingStateTableKey(table_key))
        except ValueError as exc:
            supported = ", ".join(key.value for key in PREPROCESSING_STATE_TABLE_KEYS)
            raise DatasetBuildError(
                "dataset preprocessing stage metadata contains unknown table key: "
                f"stage={stage_key or '<empty>'!r}, field={role!r}[{index}], "
                f"table={table_key!r}, supported tables: {supported}"
            ) from exc
    return tuple(normalized)


PREPROCESSING_STAGE_REGISTRY: tuple[PreprocessingStageMetadata, ...] = (
    PreprocessingStageMetadata(
        stage_key=DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION,
        display_label=DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION,
        provenance_stage=DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION,
        operation_name=_resolve_site_sequence_operation,
        serialize_parameters=_resolve_site_sequence_parameters,
        consumed_input_tables=(PreprocessingStateTableKey.DATASET_SITE_METADATA,),
        produced_output_tables=(PreprocessingStateTableKey.DATASET_SITE_METADATA,),
        stage_factory=SiteSequenceResolutionStage,
        backend="phospy.sequences",
        include_when=_include_when_site_sequence_enabled,
        diagnostics_metadata={
            "diagnostics_schema_version": 1,
            "known_diagnostics_fields": _SITE_SEQUENCE_DIAGNOSTICS_FIELDS,
            "known_row_diagnostic_fields": _SITE_SEQUENCE_ROW_DIAGNOSTIC_FIELDS,
        },
    ),
    PreprocessingStageMetadata(
        stage_key=DATASET_PREPROCESSING_STAGE_MISSING_DATA,
        display_label=DATASET_PREPROCESSING_STAGE_MISSING_DATA,
        provenance_stage=DATASET_PREPROCESSING_STAGE_MISSING_DATA,
        operation_name=_resolve_missing_data_operation,
        serialize_parameters=_resolve_missing_data_parameters,
        consumed_input_tables=(
            PreprocessingStateTableKey.DATASET_PHOSPHO,
            PreprocessingStateTableKey.DATASET_SITE_METADATA,
        ),
        produced_output_tables=(
            PreprocessingStateTableKey.DATASET_PHOSPHO,
            PreprocessingStateTableKey.DATASET_SITE_METADATA,
            PreprocessingStateTableKey.REPORT_ROW_AUDIT,
        ),
        stage_factory=MissingDataStage,
        backend="pandas",
        diagnostics_metadata={
            "diagnostics_schema_version": MISSING_DATA_DIAGNOSTICS_SCHEMA_VERSION_V1,
            "known_diagnostics_fields": tuple(
                sorted(V1_KNOWN_MISSING_DATA_DIAGNOSTICS_FIELDS)
            ),
        },
    ),
    PreprocessingStageMetadata(
        stage_key=DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM,
        display_label=DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM,
        provenance_stage=DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM,
        operation_name=_resolve_intensity_transform_operation,
        serialize_parameters=_resolve_intensity_transform_parameters,
        consumed_input_tables=(
            PreprocessingStateTableKey.DATASET_PHOSPHO,
            PreprocessingStateTableKey.DATASET_TOTAL,
        ),
        produced_output_tables=(
            PreprocessingStateTableKey.DATASET_PHOSPHO,
            PreprocessingStateTableKey.DATASET_TOTAL,
        ),
        stage_factory=IntensityTransformStage,
        backend="numpy",
        diagnostics_metadata={
            "known_diagnostics_fields": _INTENSITY_TRANSFORM_DIAGNOSTICS_FIELDS
        },
    ),
    PreprocessingStageMetadata(
        stage_key=DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION,
        display_label=DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION,
        provenance_stage=DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION,
        operation_name=_resolve_total_protein_operation,
        serialize_parameters=_resolve_total_protein_parameters,
        consumed_input_tables=(
            PreprocessingStateTableKey.DATASET_PHOSPHO,
            PreprocessingStateTableKey.DATASET_TOTAL,
            PreprocessingStateTableKey.DATASET_SITE_METADATA,
        ),
        produced_output_tables=(PreprocessingStateTableKey.DATASET_PHOSPHO,),
        stage_factory=TotalProteinCorrectionStage,
        backend="pandas",
        diagnostics_metadata={
            "diagnostics_schema_version": (
                TOTAL_PROTEIN_CORRECTION_DIAGNOSTICS_SCHEMA_VERSION_V1
            ),
            "known_diagnostics_fields": tuple(
                sorted(V1_KNOWN_TOTAL_PROTEIN_DIAGNOSTICS_FIELDS)
            ),
        },
    ),
    PreprocessingStageMetadata(
        stage_key=DATASET_PREPROCESSING_STAGE_SITE_MATRIX,
        display_label=DATASET_PREPROCESSING_STAGE_SITE_MATRIX,
        provenance_stage=DATASET_PREPROCESSING_STAGE_SITE_MATRIX,
        operation_name=_resolve_site_matrix_operation,
        serialize_parameters=_resolve_site_matrix_parameters,
        consumed_input_tables=(
            PreprocessingStateTableKey.DATASET_PHOSPHO,
            PreprocessingStateTableKey.DATASET_SITE_METADATA,
        ),
        produced_output_tables=(
            PreprocessingStateTableKey.DATASET_PHOSPHO,
            PreprocessingStateTableKey.DATASET_SITE_METADATA,
            PreprocessingStateTableKey.REPORT_DUPLICATE_SITE_RESOLUTION,
            PreprocessingStateTableKey.REPORT_METADATA_CONFLICTS,
            PreprocessingStateTableKey.REPORT_ROW_AUDIT,
        ),
        stage_factory=SiteMatrixStage,
        backend="pandas",
        diagnostics_metadata={
            "known_diagnostics_fields": _SITE_MATRIX_DIAGNOSTICS_FIELDS
        },
    ),
    PreprocessingStageMetadata(
        stage_key=DATASET_PREPROCESSING_STAGE_NORMALISATION,
        display_label=DATASET_PREPROCESSING_STAGE_NORMALISATION,
        provenance_stage=DATASET_PREPROCESSING_STAGE_NORMALISATION,
        operation_name=_resolve_normalisation_operation,
        serialize_parameters=_resolve_normalisation_parameters,
        consumed_input_tables=(PreprocessingStateTableKey.DATASET_PHOSPHO,),
        produced_output_tables=(PreprocessingStateTableKey.DATASET_PHOSPHO,),
        stage_factory=NormalisationStage,
        backend="numpy",
        diagnostics_metadata={
            "known_diagnostics_fields": _NORMALISATION_DIAGNOSTICS_FIELDS
        },
    ),
    PreprocessingStageMetadata(
        stage_key=DATASET_PREPROCESSING_STAGE_COMPARISONS,
        display_label=DATASET_PREPROCESSING_STAGE_COMPARISONS,
        provenance_stage=DATASET_PREPROCESSING_STAGE_COMPARISONS,
        operation_name=_resolve_comparisons_operation,
        serialize_parameters=_resolve_comparisons_parameters,
        consumed_input_tables=(
            PreprocessingStateTableKey.DATASET_PHOSPHO,
            PreprocessingStateTableKey.DATASET_SAMPLE_METADATA,
        ),
        produced_output_tables=(
            PreprocessingStateTableKey.DATASET_COMPARISONS,
            PreprocessingStateTableKey.REPORT_COMPARISON_GROUP_STATS,
            PreprocessingStateTableKey.REPORT_COMPARISON_PAIR_STATS,
        ),
        stage_factory=ComparisonsStage,
        backend="pandas",
        diagnostics_metadata={
            "known_diagnostics_fields": _COMPARISONS_DIAGNOSTICS_FIELDS
        },
    ),
)


def _build_stage_metadata_by_key(
    registry: Sequence[PreprocessingStageMetadata],
    *,
    context: str,
) -> dict[str, PreprocessingStageMetadata]:
    by_key: dict[str, PreprocessingStageMetadata] = {}
    for metadata in registry:
        stage_key = metadata.stage_key.strip()
        if not stage_key:
            raise DatasetBuildError(
                f"{context} includes a stage with an empty stage key"
            )
        if not metadata.display_label.strip():
            raise DatasetBuildError(
                f"{context} contains stage {stage_key!r} with empty display label"
            )
        if not callable(metadata.operation_name):
            raise DatasetBuildError(
                f"{context} contains stage {stage_key!r} without operation resolver"
            )
        if not callable(metadata.serialize_parameters):
            raise DatasetBuildError(
                f"{context} contains stage {stage_key!r} without parameter serializer"
            )
        object.__setattr__(
            metadata,
            "consumed_input_tables",
            _normalize_stage_table_keys(
                stage_key=stage_key,
                table_keys=metadata.consumed_input_tables,
                role="consumed_input_tables",
            ),
        )
        object.__setattr__(
            metadata,
            "produced_output_tables",
            _normalize_stage_table_keys(
                stage_key=stage_key,
                table_keys=metadata.produced_output_tables,
                role="produced_output_tables",
            ),
        )
        if not isinstance(metadata.diagnostics_metadata, Mapping):
            raise DatasetBuildError(
                f"{context} contains stage {stage_key!r} with invalid diagnostics metadata"
            )
        if stage_key in by_key:
            raise DatasetBuildError(
                f"{context} contains duplicate stage key {stage_key!r}"
            )
        by_key[stage_key] = metadata
    return by_key


_STAGE_METADATA_BY_KEY = _build_stage_metadata_by_key(
    PREPROCESSING_STAGE_REGISTRY,
    context="dataset preprocessing stage registry",
)


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


def resolve_registered_preprocessing_stages(
    overrides: tuple[PreprocessingStageMetadata, ...] | None = None,
) -> tuple[PreprocessingStageMetadata, ...]:
    """Return ordered registry entries with optional metadata overrides applied."""

    if overrides is None:
        return PREPROCESSING_STAGE_REGISTRY

    override_by_key = _build_stage_metadata_by_key(
        overrides,
        context="dataset preprocessing stage metadata overrides",
    )
    resolved: list[PreprocessingStageMetadata] = []
    for metadata in PREPROCESSING_STAGE_REGISTRY:
        replacement = override_by_key.pop(metadata.stage_key, None)
        resolved.append(metadata if replacement is None else replacement)
    for metadata in override_by_key.values():
        resolved.append(metadata)
    return tuple(resolved)


def merge_preprocessing_stage_metadata(
    overrides: tuple[PreprocessingStageMetadata, ...] | None,
) -> dict[str, PreprocessingStageMetadata]:
    """Return default metadata map merged with optional custom overrides."""

    merged = dict(_STAGE_METADATA_BY_KEY)
    if overrides is None:
        return merged
    override_by_key = _build_stage_metadata_by_key(
        overrides,
        context="dataset preprocessing stage metadata overrides",
    )
    merged.update(override_by_key)
    return merged


def build_registered_preprocessing_stage_instances(
    metadata_registry: Sequence[PreprocessingStageMetadata],
) -> tuple[PreprocessingStage, ...]:
    """Construct stage instances from registry-owned stage factories."""

    instances: list[PreprocessingStage] = []
    for metadata in metadata_registry:
        factory = metadata.stage_factory
        if factory is None:
            continue
        instances.append(factory())
    return tuple(instances)


__all__ = [
    "PreprocessingStageMetadata",
    "get_preprocessing_stage_metadata",
    "build_registered_preprocessing_stage_instances",
    "list_registered_preprocessing_stages",
    "merge_preprocessing_stage_metadata",
    "resolve_registered_preprocessing_stages",
    "resolve_builder_provenance_stage_order",
]
