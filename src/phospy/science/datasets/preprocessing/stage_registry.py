"""Central composition for preprocessing stage contracts."""

from __future__ import annotations

from collections.abc import Mapping, Sequence

from phospy.errors.build import DatasetBuildError
from phospy.science.datasets.preprocessing.models import (
    PreprocessingPlan,
    PreprocessingStage,
)
from phospy.science.datasets.preprocessing.stage_contract import (
    PreprocessingStageContract,
    PreprocessingStageFactoryContext,
)
from phospy.science.datasets.preprocessing.stages import (
    BATCH_CORRECTION_STAGE_CONTRACT,
    COMPARISONS_STAGE_CONTRACT,
    GROUP_COVERAGE_FILTER_STAGE_CONTRACT,
    INTENSITY_TRANSFORM_STAGE_CONTRACT,
    LOCALISATION_CONFIDENCE_STAGE_CONTRACT,
    MISSING_DATA_STAGE_CONTRACT,
    NORMALISATION_STAGE_CONTRACT,
    SITE_MATRIX_STAGE_CONTRACT,
    SITE_SEQUENCE_RESOLUTION_STAGE_CONTRACT,
    TOTAL_PROTEIN_CORRECTION_STAGE_CONTRACT,
)

# Backward-compatible alias for existing imports/tests.
PreprocessingStageMetadata = PreprocessingStageContract


PREPROCESSING_STAGE_REGISTRY: tuple[PreprocessingStageContract, ...] = (
    SITE_SEQUENCE_RESOLUTION_STAGE_CONTRACT,
    LOCALISATION_CONFIDENCE_STAGE_CONTRACT,
    GROUP_COVERAGE_FILTER_STAGE_CONTRACT,
    MISSING_DATA_STAGE_CONTRACT,
    INTENSITY_TRANSFORM_STAGE_CONTRACT,
    BATCH_CORRECTION_STAGE_CONTRACT,
    TOTAL_PROTEIN_CORRECTION_STAGE_CONTRACT,
    SITE_MATRIX_STAGE_CONTRACT,
    NORMALISATION_STAGE_CONTRACT,
    COMPARISONS_STAGE_CONTRACT,
)


def _build_stage_metadata_by_key(
    registry: Sequence[PreprocessingStageContract],
    *,
    context: str,
) -> dict[str, PreprocessingStageContract]:
    by_key: dict[str, PreprocessingStageContract] = {}
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
        if not callable(metadata.validate_plan):
            raise DatasetBuildError(
                f"{context} contains stage {stage_key!r} without plan validator"
            )
        if not callable(metadata.resolve_random_seed):
            raise DatasetBuildError(
                f"{context} contains stage {stage_key!r} without random-seed resolver"
            )
        if not callable(metadata.resolve_determinism_kind):
            raise DatasetBuildError(
                f"{context} contains stage {stage_key!r} without determinism resolver"
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


def get_preprocessing_stage_metadata(stage_key: str) -> PreprocessingStageContract:
    """Return stage contract metadata for ``stage_key``."""

    metadata = _STAGE_METADATA_BY_KEY.get(stage_key)
    if metadata is not None:
        return metadata
    supported = ", ".join(sorted(_STAGE_METADATA_BY_KEY))
    raise DatasetBuildError(
        "dataset preprocessing stage metadata is not registered for stage "
        f"{stage_key!r}; supported stages: {supported}"
    )


def list_registered_preprocessing_stages() -> tuple[PreprocessingStageContract, ...]:
    """Return the shared preprocessing stage contract registry."""

    return PREPROCESSING_STAGE_REGISTRY


def resolve_registered_preprocessing_stages(
    overrides: tuple[PreprocessingStageContract, ...] | None = None,
) -> tuple[PreprocessingStageContract, ...]:
    """Return ordered registry entries with optional metadata overrides applied."""

    if overrides is None:
        return PREPROCESSING_STAGE_REGISTRY

    override_by_key = _build_stage_metadata_by_key(
        overrides,
        context="dataset preprocessing stage metadata overrides",
    )
    resolved: list[PreprocessingStageContract] = []
    for metadata in PREPROCESSING_STAGE_REGISTRY:
        replacement = override_by_key.pop(metadata.stage_key, None)
        resolved.append(metadata if replacement is None else replacement)
    for metadata in override_by_key.values():
        resolved.append(metadata)
    return tuple(resolved)


def merge_preprocessing_stage_metadata(
    overrides: tuple[PreprocessingStageContract, ...] | None,
) -> dict[str, PreprocessingStageContract]:
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
    metadata_registry: Sequence[PreprocessingStageContract],
    *,
    context: PreprocessingStageFactoryContext | None = None,
) -> tuple[PreprocessingStage, ...]:
    """Construct stage instances from registry-owned stage factories."""

    factory_context = context or PreprocessingStageFactoryContext()
    instances: list[PreprocessingStage] = []
    for metadata in metadata_registry:
        factory = metadata.stage_factory
        if not callable(factory):
            continue
        instances.append(factory(factory_context))
    return tuple(instances)


def resolve_builder_provenance_stage_order(
    plan: PreprocessingPlan,
) -> tuple[PreprocessingStageContract, ...]:
    """Return stage metadata sequence used by builder provenance tables."""

    return tuple(
        metadata
        for metadata in PREPROCESSING_STAGE_REGISTRY
        if metadata.include_in_builder_provenance and metadata.include_when(plan)
    )


__all__ = [
    "PreprocessingStageMetadata",
    "build_registered_preprocessing_stage_instances",
    "get_preprocessing_stage_metadata",
    "list_registered_preprocessing_stages",
    "merge_preprocessing_stage_metadata",
    "resolve_builder_provenance_stage_order",
    "resolve_registered_preprocessing_stages",
]
