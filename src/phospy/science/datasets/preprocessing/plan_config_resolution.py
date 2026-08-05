"""Public preprocessing config resolution for internal plan construction."""

from __future__ import annotations

from dataclasses import dataclass
from typing import cast

from phospy.science.configs.dataset import DatasetPreprocessingConfig
from phospy.science.configs.preprocessing import (
    CorrectionMissingnessPolicy,
    DatasetBatchCorrectionConfig,
    DatasetSiteSequenceConflictPolicy,
    InternalBatchCorrectionRequest,
    SpsRuvBatchCorrectionConfig,
)
from phospy.science.datasets.preprocessing.policy_models import (
    ComparisonBuildingPolicy,
    IntensityTransformPolicy,
    LocalisationEligibilityMode,
    MissingDataPolicy,
    NormalisationPolicy,
    SiteMatrixDuplicateSitePolicy,
    SiteMatrixMissingDataPolicy,
    SiteMatrixPolicy,
    SiteSequenceConflictPolicy,
    SiteSequenceResolutionMode,
    TotalProteinCorrectionPolicy,
)


@dataclass(frozen=True, slots=True)
class PreprocessingConfigPolicyResolution:
    """Typed policy facts resolved from public preprocessing config."""

    site_sequence_resolution_enabled: bool
    intensity_transform_policy: IntensityTransformPolicy
    normalisation_policy: NormalisationPolicy
    site_matrix_policy: SiteMatrixPolicy
    site_matrix_duplicate_site_policy: SiteMatrixDuplicateSitePolicy
    site_matrix_missing_data_policy: SiteMatrixMissingDataPolicy
    comparison_building_policy: ComparisonBuildingPolicy
    site_sequence_resolution_mode: SiteSequenceResolutionMode
    localisation_mode: LocalisationEligibilityMode
    missing_data_policy: MissingDataPolicy
    batch_correction_method: str
    batch_correction_condition_column: str
    batch_correction_condition_columns: tuple[str, ...]
    batch_correction_replicate_column: str | None
    batch_correction_control_site_set: object | None
    batch_correction_missingness_policy: CorrectionMissingnessPolicy | None
    batch_correction_internal_request: InternalBatchCorrectionRequest | None
    batch_correction_preserve_condition_effects: bool
    total_correction_policy: TotalProteinCorrectionPolicy


class PreprocessingConfigPolicyResolver:
    """Resolve public config policy values before stage-order planning."""

    def run(
        self,
        config: DatasetPreprocessingConfig,
    ) -> PreprocessingConfigPolicyResolution:
        site_sequence_resolution_enabled = (
            config.site_sequence_resolution.fasta_path is not None
        )
        intensity_transform_policy = IntensityTransformPolicy.parse(
            config.intensity_transform.policy,
            field_name="preprocessing_config.intensity_transform.policy",
        )
        normalisation_policy = NormalisationPolicy.parse(
            config.normalisation.policy,
            field_name="preprocessing_config.normalisation.policy",
        )
        site_matrix_policy = SiteMatrixPolicy.parse(
            config.site_matrix.policy,
            field_name="preprocessing_config.site_matrix.policy",
        )
        site_matrix_duplicate_site_policy = SiteMatrixDuplicateSitePolicy.parse(
            config.site_matrix.duplicate_site_policy,
            field_name="preprocessing_config.site_matrix.duplicate_site_policy",
        )
        site_matrix_missing_data_policy = SiteMatrixMissingDataPolicy.parse(
            config.site_matrix.missing_data_policy,
            field_name="preprocessing_config.site_matrix.missing_data_policy",
        )
        comparison_building_policy = ComparisonBuildingPolicy.parse(
            config.comparisons.policy,
            field_name="preprocessing_config.comparisons.policy",
        )
        site_sequence_resolution_mode = SiteSequenceResolutionMode.parse(
            config.site_sequence_resolution.mode,
            field_name="preprocessing_config.site_sequence_resolution.mode",
        )
        localisation_mode = LocalisationEligibilityMode.parse(
            config.localisation.mode,
            field_name="preprocessing_config.localisation.mode",
        )
        missing_data_policy = MissingDataPolicy.parse(
            config.missing_data.policy,
            field_name="preprocessing_config.missing_data.policy",
        )
        return PreprocessingConfigPolicyResolution(
            site_sequence_resolution_enabled=site_sequence_resolution_enabled,
            intensity_transform_policy=intensity_transform_policy,
            normalisation_policy=normalisation_policy,
            site_matrix_policy=site_matrix_policy,
            site_matrix_duplicate_site_policy=site_matrix_duplicate_site_policy,
            site_matrix_missing_data_policy=site_matrix_missing_data_policy,
            comparison_building_policy=comparison_building_policy,
            site_sequence_resolution_mode=site_sequence_resolution_mode,
            localisation_mode=localisation_mode,
            missing_data_policy=missing_data_policy,
            batch_correction_method=str(config.batch_correction.method).strip(),
            batch_correction_condition_column=_resolve_batch_condition_column(config),
            batch_correction_condition_columns=_resolve_batch_condition_columns(config),
            batch_correction_replicate_column=_resolve_batch_replicate_column(config),
            batch_correction_control_site_set=_resolve_batch_control_site_set(config),
            batch_correction_missingness_policy=_resolve_batch_missingness_policy(
                config
            ),
            batch_correction_internal_request=_resolve_batch_internal_request(config),
            batch_correction_preserve_condition_effects=(
                _resolve_batch_preserve_condition_effects(config)
            ),
            total_correction_policy=TotalProteinCorrectionPolicy.parse(
                config.total_protein_correction.policy,
                field_name="preprocessing_config.total_protein_correction.policy",
            ),
        )


def resolve_site_sequence_resolution_conflict_policy(
    *,
    mode: SiteSequenceResolutionMode,
    conflict_policy: DatasetSiteSequenceConflictPolicy | None,
) -> SiteSequenceConflictPolicy:
    if conflict_policy is not None:
        return SiteSequenceConflictPolicy.parse(
            conflict_policy,
            field_name="preprocessing_config.site_sequence_resolution.conflict_policy",
        )
    if mode is SiteSequenceResolutionMode.REPLACE_EXISTING:
        return SiteSequenceConflictPolicy.REPLACE_EXISTING
    return SiteSequenceConflictPolicy.PRESERVE_EXISTING


def _sps_ruv_config(
    config: DatasetPreprocessingConfig,
) -> SpsRuvBatchCorrectionConfig | None:
    batch_correction = config.batch_correction
    if isinstance(batch_correction, SpsRuvBatchCorrectionConfig):
        return batch_correction
    return None


def _resolve_batch_condition_column(config: DatasetPreprocessingConfig) -> str:
    sps_config = _sps_ruv_config(config)
    if sps_config is not None:
        return sps_config.condition_columns[0]
    batch_correction = cast(DatasetBatchCorrectionConfig, config.batch_correction)
    return batch_correction.condition_column


def _resolve_batch_condition_columns(
    config: DatasetPreprocessingConfig,
) -> tuple[str, ...]:
    sps_config = _sps_ruv_config(config)
    if sps_config is not None:
        return sps_config.condition_columns
    batch_correction = cast(DatasetBatchCorrectionConfig, config.batch_correction)
    return (batch_correction.condition_column,)


def _resolve_batch_replicate_column(config: DatasetPreprocessingConfig) -> str | None:
    sps_config = _sps_ruv_config(config)
    if sps_config is not None:
        return sps_config.replicate_column
    return None


def _resolve_batch_internal_request(
    config: DatasetPreprocessingConfig,
) -> InternalBatchCorrectionRequest | None:
    sps_config = _sps_ruv_config(config)
    if sps_config is not None:
        return sps_config.to_internal_request()
    return None


def _resolve_batch_control_site_set(
    config: DatasetPreprocessingConfig,
) -> object | None:
    sps_config = _sps_ruv_config(config)
    if sps_config is not None:
        return sps_config.control_site_set
    return None


def _resolve_batch_missingness_policy(
    config: DatasetPreprocessingConfig,
) -> CorrectionMissingnessPolicy | None:
    sps_config = _sps_ruv_config(config)
    if sps_config is not None:
        return sps_config.missingness_policy
    return None


def _resolve_batch_preserve_condition_effects(
    config: DatasetPreprocessingConfig,
) -> bool:
    if _sps_ruv_config(config) is not None:
        return True
    batch_correction = cast(DatasetBatchCorrectionConfig, config.batch_correction)
    return cast(bool, batch_correction.preserve_condition_effects)


__all__ = [
    "PreprocessingConfigPolicyResolution",
    "PreprocessingConfigPolicyResolver",
    "resolve_site_sequence_resolution_conflict_policy",
]
