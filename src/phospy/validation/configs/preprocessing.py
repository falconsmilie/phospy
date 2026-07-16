"""Compatibility route for preprocessing public-config validation helpers."""

from __future__ import annotations

from phospy.contracts.configs.preprocessing._validation import (
    UNSUPPORTED_RUV_III_STYLE_METHOD_MESSAGE,
    reject_unsupported_ruv_iii_style_method,
    validate_batch_correction_config,
    validate_comparison_building_config,
    validate_group_coverage_filter_config,
    validate_intensity_transform_config,
    validate_internal_batch_correction_request,
    validate_localisation_config,
    validate_missing_data_config,
    validate_normalisation_config,
    validate_preprocessing_section_type,
    validate_protein_aware_preparation_config,
    validate_protein_aware_sample_alignment_config,
    validate_site_matrix_config,
    validate_site_sequence_resolution_config,
    validate_total_protein_correction_config,
    validate_total_protein_correction_identity_config,
)

__all__ = [
    "UNSUPPORTED_RUV_III_STYLE_METHOD_MESSAGE",
    "reject_unsupported_ruv_iii_style_method",
    "validate_batch_correction_config",
    "validate_comparison_building_config",
    "validate_group_coverage_filter_config",
    "validate_intensity_transform_config",
    "validate_internal_batch_correction_request",
    "validate_localisation_config",
    "validate_missing_data_config",
    "validate_normalisation_config",
    "validate_preprocessing_section_type",
    "validate_protein_aware_preparation_config",
    "validate_protein_aware_sample_alignment_config",
    "validate_site_matrix_config",
    "validate_site_sequence_resolution_config",
    "validate_total_protein_correction_config",
    "validate_total_protein_correction_identity_config",
]
