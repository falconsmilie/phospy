"""Policy construction and validation for dataset evidence resolution."""

from __future__ import annotations

import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.science.evidence.dataset_resolution.models import (
    _INPUT_MEANING_BY_SCALE,
    _OUTPUT_MEANING_BY_SCALE,
    DATASET_MULTI_SITE_POLICY_EXCLUDE_FROM_SEQUENCE_SCORING,
    DATASET_MULTI_SITE_POLICY_REJECT,
    DATASET_MULTI_SITE_POLICY_SPLIT,
    DATASET_PEPTIDE_TO_SITE_AGGREGATION_POLICY_LINEAR_ALLOCATED_MEAN_V1,
    SUPPORTED_DATASET_MULTI_SITE_POLICIES,
    PeptideEvidenceInputQuantitativeMeaning,
    PeptideToSiteAggregationPolicy,
    PeptideToSiteAllocationDomain,
    _fractional_mapping_mask,
    _mapping_fractions,
    _normalize_input_quantitative_meaning,
    _normalize_intensity_scale_kind,
    _raise_fractional_allocation_for_non_linear_input,
)
from phospy.science.evidence.multi_site import (
    MULTI_SITE_POLICY_ERROR,
    MULTI_SITE_POLICY_EXCLUDE_FROM_SEQUENCE_SCORING,
    MULTI_SITE_POLICY_SPLIT_EQUAL_WEIGHT,
    MultiSiteHandlingConfig,
)
from phospy.science.transformations.models import IntensityScaleKind

_REMOVED_DATASET_MULTI_SITE_POLICY_KEEP_JOINT = "keep_joint"
_POLICY_TO_MULTI_SITE_HANDLING_POLICY: dict[str, str] = {
    DATASET_MULTI_SITE_POLICY_REJECT: MULTI_SITE_POLICY_ERROR,
    DATASET_MULTI_SITE_POLICY_EXCLUDE_FROM_SEQUENCE_SCORING: (
        MULTI_SITE_POLICY_EXCLUDE_FROM_SEQUENCE_SCORING
    ),
    DATASET_MULTI_SITE_POLICY_SPLIT: MULTI_SITE_POLICY_SPLIT_EQUAL_WEIGHT,
}


def build_peptide_to_site_aggregation_policy(
    *,
    input_intensity_scale: IntensityScaleKind | str,
    input_quantitative_meaning: PeptideEvidenceInputQuantitativeMeaning
    | str
    | None = None,
    mapping_rows: pd.DataFrame,
) -> PeptideToSiteAggregationPolicy:
    """Build the current typed sample-intensity peptide-to-site estimand contract."""

    input_scale = _normalize_intensity_scale_kind(
        input_intensity_scale,
        field_name="peptide_to_site_aggregation_policy.input_intensity_scale",
    )
    resolved_input_meaning = (
        _INPUT_MEANING_BY_SCALE[input_scale]
        if input_quantitative_meaning is None
        else _normalize_input_quantitative_meaning(input_quantitative_meaning)
    )
    mapping_fractions = _mapping_fractions(mapping_rows)
    fractional_mask = _fractional_mapping_mask(mapping_fractions)
    fractional_mapping_present = bool(fractional_mask.any())
    if fractional_mapping_present and input_scale is not IntensityScaleKind.LINEAR:
        _raise_fractional_allocation_for_non_linear_input(
            rows=mapping_rows,
            input_intensity_scale=input_scale,
            fractional_mask=fractional_mask,
        )
    allocation_domain = (
        PeptideToSiteAllocationDomain.LINEAR_ABUNDANCE
        if input_scale is IntensityScaleKind.LINEAR
        else PeptideToSiteAllocationDomain.DECLARED_SCALE_UNIT_MAPPING_PASSTHROUGH
    )
    return PeptideToSiteAggregationPolicy(
        policy_id=DATASET_PEPTIDE_TO_SITE_AGGREGATION_POLICY_LINEAR_ALLOCATED_MEAN_V1,
        input_scale=input_scale,
        input_quantitative_meaning=resolved_input_meaning,
        output_scale=input_scale,
        output_quantitative_meaning=_OUTPUT_MEANING_BY_SCALE[input_scale],
        allocation_domain=allocation_domain,
        fractional_mapping_present=fractional_mapping_present,
    )


def build_multi_site_handling_config_for_dataset_policy(
    *,
    multi_site_policy: str,
) -> MultiSiteHandlingConfig:
    """Translate dataset-builder multi-site policy to evidence config."""

    validate_dataset_multi_site_policy(
        multi_site_policy,
        field_name="dataset build request multi_site_policy",
    )
    resolved_policy = _POLICY_TO_MULTI_SITE_HANDLING_POLICY[multi_site_policy]
    return MultiSiteHandlingConfig(
        statistical_modeling_policy=resolved_policy,
        kinase_sequence_scoring_policy=resolved_policy,
    )


def validate_dataset_multi_site_policy(policy: object, *, field_name: str) -> None:
    if policy == _REMOVED_DATASET_MULTI_SITE_POLICY_KEEP_JOINT:
        raise PhosPyInputError(
            f"{field_name}='keep_joint' is no longer supported for "
            "AnalysisReadyDatasetBuilder peptide-evidence requests because "
            "unresolved joint evidence cannot satisfy the strict site-level "
            "identity contract. Use multi_site_policy='split' to allocate "
            "ambiguous evidence to strict site rows, 'reject' to fail on "
            "ambiguous peptide rows, or 'exclude_from_sequence_scoring' to "
            "remove ambiguous rows from the analysis-ready build."
        )
    if (
        not isinstance(policy, str)
        or policy not in SUPPORTED_DATASET_MULTI_SITE_POLICIES
    ):
        supported = ", ".join(
            repr(value) for value in SUPPORTED_DATASET_MULTI_SITE_POLICIES
        )
        raise PhosPyInputError(f"{field_name} must be one of: {supported}")


__all__ = [
    "build_multi_site_handling_config_for_dataset_policy",
    "build_peptide_to_site_aggregation_policy",
    "validate_dataset_multi_site_policy",
]
