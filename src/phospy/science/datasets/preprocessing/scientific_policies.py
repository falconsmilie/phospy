"""Scientific policy records for dataset preprocessing provenance."""

from __future__ import annotations

from dataclasses import dataclass

from phospy.provenance.scientific_policy_models import (
    ScientificPolicyId,
    ScientificPolicyRecord,
)


@dataclass(frozen=True, slots=True)
class PreprocessingStageOrderPolicy:
    """Executable metadata policy for preprocessing stage-order behavior."""

    configured_stage_order: tuple[str, ...]
    default_stage_order: tuple[str, ...]
    supported_stage_order: tuple[str, ...]

    @property
    def record(self) -> ScientificPolicyRecord:
        return build_preprocessing_stage_order_policy(
            configured_stage_order=self.configured_stage_order,
            default_stage_order=self.default_stage_order,
            supported_stage_order=self.supported_stage_order,
        )


def build_preprocessing_stage_order_policy(
    *,
    configured_stage_order: tuple[str, ...],
    default_stage_order: tuple[str, ...],
    supported_stage_order: tuple[str, ...],
) -> ScientificPolicyRecord:
    return ScientificPolicyRecord(
        id=ScientificPolicyId.PREPROCESSING_STAGE_ORDER,
        name="Preprocessing Stage Order Policy",
        version="1",
        description=(
            "Defines the explicit stage execution order used to transform dataset "
            "inputs into analysis-ready workflow matrices."
        ),
        parameters={
            "configured_stage_order": " -> ".join(configured_stage_order),
            "configured_stage_count": int(len(configured_stage_order)),
            "default_stage_order": " -> ".join(default_stage_order),
            "default_stage_count": int(len(default_stage_order)),
            "supported_stage_order": " -> ".join(supported_stage_order),
            "supported_stage_count": int(len(supported_stage_order)),
        },
        assumptions=(
            "Stage order can change output row retention, transformed values, and "
            "derived comparison tables.",
            "Configured order must be interpreted as part of the scientific method "
            "for reproducibility.",
        ),
        output_scale="Ordered preprocessing execution plan for dataset construction.",
        quantitative_meaning="preprocessing_execution_order",
    )


DUPLICATE_SITE_RESOLUTION_AGGREGATE_MEAN_POLICY = ScientificPolicyRecord(
    id=ScientificPolicyId.DUPLICATE_SITE_RESOLUTION,
    name="duplicate_site_resolution_aggregate_mean_v1",
    version="1",
    description=(
        "Resolves duplicate site_key rows by aggregating duplicate rows "
        "with a column-wise arithmetic mean."
    ),
    parameters={
        "duplicate_site_policy": "aggregate_mean",
        "aggregation_function": "column_mean",
    },
    assumptions=(
        "All duplicate rows contribute numerically to one retained site row.",
        "Aggregation can blend peptide-context-specific measurements.",
    ),
    output_scale="Deduplicated site-matrix rows for downstream scoring.",
    quantitative_meaning="duplicate_site_resolved_matrix",
)


def build_duplicate_site_resolution_policy(
    *,
    duplicate_site_policy: str,
) -> ScientificPolicyRecord:
    if duplicate_site_policy == "aggregate_mean":
        return DUPLICATE_SITE_RESOLUTION_AGGREGATE_MEAN_POLICY
    return ScientificPolicyRecord(
        id=ScientificPolicyId.DUPLICATE_SITE_RESOLUTION,
        name=f"duplicate_site_resolution_{duplicate_site_policy}_v1",
        version="1",
        description=(
            "Resolves duplicate site_key rows according to the configured "
            "site-matrix duplicate-site policy."
        ),
        parameters={"duplicate_site_policy": str(duplicate_site_policy)},
        assumptions=(
            "Duplicate-site resolution policy changes retained rows and/or values.",
            "Resolved rows become the authoritative site matrix for downstream use.",
        ),
        output_scale="Deduplicated site-matrix rows for downstream scoring.",
        quantitative_meaning="duplicate_site_resolved_matrix",
    )


__all__ = [
    "DUPLICATE_SITE_RESOLUTION_AGGREGATE_MEAN_POLICY",
    "PreprocessingStageOrderPolicy",
    "build_duplicate_site_resolution_policy",
    "build_preprocessing_stage_order_policy",
]
