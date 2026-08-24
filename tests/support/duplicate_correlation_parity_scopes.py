from __future__ import annotations

from enum import StrEnum
from types import MappingProxyType


class DuplicateCorrelationParityScope(StrEnum):
    """Declared parity scope for each duplicate-correlation limma fixture."""

    FULL_PIPELINE = "full_pipeline"
    ESTIMATOR_AND_GLS_ONLY = "estimator_and_gls_only"


FULL_PIPELINE_FIXTURE_IDS = (
    "fixture_a_complete_pairs",
    "fixture_b_three_observation_blocks",
    "fixture_c_incomplete_unequal_blocks",
)

ESTIMATOR_AND_GLS_ONLY_FIXTURE_IDS = ("fixture_d_feature_level_failures",)

DUPLICATE_CORRELATION_PARITY_SCOPE_FIXTURE_IDS = MappingProxyType(
    {
        DuplicateCorrelationParityScope.FULL_PIPELINE: FULL_PIPELINE_FIXTURE_IDS,
        DuplicateCorrelationParityScope.ESTIMATOR_AND_GLS_ONLY: (
            ESTIMATOR_AND_GLS_ONLY_FIXTURE_IDS
        ),
    }
)

DUPLICATE_CORRELATION_FIXTURE_PARITY_SCOPES = MappingProxyType(
    {
        fixture_id: scope
        for scope, fixture_ids in (
            DUPLICATE_CORRELATION_PARITY_SCOPE_FIXTURE_IDS.items()
        )
        for fixture_id in fixture_ids
    }
)
