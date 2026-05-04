from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from phospy.activities.threshold_membership import threshold_membership_mask_array
from phospy.api.configs import (
    DatasetMissingDataConfig,
    DatasetPreprocessingConfig,
    DatasetTotalProteinCorrectionConfig,
)
from phospy.datasets.preprocessing.models import PreprocessingPlan
from phospy.errors.input import PhosPyInputError
from phospy.policy_models import (
    DownstreamScoreSource,
    MissingDataPolicy,
    ThresholdMode,
    TotalProteinCorrectionPolicy,
)
from phospy.prediction.candidates import build_candidate_substrate_list
from phospy.prediction.scoring import select_downstream_score_matrix


def test_preprocessing_plan_converts_public_strings_to_internal_policy_enums() -> None:
    plan = PreprocessingPlan.from_config(
        DatasetPreprocessingConfig(
            missing_data=DatasetMissingDataConfig(policy="forbid"),
            total_protein_correction=DatasetTotalProteinCorrectionConfig(policy="none"),
        )
    )

    assert isinstance(plan.missing_data_policy, MissingDataPolicy)
    assert plan.missing_data_policy is MissingDataPolicy.FORBID
    assert isinstance(
        plan.total_protein_correction_policy, TotalProteinCorrectionPolicy
    )
    assert plan.total_protein_correction_policy is TotalProteinCorrectionPolicy.NONE


def test_preprocessing_plan_rejects_unknown_policy_strings_without_defaulting() -> None:
    with pytest.raises(PhosPyInputError, match="must be one of:"):
        PreprocessingPlan(
            missing_data_policy="unknown_policy",  # type: ignore[arg-type]
            total_protein_correction_policy=TotalProteinCorrectionPolicy.NONE,
        )

    with pytest.raises(PhosPyInputError, match="must be one of:"):
        PreprocessingPlan(
            missing_data_policy=MissingDataPolicy.FORBID,
            total_protein_correction_policy="unknown_policy",  # type: ignore[arg-type]
        )


def test_downstream_score_selection_returns_enum_source() -> None:
    profile = pd.DataFrame({"K1": [0.1]}, index=["S1"])
    combined = pd.DataFrame({"K1": [0.2]}, index=["S1"])

    selected, source = select_downstream_score_matrix(
        profile_scores=profile,
        rank_weighted_fusion_scores=combined,
    )

    assert selected is combined
    assert source is DownstreamScoreSource.RANK_WEIGHTED_FUSION_SCORES


def test_threshold_mode_drives_membership_logic_with_enum_or_string_input() -> None:
    scores = np.array([0.5, 0.4, np.nan])
    strict = threshold_membership_mask_array(
        scores,
        threshold=0.5,
        threshold_mode=ThresholdMode.GREATER_THAN,
    )
    inclusive = threshold_membership_mask_array(
        scores,
        threshold=0.5,
        threshold_mode="score >= threshold",
    )

    assert strict.tolist() == [False, False, False]
    assert inclusive.tolist() == [True, False, False]


def test_candidate_selection_rejects_unknown_threshold_mode() -> None:
    scores = pd.DataFrame({"K1": [0.7]}, index=["S1"])
    with pytest.raises(PhosPyInputError, match="must be one of:"):
        build_candidate_substrate_list(
            scores=scores,
            top=1,
            score_threshold=0.0,
            inclusion=1,
            threshold_mode="invalid",  # type: ignore[arg-type]
        )
