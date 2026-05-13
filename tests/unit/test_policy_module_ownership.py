from __future__ import annotations

import importlib

import pytest

from phospy.errors.input import PhosPyInputError
from phospy.policies import PolicyEnum
from phospy.science.datasets.preprocessing.policy_models import (
    MissingDataPolicy,
    SiteMatrixPolicy,
    TotalProteinCorrectionPolicy,
)
from phospy.science.differential.policy_models import TechnicalReplicatePolicy
from phospy.science.scoring.policy_models import DownstreamScoreSource, ThresholdMode


def test_preprocessing_policies_are_owned_by_preprocessing_module() -> None:
    assert (
        MissingDataPolicy.__module__
        == "phospy.science.datasets.preprocessing.policy_models"
    )
    assert (
        SiteMatrixPolicy.__module__
        == "phospy.science.datasets.preprocessing.policy_models"
    )
    assert (
        TotalProteinCorrectionPolicy.__module__
        == "phospy.science.datasets.preprocessing.policy_models"
    )


def test_differential_policy_is_owned_by_differential_module() -> None:
    assert (
        TechnicalReplicatePolicy.__module__
        == "phospy.science.differential.policy_models"
    )


def test_scoring_policies_are_owned_by_scoring_module() -> None:
    assert ThresholdMode.__module__ == "phospy.science.scoring.policy_models"
    assert DownstreamScoreSource.__module__ == "phospy.science.scoring.policy_models"


def test_policy_base_is_owned_by_policies_infrastructure_module() -> None:
    assert PolicyEnum.__module__ == "phospy.policies.policy_base"


def test_root_policy_models_module_is_removed() -> None:
    with pytest.raises(ModuleNotFoundError):
        importlib.import_module("phospy.policy_models")


def test_policy_parse_behaviour_is_preserved_for_owner_modules() -> None:
    assert MissingDataPolicy.parse("forbid", field_name="x") is MissingDataPolicy.FORBID
    assert (
        TechnicalReplicatePolicy.parse("mean", field_name="x")
        is TechnicalReplicatePolicy.MEAN
    )
    assert (
        ThresholdMode.parse("score >= threshold", field_name="x")
        is ThresholdMode.GREATER_THAN_OR_EQUAL
    )


@pytest.mark.parametrize(
    ("enum_type", "value"),
    [
        (MissingDataPolicy, "not_valid"),
        (TechnicalReplicatePolicy, "not_valid"),
        (ThresholdMode, "not_valid"),
    ],
)
def test_policy_parse_invalid_values_raise_input_error(
    enum_type: type[PolicyEnum],
    value: str,
) -> None:
    with pytest.raises(PhosPyInputError):
        enum_type.parse(value, field_name="x")
