from __future__ import annotations

import importlib

import pytest

from phospy.differential.policies import TechnicalReplicatePolicy
from phospy.errors.input import PhosPyInputError
from phospy.kinase.policies import ThresholdMode
from phospy.kinase.result_boundary import DownstreamScoreSource
from phospy.kinase.scoring.policies import ProfileSelfInclusionPolicy
from phospy.policies import PolicyEnum
from phospy.preparation_records.policies import (
    MissingDataPolicy,
    TotalProteinCorrectionPolicy,
)
from phospy.preprocessing.contracts.policies import (
    SiteMatrixPolicy,
)


def test_preparation_record_policies_have_direct_owner_identity() -> None:
    assert MissingDataPolicy.__module__ == "phospy.preparation_records.policies"
    assert SiteMatrixPolicy.__module__ == "phospy.preprocessing.contracts.policies"
    assert (
        TotalProteinCorrectionPolicy.__module__ == "phospy.preparation_records.policies"
    )


def test_differential_policy_is_owned_by_differential_module() -> None:
    assert TechnicalReplicatePolicy.__module__ == "phospy.differential.policies"


def test_scoring_policies_have_their_direct_owner_identity() -> None:
    assert ThresholdMode.__module__ == "phospy.kinase.policies"
    assert DownstreamScoreSource.__module__ == "phospy.kinase.result_boundary"
    assert ProfileSelfInclusionPolicy.__module__ == "phospy.kinase.scoring.policies"


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
    assert (
        ProfileSelfInclusionPolicy.parse("allow", field_name="x")
        is ProfileSelfInclusionPolicy.ALLOW
    )


@pytest.mark.parametrize(
    ("enum_type", "value"),
    [
        (MissingDataPolicy, "not_valid"),
        (TechnicalReplicatePolicy, "not_valid"),
        (ThresholdMode, "not_valid"),
        (ProfileSelfInclusionPolicy, "not_valid"),
    ],
)
def test_policy_parse_invalid_values_raise_input_error(
    enum_type: type[PolicyEnum],
    value: str,
) -> None:
    with pytest.raises(PhosPyInputError):
        enum_type.parse(value, field_name="x")
