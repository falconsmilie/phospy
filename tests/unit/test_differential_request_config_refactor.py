from __future__ import annotations

from typing import get_type_hints

import phospy.api.requests as request_models
from phospy.api import (
    DifferentialAnalysisConfig,
    DifferentialAnalysisRequest,
    MultipleTestingConfig,
    TechnicalReplicatePolicy,
)
from phospy.api.configs import DifferentialAnalysisConfig as ConfigFromConfigs
from phospy.differential.policy_models import (
    TechnicalReplicatePolicy as PolicyFromDifferentialPolicyModels,
)


def test_differential_request_is_thin_shape_without_post_init() -> None:
    assert "__post_init__" not in DifferentialAnalysisRequest.__dict__
    hints = get_type_hints(DifferentialAnalysisRequest)
    assert set(hints) == {"dataset", "design", "contrasts", "config"}
    assert hints["config"] is DifferentialAnalysisConfig


def test_differential_request_construction_does_not_coerce_inputs() -> None:
    contrasts = ["not-a-contrast"]
    config = DifferentialAnalysisConfig(
        technical_replicate_policy="mean",  # type: ignore[arg-type]
        allow_design_subset=1,  # type: ignore[arg-type]
        minimum_condition_replicates=0,
    )
    request = DifferentialAnalysisRequest(
        dataset=object(),  # type: ignore[arg-type]
        design=object(),  # type: ignore[arg-type]
        contrasts=contrasts,  # type: ignore[arg-type]
        config=config,
    )
    assert request.contrasts is contrasts
    assert request.config is config


def test_api_requests_no_longer_owns_differential_policy_or_testing_config() -> None:
    assert not hasattr(request_models, "TechnicalReplicatePolicy")
    assert not hasattr(request_models, "MultipleTestingConfig")


def test_technical_replicate_policy_is_owned_by_differential_policy_models() -> None:
    assert TechnicalReplicatePolicy is PolicyFromDifferentialPolicyModels


def test_differential_config_is_exported_from_api_configs() -> None:
    assert ConfigFromConfigs is DifferentialAnalysisConfig


def test_differential_config_and_policy_are_exported_from_api() -> None:
    assert isinstance(DifferentialAnalysisConfig(), DifferentialAnalysisConfig)
    assert isinstance(MultipleTestingConfig(), MultipleTestingConfig)
    assert TechnicalReplicatePolicy.REJECT.value == "reject"
