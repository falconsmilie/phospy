from __future__ import annotations

from dataclasses import asdict

from phospy.api import (
    PredictionRunConfig,
    SignalomeRunConfig,
    SimpleKinaseWorkflowConfigSnapshot,
)
from phospy.prediction import KinaseProfilePolicy
from phospy.signalomes import SignalomeModuleSelectionPolicy


def test_prediction_run_config_constructor_normalizes_profile_policy() -> None:
    config = PredictionRunConfig(
        profile_policy={"missing_value_strategy": "median_skipna"}
    )

    assert isinstance(config.profile_policy, KinaseProfilePolicy)
    assert config.profile_policy.missing_value_strategy == "median_skipna"


def test_prediction_run_config_from_value_normalizes_profile_policy() -> None:
    config = PredictionRunConfig.from_value(
        {"profile_policy": {"missing_value_strategy": "median_skipna"}}
    )

    assert isinstance(config.profile_policy, KinaseProfilePolicy)
    assert config.profile_policy.missing_value_strategy == "median_skipna"


def test_prediction_config_snapshot_serializes_normalized_profile_policy() -> None:
    snapshot = SimpleKinaseWorkflowConfigSnapshot.from_workflow_inputs(
        prediction_config={
            "profile_policy": {"missing_value_strategy": "median_skipna"}
        }
    )

    assert snapshot.prediction_config["profile_policy"] == {
        "missing_value_strategy": "median_skipna"
    }


def test_signalome_run_config_constructor_normalizes_module_selection_policy() -> None:
    config = SignalomeRunConfig(
        module_selection_policy={
            "strategy": "single_module",
            "primary_threshold": 0.4,
            "fallback_threshold": 0.2,
            "max_clusters": 6,
        }
    )

    assert isinstance(config.module_selection_policy, SignalomeModuleSelectionPolicy)
    assert config.module_selection_policy.strategy == "single_module"


def test_signalome_run_config_from_value_normalizes_module_selection_policy() -> None:
    config = SignalomeRunConfig.from_value(
        {
            "module_selection_policy": {
                "strategy": "single_module",
                "primary_threshold": 0.4,
                "fallback_threshold": 0.2,
                "max_clusters": 6,
            }
        }
    )

    assert isinstance(config.module_selection_policy, SignalomeModuleSelectionPolicy)
    assert config.module_selection_policy.strategy == "single_module"


def test_signalome_run_config_serialization_uses_normalized_policy_payload() -> None:
    config = SignalomeRunConfig(
        module_selection_policy={
            "strategy": "single_module",
            "primary_threshold": 0.4,
            "fallback_threshold": 0.2,
            "max_clusters": 6,
        }
    )

    assert asdict(config)["module_selection_policy"] == {
        "strategy": "single_module",
        "primary_threshold": 0.4,
        "fallback_threshold": 0.2,
        "max_clusters": 6,
    }
