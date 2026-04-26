from __future__ import annotations

import pytest

from phospy.errors.input import PhosPyInputError
from phospy.io.bundles.kinase import KinaseWorkflowConfigSnapshot


def test_kinase_snapshot_payload_round_trip_preserves_fields() -> None:
    payload = {
        "scoring_config": {
            "min_substrates": 3,
            "include_diagnostic_scoring_tables": False,
            "profile_missing_value_strategy": "median_skipna",
        },
        "prediction_config": {
            "top_k": 5,
            "deterministic_max_selected_kinases": 7,
            "adaptive_ensemble_runs": 7,
            "mode": "adaptive_ensemble",
            "adaptive_policy": "r_parity",
            "n_iterations": 3,
            "random_state": 11,
        },
        "activity_config": {
            "enabled": True,
            "threshold": 0.4,
            "min_substrates": 2,
            "top_n_substrates": 8,
        },
    }

    snapshot = KinaseWorkflowConfigSnapshot.from_payload(payload)
    assert snapshot.to_payload() == payload


def test_kinase_snapshot_requires_scoring_config_object() -> None:
    with pytest.raises(PhosPyInputError, match="config snapshot.scoring_config"):
        KinaseWorkflowConfigSnapshot.from_payload(
            {
                "scoring_config": "invalid",
                "prediction_config": {"top_k": 2, "ensemble_size": 2},
                "activity_config": None,
            }
        )


def test_kinase_snapshot_compatibility_payload_defaults_diagnostic_tables_to_true() -> (
    None
):
    snapshot = KinaseWorkflowConfigSnapshot.from_payload(
        {
            "scoring_config": {"min_substrates": 2},
            "prediction_config": {"top_k": 2, "ensemble_size": 2},
            "activity_config": None,
        }
    )
    assert snapshot.scoring_config.include_diagnostic_scoring_tables is True
    assert snapshot.scoring_config.profile_missing_value_strategy == "strict"
    assert snapshot.prediction_config.mode == "deterministic_ranking"
    assert snapshot.prediction_config.adaptive_policy == "stable"
    assert snapshot.prediction_config.n_iterations == 5
    assert snapshot.prediction_config.random_state is None


def test_kinase_snapshot_rejects_legacy_and_explicit_prediction_size_conflicts() -> (
    None
):
    with pytest.raises(PhosPyInputError, match="cannot be combined"):
        KinaseWorkflowConfigSnapshot.from_payload(
            {
                "scoring_config": {"min_substrates": 2},
                "prediction_config": {
                    "top_k": 2,
                    "ensemble_size": 2,
                    "adaptive_ensemble_runs": 4,
                },
                "activity_config": None,
            }
        )
