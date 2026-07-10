from __future__ import annotations

import pytest

from phospy.api.configs import KinaseAttritionPolicy
from phospy.errors.input import PhosPyInputError
from phospy.io.bundles.kinase import KinaseWorkflowConfigSnapshot


def test_kinase_snapshot_payload_round_trip_preserves_fields() -> None:
    payload = {
        "scoring_config": {
            "min_substrates": 3,
            "include_diagnostic_scoring_tables": False,
            "include_substrate_contributions": True,
            "profile_missing_value_strategy": "median_skipna",
            "profile_self_inclusion_policy": "leave_one_out",
            "attrition_policy": {
                "minimum_reference_overlap_fraction": 0.2,
                "minimum_sequence_supported_fraction": 0.4,
                "minimum_scored_fraction": 0.6,
                "on_violation": "error",
            },
            "allow_mixed_total_protein_quantitative_meaning": True,
            "reference_context_compatibility_policy": "require_known_match",
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
            "method": "simplified_weighted_substrate_activity",
            "threshold": 0.4,
            "min_substrates": 2,
            "top_n_substrates": 8,
            "ksea_min_substrates": 5,
            "ksea_evidence_threshold": None,
            "ksea_p_value_method": "normal_approximation",
            "ksea_adjust_p_values": True,
        },
    }

    snapshot = KinaseWorkflowConfigSnapshot.from_payload(payload)
    assert snapshot.scoring_config.attrition_policy == KinaseAttritionPolicy(
        minimum_reference_overlap_fraction=0.2,
        minimum_sequence_supported_fraction=0.4,
        minimum_scored_fraction=0.6,
        on_violation="error",
    )
    assert snapshot.to_payload() == payload


def test_kinase_snapshot_requires_scoring_config_object() -> None:
    with pytest.raises(PhosPyInputError, match="config snapshot.scoring_config"):
        KinaseWorkflowConfigSnapshot.from_payload(
            {
                "scoring_config": "invalid",
                "prediction_config": {"top_k": 2},
                "activity_config": None,
            }
        )


def test_kinase_snapshot_rejects_partial_payload_without_required_fields() -> None:
    with pytest.raises(
        PhosPyInputError,
        match=(
            "config snapshot.scoring_config is missing required field\\(s\\): "
            "include_diagnostic_scoring_tables, profile_missing_value_strategy"
        ),
    ):
        KinaseWorkflowConfigSnapshot.from_payload(
            {
                "scoring_config": {"min_substrates": 2},
                "prediction_config": {"top_k": 2},
                "activity_config": None,
            }
        )


def test_kinase_snapshot_rejects_removed_ensemble_size_alias() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="contains unsupported field\\(s\\): ensemble_size",
    ):
        KinaseWorkflowConfigSnapshot.from_payload(
            {
                "scoring_config": {
                    "min_substrates": 2,
                    "include_diagnostic_scoring_tables": True,
                    "profile_missing_value_strategy": "strict",
                },
                "prediction_config": {
                    "top_k": 2,
                    "deterministic_max_selected_kinases": 5,
                    "adaptive_ensemble_runs": 5,
                    "mode": "deterministic_ranking",
                    "adaptive_policy": "stable",
                    "n_iterations": 5,
                    "random_state": None,
                    "ensemble_size": 2,
                },
                "activity_config": None,
            }
        )
