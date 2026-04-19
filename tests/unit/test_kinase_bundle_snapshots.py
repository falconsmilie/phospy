from __future__ import annotations

import pytest

from phospy.errors.input import PhosPyInputError
from phospy.io.bundles.kinase import KinaseWorkflowConfigSnapshot


def test_kinase_snapshot_payload_round_trip_preserves_fields() -> None:
    payload = {
        "scoring_config": {"min_substrates": 3},
        "prediction_config": {"top_k": 5, "ensemble_size": 7},
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
