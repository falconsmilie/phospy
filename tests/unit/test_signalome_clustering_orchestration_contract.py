from __future__ import annotations

from pathlib import Path
from typing import get_type_hints

import phospy.science.signalomes.clustering.exact_python as exact_backend_exports
import phospy.science.signalomes.clustering.orchestration as orchestration_exports
from phospy.science.signalomes.clustering.models import SignalomeClusteringEngineRequest
from phospy.science.signalomes.clustering.policies import (
    SignalomeCandidateScoringPolicy,
)


def test_candidate_scoring_policy_boundary_is_typed() -> None:
    hints = get_type_hints(SignalomeClusteringEngineRequest)
    assert hints["candidate_scoring_policy"] == (SignalomeCandidateScoringPolicy | None)


def test_candidate_scoring_policy_boundary_does_not_use_arg_type_ignore() -> None:
    root = Path(__file__).resolve().parents[2]
    adapter_source = (
        root
        / "src"
        / "phospy"
        / "science"
        / "signalomes"
        / "clustering"
        / "tree_engine_adapter.py"
    ).read_text(encoding="utf-8")
    orchestration_source = (
        root
        / "src"
        / "phospy"
        / "science"
        / "signalomes"
        / "clustering"
        / "orchestration.py"
    ).read_text(encoding="utf-8")
    assert "candidate_scoring_policy" in adapter_source
    assert "type: ignore[arg-type]" not in adapter_source
    assert "type: ignore[arg-type]" not in orchestration_source


def test_orchestration_and_exact_backend_all_do_not_publish_internal_names() -> None:
    assert all(not symbol.startswith("_") for symbol in orchestration_exports.__all__)
    assert all(not symbol.startswith("_") for symbol in exact_backend_exports.__all__)
