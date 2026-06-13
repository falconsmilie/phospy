from __future__ import annotations

from phospy.science.differential.models import DifferentialAnalysisResult


def test_enrichment_workflow_does_not_extend_differential_result_model() -> None:
    assert not hasattr(DifferentialAnalysisResult, "enrichment")
    assert not hasattr(DifferentialAnalysisResult, "enrich")
    assert not hasattr(DifferentialAnalysisResult, "run_enrichment")
