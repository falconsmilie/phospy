from __future__ import annotations

from pathlib import Path

from phospy.activities import KinaseActivityAnalyzer, KinaseActivityResult


def test_activity_analysis_and_result_models_live_under_activities_package() -> None:
    assert KinaseActivityAnalyzer.__module__ == "phospy.activities.analysis"
    assert KinaseActivityResult.__module__ == "phospy.activities.results"


def test_legacy_activity_analyzer_module_has_been_removed() -> None:
    package_root = Path(__file__).resolve().parents[1] / "src" / "phospy"
    assert not (package_root / "activities" / "analyzer.py").exists()
    assert (package_root / "activities" / "analysis.py").exists()
    assert (package_root / "activities" / "results.py").exists()
