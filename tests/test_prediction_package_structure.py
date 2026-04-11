from __future__ import annotations

from pathlib import Path


def test_prediction_package_contains_explicit_engine_result_and_scoring_modules() -> (
    None
):
    package_root = Path(__file__).resolve().parents[1] / "src" / "phospy" / "prediction"

    assert (package_root / "engines.py").exists()
    assert (package_root / "results.py").exists()
    assert (package_root / "scoring.py").exists()


def test_legacy_prediction_modules_and_root_scoring_module_have_been_removed() -> None:
    package_root = Path(__file__).resolve().parents[1] / "src" / "phospy"
    prediction_root = package_root / "prediction"

    assert not (package_root / "scoring.py").exists()
    assert not (prediction_root / "models.py").exists()
    assert not (prediction_root / "service.py").exists()
    assert not (prediction_root / "workflows.py").exists()
