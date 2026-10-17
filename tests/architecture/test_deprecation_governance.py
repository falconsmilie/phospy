from __future__ import annotations

import json
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]


def test_runtime_deprecation_registry_and_warning_filter_are_absent() -> None:
    assert not (PROJECT_ROOT / "src" / "phospy" / "_deprecations.py").exists()
    assert not (PROJECT_ROOT / "src" / "phospy" / "api" / "_compat.py").exists()
    pyproject = (PROJECT_ROOT / "pyproject.toml").read_text(encoding="utf-8")
    assert "PhosPyDeprecationWarning" not in pyproject
    assert "phospy._deprecations" not in pyproject
    baseline = json.loads(
        (PROJECT_ROOT / "architecture/public-api-baseline.json").read_text(
            encoding="utf-8"
        )
    )
    assert baseline["compatibility_registry_present"] is False
    assert "deprecations" not in baseline
