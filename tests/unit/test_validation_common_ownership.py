from __future__ import annotations

from pathlib import Path


def test_obsolete_validation_common_forwarding_routes_are_absent() -> None:
    common = (
        Path(__file__).resolve().parents[2] / "src" / "phospy" / "validation" / "common"
    )

    assert not (common / "dataframes.py").exists()
    assert not (common / "missing_values.py").exists()
    assert not (common / "numeric_frames.py").exists()
