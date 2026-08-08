from __future__ import annotations

from pathlib import Path

from tools.testing.pyright_strict_coverage import pyright_strict_coverage

ROOT = Path(__file__).resolve().parents[2]


def test_pyright_strict_coverage_reports_first_wave_increase() -> None:
    coverage = pyright_strict_coverage(ROOT)

    assert coverage.included_source_files > 0
    assert coverage.strict_source_files >= 95
    assert coverage.strict_source_percent >= 15.0
