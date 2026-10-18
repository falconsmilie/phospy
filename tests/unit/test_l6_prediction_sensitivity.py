from __future__ import annotations

import hashlib
import inspect
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from tests.support.l6_prediction_sensitivity import (
    REFERENCE_CANDIDATES_PATH,
    REFERENCE_MATRIX_PATH,
    apply_boundary_perturbation,
    build_l6_prediction_sensitivity_report,
    derive_l6_prediction_sensitivity_budget,
    render_l6_prediction_sensitivity_report,
)

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT
    / "tests"
    / "fixtures"
    / "rewrite_parity"
    / "r_reference_l6_prediction"
    / "sensitivity_report.json"
)


def _reference() -> pd.DataFrame:
    frame = pd.read_csv(REFERENCE_MATRIX_PATH, index_col=0).astype(float)
    frame.columns = frame.columns.str.strip().str.upper()
    return frame


def _ranked_selected_rows(frame: pd.DataFrame) -> tuple[tuple[str, str], ...]:
    candidates = pd.read_csv(REFERENCE_CANDIDATES_PATH).astype(str)
    candidates.loc[:, "kinase"] = candidates.loc[:, "kinase"].str.strip().str.upper()
    return tuple(
        (str(kinase), str(site_id))
        for kinase, group in candidates.groupby("kinase", sort=False)
        for site_id in frame.loc[group.loc[:, "site_id"].tolist(), kinase]
        .sort_values(ascending=False, kind="mergesort")
        .index
    )


def test_l6_sensitivity_budget_is_reference_only_deterministic_and_current(
    monkeypatch,
) -> None:
    assert inspect.signature(derive_l6_prediction_sensitivity_budget).parameters == {}
    first = build_l6_prediction_sensitivity_report()
    monkeypatch.setenv("PHOSPY_L6_PREDICTION_LANE", "unrelated-platform-observation")
    monkeypatch.setenv("PHOSPY_L6_PREDICTION_EVIDENCE_PATH", "unrelated.json")
    second = build_l6_prediction_sensitivity_report()

    assert first == second
    assert (
        REPORT.read_text(encoding="utf-8") == render_l6_prediction_sensitivity_report()
    )
    assert (
        first["reference"]["prediction_matrix"]["sha256"]
        == hashlib.sha256(REFERENCE_MATRIX_PATH.read_bytes()).hexdigest()
    )
    assert "not inputs" in first["derivation_independence"]


def test_l6_sensitivity_report_is_byte_stable_across_processes() -> None:
    command = [
        sys.executable,
        "-c",
        (
            "import sys; "
            "from tests.support.l6_prediction_sensitivity import "
            "render_l6_prediction_sensitivity_report; "
            "sys.stdout.buffer.write("
            "render_l6_prediction_sensitivity_report().encode('utf-8'))"
        ),
    ]

    first = subprocess.run(
        command,
        cwd=ROOT,
        check=True,
        capture_output=True,
    ).stdout
    second = subprocess.run(
        command,
        cwd=ROOT,
        check=True,
        capture_output=True,
    ).stdout

    assert first == second
    assert first.decode("utf-8") == REPORT.read_text(encoding="utf-8")


def test_permitted_perturbation_retains_every_protected_decision() -> None:
    reference = _reference()
    budget = derive_l6_prediction_sensitivity_budget()
    perturbed = apply_boundary_perturbation(
        reference,
        amplitude=budget.max_abs_diff_ceiling * (1.0 - 1e-9),
    )
    delta = (perturbed - reference).abs().to_numpy(dtype=float)

    assert _ranked_selected_rows(perturbed) == _ranked_selected_rows(reference)
    assert all(
        perturbed.loc[:, kinase]
        .sort_values(ascending=False, kind="mergesort")
        .index.equals(
            reference.loc[:, kinase]
            .sort_values(ascending=False, kind="mergesort")
            .index
        )
        for kinase in reference.columns
    )
    assert perturbed.index.equals(reference.index)
    assert perturbed.columns.equals(reference.columns)
    assert np.array_equal(np.isfinite(perturbed), np.isfinite(reference))
    assert float(delta.mean()) <= budget.mean_abs_diff_ceiling
    assert float(delta.max()) <= budget.max_abs_diff_ceiling
    report = build_l6_prediction_sensitivity_report()
    assert report["protected_decisions"]["structural"] == [
        "class assignment",
        "finite and missing-value masks",
        "site and kinase attrition",
        "provenance records",
        "diagnostic records",
    ]


def test_beyond_safe_region_changes_candidate_decision_before_correlation_fails() -> (
    None
):
    reference = _reference()
    budget = derive_l6_prediction_sensitivity_budget()
    amplitude = np.nextafter(budget.boundary.max_abs_failure_boundary, np.inf)
    perturbed = apply_boundary_perturbation(reference, amplitude=amplitude)
    delta = (perturbed - reference).abs().to_numpy(dtype=float)
    correlation = float(
        np.corrcoef(
            reference.to_numpy(dtype=float).ravel(),
            perturbed.to_numpy(dtype=float).ravel(),
        )[0, 1]
    )

    assert _ranked_selected_rows(perturbed) != _ranked_selected_rows(reference)
    assert float(delta.mean()) > budget.mean_abs_diff_ceiling
    assert float(delta.max()) > budget.max_abs_diff_ceiling
    assert correlation >= 0.999999


def test_candidate_measurements_require_a_prederived_budget() -> None:
    source = inspect.getsource(
        __import__(
            "tests.parity.test_l6_prediction_parity",
            fromlist=[
                "test_l6_prediction_matrix_numeric_parity_matches_promoted_reference_surfaces"
            ],
        ).test_l6_prediction_matrix_numeric_parity_matches_promoted_reference_surfaces
    )
    assert source.index("derive_l6_prediction_sensitivity_budget()") < source.index(
        "collect_l6_prediction_parity_metrics()"
    )
