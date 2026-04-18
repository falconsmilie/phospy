from __future__ import annotations

from pathlib import Path

import pandas as pd
import pandas.testing as pdt
import pytest

from legacy_archive.phospy_legacy.activities.scoring import build_kinase_target_table
from phospy.activities.scoring import compute_activity_from_inputs
from phospy.validation.workflows.activity import KinaseActivityInputValidator

ROOT = Path(__file__).resolve().parents[2]
R_REFERENCE_L6 = ROOT / "tests_legacy" / "fixtures" / "r_reference_l6"

pytestmark = pytest.mark.parity


def _activity_result():
    pred_mat = pd.read_csv(R_REFERENCE_L6 / "predMat.csv", index_col=0)
    phospho_matrix = pd.read_csv(R_REFERENCE_L6 / "l6_phospho_matrix.csv", index_col=0)
    inputs = KinaseActivityInputValidator().run(
        pred_mat=pred_mat,
        phospho_matrix=phospho_matrix,
        threshold=0.6,
        min_substrates=3,
        top_n_substrates=20,
    )
    return compute_activity_from_inputs(inputs), pred_mat


def test_weighted_activity_matches_legacy_reference_fixture() -> None:
    result, _ = _activity_result()
    expected = pd.read_csv(R_REFERENCE_L6 / "kinase_activity_matrix.csv", index_col=0)
    expected.index.name = "kinase"
    pdt.assert_frame_equal(result.weighted_activity.sort_index(), expected.sort_index())


def test_ksea_outputs_match_legacy_reference_fixture() -> None:
    result, _ = _activity_result()
    expected_scores = pd.read_csv(R_REFERENCE_L6 / "ksea_scores.csv", index_col=0)
    expected_scores.index.name = "kinase"
    expected_counts = pd.read_csv(R_REFERENCE_L6 / "ksea_counts.csv", index_col=0).iloc[
        :, 0
    ]
    expected_counts.index.name = "kinase"
    expected_counts.name = "n_substrates"

    pdt.assert_frame_equal(
        result.ksea_scores.sort_index(), expected_scores.sort_index()
    )
    pdt.assert_series_equal(
        result.ksea_counts.sort_index(), expected_counts.sort_index()
    )


def test_target_outputs_match_legacy_reference_fixture_and_kernel() -> None:
    result, pred_mat = _activity_result()
    expected_counts = pd.read_csv(
        R_REFERENCE_L6 / "kinase_target_counts.csv",
        index_col=0,
    ).iloc[:, 0]
    expected_counts.index.name = "kinase"
    expected_counts.name = "n_targets"
    pdt.assert_series_equal(
        result.target_counts.sort_index(),
        expected_counts.sort_index(),
    )

    legacy_target_table = build_kinase_target_table(pred_mat=pred_mat, threshold=0.6)
    pdt.assert_frame_equal(
        result.target_table.reset_index(drop=True),
        legacy_target_table.reset_index(drop=True),
    )
