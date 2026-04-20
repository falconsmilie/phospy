from __future__ import annotations

from pathlib import Path

import pandas as pd
import pandas.testing as pdt
import pytest

from legacy_archive.phospy_legacy.activities.scoring import (
    build_kinase_target_table as legacy_build_kinase_target_table,
)
from legacy_archive.phospy_legacy.activities.scoring import (
    compute_activity_from_inputs as legacy_compute_activity_from_inputs,
)
from legacy_archive.phospy_legacy.activities.scoring import (
    count_predicted_targets as legacy_count_predicted_targets,
)
from legacy_archive.phospy_legacy.validation.requests.analysis import (
    validate_analysis_request as legacy_validate_analysis_request,
)
from phospy.activities.scoring import compute_activity_from_inputs
from phospy.validation.workflows.activity import KinaseActivityInputValidator
from tests.support.parity_reporting import (
    format_shape,
    record_parity_metrics,
)

ROOT = Path(__file__).resolve().parents[2]
R_REFERENCE_L6 = ROOT / "tests" / "fixtures" / "rewrite_parity" / "r_reference_l6"
R_REFERENCE_PROVENANCE = R_REFERENCE_L6 / "PROVENANCE.md"
ACTIVITY_PARITY_FIXTURE_FILES = (
    "l6_phospho_matrix.csv",
    "native_profile_scores.csv",
    "predMat.csv",
    "kinase_activity_matrix.csv",
    "ksea_scores.csv",
    "ksea_counts.csv",
    "kinase_target_counts.csv",
    "kinase_target_table.csv",
)

pytestmark = [pytest.mark.parity, pytest.mark.activity_parity]


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


def _legacy_activity_result(
    *,
    pred_mat: pd.DataFrame,
    phospho_matrix: pd.DataFrame,
    threshold: float,
    min_substrates: int,
    top_n_substrates: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.Series, pd.Series, pd.DataFrame]:
    validated = legacy_validate_analysis_request(
        pred_mat=pred_mat,
        phospho_matrix=phospho_matrix,
        threshold=threshold,
        min_substrates=min_substrates,
        top_n_substrates=top_n_substrates,
    )
    weighted_activity, ksea_scores, ksea_counts = legacy_compute_activity_from_inputs(
        validated
    )
    weighted_activity.index.name = "kinase"
    ksea_scores.index.name = "kinase"
    target_counts = legacy_count_predicted_targets(
        pred_mat=validated.pred_mat,
        threshold=validated.threshold,
    )
    target_table = legacy_build_kinase_target_table(
        pred_mat=validated.pred_mat,
        threshold=validated.threshold,
    )
    return weighted_activity, ksea_scores, ksea_counts, target_counts, target_table


def test_activity_parity_fixture_set_is_present_readable_and_provenanced(
    request: pytest.FixtureRequest,
) -> None:
    for file_name in ACTIVITY_PARITY_FIXTURE_FILES:
        fixture_path = R_REFERENCE_L6 / file_name
        assert fixture_path.is_file(), f"missing activity parity fixture: {file_name}"
        assert fixture_path.stat().st_size > 0, f"empty activity fixture: {file_name}"
        loaded = pd.read_csv(fixture_path)
        assert not loaded.empty, f"unreadable or empty activity fixture: {file_name}"

    assert R_REFERENCE_PROVENANCE.is_file()
    provenance_text = R_REFERENCE_PROVENANCE.read_text(encoding="utf-8")
    assert "Promoted from" in provenance_text
    assert "Fixture ownership" in provenance_text
    for file_name in ACTIVITY_PARITY_FIXTURE_FILES:
        assert f"`{file_name}`" in provenance_text
    record_parity_metrics(
        request.config,
        family="activity_stage",
        metrics=[
            ("fixture files checked", len(ACTIVITY_PARITY_FIXTURE_FILES)),
        ],
        notes=(f"fixture provenance: {R_REFERENCE_PROVENANCE.as_posix()}",),
    )


def test_weighted_activity_matches_legacy_reference_fixture(
    request: pytest.FixtureRequest,
) -> None:
    result, _ = _activity_result()
    expected = pd.read_csv(R_REFERENCE_L6 / "kinase_activity_matrix.csv", index_col=0)
    expected.index.name = "kinase"
    pdt.assert_frame_equal(result.weighted_activity.sort_index(), expected.sort_index())
    aligned = result.weighted_activity.sort_index() - expected.sort_index()
    absolute_delta = aligned.abs()
    record_parity_metrics(
        request.config,
        family="activity_stage",
        metrics=[
            ("weighted activity shape", format_shape(*result.weighted_activity.shape)),
            (
                "weighted activity mean abs diff",
                float(absolute_delta.to_numpy().mean()),
            ),
            ("weighted activity max abs diff", float(absolute_delta.to_numpy().max())),
        ],
    )


def test_ksea_outputs_match_legacy_reference_fixture(
    request: pytest.FixtureRequest,
) -> None:
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
    score_delta = (result.ksea_scores.sort_index() - expected_scores.sort_index()).abs()
    count_delta = (result.ksea_counts.sort_index() - expected_counts.sort_index()).abs()
    record_parity_metrics(
        request.config,
        family="activity_stage",
        metrics=[
            ("ksea score shape", format_shape(*result.ksea_scores.shape)),
            ("ksea count kinases", int(result.ksea_counts.shape[0])),
            ("ksea total substrate count", int(result.ksea_counts.sum())),
            ("ksea score mean abs diff", float(score_delta.to_numpy().mean())),
            ("ksea score max abs diff", float(score_delta.to_numpy().max())),
            ("ksea count max abs diff", float(count_delta.max())),
        ],
    )


def test_target_outputs_match_legacy_reference_fixture_and_kernel(
    request: pytest.FixtureRequest,
) -> None:
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

    expected_target_table = pd.read_csv(R_REFERENCE_L6 / "kinase_target_table.csv")
    pdt.assert_frame_equal(
        result.target_table.reset_index(drop=True),
        expected_target_table.reset_index(drop=True),
    )

    legacy_target_table = legacy_build_kinase_target_table(
        pred_mat=pred_mat, threshold=0.6
    )
    pdt.assert_frame_equal(
        result.target_table.reset_index(drop=True),
        legacy_target_table.reset_index(drop=True),
    )
    target_count_delta = (
        result.target_counts.sort_index() - expected_counts.sort_index()
    ).abs()
    record_parity_metrics(
        request.config,
        family="activity_stage",
        metrics=[
            ("target count kinases", int(result.target_counts.shape[0])),
            ("target total count", int(result.target_counts.sum())),
            ("target count max abs diff", float(target_count_delta.max())),
            ("target table row count", int(result.target_table.shape[0])),
        ],
    )


def test_activity_kernels_match_legacy_on_ties_missing_values_and_partial_overlap() -> (
    None
):
    pred_mat = pd.DataFrame(
        {
            "MAP2K6": [0.9, 0.9, 0.85, 0.2, 0.0],
            "AKT1": [0.7, 0.61, 0.61, 0.59, 0.0],
        },
        index=["SITE_1", "SITE_2", "SITE_3", "SITE_4", "SITE_5"],
    )
    phospho_matrix = pd.DataFrame(
        {
            "sample_a": [10.0, 1.0, float("nan"), 5.0, -2.0, 99.0],
            "sample_b": [20.0, 2.0, 3.0, float("nan"), -1.0, 88.0],
        },
        index=["SITE_1", "SITE_2", "SITE_3", "SITE_4", "SITE_5", "SITE_EXTRA"],
    )
    threshold = 0.6
    min_substrates = 2
    top_n_substrates = 3

    rewrite_inputs = KinaseActivityInputValidator().run(
        pred_mat=pred_mat,
        phospho_matrix=phospho_matrix,
        threshold=threshold,
        min_substrates=min_substrates,
        top_n_substrates=top_n_substrates,
    )
    rewrite_result = compute_activity_from_inputs(rewrite_inputs)
    (
        legacy_weighted,
        legacy_ksea_scores,
        legacy_ksea_counts,
        legacy_target_counts,
        legacy_target_table,
    ) = _legacy_activity_result(
        pred_mat=pred_mat,
        phospho_matrix=phospho_matrix,
        threshold=threshold,
        min_substrates=min_substrates,
        top_n_substrates=top_n_substrates,
    )

    pdt.assert_frame_equal(
        rewrite_result.weighted_activity.sort_index(),
        legacy_weighted.sort_index(),
    )
    pdt.assert_frame_equal(
        rewrite_result.ksea_scores.sort_index(),
        legacy_ksea_scores.sort_index(),
    )
    pdt.assert_series_equal(
        rewrite_result.ksea_counts.sort_index(),
        legacy_ksea_counts.sort_index(),
    )
    pdt.assert_series_equal(
        rewrite_result.target_counts.sort_index(),
        legacy_target_counts.sort_index(),
    )
    pdt.assert_frame_equal(
        rewrite_result.target_table.reset_index(drop=True),
        legacy_target_table.reset_index(drop=True),
    )
