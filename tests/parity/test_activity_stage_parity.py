from __future__ import annotations

import pandas as pd
import pandas.testing as pdt
import pytest

from phospy.activities.scoring import compute_activity_from_inputs
from phospy.validation.workflows.activity import KinaseActivityInputValidator
from tests.support.parity_reporting import (
    format_shape,
    record_parity_metrics,
)
from tests.support.rewrite_fixture_data import (
    ACTIVITY_PARITY_FIXTURE_FILES,
    ACTIVITY_REFERENCE_PROVENANCE,
    activity_parity_fixture_paths,
    load_activity_reference_ksea_counts,
    load_activity_reference_ksea_scores,
    load_activity_reference_predmat,
    load_activity_reference_provenance_text,
    load_activity_reference_target_counts,
    load_activity_reference_target_table,
    load_activity_reference_weighted_activity,
    load_rat_l6_phospho,
)

pytestmark = [pytest.mark.parity, pytest.mark.activity_parity]


def _activity_result():
    pred_mat = load_activity_reference_predmat()
    phospho_matrix = load_rat_l6_phospho()
    inputs = KinaseActivityInputValidator().run(
        pred_mat=pred_mat,
        phospho_matrix=phospho_matrix,
        threshold=0.6,
        min_substrates=3,
        top_n_substrates=20,
    )
    return compute_activity_from_inputs(inputs)


def test_activity_parity_fixture_set_is_present_readable_and_provenanced(
    request: pytest.FixtureRequest,
) -> None:
    for fixture_path in activity_parity_fixture_paths():
        assert fixture_path.is_file(), (
            f"missing activity parity fixture: {fixture_path.name}"
        )
        assert fixture_path.stat().st_size > 0, (
            f"empty activity fixture: {fixture_path.name}"
        )
        loaded = pd.read_csv(fixture_path)
        assert not loaded.empty, (
            f"unreadable or empty activity fixture: {fixture_path.name}"
        )

    assert ACTIVITY_REFERENCE_PROVENANCE.is_file()
    provenance_text = load_activity_reference_provenance_text()
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
        notes=(f"fixture provenance: {ACTIVITY_REFERENCE_PROVENANCE.as_posix()}",),
    )


def test_weighted_activity_matches_rewrite_reference_fixture(
    request: pytest.FixtureRequest,
) -> None:
    result = _activity_result()
    expected = load_activity_reference_weighted_activity()
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


def test_ksea_outputs_match_rewrite_reference_fixture(
    request: pytest.FixtureRequest,
) -> None:
    result = _activity_result()
    expected_scores = load_activity_reference_ksea_scores()
    expected_counts = load_activity_reference_ksea_counts()

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


def test_target_outputs_match_rewrite_reference_fixture(
    request: pytest.FixtureRequest,
) -> None:
    result = _activity_result()
    expected_counts = load_activity_reference_target_counts()
    pdt.assert_series_equal(
        result.target_counts.sort_index(),
        expected_counts.sort_index(),
    )

    expected_target_table = load_activity_reference_target_table()
    pdt.assert_frame_equal(
        result.target_table.reset_index(drop=True),
        expected_target_table.reset_index(drop=True),
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
