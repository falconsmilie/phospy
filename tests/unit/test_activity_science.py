from __future__ import annotations

import pandas as pd
import pandas.testing as pdt
import pytest

from phospy.activities.methods.ksea_zscore import (
    KSEA_STATUS_COMPUTED,
    KSEA_STATUS_INSUFFICIENT_SUBSTRATES,
    KSEA_STATUS_ZERO_BACKGROUND_VARIANCE,
    KseaZScoreActivityMethod,
)
from phospy.activities.models import KinaseActivityInputs, PredMatOverlapSummary
from phospy.activities.scoring import (
    SimplifiedWeightedSubstrateActivityPolicy,
    compute_activity_from_inputs,
)
from phospy.activities.threshold_membership import (
    THRESHOLD_MEMBERSHIP_DESCRIPTION,
    THRESHOLD_MEMBERSHIP_OPERATOR,
    THRESHOLD_MEMBERSHIP_RULE,
    ActivityThresholdMembershipDiagnostics,
    threshold_membership_mask_array,
)
from phospy.errors.workflows import WorkflowBoundaryError
from phospy.scientific_policies import ScientificPolicyId


def _inputs(
    *,
    pred_mat: pd.DataFrame,
    phospho_matrix: pd.DataFrame,
    threshold: float,
    min_substrates: int,
    top_n_substrates: int,
) -> KinaseActivityInputs:
    overlap_count = int(pred_mat.index.intersection(phospho_matrix.index).size)
    return KinaseActivityInputs(
        pred_mat=pred_mat,
        phospho_matrix=phospho_matrix,
        threshold=threshold,
        min_substrates=min_substrates,
        top_n_substrates=top_n_substrates,
        overlap_summary=PredMatOverlapSummary(
            overlap_count=overlap_count,
            pred_mat_rows=int(pred_mat.index.size),
            phospho_rows=int(phospho_matrix.index.size),
        ),
    )


def _ksea_result(
    *,
    pred_mat: pd.DataFrame,
    phospho_matrix: pd.DataFrame,
    evidence_threshold: float = 0.5,
    min_substrates: int = 2,
    adjust_p_values: bool = True,
):
    return KseaZScoreActivityMethod(
        evidence_threshold=evidence_threshold,
        min_substrates=min_substrates,
        adjust_p_values=adjust_p_values,
    ).run(
        _inputs(
            pred_mat=pred_mat,
            phospho_matrix=phospho_matrix,
            threshold=evidence_threshold,
            min_substrates=min_substrates,
            top_n_substrates=1,
        )
    )


def test_ksea_basic_zscore_calculation_matches_hand_computed_values() -> None:
    pred_mat = pd.DataFrame(
        {"K1": [0.9, 0.8, 0.1, 0.2]},
        index=["S1;S1;", "S2;S2;", "S3;S3;", "S4;S4;"],
    )
    phospho = pd.DataFrame({"c1": [1.0, 2.0, 3.0, 4.0]}, index=pred_mat.index.copy())

    result = _ksea_result(
        pred_mat=pred_mat,
        phospho_matrix=phospho,
        evidence_threshold=0.5,
        min_substrates=2,
    )

    assert result.weighted_activity.at["K1", "c1"] == pytest.approx(-1.0954451150103324)
    stats = result.statistics_table
    assert stats is not None
    row = stats.iloc[0]
    assert row["computability_status"] == KSEA_STATUS_COMPUTED
    assert row["p_value"] == pytest.approx(0.27332167829229814)
    assert row["n_substrates"] == 2
    assert row["n_background_sites"] == 4


def test_ksea_activity_scores_exposes_primary_zscore_matrix() -> None:
    pred_mat = pd.DataFrame(
        {"K1": [0.9, 0.8, 0.1, 0.2]},
        index=["S1;S1;", "S2;S2;", "S3;S3;", "S4;S4;"],
    )
    phospho = pd.DataFrame({"c1": [1.0, 2.0, 3.0, 4.0]}, index=pred_mat.index.copy())

    result = _ksea_result(
        pred_mat=pred_mat,
        phospho_matrix=phospho,
        evidence_threshold=0.5,
        min_substrates=2,
    )

    assert result.activity_method.activity_method_id == "ksea_zscore_v1"
    assert result.activity_scores.at["K1", "c1"] == pytest.approx(-1.0954451150103324)
    pdt.assert_frame_equal(result.activity_scores, result.weighted_activity)
    pdt.assert_frame_equal(result.to_dataframe(), result.activity_scores)


def test_ksea_computes_each_kinase_condition_pair_independently() -> None:
    pred_mat = pd.DataFrame(
        {
            "K1": [0.9, 0.9, 0.1],
            "K2": [0.1, 0.8, 0.8],
        },
        index=["S1;S1;", "S2;S2;", "S3;S3;"],
    )
    phospho = pd.DataFrame(
        {
            "c1": [1.0, 2.0, 3.0],
            "c2": [3.0, 2.0, 1.0],
        },
        index=pred_mat.index.copy(),
    )

    result = _ksea_result(
        pred_mat=pred_mat,
        phospho_matrix=phospho,
        evidence_threshold=0.5,
        min_substrates=2,
    )

    assert result.weighted_activity.at["K1", "c1"] == pytest.approx(-0.7071067811865476)
    assert result.weighted_activity.at["K2", "c1"] == pytest.approx(0.7071067811865476)
    assert result.weighted_activity.at["K1", "c2"] == pytest.approx(0.7071067811865476)
    assert result.weighted_activity.at["K2", "c2"] == pytest.approx(-0.7071067811865476)


def test_ksea_reports_insufficient_substrates_without_dropping_pairs() -> None:
    pred_mat = pd.DataFrame(
        {"K1": [0.9, 0.1]},
        index=["S1;S1;", "S2;S2;"],
    )
    phospho = pd.DataFrame({"c1": [1.0, 2.0]}, index=pred_mat.index.copy())

    result = _ksea_result(
        pred_mat=pred_mat,
        phospho_matrix=phospho,
        evidence_threshold=0.5,
        min_substrates=2,
    )

    assert pd.isna(result.weighted_activity.at["K1", "c1"])
    stats = result.statistics_table
    assert stats is not None
    assert int(stats.shape[0]) == 1
    assert stats.at[0, "computability_status"] == KSEA_STATUS_INSUFFICIENT_SUBSTRATES


def test_ksea_evidence_threshold_is_inclusive_and_ignores_missing_values() -> None:
    pred_mat = pd.DataFrame(
        {"K1": [0.5, 0.49, float("nan")]},
        index=["S1;S1;", "S2;S2;", "S3;S3;"],
    )
    phospho = pd.DataFrame({"c1": [1.0, 2.0, 3.0]}, index=pred_mat.index.copy())

    result = _ksea_result(
        pred_mat=pred_mat,
        phospho_matrix=phospho,
        evidence_threshold=0.5,
        min_substrates=1,
    )

    assert result.target_counts.to_dict() == {"K1": 1}
    stats = result.statistics_table
    assert stats is not None
    assert stats.at[0, "n_substrates"] == 1
    assert stats.at[0, "computability_status"] == KSEA_STATUS_COMPUTED


def test_activity_threshold_membership_policy_is_explicit_and_inclusive() -> None:
    assert THRESHOLD_MEMBERSHIP_RULE == "score >= threshold"
    assert THRESHOLD_MEMBERSHIP_OPERATOR == ">="
    assert (
        THRESHOLD_MEMBERSHIP_DESCRIPTION
        == "scores greater than or equal to the threshold are included"
    )


def test_activity_threshold_membership_boundary_below_equal_above_is_centralised() -> (
    None
):
    mask = threshold_membership_mask_array(
        pd.Series([0.49, 0.5, 0.51], dtype=float).to_numpy(dtype=float, copy=False),
        threshold=0.5,
    )
    assert mask.tolist() == [False, True, True]


def test_activity_threshold_membership_diagnostics_from_payload_parses_numeric_values() -> (
    None
):
    diagnostics = ActivityThresholdMembershipDiagnostics.from_payload(
        {
            "threshold_parameter": "threshold",
            "threshold_value": "0.5",
            "operator": ">=",
            "rule": "score >= threshold",
            "description": "scores greater than or equal to the threshold are included",
        }
    )
    assert diagnostics.threshold_value == pytest.approx(0.5)


def test_activity_threshold_membership_diagnostics_from_payload_preserves_float_coercion() -> (
    None
):
    diagnostics = ActivityThresholdMembershipDiagnostics.from_payload(
        {
            "threshold_parameter": "threshold",
            "threshold_value": True,
            "operator": ">=",
            "rule": "score >= threshold",
            "description": "scores greater than or equal to the threshold are included",
        }
    )
    assert diagnostics.threshold_value == pytest.approx(1.0)


def test_ksea_diagnostics_report_threshold_operator_and_description() -> None:
    pred_mat = pd.DataFrame(
        {"K1": [0.49, 0.5, 0.51]},
        index=["S1;S1;", "S2;S2;", "S3;S3;"],
    )
    phospho = pd.DataFrame({"c1": [1.0, 2.0, 3.0]}, index=pred_mat.index.copy())

    result = _ksea_result(
        pred_mat=pred_mat,
        phospho_matrix=phospho,
        evidence_threshold=0.5,
        min_substrates=1,
    )

    assert result.threshold_membership_diagnostics is not None
    assert result.threshold_membership_diagnostics.threshold_parameter == (
        "evidence_threshold"
    )
    assert result.threshold_membership_diagnostics.threshold_value == pytest.approx(0.5)
    assert result.threshold_membership_diagnostics.operator == (
        THRESHOLD_MEMBERSHIP_OPERATOR
    )

    stats = result.statistics_table
    assert stats is not None
    assert stats.at[0, "evidence_threshold_operator"] == THRESHOLD_MEMBERSHIP_OPERATOR
    assert (
        stats.at[0, "evidence_threshold_description"]
        == THRESHOLD_MEMBERSHIP_DESCRIPTION
    )


def test_weighted_diagnostics_report_threshold_operator_and_description() -> None:
    pred_mat = pd.DataFrame(
        {"K1": [0.49, 0.5, 0.51]},
        index=["S1;S1;", "S2;S2;", "S3;S3;"],
    )
    phospho = pd.DataFrame({"c1": [1.0, 2.0, 3.0]}, index=pred_mat.index.copy())

    result = compute_activity_from_inputs(
        _inputs(
            pred_mat=pred_mat,
            phospho_matrix=phospho,
            threshold=0.5,
            min_substrates=1,
            top_n_substrates=3,
        )
    )

    assert result.threshold_membership_diagnostics is not None
    assert result.threshold_membership_diagnostics.threshold_parameter == "threshold"
    assert result.threshold_membership_diagnostics.threshold_value == pytest.approx(0.5)
    assert result.threshold_membership_diagnostics.operator == (
        THRESHOLD_MEMBERSHIP_OPERATOR
    )
    assert (
        result.threshold_membership_diagnostics.description
        == THRESHOLD_MEMBERSHIP_DESCRIPTION
    )


def test_weighted_and_ksea_share_boundary_threshold_membership_and_counts() -> None:
    pred_mat = pd.DataFrame(
        {"K1": [0.49, 0.5, 0.51]},
        index=["S1;S1;", "S2;S2;", "S3;S3;"],
    )
    phospho = pd.DataFrame({"c1": [1.0, 2.0, 3.0]}, index=pred_mat.index.copy())

    weighted = compute_activity_from_inputs(
        _inputs(
            pred_mat=pred_mat,
            phospho_matrix=phospho,
            threshold=0.5,
            min_substrates=1,
            top_n_substrates=3,
        )
    )
    ksea = _ksea_result(
        pred_mat=pred_mat,
        phospho_matrix=phospho,
        evidence_threshold=0.5,
        min_substrates=1,
    )

    expected_sites = {"S2;S2;", "S3;S3;"}
    assert weighted.thresholded_substrate_counts.to_dict() == {"K1": 2}
    assert weighted.target_counts.to_dict() == {"K1": 2}
    assert ksea.thresholded_substrate_counts.to_dict() == {"K1": 2}
    assert ksea.target_counts.to_dict() == {"K1": 2}

    weighted_sites = set(weighted.target_table.loc[:, "site_id"].astype(str))
    ksea_sites = set(ksea.target_table.loc[:, "site_id"].astype(str))
    assert weighted_sites == expected_sites
    assert ksea_sites == expected_sites
    assert weighted_sites == ksea_sites
    assert "S1;S1;" not in weighted_sites

    assert weighted.thresholded_substrate_mean_activity.at["K1", "c1"] == pytest.approx(
        2.5
    )
    assert ksea.activity_substrate_counts is not None
    assert ksea.activity_substrate_counts.at["K1", "c1"] == 2
    stats = ksea.statistics_table
    assert stats is not None
    assert stats.at[0, "n_substrates"] == 2
    assert stats.at[0, "evidence_threshold_operator"] == THRESHOLD_MEMBERSHIP_OPERATOR
    assert (
        stats.at[0, "evidence_threshold_description"]
        == THRESHOLD_MEMBERSHIP_DESCRIPTION
    )


def test_ksea_reports_zero_background_variance_as_not_computable() -> None:
    pred_mat = pd.DataFrame(
        {"K1": [0.8, 0.8, 0.8]},
        index=["S1;S1;", "S2;S2;", "S3;S3;"],
    )
    phospho = pd.DataFrame({"c1": [5.0, 5.0, 5.0]}, index=pred_mat.index.copy())

    result = _ksea_result(
        pred_mat=pred_mat,
        phospho_matrix=phospho,
        evidence_threshold=0.5,
        min_substrates=1,
    )

    stats = result.statistics_table
    assert stats is not None
    assert stats.at[0, "computability_status"] == KSEA_STATUS_ZERO_BACKGROUND_VARIANCE
    assert pd.isna(stats.at[0, "z_score"])


def test_ksea_excludes_non_finite_phosphosite_values_per_condition() -> None:
    pred_mat = pd.DataFrame(
        {"K1": [0.9, 0.9, 0.9]},
        index=["S1;S1;", "S2;S2;", "S3;S3;"],
    )
    phospho = pd.DataFrame(
        {
            "c1": [1.0, float("nan"), 3.0],
            "c2": [float("nan"), float("nan"), 2.0],
        },
        index=pred_mat.index.copy(),
    )

    result = _ksea_result(
        pred_mat=pred_mat,
        phospho_matrix=phospho,
        evidence_threshold=0.5,
        min_substrates=2,
    )

    stats = result.statistics_table
    assert stats is not None
    assert result.activity_substrate_counts is not None
    c1 = stats.loc[stats["condition"] == "c1"].iloc[0]
    c2 = stats.loc[stats["condition"] == "c2"].iloc[0]
    assert c1["n_background_sites"] == 2
    assert c1["n_substrates"] == 2
    assert c1["computability_status"] == KSEA_STATUS_COMPUTED
    assert c2["n_background_sites"] == 1
    assert c2["n_substrates"] == 1
    assert c2["computability_status"] == KSEA_STATUS_INSUFFICIENT_SUBSTRATES
    assert result.activity_substrate_counts.at["K1", "c1"] == 2
    assert result.activity_substrate_counts.at["K1", "c2"] == 1
    assert result.thresholded_substrate_counts.to_dict() == {"K1": 3}
    assert result.target_counts.to_dict() == {"K1": 3}
    assert (
        result.count_field_semantics["thresholded_substrate_counts"]
        == "global post-threshold evidence membership count before "
        "condition-specific finite-value filtering"
    )


def test_ksea_p_value_uses_two_sided_normal_approximation() -> None:
    pred_mat = pd.DataFrame(
        {"K1": [0.9, 0.8, 0.1, 0.2]},
        index=["S1;S1;", "S2;S2;", "S3;S3;", "S4;S4;"],
    )
    phospho = pd.DataFrame({"c1": [1.0, 2.0, 3.0, 4.0]}, index=pred_mat.index.copy())

    result = _ksea_result(
        pred_mat=pred_mat,
        phospho_matrix=phospho,
        evidence_threshold=0.5,
        min_substrates=2,
    )
    stats = result.statistics_table
    assert stats is not None
    assert stats.at[0, "p_value"] == pytest.approx(0.27332167829229814)


def test_ksea_q_values_are_benjamini_hochberg_adjusted_per_condition() -> None:
    pred_mat = pd.DataFrame(
        {
            "K1": [0.9, 0.9, 0.1, 0.1],
            "K2": [0.1, 0.1, 0.9, 0.9],
            "K3": [0.9, 0.1, 0.1, 0.9],
        },
        index=["S1;S1;", "S2;S2;", "S3;S3;", "S4;S4;"],
    )
    phospho = pd.DataFrame({"c1": [1.0, 2.0, 3.0, 4.0]}, index=pred_mat.index.copy())

    result = _ksea_result(
        pred_mat=pred_mat,
        phospho_matrix=phospho,
        evidence_threshold=0.5,
        min_substrates=2,
        adjust_p_values=True,
    )

    stats = result.statistics_table
    assert stats is not None
    c1_rows = stats.loc[stats["condition"] == "c1"].sort_values("kinase")
    q_values = c1_rows.loc[:, "q_value"].to_numpy(dtype=float)
    assert q_values[0] == pytest.approx(0.4099825174384472)
    assert q_values[1] == pytest.approx(0.4099825174384472)
    assert q_values[2] == pytest.approx(1.0)


def test_ksea_activity_substrate_counts_match_statistics_table_n_substrates() -> None:
    pred_mat = pd.DataFrame(
        {
            "K1": [0.9, 0.9, 0.1],
            "K2": [0.9, 0.1, 0.9],
        },
        index=["S1;S1;", "S2;S2;", "S3;S3;"],
    )
    phospho = pd.DataFrame(
        {
            "c1": [1.0, float("nan"), 3.0],
            "c2": [4.0, 5.0, float("nan")],
        },
        index=pred_mat.index.copy(),
    )

    result = _ksea_result(
        pred_mat=pred_mat,
        phospho_matrix=phospho,
        evidence_threshold=0.5,
        min_substrates=1,
    )

    assert result.activity_substrate_counts is not None
    stats = result.statistics_table
    assert stats is not None
    expected = (
        stats.pivot(index="kinase", columns="condition", values="n_substrates")
        .reindex(index=result.activity_substrate_counts.index)
        .reindex(columns=result.activity_substrate_counts.columns)
        .astype("int64")
    )
    expected.index.name = result.activity_substrate_counts.index.name
    expected.columns.name = result.activity_substrate_counts.columns.name
    pdt.assert_frame_equal(result.activity_substrate_counts, expected)


def test_weighted_activity_ignores_missing_values_per_sample() -> None:
    pred_mat = pd.DataFrame(
        {"PRKACA": [0.9, 0.8, 0.7]},
        index=["A;S1;", "B;S2;", "C;S3;"],
    )
    phospho_matrix = pd.DataFrame(
        {
            "phospho_corrected_1": [10.0, float("nan"), 1.0],
            "phospho_corrected_2": [20.0, 6.0, float("nan")],
        },
        index=["A;S1;", "B;S2;", "C;S3;"],
    )

    result = compute_activity_from_inputs(
        _inputs(
            pred_mat=pred_mat,
            phospho_matrix=phospho_matrix,
            threshold=0.6,
            min_substrates=3,
            top_n_substrates=3,
        )
    )

    assert result.weighted_activity.at[
        "PRKACA", "phospho_corrected_1"
    ] == pytest.approx(6.0625)
    assert result.weighted_activity.at[
        "PRKACA", "phospho_corrected_2"
    ] == pytest.approx((20 * 0.9 + 6 * 0.8) / (0.9 + 0.8))


def test_weighted_activity_scores_exposes_primary_weighted_matrix() -> None:
    pred_mat = pd.DataFrame(
        {"PRKACA": [0.9, 0.8, 0.7]},
        index=["A;S1;", "B;S2;", "C;S3;"],
    )
    phospho_matrix = pd.DataFrame(
        {
            "phospho_corrected_1": [10.0, float("nan"), 1.0],
            "phospho_corrected_2": [20.0, 6.0, float("nan")],
        },
        index=pred_mat.index.copy(),
    )

    result = compute_activity_from_inputs(
        _inputs(
            pred_mat=pred_mat,
            phospho_matrix=phospho_matrix,
            threshold=0.6,
            min_substrates=3,
            top_n_substrates=3,
        )
    )

    assert result.activity_method.activity_method_id == (
        "simplified_weighted_substrate_activity_v1"
    )
    assert result.activity_scores.at["PRKACA", "phospho_corrected_1"] == pytest.approx(
        6.0625
    )
    assert result.activity_scores.at["PRKACA", "phospho_corrected_2"] == pytest.approx(
        (20 * 0.9 + 6 * 0.8) / (0.9 + 0.8)
    )
    pdt.assert_frame_equal(result.activity_scores, result.weighted_activity)
    pdt.assert_frame_equal(result.to_dataframe(), result.activity_scores)


def test_thresholded_substrate_mean_activity_respects_threshold_and_min_substrates() -> (
    None
):
    pred_mat = pd.DataFrame(
        {
            "MAP2K6": [0.9, 0.8, 0.2],
            "AKT1": [0.95, 0.7, 0.61],
        },
        index=["A;S1;", "B;S2;", "C;S3;"],
    )
    phospho_matrix = pd.DataFrame(
        {
            "sample_a": [1.0, 2.0, 3.0],
            "sample_b": [2.0, 4.0, 6.0],
        },
        index=pred_mat.index.copy(),
    )

    result = compute_activity_from_inputs(
        _inputs(
            pred_mat=pred_mat,
            phospho_matrix=phospho_matrix,
            threshold=0.6,
            min_substrates=2,
            top_n_substrates=3,
        )
    )

    assert result.thresholded_substrate_counts.to_dict() == {"AKT1": 3, "MAP2K6": 2}
    assert result.activity_substrate_counts is None
    assert result.thresholded_substrate_mean_activity.at[
        "MAP2K6", "sample_a"
    ] == pytest.approx(1.5)
    assert result.thresholded_substrate_mean_activity.at[
        "AKT1", "sample_b"
    ] == pytest.approx(4.0)


def test_thresholded_substrate_mean_activity_ignores_missing_values_per_sample() -> (
    None
):
    pred_mat = pd.DataFrame(
        {"PRKACA": [0.9, 0.8, 0.7]},
        index=["A;S1;", "B;S2;", "C;S3;"],
    )
    phospho_matrix = pd.DataFrame(
        {
            "phospho_corrected_1": [10.0, float("nan"), 1.0],
            "phospho_corrected_2": [20.0, 6.0, float("nan")],
        },
        index=pred_mat.index.copy(),
    )

    result = compute_activity_from_inputs(
        _inputs(
            pred_mat=pred_mat,
            phospho_matrix=phospho_matrix,
            threshold=0.6,
            min_substrates=3,
            top_n_substrates=3,
        )
    )

    assert result.thresholded_substrate_mean_activity.at[
        "PRKACA", "phospho_corrected_1"
    ] == pytest.approx(5.5)
    assert result.thresholded_substrate_mean_activity.at[
        "PRKACA", "phospho_corrected_2"
    ] == pytest.approx(13.0)
    assert result.thresholded_substrate_counts.to_dict() == {"PRKACA": 3}


def test_top_n_substrate_selection_is_deterministic_for_ties() -> None:
    pred_mat = pd.DataFrame(
        {"MAP2K6": [0.9, 0.9, 0.2]},
        index=["P1;S1;", "P2;S2;", "P3;S3;"],
    )
    phospho_matrix = pd.DataFrame(
        {"sample_a": [10.0, 1.0, 100.0]},
        index=pred_mat.index.copy(),
    )

    result = compute_activity_from_inputs(
        _inputs(
            pred_mat=pred_mat,
            phospho_matrix=phospho_matrix,
            threshold=0.6,
            min_substrates=2,
            top_n_substrates=2,
        )
    )

    assert result.weighted_activity.at["MAP2K6", "sample_a"] == pytest.approx(5.5)


def test_target_count_and_target_table_outputs_are_consistent() -> None:
    pred_mat = pd.DataFrame(
        {
            "MAP2K6": [0.9, 0.5, 0.0],
            "AKT1": [0.4, 0.2, 0.1],
        },
        index=["MAPK14;Y182;", "GSK3B;S9;", "AKT1;T308;"],
    )
    phospho_matrix = pd.DataFrame(
        {
            "sample_a": [1.0, 2.0, 3.0],
            "sample_b": [2.0, 4.0, 6.0],
        },
        index=pred_mat.index.copy(),
    )

    result = compute_activity_from_inputs(
        _inputs(
            pred_mat=pred_mat,
            phospho_matrix=phospho_matrix,
            threshold=0.3,
            min_substrates=2,
            top_n_substrates=2,
        )
    )

    assert result.target_counts.to_dict() == {"MAP2K6": 2, "AKT1": 1}
    assert set(result.target_table.columns) == {"site_id", "kinase", "score"}
    assert int(result.target_table.shape[0]) == 3


def test_activity_stage_raises_when_all_candidates_are_filtered() -> None:
    pred_mat = pd.DataFrame(
        {"MAP2K6": [0.8, 0.7]},
        index=["MAPK14;Y182;", "GSK3B;S9;"],
    )
    phospho_matrix = pd.DataFrame(
        {"sample_a": [1.0, 2.0], "sample_b": [2.0, 4.0]},
        index=pred_mat.index.copy(),
    )

    with pytest.raises(
        WorkflowBoundaryError, match="seam=kinase.activity.valid_candidates"
    ):
        compute_activity_from_inputs(
            _inputs(
                pred_mat=pred_mat,
                phospho_matrix=phospho_matrix,
                threshold=0.95,
                min_substrates=3,
                top_n_substrates=2,
            )
        )


def test_activity_policy_metadata_captures_runtime_parameters() -> None:
    policy = SimplifiedWeightedSubstrateActivityPolicy(
        threshold=0.6,
        min_substrates=3,
        top_n_substrates=20,
    )
    record = policy.record

    assert record.id == ScientificPolicyId.SIMPLIFIED_WEIGHTED_SUBSTRATE_ACTIVITY
    assert record.parameters["threshold"] == pytest.approx(0.6)
    assert record.parameters["min_substrates"] == 3
    assert record.parameters["top_n_substrates"] == 20


def test_activity_result_exposes_explicit_method_identity_without_changing_scores() -> (
    None
):
    pred_mat = pd.DataFrame(
        {"PRKACA": [0.9, 0.8, 0.7]},
        index=["A;S1;", "B;S2;", "C;S3;"],
    )
    phospho_matrix = pd.DataFrame(
        {
            "phospho_corrected_1": [10.0, float("nan"), 1.0],
            "phospho_corrected_2": [20.0, 6.0, float("nan")],
        },
        index=pred_mat.index.copy(),
    )

    result = compute_activity_from_inputs(
        _inputs(
            pred_mat=pred_mat,
            phospho_matrix=phospho_matrix,
            threshold=0.6,
            min_substrates=3,
            top_n_substrates=3,
        )
    )

    assert result.activity_method.activity_method_id == (
        "simplified_weighted_substrate_activity_v1"
    )
    assert result.activity_method.activity_method_family == (
        "heuristic_weighted_substrate_score"
    )
    assert result.activity_method.activity_method_label == (
        "simplified weighted substrate activity"
    )
    assert result.activity_method.is_ksea is False
    assert result.activity_method.is_phosr_kinase_activity_equivalent is False
    assert result.weighted_activity.at[
        "PRKACA", "phospho_corrected_1"
    ] == pytest.approx(6.0625)
    assert result.weighted_activity.at[
        "PRKACA", "phospho_corrected_2"
    ] == pytest.approx((20 * 0.9 + 6 * 0.8) / (0.9 + 0.8))
