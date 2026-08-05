from __future__ import annotations

from io import StringIO

import numpy as np
import pandas as pd
import pandas.testing as pdt
import pytest

from phospy.errors.validation import PhosPyValidationError
from phospy.errors.workflows import WorkflowBoundaryError
from phospy.provenance.hashing import fingerprint_table_normalized_axes
from phospy.provenance.scientific_policy_models import ScientificPolicyId
from phospy.science.activities.membership import (
    ACTIVITY_MEMBERSHIP_SELECTION_POLICY_VERSION,
    ACTIVITY_MEMBERSHIP_SOURCE_FIXED_EXTERNAL_REFERENCE,
    KSEA_MEMBERSHIP_ELIGIBLE_REASON,
    KSEA_MEMBERSHIP_MISSING_PROVENANCE_REASON,
    ActivityMembershipSelection,
)
from phospy.science.activities.methods.ksea_zscore import (
    KSEA_STATUS_COMPUTED,
    KSEA_STATUS_INSUFFICIENT_SUBSTRATES,
    KSEA_STATUS_ZERO_BACKGROUND_VARIANCE,
    KseaZScoreActivityMethod,
)
from phospy.science.activities.methods.ssgsea_substrate_enrichment import (
    SSGSEA_Q_VALUE_METHOD_BENJAMINI_HOCHBERG,
    SSGSEA_SIGNIFICANCE_STATUS_AVAILABLE,
    SSGSEA_SIGNIFICANCE_STATUS_P_VALUE_AVAILABLE_Q_VALUE_DISABLED,
    SSGSEA_SIGNIFICANCE_STATUS_UNAVAILABLE_NO_PERMUTATIONS,
    SSGSEA_STATUS_COMPUTED,
    SSGSEA_STATUS_INSUFFICIENT_SUBSTRATES,
    SsgseaSubstrateEnrichmentActivityMethod,
    _derive_ssgsea_permutation_seed,
)
from phospy.science.activities.models import (
    ActivityMethodSummary,
    KinaseActivityInputs,
    KinaseActivityResult,
    KseaZScoreActivityDiagnostics,
    PredMatOverlapSummary,
    SsgseaSubstrateEnrichmentActivityDiagnostics,
    WeightedSubstrateActivityDiagnostics,
)
from phospy.science.activities.scientific_policies import (
    SSGSEA_PERMUTATION_RNG_SEED_POLICY,
    SSGSEA_PERMUTATION_RNG_SEED_POLICY_VERSION,
    SSGSEA_SUBSTRATE_ENRICHMENT_ACTIVITY_POLICY_VERSION,
)
from phospy.science.activities.scoring import (
    SimplifiedWeightedSubstrateActivityPolicy,
    compute_activity_from_inputs,
)
from phospy.science.activities.semantics import (
    ActivityAggregationMetadata,
    ActivityAggregationRecord,
    ActivityInputMatrix,
    ActivityProfileAxis,
    ActivityProfileMetadata,
    ActivityQuantitativeSemantics,
)
from phospy.science.activities.statistics import (
    benjamini_hochberg_q_values,
    two_sided_normal_p_value,
)
from phospy.science.activities.threshold_membership import (
    THRESHOLD_MEMBERSHIP_DESCRIPTION,
    THRESHOLD_MEMBERSHIP_OPERATOR,
    THRESHOLD_MEMBERSHIP_RULE,
    ActivityThresholdMembershipDiagnostics,
    threshold_membership_mask_array,
)
from tests.support.site_keys import site_key_index_from_display_ids


def _site_key_index(display_ids: list[str]) -> pd.Index:
    return site_key_index_from_display_ids(
        display_ids,
        protein_namespace="gene_symbol",
    )


def _with_site_key_index(frame: pd.DataFrame) -> pd.DataFrame:
    labels = frame.index.astype(str).tolist()
    if all(label.startswith("phospy:v1|") for label in labels):
        return frame
    converted = frame.copy(deep=True)
    converted.index = _site_key_index(labels)
    return converted


def _inputs(
    *,
    pred_mat: pd.DataFrame,
    phospho_matrix: pd.DataFrame,
    threshold: float,
    min_substrates: int,
    top_n_substrates: int,
    activity_input: ActivityInputMatrix | None = None,
    membership_selection: ActivityMembershipSelection | None = None,
) -> KinaseActivityInputs:
    pred_mat = _with_site_key_index(pred_mat)
    phospho_matrix = _with_site_key_index(phospho_matrix)
    if activity_input is None:
        activity_input = ActivityInputMatrix.sample_level_abundance(
            phospho_matrix,
            _assume_owned=True,
        )
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
        activity_input=activity_input,
        membership_selection=membership_selection,
    )


def _ksea_result(
    *,
    pred_mat: pd.DataFrame,
    phospho_matrix: pd.DataFrame,
    evidence_threshold: float = 0.5,
    min_substrates: int = 2,
    adjust_p_values: bool = True,
):
    keyed_pred_mat = _with_site_key_index(pred_mat)
    effect_matrix = _with_site_key_index(phospho_matrix)
    return KseaZScoreActivityMethod(
        evidence_threshold=evidence_threshold,
        min_substrates=min_substrates,
        adjust_p_values=adjust_p_values,
    ).run(
        _inputs(
            pred_mat=keyed_pred_mat,
            phospho_matrix=effect_matrix,
            threshold=evidence_threshold,
            min_substrates=min_substrates,
            top_n_substrates=1,
            activity_input=ActivityInputMatrix.standardised_effect(
                effect_matrix,
                _assume_owned=True,
            ),
            membership_selection=_eligible_fixed_membership_selection(
                keyed_pred_mat,
                threshold=evidence_threshold,
            ),
        )
    )


def _eligible_fixed_membership_selection(
    pred_mat: pd.DataFrame,
    *,
    threshold: float,
) -> ActivityMembershipSelection:
    selected_substrates = pred_mat.loc[
        pred_mat.ge(float(threshold)).any(axis=1)
    ].index.astype(str)
    return ActivityMembershipSelection(
        source_category=ACTIVITY_MEMBERSHIP_SOURCE_FIXED_EXTERNAL_REFERENCE,
        selection_method="unit_test_fixed_membership",
        selection_method_version=ACTIVITY_MEMBERSHIP_SELECTION_POLICY_VERSION,
        score_source="fixed_external_reference_membership_scores",
        threshold_top_k_policy={
            "evidence_threshold": float(threshold),
            "evidence_threshold_operator": ">=",
            "top_k": None,
        },
        source_reference_fingerprints=(
            fingerprint_table_normalized_axes(
                pred_mat,
                name="tests.fixed_external_membership_scores",
            ),
        ),
        quantitative_dataset_fingerprint=None,
        consumed_tested_matrix=False,
        selected_kinase_universe=pred_mat.columns.astype(str).tolist(),
        selected_substrate_universe=selected_substrates.tolist(),
        inferential_eligible=True,
        inferential_eligibility_reason=KSEA_MEMBERSHIP_ELIGIBLE_REASON,
    )


def _membership(kinase_to_sites: dict[str, list[str]]) -> pd.DataFrame:
    rows: list[dict[str, str]] = []
    for kinase, display_ids in kinase_to_sites.items():
        for site_key in _site_key_index(display_ids).astype(str).tolist():
            rows.append({"kinase": kinase, "substrate_site": site_key})
    return pd.DataFrame.from_records(rows)


def _statistics_table(
    profile_ids: list[str],
    *,
    include_condition: bool = False,
    condition_values: list[str] | None = None,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for position, profile_id in enumerate(profile_ids):
        row: dict[str, object] = {
            "kinase": "K1",
            "profile_id": profile_id,
            "z_score": 1.0,
            "p_value": 0.05,
            "q_value": 0.05,
            "n_substrates": 2,
            "n_background_sites": 4,
            "evidence_threshold": 0.5,
            "evidence_threshold_operator": ">=",
            "evidence_threshold_description": "unit-test threshold",
            "min_substrates": 1,
            "computability_status": "computed",
            "reason": "",
        }
        if include_condition:
            condition_id = (
                profile_id if condition_values is None else condition_values[position]
            )
            row["condition"] = condition_id
        rows.append(row)
    columns = [
        "kinase",
        "profile_id",
        "z_score",
        "p_value",
        "q_value",
        "n_substrates",
        "n_background_sites",
        "evidence_threshold",
        "evidence_threshold_operator",
        "evidence_threshold_description",
        "min_substrates",
        "computability_status",
        "reason",
    ]
    if include_condition:
        columns.insert(2, "condition")
    return pd.DataFrame.from_records(rows, columns=columns)


def _activity_result_from_statistics_table(
    statistics_table: pd.DataFrame,
    *,
    activity_input: ActivityInputMatrix,
) -> KinaseActivityResult:
    kinase_index = pd.Index(["K1"], name="kinase")
    activity_matrix = pd.DataFrame(
        [[1.0 for _ in activity_input.profile_metadata.profile_ids]],
        index=kinase_index,
        columns=pd.Index(activity_input.profile_metadata.profile_ids),
        dtype=float,
    )
    substrate_count_matrix = pd.DataFrame(
        [[2 for _ in activity_input.profile_metadata.profile_ids]],
        index=kinase_index,
        columns=pd.Index(activity_input.profile_metadata.profile_ids),
        dtype="int64",
    )
    return KinaseActivityResult(
        activity_matrix=activity_matrix,
        substrate_count_matrix=substrate_count_matrix,
        statistics_table=statistics_table,
        input_semantics=activity_input.semantics,
        profile_metadata=activity_input.profile_metadata,
    )


def _ssgsea_result(
    *,
    effect_matrix: pd.DataFrame,
    kinase_to_sites: dict[str, list[str]],
    min_substrates: int = 2,
    ranking_direction: str = "descending",
    permutation_count: int = 0,
    random_seed: int | None = 0,
    adjust_p_values: bool = True,
):
    return SsgseaSubstrateEnrichmentActivityMethod(
        min_substrates=min_substrates,
        ranking_direction=ranking_direction,
        permutation_count=permutation_count,
        random_seed=random_seed,
        adjust_p_values=adjust_p_values,
    ).run(
        activity_input=ActivityInputMatrix.standardised_effect(
            _with_site_key_index(effect_matrix),
        ),
        kinase_substrate_membership=_membership(kinase_to_sites),
    )


def _sort_named_frame(frame: pd.DataFrame) -> pd.DataFrame:
    return frame.sort_index(axis=0).sort_index(axis=1)


def _assert_optional_named_frame_equal(
    left: pd.DataFrame | None,
    right: pd.DataFrame | None,
    *,
    index_labels: list[str] | None = None,
    column_labels: list[str] | None = None,
) -> None:
    if left is None:
        assert right is None
        return
    assert right is not None
    left_frame = left
    right_frame = right
    if index_labels is not None:
        left_frame = left_frame.loc[index_labels, :]
        right_frame = right_frame.loc[index_labels, :]
    if column_labels is not None:
        left_frame = left_frame.loc[:, column_labels]
        right_frame = right_frame.loc[:, column_labels]
    pdt.assert_frame_equal(
        _sort_named_frame(left_frame),
        _sort_named_frame(right_frame),
    )


def _assert_named_ssgsea_result_equal(
    left,
    right,
    *,
    index_labels: list[str] | None = None,
    column_labels: list[str] | None = None,
    compare_q_values: bool = True,
) -> None:
    _assert_optional_named_frame_equal(
        left.activity_matrix,
        right.activity_matrix,
        index_labels=index_labels,
        column_labels=column_labels,
    )
    _assert_optional_named_frame_equal(
        left.substrate_count_matrix,
        right.substrate_count_matrix,
        index_labels=index_labels,
        column_labels=column_labels,
    )
    _assert_optional_named_frame_equal(
        left.p_value_matrix,
        right.p_value_matrix,
        index_labels=index_labels,
        column_labels=column_labels,
    )
    if compare_q_values:
        _assert_optional_named_frame_equal(
            left.q_value_matrix,
            right.q_value_matrix,
            index_labels=index_labels,
            column_labels=column_labels,
        )

    left_stats = left.statistics_table
    right_stats = right.statistics_table
    assert left_stats is not None
    assert right_stats is not None
    stat_columns = [
        "enrichment_score",
        "p_value",
        "q_value",
        "significance_status",
        "n_substrates",
        "n_background_sites",
        "permutation_count",
        "random_seed",
        "computability_status",
    ]
    if not compare_q_values:
        stat_columns.remove("q_value")
    left_named = left_stats.set_index(["kinase", "profile_id"]).sort_index()
    right_named = right_stats.set_index(["kinase", "profile_id"]).sort_index()
    if index_labels is not None or column_labels is not None:
        kinase_labels = (
            set(index_labels)
            if index_labels is not None
            else {str(value[0]) for value in left_named.index}
        )
        profile_labels = (
            set(column_labels)
            if column_labels is not None
            else {str(value[1]) for value in left_named.index}
        )
        selected_index = [
            value
            for value in left_named.index
            if str(value[0]) in kinase_labels and str(value[1]) in profile_labels
        ]
        left_named = left_named.loc[selected_index, :]
        right_named = right_named.loc[selected_index, :]
    pdt.assert_frame_equal(
        left_named.loc[:, stat_columns],
        right_named.loc[:, stat_columns],
    )


def test_ssgsea_deterministic_rank_score_on_synthetic_data() -> None:
    effect_matrix = pd.DataFrame(
        {"c1": [4.0, 3.0, 2.0, 1.0]},
        index=["S1;S1;", "S2;S2;", "S3;S3;", "S4;S4;"],
    )

    result = _ssgsea_result(
        effect_matrix=effect_matrix,
        kinase_to_sites={
            "K_TOP": ["S1;S1;", "S2;S2;"],
            "K_BOTTOM": ["S3;S3;", "S4;S4;"],
        },
    )

    assert result.activity_method.activity_method_id == (
        "ssgsea_substrate_enrichment_activity_v1"
    )
    assert result.input_semantics.profile_axis is ActivityProfileAxis.EFFECT
    assert result.profile_metadata.profile_ids == ("c1",)
    assert result.activity_matrix.at["K_TOP", "c1"] == pytest.approx(0.5)
    assert result.activity_matrix.at["K_BOTTOM", "c1"] == pytest.approx(-0.5)
    assert result.substrate_count_matrix.at["K_TOP", "c1"] == 2
    assert result.substrate_count_matrix.at["K_BOTTOM", "c1"] == 2
    stats = result.statistics_table
    assert stats is not None
    assert "profile_id" in stats.columns
    assert "condition" not in stats.columns
    assert set(stats["profile_id"]) == {"c1"}
    top = stats.loc[stats["kinase"] == "K_TOP"].iloc[0]
    assert top["computability_status"] == SSGSEA_STATUS_COMPUTED
    assert top["significance_status"] == (
        SSGSEA_SIGNIFICANCE_STATUS_UNAVAILABLE_NO_PERMUTATIONS
    )
    assert top["enrichment_score"] == pytest.approx(0.5)
    assert pd.isna(top["p_value"])
    assert pd.isna(top["q_value"])
    assert top["ranking_direction"] == "descending"
    assert result.p_value_matrix is None
    assert result.q_value_matrix is None


def test_ssgsea_minimum_substrate_filtering_retains_diagnostic_pair() -> None:
    effect_matrix = pd.DataFrame(
        {"c1": [4.0, 3.0, 2.0, 1.0]},
        index=["S1;S1;", "S2;S2;", "S3;S3;", "S4;S4;"],
    )

    result = _ssgsea_result(
        effect_matrix=effect_matrix,
        kinase_to_sites={"K1": ["S1;S1;", "S2;S2;"]},
        min_substrates=3,
    )

    assert pd.isna(result.activity_matrix.at["K1", "c1"])
    assert result.substrate_count_matrix.at["K1", "c1"] == 2
    stats = result.statistics_table
    assert stats is not None
    assert int(stats.shape[0]) == 1
    assert stats.at[0, "computability_status"] == (
        SSGSEA_STATUS_INSUFFICIENT_SUBSTRATES
    )
    assert result.method_summary is not None
    assert result.method_summary.kinase_profile_pairs_insufficient_substrates == 1


def test_ssgsea_effect_statistics_use_neutral_profile_identifiers() -> None:
    effect_matrix = pd.DataFrame(
        {"effect_a": [4.0, 3.0, 2.0, 1.0]},
        index=["S1;S1;", "S2;S2;", "S3;S3;", "S4;S4;"],
    )

    result = _ssgsea_result(
        effect_matrix=effect_matrix,
        kinase_to_sites={"K1": ["S1;S1;", "S2;S2;"]},
    )

    assert result.input_semantics.profile_axis is ActivityProfileAxis.EFFECT
    assert result.profile_metadata.profile_ids == ("effect_a",)
    stats = result.statistics_table
    assert stats is not None
    assert "profile_id" in stats.columns
    assert "condition" not in stats.columns
    assert set(stats["profile_id"]) == {"effect_a"}


def test_ssgsea_p_value_adjustment_is_bh_per_profile_when_enabled() -> None:
    effect_matrix = pd.DataFrame(
        {"c1": [6.0, 5.0, 4.0, 3.0, 2.0, 1.0]},
        index=["S1;S1;", "S2;S2;", "S3;S3;", "S4;S4;", "S5;S5;", "S6;S6;"],
    )

    result = _ssgsea_result(
        effect_matrix=effect_matrix,
        kinase_to_sites={
            "K_TOP": ["S1;S1;", "S2;S2;"],
            "K_MID": ["S3;S3;", "S4;S4;"],
            "K_BOTTOM": ["S5;S5;", "S6;S6;"],
        },
        permutation_count=25,
        random_seed=17,
        adjust_p_values=True,
    )

    assert result.p_value_matrix is not None
    assert result.q_value_matrix is not None
    stats = result.statistics_table
    assert stats is not None
    assert set(stats["significance_status"]) == {SSGSEA_SIGNIFICANCE_STATUS_AVAILABLE}
    computed = stats.loc[
        stats["computability_status"] == SSGSEA_STATUS_COMPUTED,
        :,
    ].copy()
    expected_q = benjamini_hochberg_q_values(computed.loc[:, "p_value"].astype(float))
    pdt.assert_series_equal(
        computed.loc[:, "q_value"].astype(float),
        expected_q,
        check_names=False,
    )
    for kinase, q_value in expected_q.items():
        kinase_name = str(computed.at[kinase, "kinase"])
        assert result.q_value_matrix.at[kinase_name, "c1"] == pytest.approx(q_value)


def test_ssgsea_seed_reproducibility_for_permutation_p_values() -> None:
    effect_matrix = pd.DataFrame(
        {"c1": [6.0, 5.0, 4.0, 3.0, 2.0, 1.0]},
        index=["S1;S1;", "S2;S2;", "S3;S3;", "S4;S4;", "S5;S5;", "S6;S6;"],
    )
    kwargs = {
        "effect_matrix": effect_matrix,
        "kinase_to_sites": {
            "K_TOP": ["S1;S1;", "S2;S2;"],
            "K_BOTTOM": ["S5;S5;", "S6;S6;"],
        },
        "permutation_count": 30,
        "random_seed": 41,
    }

    first = _ssgsea_result(**kwargs)
    second = _ssgsea_result(**kwargs)

    assert first.p_value_matrix is not None
    assert second.p_value_matrix is not None
    assert first.q_value_matrix is not None
    assert second.q_value_matrix is not None
    pdt.assert_frame_equal(first.activity_matrix, second.activity_matrix)
    pdt.assert_frame_equal(first.p_value_matrix, second.p_value_matrix)
    pdt.assert_frame_equal(first.q_value_matrix, second.q_value_matrix)


def test_ssgsea_permutation_results_are_kinase_order_invariant() -> None:
    effect_matrix = pd.DataFrame(
        {
            "c1": [8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0],
            "c2": [1.0, 4.0, 2.0, 8.0, 3.0, 7.0, 5.0, 6.0],
        },
        index=[
            "S1;S1;",
            "S2;S2;",
            "S3;S3;",
            "S4;S4;",
            "S5;S5;",
            "S6;S6;",
            "S7;S7;",
            "S8;S8;",
        ],
    )
    kinase_to_sites = {
        "K_TOP": ["S1;S1;", "S2;S2;", "S4;S4;"],
        "K_MID": ["S3;S3;", "S5;S5;", "S7;S7;"],
        "K_BOTTOM": ["S6;S6;", "S7;S7;", "S8;S8;"],
    }

    first = _ssgsea_result(
        effect_matrix=effect_matrix,
        kinase_to_sites=kinase_to_sites,
        permutation_count=80,
        random_seed=37,
    )
    second = _ssgsea_result(
        effect_matrix=effect_matrix,
        kinase_to_sites=dict(reversed(list(kinase_to_sites.items()))),
        permutation_count=80,
        random_seed=37,
    )

    _assert_named_ssgsea_result_equal(first, second)


def test_ssgsea_permutation_results_are_profile_order_invariant() -> None:
    effect_matrix = pd.DataFrame(
        {
            "c1": [8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0],
            "c2": [1.0, 4.0, 2.0, 8.0, 3.0, 7.0, 5.0, 6.0],
        },
        index=[
            "S1;S1;",
            "S2;S2;",
            "S3;S3;",
            "S4;S4;",
            "S5;S5;",
            "S6;S6;",
            "S7;S7;",
            "S8;S8;",
        ],
    )
    kinase_to_sites = {
        "K_TOP": ["S1;S1;", "S2;S2;", "S4;S4;"],
        "K_MID": ["S3;S3;", "S5;S5;", "S7;S7;"],
        "K_BOTTOM": ["S6;S6;", "S7;S7;", "S8;S8;"],
    }

    first = _ssgsea_result(
        effect_matrix=effect_matrix,
        kinase_to_sites=kinase_to_sites,
        permutation_count=80,
        random_seed=37,
    )
    second = _ssgsea_result(
        effect_matrix=effect_matrix.loc[:, ["c2", "c1"]],
        kinase_to_sites=kinase_to_sites,
        permutation_count=80,
        random_seed=37,
    )

    _assert_named_ssgsea_result_equal(first, second)


def test_ssgsea_permutation_results_ignore_unrelated_kinase_insertion() -> None:
    effect_matrix = pd.DataFrame(
        {
            "c1": [8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0],
            "c2": [1.0, 4.0, 2.0, 8.0, 3.0, 7.0, 5.0, 6.0],
        },
        index=[
            "S1;S1;",
            "S2;S2;",
            "S3;S3;",
            "S4;S4;",
            "S5;S5;",
            "S6;S6;",
            "S7;S7;",
            "S8;S8;",
        ],
    )
    base_kinase_to_sites = {
        "K_TOP": ["S1;S1;", "S2;S2;", "S4;S4;"],
        "K_BOTTOM": ["S6;S6;", "S7;S7;", "S8;S8;"],
    }
    extended_kinase_to_sites = {
        "K_UNRELATED": ["S2;S2;", "S5;S5;", "S8;S8;"],
        **base_kinase_to_sites,
    }

    first = _ssgsea_result(
        effect_matrix=effect_matrix,
        kinase_to_sites=base_kinase_to_sites,
        permutation_count=80,
        random_seed=37,
        adjust_p_values=False,
    )
    second = _ssgsea_result(
        effect_matrix=effect_matrix,
        kinase_to_sites=extended_kinase_to_sites,
        permutation_count=80,
        random_seed=37,
        adjust_p_values=False,
    )

    _assert_named_ssgsea_result_equal(
        first,
        second,
        index_labels=["K_TOP", "K_BOTTOM"],
        compare_q_values=False,
    )
    assert first.statistics_table is not None
    assert set(first.statistics_table["significance_status"]) == {
        SSGSEA_SIGNIFICANCE_STATUS_P_VALUE_AVAILABLE_Q_VALUE_DISABLED
    }


def test_ssgsea_permutation_seed_derivation_changes_with_global_seed() -> None:
    first_seed = _derive_ssgsea_permutation_seed(
        random_seed=37,
        profile_id="c1",
        kinase_name="K_TOP",
    )
    second_seed = _derive_ssgsea_permutation_seed(
        random_seed=38,
        profile_id="c1",
        kinase_name="K_TOP",
    )

    assert first_seed != second_seed


def test_ssgsea_permutation_results_survive_input_serialization_round_trip() -> None:
    effect_matrix = _with_site_key_index(
        pd.DataFrame(
            {
                "c1": [8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0],
                "c2": [1.0, 4.0, 2.0, 8.0, 3.0, 7.0, 5.0, 6.0],
            },
            index=[
                "S1;S1;",
                "S2;S2;",
                "S3;S3;",
                "S4;S4;",
                "S5;S5;",
                "S6;S6;",
                "S7;S7;",
                "S8;S8;",
            ],
        )
    )
    membership = _membership(
        {
            "K_TOP": ["S1;S1;", "S2;S2;", "S4;S4;"],
            "K_MID": ["S3;S3;", "S5;S5;", "S7;S7;"],
            "K_BOTTOM": ["S6;S6;", "S7;S7;", "S8;S8;"],
        }
    )
    method = SsgseaSubstrateEnrichmentActivityMethod(
        min_substrates=2,
        permutation_count=80,
        random_seed=37,
    )

    first = method.run(
        activity_input=ActivityInputMatrix.standardised_effect(effect_matrix),
        kinase_substrate_membership=membership,
    )
    restored_effect_matrix = pd.read_json(
        StringIO(effect_matrix.to_json(orient="split")),
        orient="split",
    )
    restored_membership = pd.read_json(
        StringIO(membership.to_json(orient="records")),
        orient="records",
    )
    second = method.run(
        activity_input=ActivityInputMatrix.standardised_effect(restored_effect_matrix),
        kinase_substrate_membership=restored_membership,
    )

    _assert_named_ssgsea_result_equal(first, second)


def test_ssgsea_result_populates_contract_and_policy_provenance() -> None:
    effect_matrix = pd.DataFrame(
        {"c1": [4.0, 3.0, 2.0, 1.0]},
        index=["S1;S1;", "S2;S2;", "S3;S3;", "S4;S4;"],
    )

    result = _ssgsea_result(
        effect_matrix=effect_matrix,
        kinase_to_sites={"K1": ["S1;S1;", "S2;S2;"]},
        permutation_count=5,
        random_seed=3,
    )

    pdt.assert_frame_equal(result.to_dataframe(), result.activity_matrix)
    assert result.activity_substrate_counts is not None
    pdt.assert_frame_equal(
        result.substrate_count_matrix,
        result.activity_substrate_counts,
    )
    assert isinstance(
        result.method_diagnostics,
        SsgseaSubstrateEnrichmentActivityDiagnostics,
    )
    assert result.policy_provenance
    policy = result.policy_provenance[0]
    assert policy.id == ScientificPolicyId.SSGSEA_SUBSTRATE_ENRICHMENT_ACTIVITY
    payload = policy.to_payload()
    assert payload["id"] == "ssgsea_substrate_enrichment_activity_v1"
    assert payload["version"] == SSGSEA_SUBSTRATE_ENRICHMENT_ACTIVITY_POLICY_VERSION
    parameters = payload["parameters"]
    assert isinstance(parameters, dict)
    assert parameters["method_version"] == (
        SSGSEA_SUBSTRATE_ENRICHMENT_ACTIVITY_POLICY_VERSION
    )
    assert parameters["random_seed"] == 3
    assert parameters["permutation_rng_seed_policy"] == (
        SSGSEA_PERMUTATION_RNG_SEED_POLICY
    )
    assert parameters["permutation_rng_seed_policy_version"] == (
        SSGSEA_PERMUTATION_RNG_SEED_POLICY_VERSION
    )
    assert parameters["q_value_method"] == SSGSEA_Q_VALUE_METHOD_BENJAMINI_HOCHBERG


def test_ssgsea_activity_input_requires_contrast_or_effect_semantics() -> None:
    effect_matrix = _with_site_key_index(
        pd.DataFrame(
            {"contrast_a": [4.0, 3.0, 2.0, 1.0]},
            index=["S1;S1;", "S2;S2;", "S3;S3;", "S4;S4;"],
        )
    )
    membership = _membership({"K1": ["S1;S1;", "S2;S2;"]})
    method = SsgseaSubstrateEnrichmentActivityMethod(
        min_substrates=2,
        permutation_count=0,
    )

    with pytest.raises(
        WorkflowBoundaryError,
        match="requires explicit contrast/effect input",
    ):
        method.run(
            activity_input=ActivityInputMatrix.sample_level_abundance(effect_matrix),
            kinase_substrate_membership=membership,
        )

    result = method.run(
        activity_input=ActivityInputMatrix.contrast_log_fold_change(effect_matrix),
        kinase_substrate_membership=membership,
    )

    assert result.input_semantics.profile_axis is ActivityProfileAxis.CONTRAST
    assert (
        result.input_semantics.quantitative_semantics
        is ActivityQuantitativeSemantics.CONTRAST_LOG_FOLD_CHANGE
    )
    assert result.profile_metadata.contrast_ids == ("contrast_a",)
    assert result.activity_matrix.columns.name == "profile_id"
    stats = result.statistics_table
    assert stats is not None
    assert "profile_id" in stats.columns
    assert "condition" not in stats.columns
    assert set(stats["profile_id"]) == {"contrast_a"}
    assert all(
        "condition" not in value for value in result.count_field_semantics.values()
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

    assert result.activity_matrix.at["K1", "c1"] == pytest.approx(-1.0954451150103324)
    stats = result.statistics_table
    assert stats is not None
    row = stats.iloc[0]
    assert row["computability_status"] == KSEA_STATUS_COMPUTED
    assert row["p_value"] == pytest.approx(0.27332167829229814)
    assert row["n_substrates"] == 2
    assert row["n_background_sites"] == 4


def test_activity_scores_compatibility_alias_matches_activity_matrix() -> None:
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
    assert result.activity_matrix.at["K1", "c1"] == pytest.approx(-1.0954451150103324)
    with pytest.warns(
        DeprecationWarning,
        match="KinaseActivityResult.activity_scores.*activity_matrix",
    ):
        activity_scores = result.activity_scores
    pdt.assert_frame_equal(activity_scores, result.activity_matrix)


def test_ksea_result_populates_extensible_activity_contract() -> None:
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

    assert result.input_semantics.profile_axis is ActivityProfileAxis.EFFECT
    assert result.profile_metadata.profile_ids == ("c1",)
    pdt.assert_frame_equal(result.to_dataframe(), result.activity_matrix)
    assert result.p_value_matrix is not None
    assert result.q_value_matrix is not None
    assert result.confidence_interval_low is None
    assert result.confidence_interval_high is None
    assert result.p_value_matrix.at["K1", "c1"] == pytest.approx(0.27332167829229814)
    assert result.q_value_matrix.at["K1", "c1"] == pytest.approx(0.27332167829229814)
    assert result.activity_substrate_counts is not None
    pdt.assert_frame_equal(
        result.substrate_count_matrix, result.activity_substrate_counts
    )
    assert isinstance(result.method_diagnostics, KseaZScoreActivityDiagnostics)
    assert result.method_diagnostics.statistics_table is not None
    stats = result.statistics_table
    assert stats is not None
    assert "profile_id" in stats.columns
    assert "condition" not in stats.columns
    assert set(stats["profile_id"]) == {"c1"}
    assert all(
        "condition" not in value for value in result.count_field_semantics.values()
    )
    assert result.policy_provenance
    policy = result.policy_provenance[0]
    assert policy.id == ScientificPolicyId.KSEA_ZSCORE_ACTIVITY
    assert policy.to_payload()["id"] == "ksea_zscore_activity_v1"


def test_ksea_sample_statistics_use_profile_ids_and_adjust_p_values_per_profile() -> (
    None
):
    pred_mat = _with_site_key_index(
        pd.DataFrame(
            {
                "K1": [0.9, 0.9, 0.1, 0.1],
                "K2": [0.1, 0.1, 0.9, 0.9],
                "K3": [0.9, 0.1, 0.1, 0.9],
            },
            index=["S1;S1;", "S2;S2;", "S3;S3;", "S4;S4;"],
        )
    )
    phospho = _with_site_key_index(
        pd.DataFrame(
            {
                "sample_a": [1.0, 2.0, 3.0, 4.0],
                "sample_b": [4.0, 1.0, 2.0, 3.0],
            },
            index=["S1;S1;", "S2;S2;", "S3;S3;", "S4;S4;"],
        )
    )
    activity_input = ActivityInputMatrix.sample_level_abundance(
        phospho,
        _assume_owned=True,
    )

    result = KseaZScoreActivityMethod(
        evidence_threshold=0.5,
        min_substrates=2,
        adjust_p_values=True,
    ).run(
        KinaseActivityInputs(
            pred_mat=pred_mat,
            phospho_matrix=phospho,
            threshold=0.5,
            min_substrates=2,
            top_n_substrates=1,
            overlap_summary=PredMatOverlapSummary(
                overlap_count=int(pred_mat.index.intersection(phospho.index).size),
                pred_mat_rows=int(pred_mat.index.size),
                phospho_rows=int(phospho.index.size),
            ),
            activity_input=activity_input,
            membership_selection=_eligible_fixed_membership_selection(
                pred_mat,
                threshold=0.5,
            ),
        )
    )

    assert result.input_semantics.profile_axis is ActivityProfileAxis.SAMPLE
    assert result.profile_metadata.sample_ids == ("sample_a", "sample_b")
    stats = result.statistics_table
    assert stats is not None
    assert "profile_id" in stats.columns
    assert "condition" not in stats.columns
    assert set(stats["profile_id"]) == set(result.profile_metadata.sample_ids)
    for profile_id in result.profile_metadata.sample_ids:
        computed = stats.loc[
            (stats["profile_id"] == profile_id)
            & (stats["computability_status"] == KSEA_STATUS_COMPUTED),
            :,
        ]
        expected_q_values = benjamini_hochberg_q_values(
            computed.loc[:, "p_value"].astype(float)
        )
        pdt.assert_series_equal(
            computed.loc[:, "q_value"].astype(float),
            expected_q_values,
            check_names=False,
        )


def test_ksea_computes_each_kinase_profile_pair_independently() -> None:
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

    assert result.activity_matrix.at["K1", "c1"] == pytest.approx(-0.7071067811865476)
    assert result.activity_matrix.at["K2", "c1"] == pytest.approx(0.7071067811865476)
    assert result.activity_matrix.at["K1", "c2"] == pytest.approx(0.7071067811865476)
    assert result.activity_matrix.at["K2", "c2"] == pytest.approx(-0.7071067811865476)


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

    assert pd.isna(result.activity_matrix.at["K1", "c1"])
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

    expected_sites = set(_site_key_index(["S2;S2;", "S3;S3;"]).astype(str))
    assert weighted.thresholded_substrate_counts.to_dict() == {"K1": 2}
    assert weighted.target_counts.to_dict() == {"K1": 2}
    assert ksea.thresholded_substrate_counts.to_dict() == {"K1": 2}
    assert ksea.target_counts.to_dict() == {"K1": 2}

    weighted_sites = set(weighted.target_table.loc[:, "site_id"].astype(str))
    ksea_sites = set(ksea.target_table.loc[:, "site_id"].astype(str))
    assert weighted_sites == expected_sites
    assert ksea_sites == expected_sites
    assert weighted_sites == ksea_sites
    assert _site_key_index(["S1;S1;"])[0] not in weighted_sites

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


def test_ksea_excludes_non_finite_phosphosite_values_per_profile() -> None:
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
    c1 = stats.loc[stats["profile_id"] == "c1"].iloc[0]
    c2 = stats.loc[stats["profile_id"] == "c2"].iloc[0]
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
        "profile-specific finite-value filtering"
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


def test_ksea_fixed_external_membership_policy_allows_ordinary_p_q_values() -> None:
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

    assert result.membership_selection is not None
    assert result.membership_selection.source_category == (
        ACTIVITY_MEMBERSHIP_SOURCE_FIXED_EXTERNAL_REFERENCE
    )
    assert result.membership_selection.consumed_tested_matrix is False
    assert result.membership_selection.inferential_eligible is True
    assert result.p_value_matrix is not None
    assert result.q_value_matrix is not None
    stats = result.statistics_table
    assert stats is not None
    assert stats.at[0, "p_value"] == pytest.approx(0.27332167829229814)
    assert stats.at[0, "q_value"] == pytest.approx(0.27332167829229814)


def test_ksea_missing_membership_provenance_downgrades_to_descriptive_z_scores() -> (
    None
):
    pred_mat = pd.DataFrame(
        {"K1": [0.9, 0.8, 0.1, 0.2]},
        index=["S1;S1;", "S2;S2;", "S3;S3;", "S4;S4;"],
    )
    phospho = pd.DataFrame({"c1": [1.0, 2.0, 3.0, 4.0]}, index=pred_mat.index.copy())
    keyed_pred_mat = _with_site_key_index(pred_mat)
    effect_matrix = _with_site_key_index(phospho)

    result = KseaZScoreActivityMethod(
        evidence_threshold=0.5,
        min_substrates=2,
        adjust_p_values=True,
    ).run(
        _inputs(
            pred_mat=keyed_pred_mat,
            phospho_matrix=effect_matrix,
            threshold=0.5,
            min_substrates=2,
            top_n_substrates=1,
            activity_input=ActivityInputMatrix.standardised_effect(
                effect_matrix,
                _assume_owned=True,
            ),
        )
    )

    assert result.activity_matrix.at["K1", "c1"] == pytest.approx(-1.0954451150103324)
    assert result.p_value_matrix is None
    assert result.q_value_matrix is None
    assert result.membership_selection is not None
    assert result.membership_selection.inferential_eligible is False
    assert (
        result.membership_selection.inferential_eligibility_reason
        == KSEA_MEMBERSHIP_MISSING_PROVENANCE_REASON
    )
    stats = result.statistics_table
    assert stats is not None
    assert pd.isna(stats.at[0, "p_value"])
    assert pd.isna(stats.at[0, "q_value"])
    assert bool(stats.at[0, "inferential_eligible"]) is False
    assert (
        stats.at[0, "inferential_reason"] == KSEA_MEMBERSHIP_MISSING_PROVENANCE_REASON
    )


def test_ksea_q_values_are_benjamini_hochberg_adjusted_per_profile() -> None:
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
    c1_rows = stats.loc[stats["profile_id"] == "c1"].sort_values("kinase")
    q_values = c1_rows.loc[:, "q_value"].to_numpy(dtype=float)
    assert q_values[0] == pytest.approx(0.4099825174384472)
    assert q_values[1] == pytest.approx(0.4099825174384472)
    assert q_values[2] == pytest.approx(1.0)


def test_adaptive_same_matrix_membership_simulation_inflates_old_normal_p_values() -> (
    None
):
    rng = np.random.default_rng(17)
    replicates = 250
    site_count = 160
    selected_count = 12
    nominal_alpha = 0.05
    false_positive_count = 0

    for _ in range(replicates):
        values = rng.normal(size=site_count)
        selected_positions = np.argpartition(-values, selected_count - 1)[
            :selected_count
        ]
        mean_background = float(values.mean())
        sd_background = float(values.std(ddof=1))
        mean_selected = float(values[selected_positions].mean())
        z_score = (
            (mean_selected - mean_background)
            * np.sqrt(float(selected_count))
            / sd_background
        )
        p_value = two_sided_normal_p_value(float(z_score))
        if p_value < nominal_alpha:
            false_positive_count += 1

    observed_rate = false_positive_count / replicates
    assert observed_rate > 0.50


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
        stats.pivot(index="kinase", columns="profile_id", values="n_substrates")
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

    assert result.activity_matrix.at["PRKACA", "phospho_corrected_1"] == pytest.approx(
        6.0625
    )
    assert result.activity_matrix.at["PRKACA", "phospho_corrected_2"] == pytest.approx(
        (20 * 0.9 + 6 * 0.8) / (0.9 + 0.8)
    )


def test_weighted_activity_compatibility_alias_matches_activity_matrix() -> None:
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
    assert result.activity_matrix.at["PRKACA", "phospho_corrected_1"] == pytest.approx(
        6.0625
    )
    assert result.activity_matrix.at["PRKACA", "phospho_corrected_2"] == pytest.approx(
        (20 * 0.9 + 6 * 0.8) / (0.9 + 0.8)
    )
    with pytest.warns(
        DeprecationWarning,
        match="KinaseActivityResult.weighted_activity.*activity_matrix",
    ):
        weighted_activity = result.weighted_activity
    pdt.assert_frame_equal(weighted_activity, result.activity_matrix)


def test_weighted_result_populates_extensible_activity_contract() -> None:
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

    pdt.assert_frame_equal(result.to_dataframe(), result.activity_matrix)
    assert result.p_value_matrix is None
    assert result.q_value_matrix is None
    assert result.confidence_interval_low is None
    assert result.confidence_interval_high is None
    assert result.substrate_count_matrix.at["PRKACA", "phospho_corrected_1"] == 2
    assert result.substrate_count_matrix.at["PRKACA", "phospho_corrected_2"] == 2
    assert result.activity_substrate_counts is None
    assert isinstance(result.method_diagnostics, WeightedSubstrateActivityDiagnostics)
    assert result.policy_provenance
    policy = result.policy_provenance[0]
    assert policy.id == ScientificPolicyId.SIMPLIFIED_WEIGHTED_SUBSTRATE_ACTIVITY
    assert policy.to_payload()["id"] == "simplified_weighted_substrate_activity_v1"


def test_activity_result_contract_handles_empty_optional_diagnostics_cleanly() -> None:
    activity_matrix = pd.DataFrame(
        {"c1": [1.0]},
        index=pd.Index(["K1"], name="kinase"),
        dtype=float,
    )
    substrate_count_matrix = pd.DataFrame(
        {"c1": [2]},
        index=pd.Index(["K1"], name="kinase"),
        dtype="int64",
    )
    activity_input = ActivityInputMatrix.sample_level_abundance(activity_matrix)

    result = KinaseActivityResult(
        activity_matrix=activity_matrix,
        substrate_count_matrix=substrate_count_matrix,
        input_semantics=activity_input.semantics,
        profile_metadata=activity_input.profile_metadata,
    )

    pdt.assert_frame_equal(result.to_dataframe(), activity_matrix)
    pdt.assert_frame_equal(result.activity_matrix, activity_matrix)
    pdt.assert_frame_equal(result.substrate_count_matrix, substrate_count_matrix)
    assert result.p_value_matrix is None
    assert result.q_value_matrix is None
    assert result.confidence_interval_low is None
    assert result.confidence_interval_high is None
    assert result.thresholded_substrate_counts.empty
    assert result.target_counts.empty
    assert result.target_table.empty
    assert result.policy_provenance == ()


def test_activity_statistics_table_rejects_missing_profile_id() -> None:
    activity_matrix = pd.DataFrame(
        {"sample_a": [1.0]},
        index=pd.Index(["site_a"], name="site_id"),
        dtype=float,
    )
    activity_input = ActivityInputMatrix.sample_level_abundance(activity_matrix)
    statistics_table = _statistics_table(["sample_a"]).drop(columns=["profile_id"])

    with pytest.raises(
        PhosPyValidationError,
        match="missing required columns: profile_id",
    ):
        _activity_result_from_statistics_table(
            statistics_table,
            activity_input=activity_input,
        )


def test_activity_statistics_table_rejects_blank_profile_id() -> None:
    activity_matrix = pd.DataFrame(
        {"sample_a": [1.0]},
        index=pd.Index(["site_a"], name="site_id"),
        dtype=float,
    )
    activity_input = ActivityInputMatrix.sample_level_abundance(activity_matrix)

    with pytest.raises(
        PhosPyValidationError,
        match="profile_id must contain stripped non-empty string values",
    ):
        _activity_result_from_statistics_table(
            _statistics_table([""]),
            activity_input=activity_input,
        )


def test_activity_statistics_table_rejects_unknown_profile_id() -> None:
    activity_matrix = pd.DataFrame(
        {"sample_a": [1.0]},
        index=pd.Index(["site_a"], name="site_id"),
        dtype=float,
    )
    activity_input = ActivityInputMatrix.sample_level_abundance(activity_matrix)

    with pytest.raises(
        PhosPyValidationError,
        match="unknown_profile_ids=\\('unknown_profile',\\)",
    ):
        _activity_result_from_statistics_table(
            _statistics_table(["unknown_profile"]),
            activity_input=activity_input,
        )


def test_activity_statistics_table_rejects_condition_claim_on_sample_profiles() -> None:
    activity_matrix = pd.DataFrame(
        {"sample_a": [1.0]},
        index=pd.Index(["site_a"], name="site_id"),
        dtype=float,
    )
    activity_input = ActivityInputMatrix.sample_level_abundance(activity_matrix)

    with pytest.raises(
        PhosPyValidationError,
        match="condition is reserved for condition-summary activity results",
    ):
        _activity_result_from_statistics_table(
            _statistics_table(["sample_a"], include_condition=True),
            activity_input=activity_input,
        )


def test_activity_method_summary_accepts_legacy_payload_but_serializes_profile_fields() -> (
    None
):
    legacy_payload = {
        "kinases_evaluated": 1,
        "kinase_condition_pairs_evaluated": 2,
        "kinase_condition_pairs_computed": 3,
        "kinase_condition_pairs_insufficient_substrates": 4,
        "kinase_condition_pairs_invalid_background_variance": 5,
        "kinase_condition_pairs_no_finite_background_values": 6,
        "kinase_condition_pairs_no_finite_substrate_values": 7,
    }

    summary = ActivityMethodSummary.from_payload(legacy_payload)

    assert summary.kinase_profile_pairs_evaluated == 2
    assert summary.kinase_profile_pairs_computed == 3
    payload = summary.to_payload()
    assert "kinase_profile_pairs_evaluated" in payload
    assert "kinase_condition_pairs_evaluated" not in payload
    with pytest.warns(
        DeprecationWarning,
        match="kinase_condition_pairs_evaluated.*kinase_profile_pairs_evaluated",
    ):
        assert summary.kinase_condition_pairs_evaluated == 2


def test_activity_method_summary_rejects_conflicting_payload_aliases() -> None:
    with pytest.raises(ValueError, match="conflicts with legacy alias"):
        ActivityMethodSummary.from_payload(
            {
                "kinases_evaluated": 1,
                "kinase_profile_pairs_evaluated": 2,
                "kinase_condition_pairs_evaluated": 3,
            }
        )


def test_legacy_condition_statistics_table_adapter_is_deprecated_and_defensive() -> (
    None
):
    activity_matrix = pd.DataFrame(
        {"sample_a": [1.0]},
        index=pd.Index(["site_a"], name="site_id"),
        dtype=float,
    )
    activity_input = ActivityInputMatrix.sample_level_abundance(activity_matrix)
    result = _activity_result_from_statistics_table(
        _statistics_table(["sample_a"]),
        activity_input=activity_input,
    )

    with pytest.warns(
        DeprecationWarning,
        match="does not establish a biological condition contract",
    ):
        legacy = result.legacy_condition_statistics_table_dataframe()

    assert legacy is not None
    assert legacy["condition"].tolist() == legacy["profile_id"].tolist()
    legacy.loc[0, "condition"] = "mutated"
    stats = result.statistics_table
    assert stats is not None
    assert "condition" not in stats.columns


def test_activity_input_matrix_owns_caller_frame_and_exports_snapshots() -> None:
    frame = pd.DataFrame(
        {"sample_a": [1.0]},
        index=pd.Index(["site_a"], name="site_id"),
        dtype=float,
    )

    activity_input = ActivityInputMatrix.sample_level_abundance(frame)
    frame.at["site_a", "sample_a"] = 99.0

    assert activity_input.frame.at["site_a", "sample_a"] == pytest.approx(1.0)
    exported = activity_input.matrix
    exported.at["site_a", "sample_a"] = 42.0
    assert activity_input.frame.at["site_a", "sample_a"] == pytest.approx(1.0)


def test_condition_summary_activity_input_marks_output_profiles_as_conditions() -> None:
    pred_mat = pd.DataFrame(
        {"K1": [1.0, 1.0]},
        index=["S1;S1;", "S2;S2;"],
    )
    condition_matrix = _with_site_key_index(
        pd.DataFrame(
            {
                "treated_mean": [2.0, 4.0],
                "control_mean": [1.0, 3.0],
            },
            index=pred_mat.index.copy(),
        )
    )
    aggregation_metadata = ActivityAggregationMetadata(
        aggregation_method="mean",
        records=(
            ActivityAggregationRecord(
                profile_id="treated_mean",
                source_profile_ids=("treated_rep1", "treated_rep2"),
            ),
            ActivityAggregationRecord(
                profile_id="control_mean",
                source_profile_ids=("control_rep1", "control_rep2"),
            ),
        ),
    )
    activity_input = ActivityInputMatrix.condition_summary_abundance(
        condition_matrix,
        aggregation_metadata=aggregation_metadata,
    )

    result = compute_activity_from_inputs(
        _inputs(
            pred_mat=pred_mat,
            phospho_matrix=condition_matrix,
            threshold=0.6,
            min_substrates=2,
            top_n_substrates=2,
            activity_input=activity_input,
        )
    )

    assert result.input_semantics.profile_axis is ActivityProfileAxis.CONDITION_SUMMARY
    assert (
        result.input_semantics.quantitative_semantics
        is ActivityQuantitativeSemantics.CONDITION_SUMMARY_ABUNDANCE
    )
    assert result.profile_metadata.aggregation_metadata == aggregation_metadata
    assert result.activity_matrix.columns.name == "condition"
    assert result.substrate_count_matrix.columns.name == "condition"
    assert result.thresholded_substrate_mean_activity.columns.name == "condition"
    assert result.activity_matrix.at["K1", "treated_mean"] == pytest.approx(3.0)
    assert result.activity_matrix.at["K1", "control_mean"] == pytest.approx(2.0)


def test_condition_summary_statistics_table_may_include_matching_condition_alias() -> (
    None
):
    condition_matrix = pd.DataFrame(
        {
            "treated_mean": [2.0, 4.0],
            "control_mean": [1.0, 3.0],
        },
        index=pd.Index(["site_a", "site_b"], name="site_id"),
        dtype=float,
    )
    aggregation_metadata = ActivityAggregationMetadata(
        aggregation_method="mean",
        records=(
            ActivityAggregationRecord(
                profile_id="treated_mean",
                source_profile_ids=("treated_rep1", "treated_rep2"),
            ),
            ActivityAggregationRecord(
                profile_id="control_mean",
                source_profile_ids=("control_rep1", "control_rep2"),
            ),
        ),
    )
    activity_input = ActivityInputMatrix.condition_summary_abundance(
        condition_matrix,
        aggregation_metadata=aggregation_metadata,
    )

    result = _activity_result_from_statistics_table(
        _statistics_table(
            ["treated_mean", "control_mean"],
            include_condition=True,
        ),
        activity_input=activity_input,
    )

    assert result.profile_metadata.aggregation_metadata == aggregation_metadata
    stats = result.statistics_table
    assert stats is not None
    assert "profile_id" in stats.columns
    assert "condition" in stats.columns
    assert stats["condition"].tolist() == stats["profile_id"].tolist()
    assert (
        "condition-specific" in result.count_field_semantics["substrate_count_matrix"]
    )


def test_condition_summary_statistics_table_rejects_mismatched_condition_alias() -> (
    None
):
    condition_matrix = pd.DataFrame(
        {"treated_mean": [2.0]},
        index=pd.Index(["site_a"], name="site_id"),
        dtype=float,
    )
    aggregation_metadata = ActivityAggregationMetadata(
        aggregation_method="mean",
        records=(
            ActivityAggregationRecord(
                profile_id="treated_mean",
                source_profile_ids=("treated_rep1", "treated_rep2"),
            ),
        ),
    )
    activity_input = ActivityInputMatrix.condition_summary_abundance(
        condition_matrix,
        aggregation_metadata=aggregation_metadata,
    )

    with pytest.raises(
        PhosPyValidationError,
        match="condition must equal profile_id",
    ):
        _activity_result_from_statistics_table(
            _statistics_table(
                ["treated_mean"],
                include_condition=True,
                condition_values=["treated_label_from_elsewhere"],
            ),
            activity_input=activity_input,
        )


def test_condition_summary_profile_metadata_requires_aggregation_metadata() -> None:
    with pytest.raises(
        WorkflowBoundaryError,
        match="condition-summary activity input requires explicit "
        "ActivityAggregationMetadata",
    ):
        ActivityProfileMetadata(
            axis=ActivityProfileAxis.CONDITION_SUMMARY,
            profile_ids=("treated_mean",),
            condition_ids=("treated_mean",),
        )


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

    assert result.activity_matrix.at["MAP2K6", "sample_a"] == pytest.approx(5.5)


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
        "simplified weighted substrate activity-like score"
    )
    assert result.activity_method.is_ksea is False
    assert result.activity_method.is_phosr_kinase_activity_equivalent is False
    assert result.activity_matrix.at["PRKACA", "phospho_corrected_1"] == pytest.approx(
        6.0625
    )
    assert result.activity_matrix.at["PRKACA", "phospho_corrected_2"] == pytest.approx(
        (20 * 0.9 + 6 * 0.8) / (0.9 + 0.8)
    )
