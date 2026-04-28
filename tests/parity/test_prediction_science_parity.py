from __future__ import annotations

import numpy as np
import pandas as pd
import pandas.testing as pdt
import pytest

from phospy.prediction.candidates import build_candidate_substrate_list
from phospy.prediction.motif_scoring import (
    SEQUENCE_SEMANTICS_CENTRED_SEQUENCE,
    score_phosphosite_motifs,
)
from phospy.prediction.scoring import fuse_profile_and_motif_scores_by_rank_weight
from tests.support.parity_reporting import (
    format_percent,
    format_shape,
    record_parity_metrics,
)
from tests.support.rewrite_fixture_data import (
    load_fragile_support_candidate_substrates,
    load_fragile_support_motif_frequency_matrices,
    load_fragile_support_motif_scores,
    load_fragile_support_motif_scores_full,
    load_fragile_support_motif_site_sequences_full,
    load_fragile_support_motif_sizes,
    load_fragile_support_profile_scores,
    load_fragile_support_profile_sizes,
    load_fragile_support_rank_weighted_fusion_scores,
    load_fragile_support_score_fusion_weights,
)

pytestmark = pytest.mark.parity
MOTIF_NUMERIC_RTOL = 0.0
# Motif scoring is deterministic; this tolerance only absorbs binary float noise.
MOTIF_NUMERIC_ATOL = 1e-12


def _mean_column_correlation(
    observed: pd.DataFrame,
    expected: pd.DataFrame,
    *,
    method: str,
) -> float:
    common_columns = observed.columns.intersection(expected.columns)
    common_index = observed.index.intersection(expected.index)
    correlations: list[float] = []
    for column in common_columns:
        correlation = observed.loc[common_index, column].corr(
            expected.loc[common_index, column],
            method=method,
        )
        if pd.notna(correlation):
            correlations.append(float(correlation))
    if not correlations:
        return float("nan")
    return float(pd.Series(correlations).mean())


def _assert_numeric_table_parity(
    observed: pd.DataFrame,
    expected: pd.DataFrame,
    *,
    rtol: float,
    atol: float,
) -> None:
    assert observed.index.tolist() == expected.index.tolist()
    assert observed.columns.tolist() == expected.columns.tolist()
    pdt.assert_frame_equal(
        observed,
        expected,
        check_dtype=False,
        check_exact=False,
        check_names=False,
        rtol=rtol,
        atol=atol,
    )


def test_motif_scoring_matches_fragile_support_reference_points(
    request: pytest.FixtureRequest,
) -> None:
    expected_full = load_fragile_support_motif_scores_full()
    expected_subset = load_fragile_support_motif_scores()
    expected_sizes = load_fragile_support_motif_sizes().sort_index()

    observed_result = score_phosphosite_motifs(
        site_sequences=load_fragile_support_motif_site_sequences_full(),
        motif_frequency_matrices=load_fragile_support_motif_frequency_matrices(),
        motif_sizes=expected_sizes,
        site_index=expected_full.index,
        min_motif_size=1,
        sequence_semantics=SEQUENCE_SEMANTICS_CENTRED_SEQUENCE,
    )
    observed_full = observed_result.motif_scores
    valid_site_ids = [
        row.site_id
        for row in observed_result.sequence_validation.rows
        if row.status == "valid"
    ]
    excluded_site_ids = observed_result.sequence_validation.excluded_site_ids
    observed_full_valid = observed_full.loc[valid_site_ids, :]
    expected_full_valid = expected_full.loc[valid_site_ids, :]
    subset_valid_index = expected_subset.index.intersection(valid_site_ids)
    observed_subset = observed_full.loc[subset_valid_index, expected_subset.columns]
    expected_subset_valid = expected_subset.loc[subset_valid_index, :]

    pdt.assert_series_equal(
        observed_result.motif_sizes.sort_index(),
        expected_sizes,
        check_dtype=False,
        check_exact=True,
    )
    full_spearman = _mean_column_correlation(
        observed_full_valid,
        expected_full_valid,
        method="spearman",
    )
    subset_spearman = _mean_column_correlation(
        observed_subset,
        expected_subset_valid,
        method="spearman",
    )
    assert full_spearman > 0.998
    assert subset_spearman > 0.998
    finite_values = observed_full.to_numpy(dtype=float)
    finite_values = finite_values[np.isfinite(finite_values)]
    assert (finite_values >= 0.0).all()
    assert (finite_values <= 1.0).all()
    if excluded_site_ids:
        assert observed_full.loc[list(excluded_site_ids), :].isna().all().all()

    full_delta = (observed_full_valid - expected_full_valid).abs().to_numpy()
    subset_delta = (observed_subset - expected_subset_valid).abs().to_numpy()
    record_parity_metrics(
        request.config,
        family="prediction_science",
        metrics=[
            ("motif score full-table parity", "pass"),
            ("motif score full-table shape", format_shape(*observed_full.shape)),
            ("motif score subset-table shape", format_shape(*observed_subset.shape)),
            ("motif score full max abs diff", float(full_delta.max())),
            ("motif score subset max abs diff", float(subset_delta.max())),
            (
                "motif score full mean Spearman column correlation",
                format_percent(full_spearman),
            ),
            (
                "motif score subset mean Spearman column correlation",
                format_percent(subset_spearman),
            ),
            ("motif sizes compared", int(expected_sizes.shape[0])),
            (
                "motif score excluded invalid sequences",
                int(len(excluded_site_ids)),
            ),
        ],
    )


def test_rank_weighted_fusion_scoring_matches_fragile_support_reference_tables(
    request: pytest.FixtureRequest,
) -> None:
    rank_weighted_fusion_scores, score_fusion_weights = (
        fuse_profile_and_motif_scores_by_rank_weight(
            motif_scores=load_fragile_support_motif_scores(),
            profile_scores=load_fragile_support_profile_scores(),
            motif_sizes=load_fragile_support_motif_sizes(),
            profile_sizes=load_fragile_support_profile_sizes(),
            allow_profile_only_fallback=False,
        )
    )
    expected_scores = load_fragile_support_rank_weighted_fusion_scores()
    expected_weights = load_fragile_support_score_fusion_weights()

    pdt.assert_frame_equal(
        rank_weighted_fusion_scores.sort_index(axis=0).sort_index(axis=1),
        expected_scores.sort_index(axis=0).sort_index(axis=1),
    )
    pdt.assert_frame_equal(
        score_fusion_weights.sort_index(),
        expected_weights.sort_index(),
    )
    aligned_observed = rank_weighted_fusion_scores.sort_index().sort_index(axis=1)
    aligned_expected = expected_scores.sort_index().sort_index(axis=1)
    absolute_delta = (aligned_observed - aligned_expected).abs()
    record_parity_metrics(
        request.config,
        family="prediction_science",
        metrics=[
            (
                "rank-weighted fusion score table shape",
                format_shape(*rank_weighted_fusion_scores.shape),
            ),
            (
                "rank-weighted fusion score mean abs diff",
                float(absolute_delta.to_numpy().mean()),
            ),
            (
                "rank-weighted fusion score max abs diff",
                float(absolute_delta.to_numpy().max()),
            ),
            (
                "rank-weighted fusion score mean Spearman column correlation",
                format_percent(
                    _mean_column_correlation(
                        aligned_observed,
                        aligned_expected,
                        method="spearman",
                    )
                ),
            ),
        ],
    )


def test_candidate_selection_matches_fragile_support_reference_table(
    request: pytest.FixtureRequest,
) -> None:
    observed = build_candidate_substrate_list(
        scores=load_fragile_support_rank_weighted_fusion_scores(),
        top=50,
        score_threshold=0.8,
        inclusion=20,
    )
    observed_rows = [
        {"kinase": kinase, "site_id": site_id}
        for kinase, site_ids in observed.items()
        for site_id in site_ids
    ]
    observed_frame = pd.DataFrame(observed_rows).sort_values(
        ["kinase", "site_id"], kind="stable"
    )
    expected = load_fragile_support_candidate_substrates().sort_values(
        ["kinase", "site_id"], kind="stable"
    )

    pdt.assert_frame_equal(
        observed_frame.reset_index(drop=True),
        expected.reset_index(drop=True),
    )
    record_parity_metrics(
        request.config,
        family="prediction_science",
        metrics=[
            ("candidate selection row count", observed_frame.shape[0]),
            (
                "candidate kinase count",
                int(observed_frame.loc[:, "kinase"].nunique()),
            ),
        ],
    )
