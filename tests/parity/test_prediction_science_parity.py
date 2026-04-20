from __future__ import annotations

import pandas as pd
import pandas.testing as pdt
import pytest

from phospy.prediction.candidates import build_candidate_substrate_list
from phospy.prediction.motif_scoring import (
    build_motif_library_from_sequences,
    score_phosphosite_motifs,
)
from phospy.prediction.scoring import combine_profile_and_motif_scores
from tests.support.parity_reporting import (
    format_percent,
    format_shape,
    record_parity_metrics,
)
from tests.support.rewrite_fixture_data import (
    load_fragile_support_candidate_substrates,
    load_fragile_support_combined_scores,
    load_fragile_support_combined_weights,
    load_fragile_support_motif_scores,
    load_fragile_support_motif_sizes,
    load_fragile_support_profile_scores,
    load_fragile_support_profile_sizes,
)

pytestmark = pytest.mark.parity


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


def test_motif_scoring_matches_fragile_support_reference_points(
    request: pytest.FixtureRequest,
) -> None:
    sequence_frame = pd.read_csv(
        "tests/fixtures/rewrite_parity/fragile_support_reference/site_sequences.csv"
    )
    sequence_series = (
        sequence_frame.set_index("site_id").loc[:, "centralized_sequence"].astype(str)
    )
    motif_sequences = pd.read_csv(
        "tests/fixtures/rewrite_parity/fragile_support_reference/motif_sequences.csv"
    )
    motif_sequence_map = {
        str(kinase): group.loc[:, "sequence"].astype(str).tolist()
        for kinase, group in motif_sequences.groupby("kinase", sort=False)
    }
    expected = load_fragile_support_motif_scores()

    motif_matrices, motif_sizes = build_motif_library_from_sequences(
        motif_sequences=motif_sequence_map,
    )
    observed = score_phosphosite_motifs(
        site_sequences=sequence_series,
        motif_frequency_matrices=motif_matrices,
        motif_sizes=motif_sizes,
        site_index=expected.index,
        min_motif_size=1,
    ).motif_scores

    assert list(observed.index) == list(expected.index)
    assert set(observed.columns) == set(expected.columns)
    assert ((observed >= 0.0) & (observed <= 1.0)).all().all()
    record_parity_metrics(
        request.config,
        family="prediction_science",
        metrics=[
            ("motif score table shape", format_shape(*observed.shape)),
            ("sites compared", observed.shape[0]),
            ("kinases compared", observed.shape[1]),
        ],
    )


def test_combined_scoring_matches_fragile_support_reference_tables(
    request: pytest.FixtureRequest,
) -> None:
    combined_scores, combined_weights = combine_profile_and_motif_scores(
        motif_scores=load_fragile_support_motif_scores(),
        profile_scores=load_fragile_support_profile_scores(),
        motif_sizes=load_fragile_support_motif_sizes(),
        profile_sizes=load_fragile_support_profile_sizes(),
        allow_profile_only_fallback=False,
    )
    expected_scores = load_fragile_support_combined_scores()
    expected_weights = load_fragile_support_combined_weights()

    pdt.assert_frame_equal(
        combined_scores.sort_index(axis=0).sort_index(axis=1),
        expected_scores.sort_index(axis=0).sort_index(axis=1),
    )
    pdt.assert_frame_equal(
        combined_weights.sort_index(),
        expected_weights.sort_index(),
    )
    aligned_observed = combined_scores.sort_index().sort_index(axis=1)
    aligned_expected = expected_scores.sort_index().sort_index(axis=1)
    absolute_delta = (aligned_observed - aligned_expected).abs()
    record_parity_metrics(
        request.config,
        family="prediction_science",
        metrics=[
            ("combined score table shape", format_shape(*combined_scores.shape)),
            ("combined score mean abs diff", float(absolute_delta.to_numpy().mean())),
            ("combined score max abs diff", float(absolute_delta.to_numpy().max())),
            (
                "combined score mean Spearman column correlation",
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
        scores=load_fragile_support_combined_scores(),
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
