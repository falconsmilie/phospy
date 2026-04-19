from __future__ import annotations

import pandas as pd
import pytest

from phospy import KinasePredictionConfig
from phospy.prediction.candidates import build_candidate_substrate_list
from phospy.prediction.execution import run_adaptive_ensemble_prediction
from tests.support.rewrite_fixture_data import (
    load_adaptive_sampling_edge_combined_scores,
)

pytestmark = pytest.mark.parity


def test_adaptive_ensemble_outputs_match_promoted_fixture_tolerances() -> None:
    combined_scores = load_adaptive_sampling_edge_combined_scores()
    candidate_substrates = build_candidate_substrate_list(
        scores=combined_scores,
        top=4,
        score_threshold=0.8,
        inclusion=1,
    )
    observed = run_adaptive_ensemble_prediction(
        prediction_score_matrix=combined_scores,
        candidate_substrates=candidate_substrates,
        prediction_config=KinasePredictionConfig(
            top_k=4,
            ensemble_size=1,
            mode="adaptive_ensemble",
            n_iterations=2,
            random_state=18,
        ),
    )
    observed_rows = [
        {
            "kinase": str(kinase),
            "site": str(site),
            "prob_class_1": float(score),
        }
        for kinase in observed.columns.astype(str)
        for site, score in observed.loc[:, kinase].items()
    ]
    observed_frame = pd.DataFrame(observed_rows)
    expected_frame = pd.read_csv(
        "tests/fixtures/rewrite_parity/adaptive_sampling_edge/trace_final_ensemble_predictions.csv"
    ).loc[:, ["kinase", "site", "prob_class_1"]]
    expected_frame = expected_frame.drop_duplicates().reset_index(drop=True)

    merged = observed_frame.merge(
        expected_frame, on=["kinase", "site"], suffixes=("_py", "_donor")
    )
    assert not merged.empty
    assert (
        float(
            merged.loc[:, "prob_class_1_py"].corr(merged.loc[:, "prob_class_1_donor"])
        )
        >= 0.999
    )
    assert (
        float(
            (merged.loc[:, "prob_class_1_py"] - merged.loc[:, "prob_class_1_donor"])
            .abs()
            .mean()
        )
        <= 0.01
    )

    observed_top_rows: list[dict[str, object]] = []
    for kinase in observed.columns.astype(str):
        ranked = observed.loc[:, kinase].sort_values(ascending=False).head(4)
        for rank, (site, score) in enumerate(ranked.items(), start=1):
            observed_top_rows.append(
                {
                    "kinase": kinase,
                    "rank": rank,
                    "site": str(site),
                    "prob_class_1": float(score),
                }
            )
    observed_top = pd.DataFrame(observed_top_rows)
    expected_top = pd.read_csv(
        "tests/fixtures/rewrite_parity/adaptive_sampling_edge/trace_final_ensemble_top.csv"
    ).drop_duplicates(subset=["kinase", "rank", "site"])
    expected_top = expected_top.loc[:, ["kinase", "rank", "site", "prob_class_1"]]

    merged_top = observed_top.merge(
        expected_top,
        on=["kinase", "rank"],
        suffixes=("_py", "_donor"),
    )
    assert not merged_top.empty
    assert (
        int((merged_top.loc[:, "site_py"] == merged_top.loc[:, "site_donor"]).sum())
        >= 7
    )
