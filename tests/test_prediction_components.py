from __future__ import annotations

import pandas as pd

from phospy.prediction.aggregation import PredictionAggregator
from phospy.prediction.candidates import CandidateSelector
from phospy.prediction.execution import TraceRecorder


def test_candidate_selector_selects_qualifying_kinases() -> None:
    scores = pd.DataFrame(
        {
            "K1": [0.95, 0.91, 0.40],
            "K2": [0.83, 0.10, 0.82],
        },
        index=["s1", "s2", "s3"],
    )

    selector = CandidateSelector()

    result = selector.select(scores, top=3, score_threshold=0.8, inclusion=2)

    assert result == {"K1": ["s1", "s2"], "K2": ["s1", "s3"]}


def test_prediction_aggregator_initializes_expected_prediction_matrix() -> None:
    feature_mat = pd.DataFrame(
        {"K1": [0.1, 0.2], "K2": [0.3, 0.4]},
        index=["s1", "s2"],
    )

    pred_matrix = PredictionAggregator.initialize_prediction_matrix(
        feature_mat=feature_mat,
        substrate_list={"K1": ["s1"], "K2": ["s2"]},
    )

    assert pred_matrix.index.tolist() == ["s1", "s2"]
    assert pred_matrix.columns.tolist() == ["K1", "K2"]
    assert float(pred_matrix.loc["s1", "K1"]) == 0.0
    assert float(pred_matrix.loc["s2", "K2"]) == 0.0


def test_trace_recorder_create_state_traces_all_kinases_by_default() -> None:
    recorder = TraceRecorder()

    state = recorder.create_state(
        substrate_list={"K1": ["s1"], "K2": ["s2"]},
        trace_level="summary",
        debug_kinases=None,
        trace_sink=None,
    )

    assert state.traced_kinases == {"K1", "K2"}
    assert state.debug_traces == {}
