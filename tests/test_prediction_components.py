from __future__ import annotations

import pandas as pd
import pytest

from phospy.errors import InputCompatibilityError, NoCandidateKinasesError
from phospy.prediction.aggregation import PredictionAggregator
from phospy.prediction.candidates import CandidateSelector
from phospy.prediction.contracts import EnsemblePredictorContract
from phospy.prediction.engines import PredictionExecutionRunner
from phospy.prediction.execution import TraceRecorder
from phospy.validation.requests import PredictionRequest


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


def test_prediction_execution_runner_raises_domain_error_without_candidates() -> None:
    scores = pd.DataFrame({"K1": [0.2, 0.1]}, index=["s1", "s2"])
    request = PredictionRequest.validate_request(
        combined_scores=scores,
        ensemble_size=2,
        top=1,
        score_threshold=0.9,
        inclusion=1,
        n_iterations=1,
        random_state=3,
        capture_debug_trace=False,
        default_svm_mode="default",
    )

    runner = PredictionExecutionRunner(
        candidate_selector=CandidateSelector(),
        prediction_aggregator=PredictionAggregator(),
        trace_recorder=TraceRecorder(),
        ensemble_predictor=None,
    )

    with pytest.raises(
        NoCandidateKinasesError,
        match=(
            r"No candidate kinases qualified for prediction from combined_scores "
            r"using top=1, score_threshold=0\.9, and inclusion=1"
        ),
    ):
        runner.run(request)


def test_prediction_execution_runner_uses_validated_combined_scores_without_recasting(
    monkeypatch,
) -> None:
    class RecordingAggregator:
        def __init__(self) -> None:
            self.seen_feature_mat: pd.DataFrame | None = None

        def initialize_prediction_matrix(
            self,
            *,
            feature_mat: pd.DataFrame,
            substrate_list: dict[str, list[str]],
        ) -> pd.DataFrame:
            self.seen_feature_mat = feature_mat
            return pd.DataFrame(
                0.0,
                index=feature_mat.index.copy(),
                columns=list(substrate_list),
                dtype=float,
            )

        def add_kinase_scores(self, *, pred_matrix: pd.DataFrame, batch) -> None:
            pred_matrix.loc[:, batch.kinase] = batch.scores

        def finalize(
            self,
            *,
            pred_matrix: pd.DataFrame,
            substrate_list,
            request,
            trace_state,
        ):
            return type(
                "PredictionResultStub",
                (),
                {"pred_matrix": pred_matrix, "substrate_list": substrate_list},
            )()

    class RecordingTraceRecorder:
        def create_state(
            self, *, substrate_list, trace_level, debug_kinases, trace_sink
        ):
            return object()

        def flush_final(self, *, trace_state) -> None:
            return None

    class RecordingEnsemblePredictor(EnsemblePredictorContract):
        def predict_kinase(
            self,
            *,
            kinase: str,
            substrates: list[str],
            feature_mat: pd.DataFrame,
            request: PredictionRequest,
            trace_state,
            sampling_session,
        ):
            return type(
                "BatchStub",
                (),
                {
                    "kinase": kinase,
                    "scores": pd.Series(
                        0.0,
                        index=feature_mat.index.copy(),
                        dtype=float,
                    ),
                },
            )()

    scores = pd.DataFrame({"K1": ["0.95", "0.91", "0.40"]}, index=["s1", "s2", "s3"])
    request = PredictionRequest.validate_request(
        combined_scores=scores,
        ensemble_size=2,
        top=3,
        score_threshold=0.8,
        inclusion=2,
        n_iterations=1,
        random_state=3,
        capture_debug_trace=False,
        default_svm_mode="default",
    )

    def forbid_dataframe_astype(self, dtype=None, copy=None, errors="raise"):
        raise AssertionError(
            "prediction runner should not recast validated combined_scores"
        )

    monkeypatch.setattr(pd.DataFrame, "astype", forbid_dataframe_astype)
    aggregator = RecordingAggregator()
    runner = PredictionExecutionRunner(
        candidate_selector=CandidateSelector(),
        prediction_aggregator=aggregator,
        trace_recorder=RecordingTraceRecorder(),
        ensemble_predictor=RecordingEnsemblePredictor(),
    )

    result = runner.run(request)

    assert aggregator.seen_feature_mat is request.combined_scores
    assert list(result.pred_matrix.columns) == ["K1"]


def test_prediction_execution_runner_rejects_non_compliant_ensemble_predictor() -> None:
    with pytest.raises(
        InputCompatibilityError,
        match="EnsemblePredictorContract",
    ):
        PredictionExecutionRunner(
            candidate_selector=CandidateSelector(),
            prediction_aggregator=PredictionAggregator(),
            trace_recorder=TraceRecorder(),
            ensemble_predictor=object(),
        )


def test_prediction_execution_runner_passes_sampling_session_to_contract_predictor() -> (
    None
):
    class RecordingAggregator:
        def initialize_prediction_matrix(
            self,
            *,
            feature_mat: pd.DataFrame,
            substrate_list: dict[str, list[str]],
        ) -> pd.DataFrame:
            return pd.DataFrame(
                0.0,
                index=feature_mat.index.copy(),
                columns=list(substrate_list),
                dtype=float,
            )

        def add_kinase_scores(self, *, pred_matrix: pd.DataFrame, batch) -> None:
            pred_matrix.loc[:, batch.kinase] = batch.scores

        def finalize(
            self, *, pred_matrix: pd.DataFrame, substrate_list, request, trace_state
        ):
            return type(
                "PredictionResultStub",
                (),
                {"pred_matrix": pred_matrix, "substrate_list": substrate_list},
            )()

    class RecordingTraceRecorder:
        def create_state(
            self, *, substrate_list, trace_level, debug_kinases, trace_sink
        ):
            return object()

        def flush_final(self, *, trace_state) -> None:
            return None

    class RecordingEnsemblePredictor(EnsemblePredictorContract):
        def __init__(self) -> None:
            self.seen_sampling_session = None

        def predict_kinase(
            self,
            *,
            kinase: str,
            substrates: list[str],
            feature_mat: pd.DataFrame,
            request: PredictionRequest,
            trace_state,
            sampling_session,
        ):
            self.seen_sampling_session = sampling_session
            return type(
                "BatchStub",
                (),
                {
                    "kinase": kinase,
                    "scores": pd.Series(
                        0.0,
                        index=feature_mat.index.copy(),
                        dtype=float,
                    ),
                },
            )()

    scores = pd.DataFrame({"K1": [0.95, 0.91, 0.40]}, index=["s1", "s2", "s3"])
    request = PredictionRequest.validate_request(
        combined_scores=scores,
        ensemble_size=2,
        top=3,
        score_threshold=0.8,
        inclusion=2,
        n_iterations=1,
        random_state=3,
        capture_debug_trace=False,
        default_svm_mode="default",
    )

    predictor = RecordingEnsemblePredictor()
    runner = PredictionExecutionRunner(
        candidate_selector=CandidateSelector(),
        prediction_aggregator=RecordingAggregator(),
        trace_recorder=RecordingTraceRecorder(),
        ensemble_predictor=predictor,
    )

    runner.run(request)

    assert predictor.seen_sampling_session is not None
