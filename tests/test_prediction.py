from __future__ import annotations

from pathlib import Path
from typing import get_type_hints

import numpy as np
import pandas as pd
import pytest

import phospy.prediction.sampling_core as _sampling_core
from phospy import SimpleKinaseWorkflow
from phospy.api import KinaseActivityConfig, PredictionRunConfig
from phospy.errors import (
    NoCandidateKinasesError,
    RequestValidationError,
    TableSchemaError,
)
from phospy.io import load_pred_mat
from phospy.prediction import (
    EnsemblePredictorContract,
    KinasePredictionResult,
    KinasePredictor,
    KinaseScoringResult,
    PredictionSamplingTrace,
    PredMatResult,
    build_candidate_substrate_list,
    prediction_debug_trace_tables,
)
from phospy.prediction.policies import (
    PredictionSamplingRandomSource,
    resolve_prediction_sampling_policy,
)
from phospy.prediction.sampling import (
    make_prediction_random_generators as _make_prediction_random_generators,
)
from phospy.prediction.sampling import multi_ada_sampling as _multi_ada_sampling
from phospy.prediction.sampling import (
    transform_resampling_probabilities as _transform_resampling_probabilities,
)
from phospy.prediction.sampling_core import (
    _resolve_final_score_series as _resolve_final_score_series,
)
from phospy.prediction.svm import _RLikeStandardScaler
from phospy.prediction.svm import make_svm as _make_svm
from phospy.prediction.svm import require_sklearn as _require_sklearn
from phospy.prediction.svm import (
    resolve_svm_probability_random_state as _resolve_svm_probability_random_state,
)
from phospy.prediction.trace_runtime import TraceSink
from phospy.prediction.traces import DirectoryTraceSink, create_trace_sink
from phospy.validation.requests import PredictionRequest

SIMPLE_WORKFLOW_FIXTURE_DIR = (
    Path(__file__).resolve().parents[1] / "examples" / "data" / "simple_workflow"
)


def make_combined_scores() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "KINASE_A": [0.95, 0.93, 0.91, 0.89, 0.20, 0.18, 0.16, 0.14],
            "KINASE_B": [0.10, 0.12, 0.14, 0.16, 0.96, 0.94, 0.92, 0.90],
        },
        index=[f"SITE_{i}" for i in range(1, 9)],
    )


def test_build_candidate_substrate_list_applies_selection_rules() -> None:
    substrate_list = build_candidate_substrate_list(
        make_combined_scores(),
        top=4,
        score_threshold=0.9,
        inclusion=2,
    )

    assert substrate_list == {
        "KINASE_A": ["SITE_1", "SITE_2", "SITE_3"],
        "KINASE_B": ["SITE_5", "SITE_6", "SITE_7"],
    }


def test_kinase_predictor_returns_probability_matrix() -> None:
    predictor = KinasePredictor()

    result = predictor.predict(
        combined_scores=make_combined_scores(),
        ensemble_size=3,
        top=4,
        score_threshold=0.85,
        inclusion=3,
        n_iterations=2,
        random_state=5,
    )

    assert isinstance(result, KinasePredictionResult)
    assert list(result.pred_matrix.columns) == ["KINASE_A", "KINASE_B"]
    assert ((result.pred_matrix >= 0.0) & (result.pred_matrix <= 1.0)).all().all()
    assert (
        result.pred_matrix.loc[
            ["SITE_1", "SITE_2", "SITE_3", "SITE_4"], "KINASE_A"
        ].mean()
        > result.pred_matrix.loc[
            ["SITE_5", "SITE_6", "SITE_7", "SITE_8"], "KINASE_A"
        ].mean()
    )
    assert (
        result.pred_matrix.loc[
            ["SITE_5", "SITE_6", "SITE_7", "SITE_8"], "KINASE_B"
        ].mean()
        > result.pred_matrix.loc[
            ["SITE_1", "SITE_2", "SITE_3", "SITE_4"], "KINASE_B"
        ].mean()
    )


def test_prediction_result_matrix_is_detached_from_combined_scores_input() -> None:
    predictor = KinasePredictor()
    combined_scores = make_combined_scores()
    original = combined_scores.copy(deep=True)

    result = predictor.predict(
        combined_scores=combined_scores,
        ensemble_size=2,
        top=4,
        score_threshold=0.85,
        inclusion=3,
        n_iterations=2,
        random_state=5,
    )

    result.pred_matrix.loc["SITE_1", "KINASE_A"] = -999.0

    pd.testing.assert_frame_equal(combined_scores, original)


def test_prediction_result_exposes_canonical_pred_mat_result() -> None:
    predictor = KinasePredictor()

    result = predictor.predict(
        combined_scores=make_combined_scores(),
        ensemble_size=2,
        top=4,
        score_threshold=0.85,
        inclusion=3,
        n_iterations=2,
        random_state=5,
    )

    assert isinstance(result.pred_mat_result, PredMatResult)
    assert result.pred_mat_result.data_frame is result.pred_matrix
    assert (
        result.pred_mat_result.phosphosite_ids.tolist()
        == result.pred_matrix.index.tolist()
    )
    assert (
        result.pred_mat_result.kinase_names.tolist()
        == result.pred_matrix.columns.tolist()
    )
    assert not hasattr(result, "pred_mat")


def test_prediction_result_pred_mat_export_round_trips_through_loader(
    tmp_path: Path,
) -> None:
    predictor = KinasePredictor()

    result = predictor.predict(
        combined_scores=make_combined_scores(),
        ensemble_size=2,
        top=4,
        score_threshold=0.85,
        inclusion=3,
        n_iterations=2,
        random_state=5,
    )
    export_path = result.pred_mat_result.to_csv(tmp_path / "predMat.csv")

    reloaded = load_pred_mat(export_path)

    pd.testing.assert_frame_equal(
        reloaded,
        result.pred_mat_result.to_frame(copy=False),
    )


def test_predict_from_scoring_result_uses_combined_scores() -> None:
    predictor = KinasePredictor()
    scoring_result = KinaseScoringResult(
        profile_scores=make_combined_scores(),
        combined_scores=make_combined_scores(),
    )

    result = predictor.predict_from_scoring_result(
        scoring_result,
        ensemble_size=2,
        top=4,
        score_threshold=0.85,
        inclusion=3,
        n_iterations=2,
        random_state=9,
    )

    assert set(result.substrate_list) == {"KINASE_A", "KINASE_B"}


def test_predict_from_scoring_result_requires_explicit_profile_fallback() -> None:
    predictor = KinasePredictor()
    scoring_result = KinaseScoringResult(profile_scores=make_combined_scores())

    with pytest.raises(ValueError, match="combined_scores"):
        predictor.predict_from_scoring_result(scoring_result)

    result = predictor.predict_from_scoring_result(
        scoring_result,
        ensemble_size=2,
        top=4,
        score_threshold=0.85,
        inclusion=3,
        n_iterations=2,
        random_state=9,
        allow_profile_only_fallback=True,
    )

    assert list(result.pred_matrix.columns) == ["KINASE_A", "KINASE_B"]


def test_predict_raises_domain_error_when_no_kinases_pass_inclusion() -> None:
    predictor = KinasePredictor()

    with pytest.raises(
        NoCandidateKinasesError,
        match=(
            r"No candidate kinases qualified for prediction from combined_scores "
            r"using top=2, score_threshold=0\.99, and inclusion=2"
        ),
    ):
        predictor.predict(
            combined_scores=make_combined_scores(),
            ensemble_size=2,
            top=2,
            score_threshold=0.99,
            inclusion=2,
            n_iterations=2,
            random_state=3,
        )


def test_predict_uses_summary_trace_level_by_default_when_debug_is_enabled() -> None:
    predictor = KinasePredictor()

    result = predictor.predict(
        combined_scores=make_combined_scores(),
        ensemble_size=2,
        top=4,
        score_threshold=0.85,
        inclusion=3,
        n_iterations=2,
        random_state=5,
        capture_debug_trace=True,
        debug_kinases=["KINASE_A"],
        debug_top_n=3,
    )

    assert result.trace_level == "summary"
    assert result.debug_traces is not None
    assert set(result.debug_traces) == {"KINASE_A"}
    trace = result.debug_traces["KINASE_A"]
    assert trace.candidate_substrates == ["SITE_1", "SITE_2", "SITE_3", "SITE_4"]
    assert len(trace.negative_pool_sites) == 4
    assert len(trace.ensemble_traces) == 2
    ensemble_trace = trace.ensemble_traces[0]
    assert len(ensemble_trace.initial_negative_sites) == 4
    assert ensemble_trace.iterations == []
    assert ensemble_trace.final_prediction_probabilities is None
    assert ensemble_trace.final_decision_values is None
    assert len(ensemble_trace.final_top_sites) == 3


def test_predict_can_stream_full_trace_tables_to_sink(tmp_path: Path) -> None:
    predictor = KinasePredictor()

    result = predictor.predict(
        combined_scores=make_combined_scores(),
        ensemble_size=2,
        top=4,
        score_threshold=0.85,
        inclusion=3,
        n_iterations=2,
        random_state=5,
        capture_debug_trace=True,
        debug_kinases=["KINASE_A"],
        debug_top_n=3,
        trace_level="full",
        trace_sink=tmp_path / "trace_output",
    )

    assert result.trace_level == "full"
    assert result.trace_sink is not None
    assert result.debug_traces is not None
    assert (tmp_path / "trace_output" / "trace_iteration_samples.csv").exists()
    ensemble_trace = result.debug_traces["KINASE_A"].ensemble_traces[0]
    assert ensemble_trace.iterations == []
    tables = prediction_debug_trace_tables(result)
    assert not tables["trace_iteration_probabilities"].empty
    assert not tables["trace_iteration_samples"].empty
    assert not tables["trace_final_ensemble_predictions"].empty
    assert not tables["trace_final_ensemble_top"].empty


def test_predict_request_transfers_owned_runtime_trace_sink_to_result_on_success() -> (
    None
):
    predictor = KinasePredictor()
    request = PredictionRequest.validate_request(
        combined_scores=make_combined_scores(),
        ensemble_size=2,
        top=4,
        score_threshold=0.85,
        inclusion=3,
        n_iterations=2,
        random_state=5,
        debug_top_n=3,
        trace_level="full",
        default_svm_mode="default",
        capture_debug_trace=True,
        debug_kinases=["KINASE_A"],
    )

    assert request.trace_sink is None

    result = predictor.predict_request(request)

    assert result.trace_level == "full"
    assert result.trace_sink is not None
    assert result.owns_trace_sink is True
    output_dir = result.trace_sink.output_dir
    assert output_dir.exists()

    tables = prediction_debug_trace_tables(result)
    assert not tables["trace_iteration_probabilities"].empty
    assert not tables["trace_iteration_samples"].empty

    result.close()
    assert not output_dir.exists()


def test_directory_trace_sink_buffers_until_read_or_flush(tmp_path: Path) -> None:
    sink = DirectoryTraceSink(tmp_path / "trace_output")

    sink.write_rows(
        "trace_iteration_samples",
        [
            {
                "kinase": "KINASE_A",
                "ensemble": 1,
                "iteration": 1,
                "class_label": 1,
                "draw": 1,
                "site": "SITE_1",
            }
        ],
    )

    csv_path = tmp_path / "trace_output" / "trace_iteration_samples.csv"
    assert not csv_path.exists()

    table = sink.read_table("trace_iteration_samples")

    assert csv_path.exists()
    assert table.loc[:, "site"].tolist() == ["SITE_1"]


def test_directory_trace_sink_auto_flushes_after_threshold(tmp_path: Path) -> None:
    sink = DirectoryTraceSink(tmp_path / "trace_output", max_buffer_rows=2)

    sink.write_rows(
        "trace_iteration_samples",
        [
            {
                "kinase": "KINASE_A",
                "ensemble": 1,
                "iteration": 1,
                "class_label": 1,
                "draw": 1,
                "site": "SITE_1",
            }
        ],
    )

    csv_path = tmp_path / "trace_output" / "trace_iteration_samples.csv"
    assert not csv_path.exists()

    sink.write_rows(
        "trace_iteration_samples",
        [
            {
                "kinase": "KINASE_A",
                "ensemble": 1,
                "iteration": 1,
                "class_label": 2,
                "draw": 2,
                "site": "SITE_2",
            }
        ],
    )

    assert csv_path.exists()
    table = pd.read_csv(csv_path)
    assert table.loc[:, "site"].tolist() == ["SITE_1", "SITE_2"]


def test_directory_trace_sink_rejects_non_positive_auto_flush_threshold(
    tmp_path: Path,
) -> None:
    with pytest.raises(ValueError, match="max_buffer_rows"):
        DirectoryTraceSink(tmp_path / "trace_output", max_buffer_rows=0)


def test_multi_ada_sampling_uses_numpy_iteration_path_when_trace_capture_is_disabled(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    train_mat = pd.DataFrame(
        {
            "feature_1": [1.0, 0.9, 0.1, 0.0],
            "feature_2": [1.0, 0.8, 0.2, 0.1],
        },
        index=["SITE_1", "SITE_2", "SITE_3", "SITE_4"],
    )
    test_mat = train_mat.copy()
    labels = np.asarray([1, 1, 2, 2], dtype=int)
    calls = {"sampling": 0, "trace": 0}

    original_sampling_state = _sampling_core._extract_iteration_sampling_state
    original_trace_payload = _sampling_core._extract_iteration_trace_payload

    def counting_sampling_state(**kwargs):
        calls["sampling"] += 1
        return original_sampling_state(**kwargs)

    def counting_trace_payload(**kwargs):
        calls["trace"] += 1
        return original_trace_payload(**kwargs)

    monkeypatch.setattr(
        _sampling_core,
        "_extract_iteration_sampling_state",
        counting_sampling_state,
    )
    monkeypatch.setattr(
        _sampling_core,
        "_extract_iteration_trace_payload",
        counting_trace_payload,
    )

    _multi_ada_sampling(
        train_mat=train_mat,
        test_mat=test_mat,
        labels=labels,
        kernel="rbf",
        n_iterations=2,
        resampling_rng=np.random.default_rng(7),
        capture_trace=False,
        trace_level="none",
        trace_sink=None,
        kinase="KINASE_A",
        ensemble_index=1,
        initial_negative_sites=["SITE_3", "SITE_4"],
        debug_top_n=2,
        svm_mode="default",
        sampling_override=None,
    )

    assert calls == {"sampling": 2, "trace": 0}


def test_multi_ada_sampling_builds_trace_payload_only_for_full_trace(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    train_mat = pd.DataFrame(
        {
            "feature_1": [1.0, 0.9, 0.1, 0.0],
            "feature_2": [1.0, 0.8, 0.2, 0.1],
        },
        index=["SITE_1", "SITE_2", "SITE_3", "SITE_4"],
    )
    test_mat = train_mat.copy()
    labels = np.asarray([1, 1, 2, 2], dtype=int)
    calls = {"sampling": 0, "trace": 0}

    class MemoryTraceSink(TraceSink):
        def __init__(self) -> None:
            self.rows: dict[str, list[dict[str, object]]] = {}

        def write_rows(self, table_name: str, rows: list[dict[str, object]]) -> None:
            self.rows.setdefault(table_name, []).extend(rows)

        def read_table(self, table_name: str) -> pd.DataFrame:
            return pd.DataFrame(self.rows.get(table_name, []))

    original_sampling_state = _sampling_core._extract_iteration_sampling_state
    original_trace_payload = _sampling_core._extract_iteration_trace_payload

    def counting_sampling_state(**kwargs):
        calls["sampling"] += 1
        return original_sampling_state(**kwargs)

    def counting_trace_payload(**kwargs):
        calls["trace"] += 1
        return original_trace_payload(**kwargs)

    monkeypatch.setattr(
        _sampling_core,
        "_extract_iteration_sampling_state",
        counting_sampling_state,
    )
    monkeypatch.setattr(
        _sampling_core,
        "_extract_iteration_trace_payload",
        counting_trace_payload,
    )

    _multi_ada_sampling(
        train_mat=train_mat,
        test_mat=test_mat,
        labels=labels,
        kernel="rbf",
        n_iterations=2,
        resampling_rng=np.random.default_rng(7),
        capture_trace=True,
        trace_level="full",
        trace_sink=MemoryTraceSink(),
        kinase="KINASE_A",
        ensemble_index=1,
        initial_negative_sites=["SITE_3", "SITE_4"],
        debug_top_n=2,
        svm_mode="default",
        sampling_override=None,
    )

    assert calls == {"sampling": 2, "trace": 2}


def test_multi_ada_sampling_uses_precomputed_arrays_without_dataframe_extraction(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    train_mat = pd.DataFrame(
        {
            "feature_1": [1.0, 0.9, 0.1, 0.0],
            "feature_2": [1.0, 0.8, 0.2, 0.1],
        },
        index=["SITE_1", "SITE_2", "SITE_3", "SITE_4"],
    )
    test_mat = train_mat.copy()
    labels = np.asarray([1, 1, 2, 2], dtype=int)
    train_values = train_mat.to_numpy(dtype=float)
    test_values = test_mat.to_numpy(dtype=float)

    def forbid_to_numpy(self, dtype=None, copy=False, na_value=None):
        raise AssertionError(
            "multi_ada_sampling should use the supplied precomputed arrays"
        )

    monkeypatch.setattr(pd.DataFrame, "to_numpy", forbid_to_numpy)

    scores, trace = _multi_ada_sampling(
        train_mat=None,
        test_mat=None,
        labels=labels,
        kernel="rbf",
        n_iterations=2,
        resampling_rng=np.random.default_rng(7),
        capture_trace=False,
        trace_level="none",
        trace_sink=None,
        kinase="KINASE_A",
        ensemble_index=1,
        initial_negative_sites=["SITE_3", "SITE_4"],
        debug_top_n=2,
        svm_mode="default",
        sampling_override=None,
        train_values=train_values,
        train_index=train_mat.index,
        test_values=test_values,
        test_index=test_mat.index,
    )

    assert trace is None
    assert scores.index.tolist() == test_mat.index.tolist()
    assert ((scores >= 0.0) & (scores <= 1.0)).all()


def test_multi_ada_sampling_requires_trace_sink_for_full_trace() -> None:
    train_mat = pd.DataFrame(
        {
            "feature_1": [1.0, 0.9, 0.1, 0.0],
            "feature_2": [1.0, 0.8, 0.2, 0.1],
        },
        index=["SITE_1", "SITE_2", "SITE_3", "SITE_4"],
    )
    test_mat = train_mat.copy()
    labels = np.asarray([1, 1, 2, 2], dtype=int)

    with pytest.raises(ValueError, match="trace sink"):
        _multi_ada_sampling(
            train_mat=train_mat,
            test_mat=test_mat,
            labels=labels,
            kernel="rbf",
            n_iterations=2,
            resampling_rng=np.random.default_rng(7),
            capture_trace=True,
            trace_level="full",
            trace_sink=None,
            kinase="KINASE_A",
            ensemble_index=1,
            initial_negative_sites=["SITE_3", "SITE_4"],
            debug_top_n=2,
            svm_mode="default",
            sampling_override=None,
        )


def test_predict_can_capture_debug_trace_for_all_kinases() -> None:
    predictor = KinasePredictor()

    result = predictor.predict(
        combined_scores=make_combined_scores(),
        ensemble_size=1,
        top=4,
        score_threshold=0.85,
        inclusion=3,
        n_iterations=2,
        random_state=5,
        capture_debug_trace=True,
    )

    assert result.debug_traces is not None
    assert set(result.debug_traces) == {"KINASE_A", "KINASE_B"}


def test_predict_from_scoring_result_passes_debug_options_through() -> None:
    predictor = KinasePredictor()
    scoring_result = KinaseScoringResult(
        profile_scores=make_combined_scores(),
        combined_scores=make_combined_scores(),
    )

    result = predictor.predict_from_scoring_result(
        scoring_result,
        ensemble_size=1,
        top=4,
        score_threshold=0.85,
        inclusion=3,
        n_iterations=2,
        random_state=9,
        capture_debug_trace=True,
        debug_kinases=["KINASE_B"],
        debug_top_n=2,
    )

    assert result.debug_traces is not None
    assert set(result.debug_traces) == {"KINASE_B"}
    assert len(result.debug_traces["KINASE_B"].ensemble_traces[0].final_top_sites) == 2


def test_build_candidate_substrate_list_preserves_input_order_for_ties() -> None:
    scores = pd.DataFrame(
        {"KINASE_A": [0.95, 0.95, 0.90, 0.10]},
        index=["SITE_A", "SITE_B", "SITE_C", "SITE_D"],
    )

    substrate_list = build_candidate_substrate_list(
        scores,
        top=3,
        score_threshold=0.5,
        inclusion=2,
    )

    assert substrate_list == {"KINASE_A": ["SITE_A", "SITE_B", "SITE_C"]}


def test_r_like_standard_scaler_uses_sample_standard_deviation() -> None:
    values = np.array([[1.0, 2.0], [3.0, 6.0], [5.0, 10.0]])

    scaler = _RLikeStandardScaler().fit(values)
    transformed = scaler.transform(values)

    assert scaler.mean_ is not None
    assert scaler.scale_ is not None
    assert np.allclose(scaler.mean_, values.mean(axis=0))
    assert np.allclose(scaler.scale_, values.std(axis=0, ddof=1))
    assert np.allclose(transformed.mean(axis=0), [0.0, 0.0])
    assert np.allclose(transformed.std(axis=0, ddof=1), [1.0, 1.0])


def test_make_svm_default_mode_uses_sklearn_defaults_explicitly() -> None:
    StandardScaler, SVC = _require_sklearn()

    model = _make_svm(
        StandardScaler=StandardScaler,
        SVC=SVC,
        kernel="rbf",
        svm_mode="default",
    )

    assert model.steps[0][1].__class__.__name__ == "StandardScaler"
    assert model.steps[1][1].gamma == "scale"


def test_make_svm_r_parity_mode_uses_r_like_scaler_and_gamma_auto() -> None:
    StandardScaler, SVC = _require_sklearn()

    model = _make_svm(
        StandardScaler=StandardScaler,
        SVC=SVC,
        kernel="rbf",
        svm_mode="r_parity",
    )

    assert isinstance(model.steps[0][1], _RLikeStandardScaler)
    assert model.steps[1][1].gamma == "auto"


def test_predict_accepts_explicit_r_parity_mode() -> None:
    predictor = KinasePredictor(svm_mode="r_parity")

    result = predictor.predict(
        combined_scores=make_combined_scores(),
        ensemble_size=2,
        top=4,
        score_threshold=0.85,
        inclusion=3,
        n_iterations=2,
        random_state=5,
    )

    assert isinstance(result, KinasePredictionResult)
    assert list(result.pred_matrix.columns) == ["KINASE_A", "KINASE_B"]


def test_predict_rejects_unknown_svm_mode() -> None:
    with pytest.raises(ValueError, match="svm_mode"):
        KinasePredictor(svm_mode="broken")  # type: ignore[arg-type]


def test_predict_from_scoring_result_allows_svm_mode_override() -> None:
    predictor = KinasePredictor(svm_mode="default")
    scoring_result = KinaseScoringResult(
        profile_scores=make_combined_scores(),
        combined_scores=make_combined_scores(),
    )

    result = predictor.predict_from_scoring_result(
        scoring_result,
        ensemble_size=2,
        top=4,
        score_threshold=0.85,
        inclusion=3,
        n_iterations=2,
        random_state=9,
        svm_mode="r_parity",
    )

    assert list(result.pred_matrix.columns) == ["KINASE_A", "KINASE_B"]


def _write_sampling_trace_fixture(trace_dir: Path) -> Path:
    trace_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        [
            {"kinase": "KINASE_A", "ensemble": 1, "draw": 1, "site": "SITE_8"},
            {"kinase": "KINASE_A", "ensemble": 1, "draw": 2, "site": "SITE_8"},
            {"kinase": "KINASE_A", "ensemble": 1, "draw": 3, "site": "SITE_7"},
            {"kinase": "KINASE_A", "ensemble": 1, "draw": 4, "site": "SITE_6"},
        ]
    ).to_csv(trace_dir / "trace_initial_negatives.csv", index=False)
    pd.DataFrame(
        [
            {
                "kinase": "KINASE_A",
                "ensemble": 1,
                "iteration": 1,
                "class_label": 1,
                "draw": 1,
                "site": "SITE_4",
            },
            {
                "kinase": "KINASE_A",
                "ensemble": 1,
                "iteration": 1,
                "class_label": 1,
                "draw": 2,
                "site": "SITE_4",
            },
            {
                "kinase": "KINASE_A",
                "ensemble": 1,
                "iteration": 1,
                "class_label": 1,
                "draw": 3,
                "site": "SITE_2",
            },
            {
                "kinase": "KINASE_A",
                "ensemble": 1,
                "iteration": 1,
                "class_label": 1,
                "draw": 4,
                "site": "SITE_1",
            },
            {
                "kinase": "KINASE_A",
                "ensemble": 1,
                "iteration": 1,
                "class_label": 2,
                "draw": 1,
                "site": "SITE_8",
            },
            {
                "kinase": "KINASE_A",
                "ensemble": 1,
                "iteration": 1,
                "class_label": 2,
                "draw": 2,
                "site": "SITE_8",
            },
            {
                "kinase": "KINASE_A",
                "ensemble": 1,
                "iteration": 1,
                "class_label": 2,
                "draw": 3,
                "site": "SITE_7",
            },
            {
                "kinase": "KINASE_A",
                "ensemble": 1,
                "iteration": 1,
                "class_label": 2,
                "draw": 4,
                "site": "SITE_6",
            },
            {
                "kinase": "KINASE_A",
                "ensemble": 1,
                "iteration": 2,
                "class_label": 1,
                "draw": 1,
                "site": "SITE_3",
            },
            {
                "kinase": "KINASE_A",
                "ensemble": 1,
                "iteration": 2,
                "class_label": 1,
                "draw": 2,
                "site": "SITE_3",
            },
            {
                "kinase": "KINASE_A",
                "ensemble": 1,
                "iteration": 2,
                "class_label": 1,
                "draw": 3,
                "site": "SITE_2",
            },
            {
                "kinase": "KINASE_A",
                "ensemble": 1,
                "iteration": 2,
                "class_label": 1,
                "draw": 4,
                "site": "SITE_1",
            },
            {
                "kinase": "KINASE_A",
                "ensemble": 1,
                "iteration": 2,
                "class_label": 2,
                "draw": 1,
                "site": "SITE_8",
            },
            {
                "kinase": "KINASE_A",
                "ensemble": 1,
                "iteration": 2,
                "class_label": 2,
                "draw": 2,
                "site": "SITE_7",
            },
            {
                "kinase": "KINASE_A",
                "ensemble": 1,
                "iteration": 2,
                "class_label": 2,
                "draw": 3,
                "site": "SITE_7",
            },
            {
                "kinase": "KINASE_A",
                "ensemble": 1,
                "iteration": 2,
                "class_label": 2,
                "draw": 4,
                "site": "SITE_6",
            },
        ]
    ).to_csv(trace_dir / "trace_iteration_samples.csv", index=False)
    return trace_dir


def test_predict_can_replay_sampling_trace_from_directory(tmp_path: Path) -> None:
    predictor = KinasePredictor()
    trace_dir = _write_sampling_trace_fixture(tmp_path / "prediction_trace")

    result = predictor.predict(
        combined_scores=make_combined_scores(),
        ensemble_size=1,
        top=4,
        score_threshold=0.85,
        inclusion=3,
        n_iterations=2,
        random_state=5,
        capture_debug_trace=True,
        debug_kinases=["KINASE_A"],
        sampling_trace=trace_dir,
        trace_level="full",
        trace_sink=tmp_path / "trace_output_from_directory",
    )

    assert result.debug_traces is not None
    ensemble_trace = result.debug_traces["KINASE_A"].ensemble_traces[0]
    assert ensemble_trace.initial_negative_sites == [
        "SITE_8",
        "SITE_8",
        "SITE_7",
        "SITE_6",
    ]
    tables = prediction_debug_trace_tables(result)
    samples = tables["trace_iteration_samples"]
    positive_iter_1 = samples.loc[
        (samples["iteration"] == 1) & (samples["class_label"] == 1), "site"
    ].tolist()
    negative_iter_1 = samples.loc[
        (samples["iteration"] == 1) & (samples["class_label"] == 2), "site"
    ].tolist()
    positive_iter_2 = samples.loc[
        (samples["iteration"] == 2) & (samples["class_label"] == 1), "site"
    ].tolist()
    negative_iter_2 = samples.loc[
        (samples["iteration"] == 2) & (samples["class_label"] == 2), "site"
    ].tolist()
    assert positive_iter_1 == ["SITE_4", "SITE_4", "SITE_2", "SITE_1"]
    assert negative_iter_1 == ["SITE_8", "SITE_8", "SITE_7", "SITE_6"]
    assert positive_iter_2 == ["SITE_3", "SITE_3", "SITE_2", "SITE_1"]
    assert negative_iter_2 == ["SITE_8", "SITE_7", "SITE_7", "SITE_6"]


def test_predict_accepts_preloaded_sampling_trace(tmp_path: Path) -> None:
    predictor = KinasePredictor()
    trace_dir = _write_sampling_trace_fixture(tmp_path / "prediction_trace")
    sampling_trace = PredictionSamplingTrace.from_trace_directory(trace_dir)

    result = predictor.predict(
        combined_scores=make_combined_scores(),
        ensemble_size=1,
        top=4,
        score_threshold=0.85,
        inclusion=3,
        n_iterations=2,
        random_state=5,
        capture_debug_trace=True,
        debug_kinases=["KINASE_A"],
        sampling_trace=sampling_trace,
        trace_level="summary",
    )

    assert result.debug_traces is not None
    ensemble_trace = result.debug_traces["KINASE_A"].ensemble_traces[0]
    assert ensemble_trace.initial_negative_sites == [
        "SITE_8",
        "SITE_8",
        "SITE_7",
        "SITE_6",
    ]


def test_sampling_trace_directory_supports_parquet_replay(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    trace_dir = tmp_path / "prediction_trace"
    trace_dir.mkdir(parents=True, exist_ok=True)

    initial_frame = pd.DataFrame(
        [
            {"kinase": "KINASE_A", "ensemble": 1, "draw": 1, "site": "SITE_8"},
            {"kinase": "KINASE_A", "ensemble": 1, "draw": 2, "site": "SITE_7"},
        ]
    )
    samples_frame = pd.DataFrame(
        [
            {
                "kinase": "KINASE_A",
                "ensemble": 1,
                "iteration": 1,
                "class_label": 1,
                "draw": 1,
                "site": "SITE_4",
            },
            {
                "kinase": "KINASE_A",
                "ensemble": 1,
                "iteration": 1,
                "class_label": 2,
                "draw": 1,
                "site": "SITE_8",
            },
        ]
    )
    (trace_dir / "trace_initial_negatives.part-000000.parquet").touch()
    (trace_dir / "trace_iteration_samples.part-000000.parquet").touch()

    def fake_read_parquet(path: Path | str) -> pd.DataFrame:
        filename = Path(path).name
        if filename.startswith("trace_initial_negatives"):
            return initial_frame.copy()
        if filename.startswith("trace_iteration_samples"):
            return samples_frame.copy()
        raise AssertionError(f"unexpected parquet path: {filename}")

    monkeypatch.setattr(pd, "read_parquet", fake_read_parquet)

    sampling_trace = PredictionSamplingTrace.from_trace_directory(trace_dir)
    ensemble_trace = sampling_trace.get_ensemble_override("KINASE_A", 1)

    assert ensemble_trace is not None
    assert ensemble_trace.initial_negative_sites == ["SITE_8", "SITE_7"]
    assert ensemble_trace.iteration_sample_sites == {1: {1: ["SITE_4"], 2: ["SITE_8"]}}


def test_sampling_trace_directory_supports_multi_part_parquet_replay_incrementally(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    trace_dir = tmp_path / "prediction_trace"
    trace_dir.mkdir(parents=True, exist_ok=True)

    initial_part_1 = pd.DataFrame(
        [
            {"kinase": "KINASE_A", "ensemble": 1, "draw": 2, "site": "SITE_7"},
            {"kinase": "KINASE_A", "ensemble": 1, "draw": 4, "site": "SITE_5"},
        ]
    )
    initial_part_2 = pd.DataFrame(
        [
            {"kinase": "KINASE_A", "ensemble": 1, "draw": 1, "site": "SITE_8"},
            {"kinase": "KINASE_A", "ensemble": 1, "draw": 3, "site": "SITE_6"},
        ]
    )
    samples_part_1 = pd.DataFrame(
        [
            {
                "kinase": "KINASE_A",
                "ensemble": 1,
                "iteration": 1,
                "class_label": 2,
                "draw": 2,
                "site": "SITE_8",
            },
            {
                "kinase": "KINASE_A",
                "ensemble": 1,
                "iteration": 1,
                "class_label": 1,
                "draw": 1,
                "site": "SITE_4",
            },
        ]
    )
    samples_part_2 = pd.DataFrame(
        [
            {
                "kinase": "KINASE_A",
                "ensemble": 1,
                "iteration": 1,
                "class_label": 2,
                "draw": 1,
                "site": "SITE_7",
            },
        ]
    )

    initial_part_path_1 = trace_dir / "trace_initial_negatives.part-000000.parquet"
    initial_part_path_2 = trace_dir / "trace_initial_negatives.part-000001.parquet"
    samples_part_path_1 = trace_dir / "trace_iteration_samples.part-000000.parquet"
    samples_part_path_2 = trace_dir / "trace_iteration_samples.part-000001.parquet"
    for path in (
        initial_part_path_1,
        initial_part_path_2,
        samples_part_path_1,
        samples_part_path_2,
    ):
        path.touch()

    parquet_frames = {
        initial_part_path_1.name: initial_part_1,
        initial_part_path_2.name: initial_part_2,
        samples_part_path_1.name: samples_part_1,
        samples_part_path_2.name: samples_part_2,
    }
    read_order: list[str] = []

    def fake_read_parquet(path: Path | str) -> pd.DataFrame:
        filename = Path(path).name
        read_order.append(filename)
        return parquet_frames[filename].copy()

    monkeypatch.setattr(pd, "read_parquet", fake_read_parquet)

    sampling_trace = PredictionSamplingTrace.from_trace_directory(trace_dir)
    ensemble_trace = sampling_trace.get_ensemble_override("KINASE_A", 1)

    assert read_order == [
        initial_part_path_1.name,
        initial_part_path_2.name,
        samples_part_path_1.name,
        samples_part_path_2.name,
    ]
    assert ensemble_trace is not None
    assert ensemble_trace.initial_negative_sites == [
        "SITE_8",
        "SITE_7",
        "SITE_6",
        "SITE_5",
    ]
    assert ensemble_trace.iteration_sample_sites == {
        1: {1: ["SITE_4"], 2: ["SITE_7", "SITE_8"]}
    }


def test_sampling_trace_directory_reports_missing_columns_for_specific_parquet_part(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    trace_dir = tmp_path / "prediction_trace"
    trace_dir.mkdir(parents=True, exist_ok=True)

    good_part = trace_dir / "trace_initial_negatives.part-000000.parquet"
    bad_part = trace_dir / "trace_initial_negatives.part-000001.parquet"
    good_part.touch()
    bad_part.touch()

    def fake_read_parquet(path: Path | str) -> pd.DataFrame:
        filename = Path(path).name
        if filename == good_part.name:
            return pd.DataFrame(
                [{"kinase": "KINASE_A", "ensemble": 1, "draw": 1, "site": "SITE_8"}]
            )
        if filename == bad_part.name:
            return pd.DataFrame([{"kinase": "KINASE_A", "ensemble": 1, "draw": 2}])
        raise AssertionError(f"unexpected parquet path: {filename}")

    monkeypatch.setattr(pd, "read_parquet", fake_read_parquet)

    with pytest.raises(TableSchemaError, match=bad_part.name):
        PredictionSamplingTrace.from_trace_directory(trace_dir)


def test_predict_raises_for_invalid_sampling_trace_sites(tmp_path: Path) -> None:
    predictor = KinasePredictor()
    trace_dir = tmp_path / "prediction_trace"
    trace_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        [
            {"kinase": "KINASE_A", "ensemble": 1, "draw": 1, "site": "SITE_8"},
            {"kinase": "KINASE_A", "ensemble": 1, "draw": 2, "site": "SITE_7"},
            {"kinase": "KINASE_A", "ensemble": 1, "draw": 3, "site": "SITE_6"},
            {"kinase": "KINASE_A", "ensemble": 1, "draw": 4, "site": "SITE_404"},
        ]
    ).to_csv(trace_dir / "trace_initial_negatives.csv", index=False)

    with pytest.raises(ValueError, match="outside the available training rows"):
        predictor.predict(
            combined_scores=make_combined_scores(),
            ensemble_size=1,
            top=4,
            score_threshold=0.85,
            inclusion=3,
            n_iterations=2,
            random_state=5,
            sampling_trace=trace_dir,
        )


def test_transform_resampling_probabilities_flattens_default_mode() -> None:
    values = np.array([0.9, 0.1], dtype=float)

    transformed = _transform_resampling_probabilities(values, svm_mode="default")
    normalized_input = values / values.sum()
    normalized_transformed = transformed / transformed.sum()

    assert normalized_transformed[0] < normalized_input[0]
    assert normalized_transformed[1] > normalized_input[1]
    assert normalized_transformed[0] > normalized_transformed[1]


def test_transform_resampling_probabilities_keeps_r_parity_weights() -> None:
    values = np.array([0.9, 0.1], dtype=float)

    transformed = _transform_resampling_probabilities(values, svm_mode="r_parity")

    assert np.allclose(transformed, values)


def test_resolve_prediction_sampling_policy_maps_public_modes() -> None:
    default_policy = resolve_prediction_sampling_policy("default")
    r_parity_policy = resolve_prediction_sampling_policy("r_parity")

    assert default_policy.name == "default"
    assert default_policy.seed_strategy == "stable_by_kinase"
    assert default_policy.resampling_weight_mode == "default"
    assert default_policy.final_score_mode == "mean_probability"
    assert r_parity_policy.name == "r_parity"
    assert r_parity_policy.seed_strategy == "global_parity"
    assert r_parity_policy.resampling_weight_mode == "r_parity"
    assert r_parity_policy.final_score_mode == "decision_sigmoid"


def test_resolve_final_score_series_keeps_default_probability_scores() -> None:
    policy = resolve_prediction_sampling_policy("default")
    pred_df = pd.DataFrame(
        {"1": [0.2, 0.8], "2": [0.8, 0.2]},
        index=["SITE_1", "SITE_2"],
    )
    decision_values = pd.Series([0.0, 2.0], index=pred_df.index, dtype=float)

    resolved = _resolve_final_score_series(
        pred_df=pred_df,
        final_decision_values=decision_values,
        sampling_policy=policy,
    )

    pd.testing.assert_series_equal(resolved, pred_df.loc[:, "1"], check_names=False)


def test_resolve_final_score_series_uses_decision_sigmoid_in_r_parity() -> None:
    policy = resolve_prediction_sampling_policy("r_parity")
    pred_df = pd.DataFrame(
        {"1": [0.2, 0.8], "2": [0.8, 0.2]},
        index=["SITE_1", "SITE_2"],
    )
    decision_values = pd.Series([-2.0, 1.0], index=pred_df.index, dtype=float)

    resolved = _resolve_final_score_series(
        pred_df=pred_df,
        final_decision_values=decision_values,
        sampling_policy=policy,
    )

    expected = pd.Series(
        1.0 / (1.0 + np.exp(-decision_values.to_numpy(dtype=float))),
        index=pred_df.index,
        dtype=float,
    )
    pd.testing.assert_series_equal(resolved, expected, check_names=False)
    assert ((resolved >= 0.0) & (resolved <= 1.0)).all()


def test_prediction_sampling_random_source_uses_global_parity_call_order() -> None:
    policy = resolve_prediction_sampling_policy("r_parity")
    source_a = PredictionSamplingRandomSource(policy=policy, random_state=17)
    source_b = PredictionSamplingRandomSource(policy=policy, random_state=17)

    first_a = source_a.generators_for_kinase(kinase="KINASE_A")
    second_a = source_a.generators_for_kinase(kinase="KINASE_B")
    first_b = source_b.generators_for_kinase(kinase="KINASE_A")
    second_b = source_b.generators_for_kinase(kinase="KINASE_B")

    assert int(first_a[0].integers(0, 1000)) == int(first_b[0].integers(0, 1000))
    assert int(first_a[1].integers(0, 1000)) == int(first_b[1].integers(0, 1000))
    assert int(second_a[0].integers(0, 1000)) == int(second_b[0].integers(0, 1000))
    assert int(second_a[1].integers(0, 1000)) == int(second_b[1].integers(0, 1000))


def test_make_prediction_random_generators_returns_independent_streams() -> None:
    negative_rng, resampling_rng = _make_prediction_random_generators(
        np.random.default_rng(17)
    )

    first_negative_draw = int(negative_rng.integers(0, 1000))
    first_resampling_draw = int(resampling_rng.integers(0, 1000))

    assert first_negative_draw != first_resampling_draw


def test_resolve_svm_probability_random_state_returns_fixed_seed() -> None:
    assert _resolve_svm_probability_random_state() == 1


def test_kinase_predictor_is_reproducible_for_same_random_state() -> None:
    predictor = KinasePredictor()

    result_a = predictor.predict(
        combined_scores=make_combined_scores(),
        ensemble_size=3,
        top=4,
        score_threshold=0.85,
        inclusion=3,
        n_iterations=2,
        random_state=11,
    )
    result_b = predictor.predict(
        combined_scores=make_combined_scores(),
        ensemble_size=3,
        top=4,
        score_threshold=0.85,
        inclusion=3,
        n_iterations=2,
        random_state=11,
    )

    pd.testing.assert_frame_equal(result_a.pred_matrix, result_b.pred_matrix)


def test_kinase_predictor_is_invariant_to_kinase_column_order_for_fixed_seed() -> None:
    predictor = KinasePredictor()
    combined_scores = make_combined_scores()

    result_original = predictor.predict(
        combined_scores=combined_scores,
        ensemble_size=3,
        top=4,
        score_threshold=0.85,
        inclusion=3,
        n_iterations=2,
        random_state=11,
    )
    result_reordered = predictor.predict(
        combined_scores=combined_scores.loc[:, ["KINASE_B", "KINASE_A"]],
        ensemble_size=3,
        top=4,
        score_threshold=0.85,
        inclusion=3,
        n_iterations=2,
        random_state=11,
    )

    pd.testing.assert_frame_equal(
        result_original.pred_matrix,
        result_reordered.pred_matrix.reindex(
            columns=result_original.pred_matrix.columns
        ),
    )


def test_create_trace_sink_owns_and_cleans_up_temp_directory() -> None:
    sink = create_trace_sink(None, fmt="csv")
    assert isinstance(sink, DirectoryTraceSink)
    output_dir = sink.output_dir

    sink.write_rows(
        "trace_iteration_samples",
        [
            {
                "kinase": "KINASE_A",
                "ensemble": 1,
                "iteration": 1,
                "class_label": 1,
                "draw": 1,
                "site": "SITE_1",
            }
        ],
    )

    assert output_dir.exists()
    sink.close()

    assert not output_dir.exists()


def test_predict_request_closes_owned_runtime_trace_sink_on_failure() -> None:
    predictor = KinasePredictor()
    request = PredictionRequest.validate_request(
        combined_scores=make_combined_scores(),
        ensemble_size=2,
        top=4,
        score_threshold=0.85,
        inclusion=3,
        n_iterations=2,
        random_state=5,
        debug_top_n=3,
        trace_level="full",
        default_svm_mode="default",
        capture_debug_trace=True,
        debug_kinases=["KINASE_A"],
    )

    created_output_dirs: list[Path] = []
    original_create_state = predictor.trace_recorder.create_state

    def create_state_and_fail(*args, **kwargs):
        original_create_state(*args, **kwargs)
        trace_sink = kwargs["trace_sink"]
        assert trace_sink is not None
        created_output_dirs.append(trace_sink.output_dir)
        raise RuntimeError("boom")

    predictor.trace_recorder.create_state = create_state_and_fail  # type: ignore[method-assign]

    with pytest.raises(RuntimeError, match="boom"):
        predictor.predict_request(request)

    assert len(created_output_dirs) == 1
    assert not created_output_dirs[0].exists()


class _TrackingTraceSink(TraceSink):
    def __init__(self) -> None:
        self.closed = False

    def write_rows(self, table_name: str, rows: list[dict[str, object]]) -> None:
        return None

    def read_table(self, table_name: str) -> pd.DataFrame:
        return pd.DataFrame()

    def close(self) -> None:
        self.closed = True


def test_prediction_result_close_does_not_close_caller_supplied_trace_sink() -> None:
    predictor = KinasePredictor()
    sink = _TrackingTraceSink()

    result = predictor.predict(
        combined_scores=make_combined_scores(),
        ensemble_size=1,
        top=4,
        score_threshold=0.85,
        inclusion=3,
        n_iterations=2,
        random_state=5,
        capture_debug_trace=True,
        debug_kinases=["KINASE_A"],
        debug_top_n=3,
        trace_level="full",
        trace_sink=sink,
    )

    assert result.trace_sink is sink
    assert result.owns_trace_sink is False

    result.close()

    assert sink.closed is False


def test_predict_request_does_not_close_caller_supplied_trace_sink_on_failure() -> None:
    predictor = KinasePredictor()
    sink = _TrackingTraceSink()
    request = PredictionRequest.validate_request(
        combined_scores=make_combined_scores(),
        ensemble_size=2,
        top=4,
        score_threshold=0.85,
        inclusion=3,
        n_iterations=2,
        random_state=5,
        debug_top_n=3,
        trace_level="full",
        trace_sink=sink,
        default_svm_mode="default",
        capture_debug_trace=True,
        debug_kinases=["KINASE_A"],
    )

    def create_state_and_fail(*args, **kwargs):
        raise RuntimeError("boom")

    predictor.trace_recorder.create_state = create_state_and_fail  # type: ignore[method-assign]

    with pytest.raises(RuntimeError, match="boom"):
        predictor.predict_request(request)

    assert sink.closed is False


def test_prediction_result_close_is_idempotent_for_owned_trace_sink() -> None:
    predictor = KinasePredictor()

    result = predictor.predict(
        combined_scores=make_combined_scores(),
        ensemble_size=1,
        top=4,
        score_threshold=0.85,
        inclusion=3,
        n_iterations=2,
        random_state=5,
        capture_debug_trace=True,
        debug_kinases=["KINASE_A"],
        debug_top_n=3,
        trace_level="full",
        trace_sink=None,
    )

    assert result.trace_sink is not None
    assert result.owns_trace_sink is True
    output_dir = result.trace_sink.output_dir
    assert output_dir.exists()

    result.close()
    result.close()

    assert not output_dir.exists()


def test_prediction_result_context_manager_closes_owned_trace_sink() -> None:
    predictor = KinasePredictor()

    with predictor.predict(
        combined_scores=make_combined_scores(),
        ensemble_size=1,
        top=4,
        score_threshold=0.85,
        inclusion=3,
        n_iterations=2,
        random_state=5,
        capture_debug_trace=True,
        debug_kinases=["KINASE_A"],
        debug_top_n=3,
        trace_level="full",
        trace_sink=None,
    ) as result:
        assert result.trace_sink is not None
        assert result.owns_trace_sink is True
        output_dir = result.trace_sink.output_dir
        assert output_dir.exists()

    assert not output_dir.exists()
    assert result.owns_trace_sink is False


def test_prediction_result_does_not_define_destructor_cleanup() -> None:
    assert "__del__" not in KinasePredictionResult.__dict__


def test_predict_rejects_invalid_request_at_boundary() -> None:
    predictor = KinasePredictor()

    with pytest.raises(RequestValidationError, match="ensemble_size"):
        predictor.predict(
            combined_scores=make_combined_scores(),
            ensemble_size=0,
        )


def test_predict_rejects_invalid_trace_sink_format_at_boundary() -> None:
    predictor = KinasePredictor()

    with pytest.raises(RequestValidationError, match="trace_sink_format"):
        predictor.predict(
            combined_scores=make_combined_scores(),
            ensemble_size=1,
            top=4,
            score_threshold=0.85,
            inclusion=3,
            n_iterations=1,
            trace_level="full",
            trace_sink_format="json",  # type: ignore[arg-type]
        )


def test_predict_rejects_trace_sink_without_full_trace_level_at_boundary(
    tmp_path: Path,
) -> None:
    predictor = KinasePredictor()

    with pytest.raises(
        RequestValidationError,
        match="trace_sink may only be provided when trace_level='full'",
    ):
        predictor.predict(
            combined_scores=make_combined_scores(),
            ensemble_size=1,
            top=4,
            score_threshold=0.85,
            inclusion=3,
            n_iterations=1,
            trace_level="summary",
            trace_sink=tmp_path / "trace_output",
        )


def test_predict_request_uses_prevalidated_request_without_revalidating_candidate_scalars() -> (
    None
):
    predictor = KinasePredictor()
    request = predictor.request_factory.create(
        combined_scores=make_combined_scores(),
        ensemble_size=1,
        top=4,
        score_threshold=0.85,
        inclusion=3,
        n_iterations=1,
        random_state=5,
        capture_debug_trace=False,
        debug_kinases=None,
        debug_top_n=3,
        svm_mode=None,
        sampling_trace=None,
        trace_level=None,
        trace_sink=None,
        trace_sink_format="csv",
    )

    original = predictor.candidate_selector.select

    def select_once(
        combined_scores: pd.DataFrame,
        *,
        top: int,
        score_threshold: float,
        inclusion: int,
    ) -> dict[str, list[str]]:
        assert top == 4
        assert inclusion == 3
        return original(
            combined_scores,
            top=top,
            score_threshold=score_threshold,
            inclusion=inclusion,
        )

    predictor.candidate_selector.select = select_once  # type: ignore[method-assign]
    result = predictor.predict_request(request)

    assert isinstance(result, KinasePredictionResult)


def test_predict_rejects_non_numeric_combined_scores_at_boundary() -> None:
    predictor = KinasePredictor()

    with pytest.raises(RequestValidationError, match="combined_scores"):
        predictor.predict(
            combined_scores=pd.DataFrame({"KINASE_A": ["bad"]}, index=["SITE_1"]),
            ensemble_size=1,
            top=1,
            score_threshold=0.8,
            inclusion=1,
            n_iterations=1,
        )


def test_build_candidate_substrate_list_rejects_invalid_combined_scores_cleanly() -> (
    None
):
    with pytest.raises(ValueError, match="combined_scores"):
        build_candidate_substrate_list(
            pd.DataFrame({"KINASE_A": ["bad"]}, index=["SITE_1"]),
            top=1,
            score_threshold=0.8,
            inclusion=1,
        )


def test_prediction_public_api_uses_explicit_ensemble_predictor_contract() -> None:
    assert (
        get_type_hints(KinasePredictor.__init__)["ensemble_predictor"]
        == EnsemblePredictorContract | None
    )


def test_prediction_public_api_has_concrete_result_return_types() -> None:
    assert get_type_hints(KinasePredictor.predict)["return"] == KinasePredictionResult
    assert (
        get_type_hints(KinasePredictor.predict_request)["return"]
        == KinasePredictionResult
    )
    assert (
        get_type_hints(KinasePredictor.predict_from_scoring_result)["return"]
        == KinasePredictionResult
    )


def test_simple_kinase_workflow_runs_from_fixture_files_and_pins_result_invariants() -> (
    None
):
    expected_site_ids = [
        "EIF4B;S422;",
        "GSK3B;S9;",
        "TBC1D1;S231;",
        "TBC1D1;T590;",
        "TSC2;S939;",
    ]
    expected_kinases = [
        "AKT1",
        "AKT3",
        "MAP4K5",
        "PRKAA1",
        "PRKACA",
        "RPS6KA1",
        "RPS6KB1",
        "Yang.S6K",
    ]

    with SimpleKinaseWorkflow(flank_size=7).run(
        total=SIMPLE_WORKFLOW_FIXTURE_DIR / "total.tsv",
        phospho=SIMPLE_WORKFLOW_FIXTURE_DIR / "phospho.tsv",
        species="rat",
        prediction_config=PredictionRunConfig(
            min_substrates=1,
            min_motif_size=1,
            ensemble_size=2,
            top=3,
            inclusion=2,
            n_iterations=2,
            random_state=7,
        ),
        activity_config=KinaseActivityConfig(
            threshold=0.1,
            min_substrates=1,
            top_n_substrates=3,
        ),
    ) as result:
        pred_mat = result.pred_mat_result.to_frame(copy=False)

        assert (
            result.analysis_ready_dataset.phospho_matrix.index.tolist()
            == expected_site_ids
        )
        assert result.pred_mat_result.phosphosite_ids.tolist() == expected_site_ids
        assert result.pred_mat_result.kinase_names.tolist() == expected_kinases
        assert pred_mat.shape == (5, 8)
        assert pred_mat.index.tolist() == expected_site_ids
        assert pred_mat.columns.tolist() == expected_kinases
        assert result.reference_bundle.species == "rat"
        assert result.reference_bundle.source_metadata.reference == "l6_native"
        assert set(result.reference_bundle.substrate_map).issuperset(expected_kinases)
        assert result.kinase_activity_result.weighted_activity.columns.tolist() == (
            result.analysis_ready_dataset.phospho_matrix.columns.tolist()
        )
        assert result.kinase_activity_result.ksea_scores.columns.tolist() == (
            result.analysis_ready_dataset.phospho_matrix.columns.tolist()
        )
        assert (
            result.kinase_activity_result.weighted_activity.index.tolist()
            == expected_kinases
        )
        assert set(result.kinase_activity_result.ksea_counts.index.tolist()) == set(
            expected_kinases
        )
        assert set(result.kinase_activity_result.target_counts.index.tolist()) == set(
            expected_kinases
        )
        assert result.kinase_activity_result.ksea_counts.is_monotonic_decreasing
        assert result.kinase_activity_result.target_counts.is_monotonic_decreasing
        assert set(result.kinase_activity_result.target_table["site_id"]) == set(
            expected_site_ids
        )
        assert set(result.kinase_activity_result.target_table["kinase"]) == set(
            expected_kinases
        )
