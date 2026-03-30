from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from phospy.prediction import (
    KinasePredictionResult,
    KinasePredictor,
    PredictionSamplingTrace,
    build_candidate_substrate_list,
    prediction_debug_trace_tables,
)
from phospy.prediction.sampling import (
    make_prediction_random_generators as _make_prediction_random_generators,
)
from phospy.prediction.sampling import multi_ada_sampling as _multi_ada_sampling
from phospy.prediction.sampling import (
    transform_resampling_probabilities as _transform_resampling_probabilities,
)
from phospy.prediction.svm import _RLikeStandardScaler
from phospy.prediction.svm import make_svm as _make_svm
from phospy.prediction.svm import require_sklearn as _require_sklearn
from phospy.prediction.svm import (
    resolve_svm_probability_random_state as _resolve_svm_probability_random_state,
)
from phospy.prediction.traces import DirectoryTraceSink
from phospy.scoring import KinaseScoringResult


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


def test_predict_returns_empty_matrix_when_no_kinases_pass_inclusion() -> None:
    predictor = KinasePredictor()
    result = predictor.predict(
        combined_scores=make_combined_scores(),
        ensemble_size=2,
        top=2,
        score_threshold=0.99,
        inclusion=2,
        n_iterations=2,
        random_state=3,
    )

    assert result.substrate_list == {}
    assert result.pred_matrix.empty
    assert list(result.pred_matrix.index) == list(make_combined_scores().index)


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
