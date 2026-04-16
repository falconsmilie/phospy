from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from phospy.errors import (
    CustomPredictorOutputError,
    InputCompatibilityError,
    NoCandidateKinasesError,
    PhospyError,
    PhospyValidationError,
    PredictionConfigurationError,
    RequestValidationError,
    TableSchemaError,
)
from phospy.prediction import (
    KinasePredictor,
    KinaseScorer,
)
from phospy.prediction.candidates import build_candidate_substrate_list
from phospy.prediction.motif_scoring import KinaseMotifScorer, create_frequency_matrix
from phospy.prediction.traces import PredictionSamplingTrace


def test_motif_scorer_rejects_missing_motif_size_with_package_error() -> None:
    with pytest.raises(InputCompatibilityError, match="motif_sizes is missing entries"):
        KinaseMotifScorer(
            motif_frequency_matrices={
                "KINASE_A": create_frequency_matrix(["AAAAA"], flank_size=2)
            },
            motif_sizes=pd.Series(dtype=float),
            flank_size=2,
        )


def test_create_frequency_matrix_rejects_inconsistent_windows_with_package_error() -> (
    None
):
    with pytest.raises(TableSchemaError, match="same window length"):
        create_frequency_matrix(["AAAAA", "AAA"], flank_size=2)


def test_kinase_scorer_rejects_mismatched_columns_with_package_error() -> None:
    scorer = KinaseScorer(
        pd.DataFrame(
            {"sample_1": [1.0], "sample_2": [2.0], "sample_3": [3.0]},
            index=["KINASE_A"],
        )
    )
    phospho_matrix = pd.DataFrame(
        {"sample_1": [1.0], "sample_2": [2.0], "sample_x": [3.0]},
        index=["SITE_A"],
    )

    with pytest.raises(
        InputCompatibilityError, match="must match kinase profile columns"
    ):
        scorer.score_phosphosite_profiles(phospho_matrix)


def test_kinase_predictor_rejects_unknown_svm_mode_with_package_error() -> None:
    with pytest.raises(PhospyValidationError, match="svm_mode"):
        KinasePredictor(svm_mode="broken")  # type: ignore[arg-type]


def test_build_candidate_substrate_list_rejects_invalid_top_with_package_error() -> (
    None
):
    scores = pd.DataFrame({"KINASE_A": [0.95]}, index=["SITE_1"])

    with pytest.raises(PhospyValidationError, match="top must be at least 1"):
        build_candidate_substrate_list(scores, top=0)


def test_prediction_sampling_trace_rejects_missing_trace_files_with_package_error(
    tmp_path: Path,
) -> None:
    with pytest.raises(TableSchemaError, match="sampling trace directory"):
        PredictionSamplingTrace.from_trace_directory(tmp_path / "missing_trace")


def test_predict_rejects_missing_negative_pool_with_configuration_error() -> None:
    predictor = KinasePredictor()
    combined_scores = pd.DataFrame(
        {"KINASE_A": [0.95, 0.93]}, index=["SITE_1", "SITE_2"]
    )

    with pytest.raises(
        PredictionConfigurationError,
        match="No negative pool sites available to train predictor for KINASE_A",
    ):
        predictor.predict(
            combined_scores=combined_scores,
            ensemble_size=2,
            top=2,
            score_threshold=0.9,
            inclusion=1,
            n_iterations=2,
            random_state=3,
        )


def test_prediction_configuration_error_is_a_package_validation_error() -> None:
    assert issubclass(PredictionConfigurationError, PhospyValidationError)
    assert issubclass(PredictionConfigurationError, PhospyError)


def test_custom_predictor_output_error_is_a_package_validation_error() -> None:
    assert issubclass(CustomPredictorOutputError, InputCompatibilityError)
    assert issubclass(CustomPredictorOutputError, PhospyValidationError)
    assert issubclass(CustomPredictorOutputError, PhospyError)


def test_no_candidate_kinases_error_is_a_package_validation_error() -> None:
    assert issubclass(NoCandidateKinasesError, InputCompatibilityError)
    assert issubclass(NoCandidateKinasesError, PhospyValidationError)
    assert issubclass(NoCandidateKinasesError, PhospyError)


def test_predict_rejects_trace_sink_without_full_trace_level_with_package_error(
    tmp_path: Path,
) -> None:
    predictor = KinasePredictor()

    with pytest.raises(
        RequestValidationError,
        match="trace_sink may only be provided when trace_level='full'",
    ):
        predictor.predict(
            combined_scores=pd.DataFrame(
                {"KINASE_A": [0.95, 0.93]}, index=["SITE_1", "SITE_2"]
            ),
            ensemble_size=2,
            top=2,
            score_threshold=0.9,
            inclusion=1,
            n_iterations=2,
            random_state=3,
            trace_level="summary",
            trace_sink=tmp_path / "trace_output",
        )


def test_predict_from_scoring_result_requires_profile_fallback_with_package_error() -> (
    None
):
    predictor = KinasePredictor()
    scoring_result = type(
        "ScoringResultLike",
        (),
        {"combined_scores": None, "profile_scores": pd.DataFrame()},
    )()

    with pytest.raises(
        InputCompatibilityError, match="allow_profile_only_fallback=True"
    ):
        predictor.predict_from_scoring_result(scoring_result)  # type: ignore[arg-type]


def test_predict_reports_strict_threshold_candidate_shortfall_diagnostics() -> None:
    predictor = KinasePredictor()
    combined_scores = pd.DataFrame(
        {
            "K1": [0.95, 0.20, 0.10],
            "K2": [0.85, 0.84, 0.10],
        },
        index=["SITE_1", "SITE_2", "SITE_3"],
    )

    with pytest.raises(
        NoCandidateKinasesError,
        match="No candidate kinases qualified for prediction",
    ) as exc_info:
        predictor.predict(
            combined_scores=combined_scores,
            ensemble_size=1,
            top=2,
            score_threshold=0.8,
            inclusion=3,
            n_iterations=1,
            random_state=3,
        )
    message = str(exc_info.value)
    assert "Evaluated 2 kinase column(s) across 3 phosphosite row(s)." in message
    assert "Effective top window per kinase=2." in message
    assert "Kinases with at least one site above score_threshold: 2." in message
    assert "Best-support kinase had 2 qualifying site(s), below inclusion=3." in message
    assert "Near-miss kinases below inclusion: K2 (2), K1 (1)" in message


def test_profile_score_fallback_reports_considered_path_on_candidate_shortfall() -> (
    None
):
    predictor = KinasePredictor()
    scoring_result = type(
        "ScoringResultLike",
        (),
        {
            "combined_scores": None,
            "profile_scores": pd.DataFrame(
                {"KINASE_A": [0.95, 0.20]},
                index=["SITE_1", "SITE_2"],
            ),
        },
    )()

    with pytest.raises(
        NoCandidateKinasesError,
        match="No candidate kinases qualified for prediction",
    ) as exc_info:
        predictor.predict_from_scoring_result(
            scoring_result,  # type: ignore[arg-type]
            ensemble_size=1,
            top=1,
            score_threshold=0.9,
            inclusion=2,
            n_iterations=1,
            random_state=7,
            allow_profile_only_fallback=True,
        )
    assert (
        "Fallback path considered: profile_scores was used because combined_scores "
        "was unavailable."
    ) in str(exc_info.value)
