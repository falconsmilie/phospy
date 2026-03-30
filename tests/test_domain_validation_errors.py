from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from phospy.analysis import KinaseActivityAnalyzer
from phospy.motifs import KinaseMotifScorer, create_frequency_matrix
from phospy.prediction import (
    KinasePredictor,
    PredictionSamplingTrace,
    build_candidate_substrate_list,
)
from phospy.profiles import KinaseProfileBuilder
from phospy.scoring import KinaseScorer
from phospy.validation.errors import (
    InputCompatibilityError,
    PhospyValidationError,
    TableSchemaError,
)


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


def test_profile_builder_rejects_unknown_aggregation_with_package_error() -> None:
    with pytest.raises(PhospyValidationError, match="aggregation must be 'median'"):
        KinaseProfileBuilder(aggregation="mean")  # type: ignore[arg-type]


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


def test_kinase_activity_analyzer_rejects_invalid_threshold_with_package_error() -> (
    None
):
    analyzer = KinaseActivityAnalyzer(
        pred_mat=pd.DataFrame({"KINASE_A": [0.9]}, index=["SITE_1"])
    )

    with pytest.raises(
        PhospyValidationError, match="threshold must be between 0.0 and 1.0"
    ):
        analyzer.build_target_table(threshold=1.5)


def test_kinase_activity_analyzer_rejects_invalid_substrate_counts_with_package_error() -> (
    None
):
    analyzer = KinaseActivityAnalyzer(
        pred_mat=pd.DataFrame({"KINASE_A": [0.9]}, index=["SITE_1"])
    )
    phospho_matrix = pd.DataFrame({"sample_1": [1.0]}, index=["SITE_1"])

    with pytest.raises(
        PhospyValidationError, match="top_n_substrates must be at least 1"
    ):
        analyzer.compute_weighted_activity(
            phospho_matrix=phospho_matrix,
            top_n_substrates=0,
        )

    with pytest.raises(
        PhospyValidationError, match="min_substrates must be at least 1"
    ):
        analyzer.compute_ksea_scores(
            phospho_matrix=phospho_matrix,
            min_substrates=0,
        )
