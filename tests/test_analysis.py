from __future__ import annotations

import pandas as pd
import pytest

from phospy import KinaseActivityAnalyzer
from phospy.constants import KINASE_OUTPUT_FILENAMES
from phospy.validation.errors import RequestValidationError, TableSchemaError
from phospy.validation.requests import KinaseActivityRequest
from phospy.validation.tables import PredMatSchema, SiteMatrixSchema


def make_pred_mat() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "PRKACA": [0.9, 0.8, 0.7],
            "BTK": [0.2, 0.85, 0.75],
        },
        index=["PRKACA;S339;", "BTK;Y551;", "LYN;Y397;"],
    )


def make_phospho_matrix() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "phospho_corrected_1": [4.0, 4.0, 4.0],
            "phospho_corrected_2": [5.0, 5.0, 5.0],
        },
        index=["PRKACA;S339;", "BTK;Y551;", "LYN;Y397;"],
    )


def test_analyzer_load_pred_mat_validates_input(tmp_path) -> None:
    pred_mat_path = tmp_path / "predMat.csv"
    pd.DataFrame(
        {
            "PRKACA": [1.2],
            "BTK": [0.8],
        },
        index=["PRKACA;S339;"],
    ).to_csv(pred_mat_path)

    with pytest.raises(TableSchemaError, match="outside the allowed range"):
        KinaseActivityAnalyzer().load_pred_mat(pred_mat_path)


def test_analyzer_analyze_returns_expected_result() -> None:
    result = KinaseActivityAnalyzer().analyze(
        pred_mat=make_pred_mat(),
        phospho_matrix=make_phospho_matrix(),
        threshold=0.6,
        min_substrates=2,
    )

    assert set(result.weighted_activity.index) == {"PRKACA", "BTK"}
    assert set(result.ksea_scores.index) == {"PRKACA", "BTK"}
    assert result.ksea_counts.to_dict() == {"PRKACA": 3, "BTK": 2}
    assert result.target_counts.to_dict() == {"PRKACA": 3, "BTK": 2}


def test_analyzer_load_and_analyze_runs_end_to_end(tmp_path) -> None:
    pred_mat_path = tmp_path / "predMat.csv"
    make_pred_mat().to_csv(pred_mat_path)

    result = KinaseActivityAnalyzer().load_and_analyze(
        pred_mat_path=pred_mat_path,
        phospho_matrix=make_phospho_matrix(),
        threshold=0.6,
        min_substrates=2,
    )

    assert set(result.weighted_activity.index) == {"PRKACA", "BTK"}
    assert set(result.target_table["kinase"]) == {"PRKACA", "BTK"}


def test_analyzer_write_outputs_writes_expected_files(tmp_path) -> None:
    analyzer = KinaseActivityAnalyzer()
    result = analyzer.analyze(
        pred_mat=make_pred_mat(),
        phospho_matrix=make_phospho_matrix(),
        threshold=0.6,
        min_substrates=2,
    )

    outdir = tmp_path / "kinase-output"
    analyzer.write_outputs(result, outdir)

    expected_files = set(KINASE_OUTPUT_FILENAMES)
    assert expected_files == {path.name for path in outdir.iterdir()}


def test_analyzer_rejects_invalid_request_threshold() -> None:
    with pytest.raises(RequestValidationError, match="threshold"):
        KinaseActivityAnalyzer().analyze(
            pred_mat=make_pred_mat(),
            phospho_matrix=make_phospho_matrix(),
            threshold=1.5,
        )


def test_analyzer_load_pred_mat_validates_through_loader_once(
    tmp_path, monkeypatch: pytest.MonkeyPatch
) -> None:
    pred_mat_path = tmp_path / "predMat.csv"
    make_pred_mat().to_csv(pred_mat_path)

    calls: list[str] = []
    original_validate = PredMatSchema.validate

    def counting_validate(df: pd.DataFrame, *, context: str) -> pd.DataFrame:
        calls.append(context)
        return original_validate(df, context=context)

    monkeypatch.setattr(PredMatSchema, "validate", staticmethod(counting_validate))

    loaded = KinaseActivityAnalyzer().load_pred_mat(pred_mat_path)

    assert loaded.equals(make_pred_mat())
    assert calls == [f"pred_mat ({pred_mat_path})"]


def test_load_and_analyze_does_not_revalidate_loaded_prediction_matrix(
    tmp_path, monkeypatch: pytest.MonkeyPatch
) -> None:
    pred_mat_path = tmp_path / "predMat.csv"
    make_pred_mat().to_csv(pred_mat_path)

    pred_calls: list[str] = []
    matrix_calls: list[str] = []
    original_pred_validate = PredMatSchema.validate
    original_matrix_validate = SiteMatrixSchema.validate

    def counting_pred_validate(df: pd.DataFrame, *, context: str) -> pd.DataFrame:
        pred_calls.append(context)
        return original_pred_validate(df, context=context)

    def counting_matrix_validate(df: pd.DataFrame, *, context: str) -> pd.DataFrame:
        matrix_calls.append(context)
        return original_matrix_validate(df, context=context)

    monkeypatch.setattr(
        PredMatSchema,
        "validate",
        staticmethod(counting_pred_validate),
    )
    monkeypatch.setattr(
        SiteMatrixSchema,
        "validate",
        staticmethod(counting_matrix_validate),
    )

    result = KinaseActivityAnalyzer().load_and_analyze(
        pred_mat_path=pred_mat_path,
        phospho_matrix=make_phospho_matrix(),
        threshold=0.6,
        min_substrates=2,
    )

    assert set(result.weighted_activity.index) == {"PRKACA", "BTK"}
    assert pred_calls == [f"pred_mat ({pred_mat_path})"]
    assert matrix_calls == ["phospho_matrix"]


def test_analyze_request_uses_validated_boundary_request(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    pred_calls: list[str] = []
    matrix_calls: list[str] = []
    original_pred_validate = PredMatSchema.validate
    original_matrix_validate = SiteMatrixSchema.validate

    def counting_pred_validate(df: pd.DataFrame, *, context: str) -> pd.DataFrame:
        pred_calls.append(context)
        return original_pred_validate(df, context=context)

    def counting_matrix_validate(df: pd.DataFrame, *, context: str) -> pd.DataFrame:
        matrix_calls.append(context)
        return original_matrix_validate(df, context=context)

    monkeypatch.setattr(
        PredMatSchema,
        "validate",
        staticmethod(counting_pred_validate),
    )
    monkeypatch.setattr(
        SiteMatrixSchema,
        "validate",
        staticmethod(counting_matrix_validate),
    )

    analyzer = KinaseActivityAnalyzer()
    request = analyzer.validate_request(
        pred_mat=make_pred_mat(),
        phospho_matrix=make_phospho_matrix(),
        threshold=0.6,
        min_substrates=2,
        top_n_substrates=20,
    )
    result = analyzer.analyze_request(request=request)

    assert set(result.weighted_activity.index) == {"PRKACA", "BTK"}
    assert pred_calls == ["pred_mat"]
    assert matrix_calls == ["phospho_matrix"]


def test_analyze_validated_request_rejects_raw_request_objects() -> None:
    request = KinaseActivityRequest.validate_request(
        threshold=0.6,
        min_substrates=2,
        top_n_substrates=20,
    )

    with pytest.raises(TypeError, match="ValidatedAnalysisRequest"):
        KinaseActivityAnalyzer.analyze_validated_request(
            request=request  # type: ignore[arg-type]
        )


def test_analyze_request_rejects_raw_request_objects() -> None:
    request = KinaseActivityRequest.validate_request(
        threshold=0.6,
        min_substrates=2,
        top_n_substrates=20,
    )

    with pytest.raises(TypeError, match="ValidatedAnalysisRequest"):
        KinaseActivityAnalyzer.analyze_request(request=request)  # type: ignore[arg-type]
