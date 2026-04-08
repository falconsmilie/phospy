from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from phospy import KinaseActivityAnalyzer, PredMatResult
from phospy.constants import KINASE_OUTPUT_FILENAMES
from phospy.io import load_pred_mat
from phospy.validation.errors import (
    NoCandidateKinasesError,
    RequestValidationError,
    TableSchemaError,
)
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
        load_pred_mat(pred_mat_path)


def test_analyzer_run_returns_expected_result() -> None:
    result = KinaseActivityAnalyzer().run(
        pred_mat=make_pred_mat(),
        phospho_matrix=make_phospho_matrix(),
        threshold=0.6,
        min_substrates=2,
    )

    assert set(result.weighted_activity.index) == {"PRKACA", "BTK"}
    assert set(result.ksea_scores.index) == {"PRKACA", "BTK"}
    assert result.ksea_counts.to_dict() == {"PRKACA": 3, "BTK": 2}
    assert result.target_counts.to_dict() == {"PRKACA": 3, "BTK": 2}


def test_kinase_activity_result_tables_are_detached_snapshots_from_inputs() -> None:
    pred_mat = make_pred_mat()
    phospho_matrix = make_phospho_matrix()

    result = KinaseActivityAnalyzer().run(
        pred_mat=pred_mat,
        phospho_matrix=phospho_matrix,
        threshold=0.6,
        min_substrates=2,
    )

    original_pred_mat = pred_mat.copy(deep=True)
    original_phospho_matrix = phospho_matrix.copy(deep=True)

    result.weighted_activity.loc["PRKACA", "phospho_corrected_1"] = -999.0
    result.ksea_scores.loc["PRKACA", "ksea_score"] = -999.0
    result.target_table.loc[result.target_table.index[0], "kinase"] = "CHANGED"

    pd.testing.assert_frame_equal(pred_mat, original_pred_mat)
    pd.testing.assert_frame_equal(phospho_matrix, original_phospho_matrix)


def test_analyzer_run_runs_end_to_end_with_loaded_pred_mat(tmp_path) -> None:
    pred_mat_path = tmp_path / "predMat.csv"
    make_pred_mat().to_csv(pred_mat_path)

    result = KinaseActivityAnalyzer().run(
        pred_mat=load_pred_mat(pred_mat_path),
        phospho_matrix=make_phospho_matrix(),
        threshold=0.6,
        min_substrates=2,
    )

    assert set(result.weighted_activity.index) == {"PRKACA", "BTK"}
    assert set(result.target_table["kinase"]) == {"PRKACA", "BTK"}


def test_analyzer_load_pred_mat_uses_injected_loader(tmp_path) -> None:
    pred_mat_path = tmp_path / "predMat.csv"
    expected = make_pred_mat()
    calls: list[Path] = []

    def fake_loader(path: str | Path) -> pd.DataFrame:
        resolved = Path(path)
        calls.append(resolved)
        return expected

    analyzer = KinaseActivityAnalyzer(pred_mat_loader=fake_loader)

    loaded = analyzer.load_pred_mat(pred_mat_path)

    assert loaded is expected
    assert calls == [pred_mat_path]


def test_analyzer_run_uses_injected_runner() -> None:
    expected = KinaseActivityAnalyzer().run(
        pred_mat=make_pred_mat(),
        phospho_matrix=make_phospho_matrix(),
        threshold=0.6,
        min_substrates=2,
    )

    class StubRunner:
        def __init__(self) -> None:
            self.request = None

        def execute(self, request):
            self.request = request
            return expected

    runner = StubRunner()
    analyzer = KinaseActivityAnalyzer(runner=runner)

    result = analyzer.run(
        pred_mat=make_pred_mat(),
        phospho_matrix=make_phospho_matrix(),
        threshold=0.7,
        min_substrates=4,
        top_n_substrates=9,
    )

    assert result is expected
    assert runner.request is not None
    assert runner.request.request.threshold == 0.7
    assert runner.request.request.min_substrates == 4
    assert runner.request.request.top_n_substrates == 9


def test_analyzer_write_outputs_uses_injected_writer(tmp_path) -> None:
    result = KinaseActivityAnalyzer().run(
        pred_mat=make_pred_mat(),
        phospho_matrix=make_phospho_matrix(),
        threshold=0.6,
        min_substrates=2,
    )

    class StubWriter:
        def __init__(self) -> None:
            self.calls: list[tuple[object, Path]] = []

        def write(self, written_result, outdir: str | Path) -> None:
            self.calls.append((written_result, Path(outdir)))

    writer = StubWriter()
    analyzer = KinaseActivityAnalyzer(result_writer=writer)
    outdir = tmp_path / "kinase-output"

    analyzer.write_outputs(result, outdir)

    assert writer.calls == [(result, outdir)]


def test_analyzer_write_outputs_writes_expected_files(tmp_path) -> None:
    analyzer = KinaseActivityAnalyzer()
    result = analyzer.run(
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
        KinaseActivityAnalyzer().run(
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


def test_run_revalidates_public_inputs_loaded_from_disk(
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

    result = KinaseActivityAnalyzer().run(
        pred_mat=load_pred_mat(pred_mat_path),
        phospho_matrix=make_phospho_matrix(),
        threshold=0.6,
        min_substrates=2,
    )

    assert set(result.weighted_activity.index) == {"PRKACA", "BTK"}
    assert pred_calls == [f"pred_mat ({pred_mat_path})", "pred_mat"]
    assert matrix_calls == ["phospho_matrix"]


def test_run_request_uses_validated_boundary_request(
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
    request = analyzer._validate_request(
        pred_mat=make_pred_mat(),
        phospho_matrix=make_phospho_matrix(),
        threshold=0.6,
        min_substrates=2,
        top_n_substrates=20,
    )
    result = analyzer._run_request(request)

    assert set(result.weighted_activity.index) == {"PRKACA", "BTK"}
    assert pred_calls == ["pred_mat"]
    assert matrix_calls == ["phospho_matrix"]


def test_analyzer_rejects_pred_mat_without_candidate_kinases() -> None:
    analyzer = KinaseActivityAnalyzer()
    phospho_matrix = make_phospho_matrix()
    empty_pred_mat = make_pred_mat().iloc[:, 0:0]

    with pytest.raises(
        NoCandidateKinasesError,
        match=(
            "pred_mat does not contain any kinase columns because no candidate "
            "kinases qualified for prediction"
        ),
    ):
        analyzer.run(
            pred_mat=empty_pred_mat,
            phospho_matrix=phospho_matrix,
        )


def test_analyzer_passes_pred_mat_result_to_validation_boundary(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured: list[object] = []
    original_validate_request = KinaseActivityAnalyzer._validate_request.__globals__[
        "validate_analysis_request"
    ]

    def capturing_validate_request(*, pred_mat, phospho_matrix, **kwargs):
        captured.append(pred_mat)
        return original_validate_request(
            pred_mat=pred_mat,
            phospho_matrix=phospho_matrix,
            **kwargs,
        )

    monkeypatch.setitem(
        KinaseActivityAnalyzer._validate_request.__globals__,
        "validate_analysis_request",
        capturing_validate_request,
    )

    pred_mat_result = PredMatResult(make_pred_mat())
    result = KinaseActivityAnalyzer().run(
        pred_mat=pred_mat_result,
        phospho_matrix=make_phospho_matrix(),
        threshold=0.6,
        min_substrates=2,
    )

    assert set(result.weighted_activity.index) == {"PRKACA", "BTK"}
    assert captured == [pred_mat_result]
