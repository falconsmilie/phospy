from __future__ import annotations

import pandas as pd
import pytest

from phospy import KinaseActivityAnalyzer
from phospy.validation.errors import TableSchemaError


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

    expected_files = {
        "kinase_activity_matrix.csv",
        "kinase_target_counts.csv",
        "kinase_target_table.csv",
        "ksea_counts.csv",
        "ksea_scores.csv",
    }
    assert expected_files == {path.name for path in outdir.iterdir()}
