from __future__ import annotations

import pandas as pd
import pytest

from phospy.api import (
    ContrastMatrix,
    DesignMatrix,
    DifferentialAnalysis,
    DifferentialAnalysisRequest,
    EmpiricalBayesConfig,
)
from phospy.errors import PhosPyInputError


def _matrix() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "A_1": [1.0, 2.0, 1.0, 1.2, 2.2],
            "A_2": [1.1, 2.1, 1.1, 1.1, 2.0],
            "B_1": [2.1, 2.0, 1.0, 1.3, 2.1],
            "B_2": [2.0, 2.2, 0.9, 1.2, 2.3],
            "C_1": [1.0, 2.0, 3.0, 1.5, 2.4],
            "C_2": [1.1, 2.1, 3.1, 1.4, 2.2],
        },
        index=pd.Index(
            ["SITE_1", "SITE_2", "SITE_3", "SITE_4", "SITE_5"],
            name="site_id",
        ),
    )


def _design() -> DesignMatrix:
    return DesignMatrix(
        pd.DataFrame(
            {
                "A": [1.0, 1.0, 0.0, 0.0, 0.0, 0.0],
                "B": [0.0, 0.0, 1.0, 1.0, 0.0, 0.0],
                "C": [0.0, 0.0, 0.0, 0.0, 1.0, 1.0],
            },
            index=pd.Index(["A_1", "A_2", "B_1", "B_2", "C_1", "C_2"], name="sample"),
        )
    )


def _contrasts() -> ContrastMatrix:
    return ContrastMatrix(
        pd.DataFrame(
            {
                "B_vs_A": [-1.0, 1.0, 0.0],
                "C_vs_A": [-1.0, 0.0, 1.0],
            },
            index=pd.Index(["A", "B", "C"], name="coefficient"),
        )
    )


def test_differential_analysis_returns_per_contrast_moderated_tables() -> None:
    result = DifferentialAnalysis().run(
        DifferentialAnalysisRequest(
            matrix=_matrix(),
            design=_design(),
            contrasts=_contrasts(),
            empirical_bayes=EmpiricalBayesConfig(method="standard"),
        )
    )

    assert set(result.contrast_tables) == {"B_vs_A", "C_vs_A"}
    for contrast_name in ("B_vs_A", "C_vs_A"):
        table = result.table_for(contrast_name)
        assert list(table.columns) == ["logFC", "t", "P.Value", "adj.P.Val"]
        assert table.shape[0] == 5
        assert table.index.tolist() == _matrix().index.tolist()
        assert (table.loc[:, "P.Value"] >= 0.0).all()
        assert (table.loc[:, "P.Value"] <= 1.0).all()
        assert (table.loc[:, "adj.P.Val"] >= 0.0).all()
        assert (table.loc[:, "adj.P.Val"] <= 1.0).all()


def test_differential_analysis_fails_on_sample_design_mismatch() -> None:
    design = _design().to_dataframe()
    design.index = pd.Index(["X1", "X2", "X3", "X4", "X5", "X6"], name="sample")
    with pytest.raises(
        PhosPyInputError,
        match="differential.matrix.columns must match differential.design.index",
    ):
        DifferentialAnalysis().run(
            DifferentialAnalysisRequest(
                matrix=_matrix(),
                design=design,
                contrasts=_contrasts(),
            )
        )


def test_differential_analysis_fails_on_contrast_design_term_mismatch() -> None:
    contrasts = _contrasts().to_dataframe().rename(index={"A": "A_wrong"})
    with pytest.raises(
        PhosPyInputError,
        match="differential.contrasts.index must match differential.design.columns",
    ):
        DifferentialAnalysis().run(
            DifferentialAnalysisRequest(
                matrix=_matrix(),
                design=_design(),
                contrasts=contrasts,
            )
        )


def test_differential_analysis_fails_when_residual_dof_is_non_positive() -> None:
    matrix = _matrix()
    identity_design = pd.DataFrame(
        {
            "A_1": [1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            "A_2": [0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
            "B_1": [0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
            "B_2": [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
            "C_1": [0.0, 0.0, 0.0, 0.0, 1.0, 0.0],
            "C_2": [0.0, 0.0, 0.0, 0.0, 0.0, 1.0],
        },
        index=matrix.columns.copy(),
    )
    identity_contrasts = pd.DataFrame(
        {"A_1_vs_A_2": [1.0, -1.0, 0.0, 0.0, 0.0, 0.0]},
        index=identity_design.columns.copy(),
    )

    with pytest.raises(
        PhosPyInputError,
        match="residual degrees of freedom must be positive",
    ):
        DifferentialAnalysis().run(
            DifferentialAnalysisRequest(
                matrix=matrix,
                design=identity_design,
                contrasts=identity_contrasts,
            )
        )
