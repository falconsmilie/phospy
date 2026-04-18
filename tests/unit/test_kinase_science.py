from __future__ import annotations

import pandas as pd
import pytest

from phospy.errors import WorkflowStageError
from phospy.workflows.kinase.science import (
    build_kinase_profiles,
    build_prediction_outputs,
)


def test_build_kinase_profiles_excludes_kinases_below_support_floor() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 2.0, 3.0],
            "sample_b": [2.0, 4.0, 6.0],
        },
        index=pd.Index(["SITE1", "SITE2", "SITE3"], name="site_id"),
    )
    kinase_substrate_map = pd.DataFrame(
        {
            "kinase": ["K1", "K1", "K2"],
            "substrate_site": ["SITE1", "SITE2", "SITE3"],
        }
    )

    result = build_kinase_profiles(
        phospho=phospho,
        kinase_substrate_map=kinase_substrate_map,
        min_substrates=2,
    )

    assert list(result.profile_matrix.index) == ["K1"]
    assert result.quantified_substrates == {"K1": ["SITE1", "SITE2"]}
    assert result.substrate_counts.to_dict() == {"K1": 2, "K2": 1}
    assert result.profile_matrix.at["K1", "sample_a"] == pytest.approx(1.5)
    assert result.profile_matrix.at["K1", "sample_b"] == pytest.approx(3.0)


def test_build_kinase_profiles_rejects_single_substrate_floor() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0],
            "sample_b": [2.0],
        },
        index=pd.Index(["SITE1"], name="site_id"),
    )
    kinase_substrate_map = pd.DataFrame(
        {
            "kinase": ["K1"],
            "substrate_site": ["SITE1"],
        }
    )

    with pytest.raises(
        WorkflowStageError,
        match="seam=kinase.science.min_substrate_floor",
    ):
        build_kinase_profiles(
            phospho=phospho,
            kinase_substrate_map=kinase_substrate_map,
            min_substrates=1,
        )


def test_build_prediction_outputs_uses_missing_values_for_unsupported_cells() -> None:
    score_matrix = pd.DataFrame(
        {"K1": [0.9, 0.1, 0.2], "K2": [0.2, 0.8, 0.3]},
        index=pd.Index(["S1", "S2", "S3"], name="site_id"),
        dtype=float,
    )
    pred_mat, substrate_list = build_prediction_outputs(
        score_matrix=score_matrix,
        selected_kinases=pd.Index(["K1", "K2"], name="kinase"),
        candidate_substrates={"K1": ["S1"], "K2": ["S2"]},
        top_k=1,
    )

    assert pred_mat.at["S1", "K1"] == pytest.approx(0.9)
    assert pred_mat.at["S2", "K2"] == pytest.approx(0.8)
    assert pd.isna(pred_mat.at["S1", "K2"])
    assert pd.isna(pred_mat.at["S2", "K1"])
    assert pd.isna(pred_mat.at["S3", "K1"])
    assert pd.isna(pred_mat.at["S3", "K2"])
    assert int(substrate_list.shape[0]) == 2
