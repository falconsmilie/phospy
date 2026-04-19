from __future__ import annotations

import pandas as pd
import pandas.testing as pdt
import pytest

from phospy.api.configs import (
    KINASE_PROFILE_MISSING_VALUE_STRATEGY_MEDIAN_SKIPNA,
    KINASE_PROFILE_MISSING_VALUE_STRATEGY_STRICT,
)
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


def test_build_kinase_profiles_supports_strict_and_median_skipna_policies() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 3.0, 9.0],
            "sample_b": [2.0, 4.0, 8.0],
            "sample_c": [3.0, float("nan"), 7.0],
        },
        index=pd.Index(["SITE1", "SITE2", "SITE3"], name="site_id"),
    )
    kinase_substrate_map = pd.DataFrame(
        {
            "kinase": ["K1", "K1", "K2", "K2"],
            "substrate_site": ["SITE1", "SITE2", "SITE2", "SITE3"],
        }
    )

    strict = build_kinase_profiles(
        phospho=phospho,
        kinase_substrate_map=kinase_substrate_map,
        min_substrates=2,
        profile_missing_value_strategy=KINASE_PROFILE_MISSING_VALUE_STRATEGY_STRICT,
    )
    skipna = build_kinase_profiles(
        phospho=phospho,
        kinase_substrate_map=kinase_substrate_map,
        min_substrates=2,
        profile_missing_value_strategy=(
            KINASE_PROFILE_MISSING_VALUE_STRATEGY_MEDIAN_SKIPNA
        ),
    )

    assert pd.isna(strict.profile_matrix.at["K1", "sample_c"])
    assert skipna.profile_matrix.at["K1", "sample_c"] == pytest.approx(3.0)
    assert pd.isna(strict.profile_matrix.at["K2", "sample_c"])
    assert skipna.profile_matrix.at["K2", "sample_c"] == pytest.approx(7.0)
    assert strict.profile_matrix.at["K1", "sample_a"] == pytest.approx(2.0)
    assert skipna.profile_matrix.at["K1", "sample_a"] == pytest.approx(2.0)


def test_build_kinase_profiles_rejects_unknown_missing_value_strategy() -> None:
    phospho = pd.DataFrame(
        {"sample_a": [1.0, 2.0], "sample_b": [3.0, 4.0]},
        index=pd.Index(["SITE1", "SITE2"], name="site_id"),
    )
    kinase_substrate_map = pd.DataFrame(
        {"kinase": ["K1", "K1"], "substrate_site": ["SITE1", "SITE2"]}
    )

    with pytest.raises(
        WorkflowStageError,
        match="seam=kinase.science.profile_missing_value_strategy",
    ):
        build_kinase_profiles(
            phospho=phospho,
            kinase_substrate_map=kinase_substrate_map,
            min_substrates=2,
            profile_missing_value_strategy="unknown",  # type: ignore[arg-type]
        )


def test_build_prediction_outputs_uses_missing_values_for_unsupported_cells() -> None:
    score_matrix = pd.DataFrame(
        {"K1": [0.9, 0.1, 0.2], "K2": [0.2, 0.8, 0.3]},
        index=pd.Index(["S1", "S2", "S3"], name="site_id"),
        dtype=float,
    )
    pred_mat, substrate_list = build_prediction_outputs(
        prediction_score_matrix=score_matrix,
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


def test_build_prediction_outputs_matches_legacy_selection_semantics() -> None:
    score_matrix = pd.DataFrame(
        {
            "K1": [0.8, 0.4, float("nan"), 0.95, 0.2, 0.7],
            "K2": [0.1, 0.6, 0.75, float("nan"), 0.3, 0.5],
        },
        index=pd.Index(["S1", "S2", "S3", "S4", "S5", "S6"], name="site_id"),
        dtype=float,
    )
    selected_kinases = pd.Index(["K1", "K2"], name="kinase")
    candidate_substrates = {
        "K1": ["S6", "MISSING", "S1", "S3", "S4", "S2", "S5"],
        "K2": ["S1", "S5", "S3", "S2", "MISSING", "S6"],
    }
    top_k = 3

    def _legacy_build_prediction_outputs() -> tuple[pd.DataFrame, pd.DataFrame]:
        pred_mat = pd.DataFrame(
            float("nan"),
            index=score_matrix.index.copy(),
            columns=selected_kinases.copy(),
        )
        pred_mat.index.name = score_matrix.index.name
        pred_mat.columns.name = "kinase"
        substrate_rows: list[dict[str, object]] = []
        for kinase in selected_kinases:
            candidate_sites = candidate_substrates.get(str(kinase), [])
            available_sites = [
                site for site in candidate_sites if site in score_matrix.index
            ]
            if not available_sites:
                continue
            ranked_sites = (
                score_matrix.loc[available_sites, kinase]
                .astype(float)
                .dropna()
                .sort_values(ascending=False)
                .head(top_k)
            )
            if ranked_sites.empty:
                continue
            pred_mat.loc[ranked_sites.index, kinase] = ranked_sites.values
            for rank, (site_id, score) in enumerate(ranked_sites.items(), start=1):
                substrate_rows.append(
                    {
                        "kinase": str(kinase),
                        "substrate_site": site_id,
                        "score": float(score),
                        "rank": rank,
                    }
                )
        substrate_list = pd.DataFrame(
            substrate_rows,
            columns=["kinase", "substrate_site", "score", "rank"],
        )
        return pred_mat, substrate_list

    expected_pred_mat, expected_substrate_list = _legacy_build_prediction_outputs()
    observed_pred_mat, observed_substrate_list = build_prediction_outputs(
        prediction_score_matrix=score_matrix,
        selected_kinases=selected_kinases,
        candidate_substrates=candidate_substrates,
        top_k=top_k,
    )

    pdt.assert_frame_equal(observed_pred_mat, expected_pred_mat, check_dtype=False)
    pdt.assert_frame_equal(
        observed_substrate_list.reset_index(drop=True),
        expected_substrate_list.reset_index(drop=True),
        check_dtype=False,
    )
