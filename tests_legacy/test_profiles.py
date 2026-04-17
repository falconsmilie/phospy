from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from phospy.prediction.profiles import (
    KinaseProfilePolicy,
    build_kinase_substrate_profiles,
)


def make_phospho_matrix() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "s1": [1.0, 3.0, 5.0],
            "s2": [2.0, 4.0, 6.0],
            "s3": [3.0, np.nan, 9.0],
        },
        index=["SITE_1", "SITE_2", "SITE_3"],
    )


def test_build_kinase_substrate_profiles_matches_phosr_summary_rules() -> None:
    phospho_matrix = make_phospho_matrix()
    substrate_map = {
        "KINASE_A": ["SITE_1"],
        "KINASE_B": ["SITE_1", "SITE_2"],
        "KINASE_C": ["MISSING_SITE"],
        "KINASE_D": ["SITE_1", "SITE_1", "SITE_3"],
    }

    result = build_kinase_substrate_profiles(
        substrate_map=substrate_map,
        phospho_matrix=phospho_matrix,
    )

    pd.testing.assert_series_equal(
        result.profile_matrix.loc["KINASE_A"],
        phospho_matrix.loc["SITE_1"].rename("KINASE_A").astype(float),
    )
    assert float(result.profile_matrix.loc["KINASE_B", "s1"]) == pytest.approx(2.0)
    assert float(result.profile_matrix.loc["KINASE_B", "s2"]) == pytest.approx(3.0)
    assert np.isnan(result.profile_matrix.loc["KINASE_B", "s3"])
    assert result.quantified_substrates["KINASE_D"] == ["SITE_1", "SITE_3"]
    assert int(result.substrate_counts.loc["KINASE_D"]) == 2
    assert "KINASE_C" not in result.profile_matrix.index
    assert int(result.substrate_counts.loc["KINASE_C"]) == 0


def test_build_kinase_substrate_profiles_applies_minimum_substrate_filter() -> None:
    phospho_matrix = make_phospho_matrix()
    substrate_map = {
        "KINASE_A": ["SITE_1"],
        "KINASE_B": ["SITE_1", "SITE_2"],
    }

    result = build_kinase_substrate_profiles(
        substrate_map=substrate_map,
        phospho_matrix=phospho_matrix,
        min_substrates=2,
    )

    assert list(result.profile_matrix.index) == ["KINASE_B"]
    assert int(result.substrate_counts.loc["KINASE_A"]) == 1


def test_build_kinase_substrate_profiles_returns_expected_profile_rows() -> None:
    phospho_matrix = make_phospho_matrix()
    substrate_map = {"KINASE_A": ["SITE_1", "SITE_2"]}

    result = build_kinase_substrate_profiles(
        substrate_map=substrate_map,
        phospho_matrix=phospho_matrix,
    )

    assert list(result.profile_matrix.index) == ["KINASE_A"]
    assert float(result.profile_matrix.loc["KINASE_A", "s1"]) == pytest.approx(2.0)


def test_kinase_profile_result_is_detached_from_input_matrix() -> None:
    phospho_matrix = make_phospho_matrix()
    original = phospho_matrix.copy(deep=True)

    result = build_kinase_substrate_profiles(
        substrate_map={"KINASE_A": ["SITE_1", "SITE_2"]},
        phospho_matrix=phospho_matrix,
    )

    result.profile_matrix.loc["KINASE_A", "s1"] = -999.0
    result.substrate_counts.loc["KINASE_A"] = 999

    pd.testing.assert_frame_equal(phospho_matrix, original)


def test_build_kinase_substrate_profiles_uses_vectorized_median(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    phospho_matrix = make_phospho_matrix()

    def raising_apply(self, *args, **kwargs):
        msg = "DataFrame.apply should not be used for profile aggregation"
        raise AssertionError(msg)

    monkeypatch.setattr(pd.DataFrame, "apply", raising_apply)

    result = build_kinase_substrate_profiles(
        substrate_map={"KINASE_A": ["SITE_1", "SITE_2"]},
        phospho_matrix=phospho_matrix,
    )

    assert float(result.profile_matrix.loc["KINASE_A", "s1"]) == pytest.approx(2.0)


def test_build_kinase_substrate_profiles_can_skip_missing_values_when_requested() -> (
    None
):
    phospho_matrix = make_phospho_matrix()

    result = build_kinase_substrate_profiles(
        substrate_map={"KINASE_A": ["SITE_1", "SITE_2"]},
        phospho_matrix=phospho_matrix,
        policy=KinaseProfilePolicy(missing_value_strategy="median_skipna"),
    )

    assert float(result.profile_matrix.loc["KINASE_A", "s1"]) == pytest.approx(2.0)
    assert float(result.profile_matrix.loc["KINASE_A", "s2"]) == pytest.approx(3.0)
    assert float(result.profile_matrix.loc["KINASE_A", "s3"]) == pytest.approx(3.0)
