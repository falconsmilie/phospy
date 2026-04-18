from __future__ import annotations

import pandas as pd
import pytest

from phospy.errors import WorkflowStageError
from phospy.workflows.kinase.science import build_kinase_profiles


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
