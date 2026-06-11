from __future__ import annotations

import pandas as pd

from phospy.science.prediction.motif_scoring import (
    KINASE_LIBRARY_SITE_STATUS_VALID_SCORED_SITE,
    KinaseLibraryMotifMatrix,
    score_kinase_library_motifs,
)


def test_motif_scoring_does_not_cross_score_ser_thr_and_tyr_lanes() -> None:
    result = score_kinase_library_motifs(
        site_sequences={
            "P1;S10;": "AST",
            "P2;Y20;": "AYP",
        },
        matrices=[
            _matrix("KINASE_ST", "ser_thr", center_residue="S", center_score=5.0),
            _matrix("KINASE_TYR", "tyr", center_residue="Y", center_score=7.0),
        ],
        flank_size=1,
    )

    assert result.raw_scores.at["P1;S10;", "KINASE_ST"] == 5.0
    assert pd.isna(result.raw_scores.at["P1;S10;", "KINASE_TYR"])
    assert pd.isna(result.raw_scores.at["P2;Y20;", "KINASE_ST"])
    assert result.raw_scores.at["P2;Y20;", "KINASE_TYR"] == 7.0
    assert set(result.site_diagnostics.loc[:, "status"]) == {
        KINASE_LIBRARY_SITE_STATUS_VALID_SCORED_SITE
    }


def _matrix(
    kinase: str,
    residue_class: str,
    *,
    center_residue: str,
    center_score: float,
) -> KinaseLibraryMotifMatrix:
    score_table = pd.DataFrame(
        0.0,
        index=pd.Index(["A", "S", "T", "Y", "P"], name="amino_acid"),
        columns=pd.Index([-1, 0, 1], name="position"),
    )
    score_table.at[center_residue, 0] = center_score
    return KinaseLibraryMotifMatrix(
        kinase=kinase,
        residue_class=residue_class,
        score_table=score_table,
    )
