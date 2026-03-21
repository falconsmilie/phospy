from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from phosrpy.motifs import (
    KinaseMotifScorer,
    create_frequency_matrix,
    frequency_scoring,
    score_phosphosite_motifs,
)


def test_create_frequency_matrix_normalizes_counts_and_ignores_gaps() -> None:
    frequency_mat = create_frequency_matrix(["A_A", "ACA"], flank_size=1)

    assert float(frequency_mat.loc["A", "p1"]) == pytest.approx(1.0)
    assert float(frequency_mat.loc["C", "p2"]) == pytest.approx(0.5)
    assert float(frequency_mat.loc["A", "p3"]) == pytest.approx(1.0)


def test_frequency_scoring_scores_valid_amino_acids_and_ignores_invalid() -> None:
    frequency_mat = create_frequency_matrix(["ACA", "AAA"], flank_size=1)
    sequence_list = pd.Series(
        ["ACA", "AXA", np.nan],
        index=["SITE_A", "SITE_X", "SITE_NA"],
    )

    result = frequency_scoring(sequence_list=sequence_list, frequency_mat=frequency_mat)

    assert float(result.loc["SITE_A"]) == pytest.approx(2.5)
    assert float(result.loc["SITE_X"]) == pytest.approx(2.0)
    assert float(result.loc["SITE_NA"]) == pytest.approx(0.0)


def test_kinase_motif_scorer_extracts_centered_windows_and_scales_scores() -> None:
    scorer = KinaseMotifScorer.from_substrate_sequences(
        {
            "KINASE_A": ["AAAAA"],
            "KINASE_B": ["TTTTT"],
        },
        flank_size=2,
    )

    result = scorer.score_sequences(
        seqs={
            "SITE_A": "QQAAAAAYY",
            "SITE_B": "QQTTTTTYY",
        },
        min_motif_size=1,
    )

    assert list(result.motif_scores.columns) == ["KINASE_A", "KINASE_B"]
    assert float(result.motif_scores.loc["SITE_A", "KINASE_A"]) == pytest.approx(1.0)
    assert float(result.motif_scores.loc["SITE_A", "KINASE_B"]) == pytest.approx(0.0)
    assert result.sequence_windows.loc["SITE_A"] == "AAAAA"


def test_score_phosphosite_motifs_filters_by_minimum_motif_size() -> None:
    motif_frequency_matrices = {
        "KINASE_A": create_frequency_matrix(["AAAAA"], flank_size=2),
        "KINASE_B": create_frequency_matrix(["TTTTT"], flank_size=2),
    }
    motif_sizes = pd.Series({"KINASE_A": 5, "KINASE_B": 1}, dtype=float)

    result = score_phosphosite_motifs(
        seqs={"SITE_A": "QQAAAAAYY"},
        motif_frequency_matrices=motif_frequency_matrices,
        motif_sizes=motif_sizes,
        min_motif_size=2,
        flank_size=2,
    )

    assert list(result.motif_scores.columns) == ["KINASE_A"]
    assert list(result.motif_sizes.index) == ["KINASE_A"]
