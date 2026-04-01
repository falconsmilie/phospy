from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import phospy.motifs as motifs
from phospy.motifs import (
    KinaseMotifScorer,
    create_frequency_matrix,
    frequency_scoring,
    score_phosphosite_motifs,
)
from phospy.validation.errors import TableSchemaError


def test_create_frequency_matrix_normalizes_counts_and_ignores_gaps() -> None:
    frequency_mat = create_frequency_matrix(["A_A", "ACA"], flank_size=1)

    assert float(frequency_mat.loc["A", "p1"]) == pytest.approx(1.0)
    assert float(frequency_mat.loc["C", "p2"]) == pytest.approx(0.5)
    assert float(frequency_mat.loc["A", "p3"]) == pytest.approx(1.0)


def test_create_frequency_matrix_normalizes_lowercase_sequences() -> None:
    lower = create_frequency_matrix(["aca"], flank_size=1)
    upper = create_frequency_matrix(["ACA"], flank_size=1)

    pd.testing.assert_frame_equal(lower, upper)


def test_frequency_scoring_normalizes_lowercase_sequences() -> None:
    frequency_mat = create_frequency_matrix(["ACA", "AAA"], flank_size=1)

    result = frequency_scoring(
        sequence_list=pd.Series(["ACA", "aca"], index=["SITE_UPPER", "SITE_LOWER"]),
        frequency_mat=frequency_mat,
    )

    assert float(result.loc["SITE_UPPER"]) == pytest.approx(2.5)
    assert float(result.loc["SITE_LOWER"]) == pytest.approx(2.5)


def test_frequency_scoring_rejects_invalid_amino_acid_characters() -> None:
    frequency_mat = create_frequency_matrix(["ACA", "AAA"], flank_size=1)

    with pytest.raises(TableSchemaError, match="invalid amino-acid characters"):
        frequency_scoring(
            sequence_list=pd.Series(["AXA", np.nan], index=["SITE_X", "SITE_NA"]),
            frequency_mat=frequency_mat,
        )


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


def test_kinase_motif_scorer_encodes_sequence_windows_once_per_batch(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    scorer = KinaseMotifScorer(
        motif_frequency_matrices={
            "KINASE_A": create_frequency_matrix(["AAAAA"], flank_size=2),
            "KINASE_B": create_frequency_matrix(["TTTTT"], flank_size=2),
            "KINASE_C": create_frequency_matrix(["SSSSS"], flank_size=2),
        },
        motif_sizes=pd.Series(
            {"KINASE_A": 2, "KINASE_B": 2, "KINASE_C": 2}, dtype=float
        ),
        flank_size=2,
    )

    encode_calls = 0
    original_encode = motifs._encode_sequence_positions

    def counting_encode(sequences: object, width: int) -> np.ndarray:
        nonlocal encode_calls
        encode_calls += 1
        return original_encode(sequences, width)

    monkeypatch.setattr(motifs, "_encode_sequence_positions", counting_encode)

    result = scorer.score_sequences(
        seqs={
            "SITE_A": "QQAAAAAYY",
            "SITE_B": "QQTTTTTYY",
            "SITE_C": "QQSSSSSYY",
        },
        min_motif_size=1,
    )

    assert encode_calls == 1
    assert list(result.motif_scores.columns) == ["KINASE_A", "KINASE_B", "KINASE_C"]


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
