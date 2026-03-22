from __future__ import annotations

import pandas as pd

from phospy import build_site_matrix


def test_build_site_matrix_creates_site_ids_and_deduplicates_by_mean() -> None:
    df = pd.DataFrame(
        {
            "gene_p_site": ["PRKACA_S339", "PRKACA_S339", "BTK_Y551"],
            "centralized_sequence": ["AAAAAA", "BBBBBB", "CCCCCC"],
            "phospho_corrected_1": [1.0, 10.0, 3.0],
            "phospho_corrected_2": [1.0, 10.0, 3.0],
            "phospho_corrected_3": [1.0, 10.0, 3.0],
            "phospho_corrected_4": [1.0, 10.0, 3.0],
            "phospho_corrected_5": [1.0, 10.0, 3.0],
            "phospho_corrected_6": [1.0, 10.0, 3.0],
        }
    )
    phosr_input, matrix, sequences = build_site_matrix(
        df=df,
        gene_p_site_col="gene_p_site",
        sequence_col="centralized_sequence",
        value_cols=[
            "phospho_corrected_1",
            "phospho_corrected_2",
            "phospho_corrected_3",
            "phospho_corrected_4",
            "phospho_corrected_5",
            "phospho_corrected_6",
        ],
    )
    assert "PRKACA;S339;" in matrix.index
    assert matrix.loc["PRKACA;S339;", "phospho_corrected_1"] == 10.0
    assert sequences.loc["PRKACA;S339;"] == "BBBBBB"
    assert phosr_input.shape[0] == 2
