from __future__ import annotations

import numpy as np
import pandas as pd

from phospy.api import (
    Contrast,
    DifferentialAnalysisRequest,
    DifferentialAnalysisWorkflow,
    ExperimentalDesign,
    Organism,
    SampleDesignRecord,
)
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.science.differential.models import (
    DIFFERENTIAL_RESULT_STATUS_COLUMN,
    DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN,
    DIFFERENTIAL_RESULT_STATUS_TESTED,
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_ALL_CONSTANT,
)
from phospy.science.statistics.multiple_testing import adjust_p_values
from tests.support.intensity_scale_states import (
    supported_log2_intensity_scale_state,
    supported_log2_processing_state,
)
from tests.support.site_keys import protein_site_key_index, site_key_context_columns

_GENES = ("MAPK14", "AKT1", "GSK3B")
_SITES = ("Y182", "T308", "S9")


def _site_index() -> pd.Index:
    return protein_site_key_index(
        protein_identifiers=_GENES,
        sites=_SITES,
    )


def _matrix() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "A_1": [5.0, 1.0, 2.0],
            "A_2": [5.0, 1.1, 2.2],
            "B_1": [5.0, 2.0, 1.8],
            "B_2": [5.0, 2.1, 2.0],
        },
        index=_site_index(),
    )


def _dataset() -> AnalysisReadyPhosphoDataset:
    site_index = _site_index()
    site_metadata = pd.DataFrame(
        {
            "site_key": site_index.astype(str).tolist(),
            "display_id": [
                f"{gene};{site};" for gene, site in zip(_GENES, _SITES, strict=True)
            ],
            **site_key_context_columns(site_index),
            "gene_symbol": list(_GENES),
            "site": list(_SITES),
            "site_sequence": [("A" * 15) + site[0] + ("A" * 15) for site in _SITES],
            "protein_id": list(_GENES),
        },
        index=site_index.copy(),
    )
    return AnalysisReadyPhosphoDataset(
        phospho=_matrix(),
        site_metadata=site_metadata,
        organism=Organism.RAT,
        intensity_scale_state=supported_log2_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_log2_processing_state(has_total_matrix=False),
    )


def _request() -> DifferentialAnalysisRequest:
    return DifferentialAnalysisRequest(
        dataset=_dataset(),
        design=ExperimentalDesign(
            samples=(
                SampleDesignRecord(
                    sample_id="A_1",
                    condition="A",
                    biological_replicate_id="A_1",
                ),
                SampleDesignRecord(
                    sample_id="A_2",
                    condition="A",
                    biological_replicate_id="A_2",
                ),
                SampleDesignRecord(
                    sample_id="B_1",
                    condition="B",
                    biological_replicate_id="B_1",
                ),
                SampleDesignRecord(
                    sample_id="B_2",
                    condition="B",
                    biological_replicate_id="B_2",
                ),
            )
        ),
        contrasts=(
            Contrast(
                name="B_vs_A",
                numerator_condition="B",
                denominator_condition="A",
            ),
        ),
    )


def test_differential_marks_all_constant_rows_as_withheld() -> None:
    result = DifferentialAnalysisWorkflow().run(_request())
    table = result.table_for("B_vs_A")
    constant_site = _site_index()[0]

    assert (
        table.at[constant_site, DIFFERENTIAL_RESULT_STATUS_COLUMN]
        == DIFFERENTIAL_RESULT_STATUS_WITHHELD_ALL_CONSTANT
    )
    assert table.loc[constant_site, ["logFC", "t", "P.Value", "adj.P.Val"]].isna().all()
    assert result.feature_eligibility is not None
    assert (
        result.feature_eligibility.at[constant_site, DIFFERENTIAL_RESULT_STATUS_COLUMN]
        == DIFFERENTIAL_RESULT_STATUS_WITHHELD_ALL_CONSTANT
    )


def test_differential_mixed_valid_and_constant_rows_tests_valid_rows() -> None:
    result = DifferentialAnalysisWorkflow().run(_request())
    table = result.table_for("B_vs_A")
    tested = (
        table[DIFFERENTIAL_RESULT_STATUS_COLUMN] == DIFFERENTIAL_RESULT_STATUS_TESTED
    )
    withheld = ~tested

    assert tested.tolist() == [False, True, True]
    assert (
        np.isfinite(table.loc[tested, ["logFC", "t", "P.Value", "adj.P.Val"]])
        .all()
        .all()
    )
    assert (
        table.loc[withheld, ["logFC", "t", "P.Value", "adj.P.Val"]].isna().all().all()
    )

    expected_adjusted = adjust_p_values(
        table.loc[tested, "P.Value"].to_numpy(dtype=float),
        method="benjamini_hochberg",
    )
    np.testing.assert_allclose(
        table.loc[tested, "adj.P.Val"].to_numpy(dtype=float),
        expected_adjusted,
        rtol=1e-12,
        atol=1e-12,
    )


def test_differential_withheld_rows_have_reason_codes() -> None:
    result = DifferentialAnalysisWorkflow().run(_request())
    table = result.table_for("B_vs_A")
    constant_site = _site_index()[0]
    row = table.loc[constant_site, :]

    assert row["site_key"] == constant_site
    assert row[DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN]
    assert "all-constant" in str(row[DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN])
    assert int(row["analysed_value_count"]) == 4
    assert int(row["observed_value_count"]) == 4
    assert int(row["invalid_numeric_value_count"]) == 0
    assert int(row["unique_observed_value_count"]) == 1
