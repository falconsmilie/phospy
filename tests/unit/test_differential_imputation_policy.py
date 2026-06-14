from __future__ import annotations

from dataclasses import replace

import numpy as np
import pandas as pd
import pytest

from phospy.api import (
    Contrast,
    DifferentialAnalysisConfig,
    DifferentialAnalysisRequest,
    DifferentialAnalysisWorkflow,
    ExperimentalDesign,
    Organism,
    SampleDesignRecord,
)
from phospy.errors import WorkflowValidationError
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.science.differential.models import (
    DIFFERENTIAL_RESULT_STATUS_TESTED,
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_HIGH_IMPUTATION,
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_INSUFFICIENT_OBSERVED,
)
from phospy.science.differential.multiple_testing import benjamini_hochberg
from tests.support.intensity_scale_states import (
    supported_log2_intensity_scale_state,
    supported_log2_processing_state,
)
from tests.support.site_keys import protein_site_key_index, site_key_context_columns

_GENES = ["MAPK14", "AKT1", "GSK3B", "RPS6"]
_SITES = ["Y182", "T308", "S9", "S235"]
_SAMPLES = ("A_1", "A_2", "A_3", "B_1", "B_2", "B_3")


def _site_index() -> pd.Index:
    return protein_site_key_index(
        protein_identifiers=_GENES,
        sites=_SITES,
    )


def _site_metadata(index: pd.Index) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "site_key": index.astype(str).tolist(),
            "display_id": [
                f"{gene};{site};" for gene, site in zip(_GENES, _SITES, strict=True)
            ],
            **site_key_context_columns(index),
            "gene_symbol": _GENES,
            "site": _SITES,
            "site_sequence": [("A" * 15) + site[0] + ("A" * 15) for site in _SITES],
            "protein_id": _GENES,
        },
        index=index.copy(),
    )


def _phospho(index: pd.Index) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "A_1": [10.0, 20.0, 5.0, 30.0],
            "A_2": [10.1, 20.2, 5.1, 30.2],
            "A_3": [9.9, 20.1, 5.2, 30.1],
            "B_1": [15.0, 24.0, 6.0, 29.0],
            "B_2": [15.1, 24.1, 6.2, 29.2],
            "B_3": [14.9, 24.2, 6.1, 29.1],
        },
        index=index.copy(),
    )


def _observed_mask(index: pd.Index) -> pd.DataFrame:
    mask = pd.DataFrame(True, index=index.copy(), columns=pd.Index(_SAMPLES))
    mask.loc[index[1], "B_3"] = False
    mask.loc[index[2], ["A_3", "B_3"]] = False
    mask.loc[index[3], ["B_2", "B_3"]] = False
    return mask


def _imputed_processing_state():
    processing_state = supported_log2_processing_state(has_total_matrix=False)
    return replace(
        processing_state,
        missing_data=replace(processing_state.missing_data, imputed=True),
    )


def _imputed_dataset(*, with_metadata: bool = True) -> AnalysisReadyPhosphoDataset:
    index = _site_index()
    return AnalysisReadyPhosphoDataset(
        phospho=_phospho(index),
        site_metadata=_site_metadata(index),
        imputation_observation_mask=(_observed_mask(index) if with_metadata else None),
        organism=Organism.RAT,
        intensity_scale_state=supported_log2_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=_imputed_processing_state(),
    )


def _design() -> ExperimentalDesign:
    return ExperimentalDesign(
        samples=tuple(
            SampleDesignRecord(
                sample_id=sample_id,
                condition=sample_id.split("_", maxsplit=1)[0],
                biological_replicate_id=sample_id,
            )
            for sample_id in _SAMPLES
        )
    )


def _contrasts() -> tuple[Contrast, ...]:
    return (
        Contrast(
            name="B_vs_A",
            numerator_condition="B",
            denominator_condition="A",
        ),
    )


def _withhold_request() -> DifferentialAnalysisRequest:
    return DifferentialAnalysisRequest(
        dataset=_imputed_dataset(),
        design=_design(),
        contrasts=_contrasts(),
        config=DifferentialAnalysisConfig(
            imputed_value_policy="withhold_imputed_features",
            imputed_value_max_fraction=0.20,
        ),
    )


def test_differential_imputation_policy_reject_is_default() -> None:
    with pytest.raises(WorkflowValidationError, match="imputed cells"):
        DifferentialAnalysisWorkflow().run(
            DifferentialAnalysisRequest(
                dataset=_imputed_dataset(),
                design=_design(),
                contrasts=_contrasts(),
            )
        )


def test_differential_imputation_policy_requires_metadata_for_non_reject() -> None:
    with pytest.raises(
        WorkflowValidationError,
        match="requires dataset-owned imputation observation metadata",
    ):
        DifferentialAnalysisWorkflow().run(
            DifferentialAnalysisRequest(
                dataset=_imputed_dataset(with_metadata=False),
                design=_design(),
                contrasts=_contrasts(),
                config=DifferentialAnalysisConfig(
                    imputed_value_policy="withhold_imputed_features",
                    imputed_value_max_fraction=0.20,
                ),
            )
        )


def test_differential_withhold_policy_marks_high_imputation_features() -> None:
    result = DifferentialAnalysisWorkflow().run(_withhold_request())
    table = result.table_for("B_vs_A")

    assert {
        "imputed_cell_count",
        "observed_cell_count",
        "imputed_fraction",
        "imputation_policy",
        "result_status",
    }.issubset(table.columns)
    assert table["imputed_cell_count"].tolist() == [0, 1, 2, 2]
    assert table["observed_cell_count"].tolist() == [6, 5, 4, 4]
    assert table["imputed_fraction"].tolist() == [0.0, 1 / 6, 2 / 6, 2 / 6]
    assert table["imputation_policy"].unique().tolist() == ["withhold_imputed_features"]
    assert table["result_status"].tolist() == [
        DIFFERENTIAL_RESULT_STATUS_TESTED,
        DIFFERENTIAL_RESULT_STATUS_TESTED,
        DIFFERENTIAL_RESULT_STATUS_WITHHELD_HIGH_IMPUTATION,
        DIFFERENTIAL_RESULT_STATUS_WITHHELD_INSUFFICIENT_OBSERVED,
    ]


def test_differential_withhold_policy_excludes_withheld_features_from_testing() -> None:
    result = DifferentialAnalysisWorkflow().run(_withhold_request())
    table = result.table_for("B_vs_A")
    tested = table["result_status"] == DIFFERENTIAL_RESULT_STATUS_TESTED
    withheld = ~tested

    assert (
        np.isfinite(table.loc[tested, ["logFC", "t", "P.Value", "adj.P.Val"]])
        .all()
        .all()
    )
    assert (
        table.loc[withheld, ["logFC", "t", "P.Value", "adj.P.Val"]].isna().all().all()
    )
    assert result.residual_variance_series().loc[table.index[withheld]].isna().all()

    expected_adjusted = benjamini_hochberg(
        table.loc[tested, "P.Value"].to_numpy(dtype=float)
    )
    np.testing.assert_allclose(
        table.loc[tested, "adj.P.Val"].to_numpy(dtype=float),
        expected_adjusted,
        rtol=1e-12,
        atol=1e-12,
    )


def test_differential_imputation_policy_is_recorded_in_provenance() -> None:
    result = DifferentialAnalysisWorkflow().run(_withhold_request())

    assert result.policy_provenance is not None
    missing_values = result.policy_provenance.missing_values
    assert missing_values.imputed_value_policy == "withhold_imputed_features"
    assert missing_values.imputed_value_max_fraction == pytest.approx(0.20)
    assert missing_values.imputation_metadata_required is True
    assert missing_values.adjusted_p_value_scope == (
        "benjamini_hochberg_adjustment_over_tested_features_only"
    )


def test_differential_imputation_policy_output_is_reproducible() -> None:
    first = DifferentialAnalysisWorkflow().run(_withhold_request()).table_for("B_vs_A")
    second = DifferentialAnalysisWorkflow().run(_withhold_request()).table_for("B_vs_A")

    pd.testing.assert_frame_equal(first, second)
