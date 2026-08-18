from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from phospy.advanced import DifferentialAnalysisConfig
from phospy.api import (
    Contrast,
    DifferentialAnalysisRequest,
    DifferentialAnalysisWorkflow,
    ExperimentalDesign,
    Organism,
    SampleDesignRecord,
)
from phospy.contracts.configs import PAIRED_DESIGN_POLICY_DUPLICATE_CORRELATION
from phospy.contracts.result_caveats import result_caveats_from_payloads
from phospy.errors import WorkflowValidationError
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.science.differential.models import (
    DIFFERENTIAL_RESULT_STATUS_TESTED,
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_HIGH_IMPUTATION,
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_INSUFFICIENT_OBSERVED,
)
from phospy.science.statistics.multiple_testing import adjust_p_values
from phospy.workflows.differential.caveats import (
    DIFFERENTIAL_IMPUTATION_WITHHOLDING_POLICY_CAVEAT_CODE,
)
from phospy.workflows.differential.imputation_inference import (
    DIFFERENTIAL_IMPUTATION_INFERENCE_STATUS_NO_TESTED_IMPUTED_VALUES,
    DIFFERENTIAL_IMPUTATION_INFERENCE_STATUS_RETAINED_IMPUTED_VALUES,
    DIFFERENTIAL_ROW_INFERENCE_STATUS_TESTED_FULLY_OBSERVED,
    DIFFERENTIAL_ROW_INFERENCE_STATUS_TESTED_RETAINED_IMPUTED_VALUES,
    DIFFERENTIAL_ROW_INFERENCE_STATUS_WITHHELD_NOT_TESTED,
)
from tests.support.analysis_ready_dataset_factories import (
    trusted_analysis_ready_dataset_from_tables,
)
from tests.support.intensity_scale_states import (
    supported_log2_intensity_scale_state,
    supported_log2_processing_state,
)
from tests.support.processing_state import (
    imputed_processing_state as valid_imputed_processing_state,
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
    return valid_imputed_processing_state(processing_state)


def _imputed_dataset(
    *,
    with_metadata: bool = True,
    observed_mask: pd.DataFrame | None = None,
) -> AnalysisReadyPhosphoDataset:
    index = _site_index()
    return trusted_analysis_ready_dataset_from_tables(
        phospho=_phospho(index),
        site_metadata=_site_metadata(index),
        imputation_observation_mask=(
            observed_mask
            if observed_mask is not None
            else (_observed_mask(index) if with_metadata else None)
        ),
        organism=Organism.RAT,
        intensity_scale_state=supported_log2_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=_imputed_processing_state(),
    )


def _design(sample_ids: tuple[str, ...] = _SAMPLES) -> ExperimentalDesign:
    return ExperimentalDesign(
        samples=tuple(
            SampleDesignRecord(
                sample_id=sample_id,
                condition=sample_id.split("_", maxsplit=1)[0],
                biological_replicate_id=sample_id,
            )
            for sample_id in sample_ids
        )
    )


def _paired_design() -> ExperimentalDesign:
    block_ids = {
        "A_1": "pair_1",
        "A_2": "pair_2",
        "A_3": "pair_3",
        "B_1": "pair_1",
        "B_2": "pair_2",
        "B_3": "pair_3",
    }
    return ExperimentalDesign(
        samples=tuple(
            SampleDesignRecord(
                sample_id=sample_id,
                condition=sample_id.split("_", maxsplit=1)[0],
                biological_replicate_id=sample_id,
                block_id=block_ids[sample_id],
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


def _caveat_by_code(result, code: str):
    matches = [caveat for caveat in result.caveats if caveat.code == code]
    assert len(matches) == 1
    return matches[0]


def test_differential_imputation_policy_reject_is_default() -> None:
    with pytest.raises(WorkflowValidationError, match="imputed cells"):
        DifferentialAnalysisWorkflow().run(
            DifferentialAnalysisRequest(
                dataset=_imputed_dataset(),
                design=_design(),
                contrasts=_contrasts(),
            )
        )


def test_differential_rejects_imputed_dataset_without_required_policy() -> None:
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
        "contains_imputed_cells",
        "observed_only_fit",
        "residual_df_adjusted_for_imputation",
        "inferential_status",
        "result_status",
    }.issubset(table.columns)
    assert table["imputed_cell_count"].tolist() == [0, 1, 2, 2]
    assert table["observed_cell_count"].tolist() == [6, 5, 4, 4]
    assert table["imputed_fraction"].tolist() == [0.0, 1 / 6, 2 / 6, 2 / 6]
    assert table["imputation_policy"].unique().tolist() == ["withhold_imputed_features"]
    assert table["contains_imputed_cells"].tolist() == [False, True, True, True]
    assert table["observed_only_fit"].unique().tolist() == [False]
    assert table["residual_df_adjusted_for_imputation"].unique().tolist() == [False]
    assert table["inferential_status"].tolist() == [
        DIFFERENTIAL_ROW_INFERENCE_STATUS_TESTED_FULLY_OBSERVED,
        DIFFERENTIAL_ROW_INFERENCE_STATUS_TESTED_RETAINED_IMPUTED_VALUES,
        DIFFERENTIAL_ROW_INFERENCE_STATUS_WITHHELD_NOT_TESTED,
        DIFFERENTIAL_ROW_INFERENCE_STATUS_WITHHELD_NOT_TESTED,
    ]
    assert table["result_status"].tolist() == [
        DIFFERENTIAL_RESULT_STATUS_TESTED,
        DIFFERENTIAL_RESULT_STATUS_TESTED,
        DIFFERENTIAL_RESULT_STATUS_WITHHELD_HIGH_IMPUTATION,
        DIFFERENTIAL_RESULT_STATUS_WITHHELD_INSUFFICIENT_OBSERVED,
    ]
    assert result.feature_eligibility is not None
    assert result.feature_eligibility["inferential_status"].tolist() == (
        table["inferential_status"].tolist()
    )

    caveat = _caveat_by_code(
        result,
        DIFFERENTIAL_IMPUTATION_WITHHOLDING_POLICY_CAVEAT_CODE,
    )
    assert caveat.severity == "warning"
    assert "retained tested rows" in caveat.message
    assert caveat.details["tested_feature_count"] == 2
    assert caveat.details["testable_feature_count"] == 2
    assert caveat.details["withheld_feature_count"] == 2
    assert caveat.details["tested_imputed_feature_count"] == 1
    assert caveat.details["tested_imputed_cell_count"] == 1
    assert caveat.details["observed_only_fit"] is False
    assert caveat.details["residual_df_adjusted_for_imputation"] is False
    assert caveat.details["inferential_status"] == (
        DIFFERENTIAL_IMPUTATION_INFERENCE_STATUS_RETAINED_IMPUTED_VALUES
    )
    assert caveat.details["adjusted_p_value_denominator_feature_count"] == 2


def test_differential_duplicate_correlation_withhold_policy_records_imputation_residual_df_contract() -> (
    None
):
    result = DifferentialAnalysisWorkflow().run(
        DifferentialAnalysisRequest(
            dataset=_imputed_dataset(),
            design=_paired_design(),
            contrasts=_contrasts(),
            config=DifferentialAnalysisConfig(
                paired_design_policy=PAIRED_DESIGN_POLICY_DUPLICATE_CORRELATION,
                imputed_value_policy="withhold_imputed_features",
                imputed_value_max_fraction=0.20,
            ),
        )
    )
    table = result.table_for("B_vs_A")

    assert table["observed_only_fit"].unique().tolist() == [False]
    assert table["residual_df_adjusted_for_imputation"].unique().tolist() == [False]
    assert table["inferential_status"].tolist() == [
        DIFFERENTIAL_ROW_INFERENCE_STATUS_TESTED_FULLY_OBSERVED,
        DIFFERENTIAL_ROW_INFERENCE_STATUS_TESTED_RETAINED_IMPUTED_VALUES,
        DIFFERENTIAL_ROW_INFERENCE_STATUS_WITHHELD_NOT_TESTED,
        DIFFERENTIAL_ROW_INFERENCE_STATUS_WITHHELD_NOT_TESTED,
    ]
    assert result.residual_degrees_of_freedom == 4.0
    assert result.policy_provenance is not None
    assert result.policy_provenance.missing_values.observed_only_fit is False
    assert (
        result.policy_provenance.missing_values.residual_df_adjusted_for_imputation
        is False
    )
    duplicate = result.policy_provenance.duplicate_correlation
    assert duplicate is not None
    assert duplicate.imputed_values_participated is True
    assert duplicate.imputed_feature_count == 1
    assert duplicate.imputed_cell_count == 1
    assert duplicate.design_rank == 2
    assert duplicate.gls_fit_status == "fit"
    assert result.workflow_provenance is not None
    assert (
        result.workflow_provenance["imputation_inference"][
            "residual_df_adjusted_for_imputation"
        ]
        is False
    )


def test_differential_imputation_subset_summary_uses_dataset_owned_summary_contract(
    monkeypatch,
) -> None:
    summary_calls: list[tuple[tuple[object, ...], tuple[object, ...]]] = []
    original_summary = (
        AnalysisReadyPhosphoDataset.imputation_observation_summary_dataframe
    )

    def summary_spy(
        self: AnalysisReadyPhosphoDataset,
        *,
        feature_ids,
        sample_ids,
    ):
        summary_calls.append((tuple(feature_ids), tuple(sample_ids)))
        return original_summary(
            self,
            feature_ids=feature_ids,
            sample_ids=sample_ids,
        )

    def raw_mask_guard(self: AnalysisReadyPhosphoDataset):
        raise AssertionError("differential interpreter must use dataset summaries")

    monkeypatch.setattr(
        AnalysisReadyPhosphoDataset,
        "imputation_observation_summary_dataframe",
        summary_spy,
    )
    monkeypatch.setattr(
        AnalysisReadyPhosphoDataset,
        "_borrow_imputation_observed_mask_frame",
        raw_mask_guard,
    )
    analysis_sample_ids = ("A_1", "A_2", "B_1", "B_2")

    result = DifferentialAnalysisWorkflow().run(
        DifferentialAnalysisRequest(
            dataset=_imputed_dataset(),
            design=_design(analysis_sample_ids),
            contrasts=_contrasts(),
            config=DifferentialAnalysisConfig(
                imputed_value_policy="withhold_imputed_features",
                imputed_value_max_fraction=0.20,
                allow_design_subset=True,
            ),
        )
    )
    table = result.table_for("B_vs_A")

    assert table["imputed_cell_count"].tolist() == [0, 0, 0, 1]
    assert table["observed_cell_count"].tolist() == [4, 4, 4, 3]
    caveat = _caveat_by_code(
        result,
        DIFFERENTIAL_IMPUTATION_WITHHOLDING_POLICY_CAVEAT_CODE,
    )
    assert caveat.details["tested_feature_count"] == 3
    assert caveat.details["tested_imputed_feature_count"] == 0
    assert caveat.details["tested_imputed_cell_count"] == 0
    assert caveat.details["inferential_status"] == (
        DIFFERENTIAL_IMPUTATION_INFERENCE_STATUS_NO_TESTED_IMPUTED_VALUES
    )
    assert summary_calls
    assert summary_calls[0] == (tuple(_site_index()), analysis_sample_ids)
    assert (tuple(_site_index()), ("A_1", "A_2")) in summary_calls
    assert (tuple(_site_index()), ("B_1", "B_2")) in summary_calls


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

    assert result.workflow_provenance is not None
    assert result.workflow_provenance["row_attrition_metrics"] == {
        "input_sites": 4,
        "sites_retained_for_model_fitting": 2,
        "sites_excluded_before_testing": 2,
        "sites_with_failed_model_fit": 0,
        "sites_included_in_multiple_testing_family": 2,
    }
    row_attrition = result.workflow_provenance["row_attrition"]
    assert [record["reason"] for record in row_attrition["records"]] == [
        "sites_excluded_before_testing"
    ]
    assert [record["removed_rows"] for record in row_attrition["records"]] == [2]
    assert (
        result.workflow_provenance["imputation_inference"][
            "adjusted_p_value_denominator_feature_count"
        ]
        == 2
    )


def test_differential_withhold_policy_reports_when_all_imputed_rows_are_withheld() -> (
    None
):
    result = DifferentialAnalysisWorkflow().run(
        DifferentialAnalysisRequest(
            dataset=_imputed_dataset(),
            design=_design(),
            contrasts=_contrasts(),
            config=DifferentialAnalysisConfig(
                imputed_value_policy="withhold_imputed_features",
                imputed_value_max_fraction=0.0,
            ),
        )
    )

    table = result.table_for("B_vs_A")
    caveat = _caveat_by_code(
        result,
        DIFFERENTIAL_IMPUTATION_WITHHOLDING_POLICY_CAVEAT_CODE,
    )

    assert table["result_status"].tolist() == [
        DIFFERENTIAL_RESULT_STATUS_TESTED,
        DIFFERENTIAL_RESULT_STATUS_WITHHELD_HIGH_IMPUTATION,
        DIFFERENTIAL_RESULT_STATUS_WITHHELD_HIGH_IMPUTATION,
        DIFFERENTIAL_RESULT_STATUS_WITHHELD_INSUFFICIENT_OBSERVED,
    ]
    assert caveat.severity == "warning"
    assert "no tested rows contained imputed values" in caveat.message
    assert caveat.details["tested_feature_count"] == 1
    assert caveat.details["withheld_feature_count"] == 3
    assert caveat.details["tested_imputed_feature_count"] == 0
    assert caveat.details["tested_imputed_cell_count"] == 0
    assert caveat.details["inferential_status"] == (
        DIFFERENTIAL_IMPUTATION_INFERENCE_STATUS_NO_TESTED_IMPUTED_VALUES
    )


def test_differential_withhold_policy_does_not_warn_when_no_analysed_cells_are_imputed() -> (
    None
):
    index = _site_index()
    fully_observed_mask = pd.DataFrame(
        True, index=index.copy(), columns=pd.Index(_SAMPLES)
    )
    result = DifferentialAnalysisWorkflow().run(
        DifferentialAnalysisRequest(
            dataset=_imputed_dataset(observed_mask=fully_observed_mask),
            design=_design(),
            contrasts=_contrasts(),
            config=DifferentialAnalysisConfig(
                imputed_value_policy="withhold_imputed_features",
                imputed_value_max_fraction=0.20,
            ),
        )
    )

    table = result.table_for("B_vs_A")
    caveat = _caveat_by_code(
        result,
        DIFFERENTIAL_IMPUTATION_WITHHOLDING_POLICY_CAVEAT_CODE,
    )

    assert table["imputed_cell_count"].tolist() == [0, 0, 0, 0]
    assert table["contains_imputed_cells"].tolist() == [False, False, False, False]
    assert caveat.severity == "info"
    assert "no analysed rows contained imputed values" in caveat.message
    assert caveat.details["tested_imputed_feature_count"] == 0
    assert caveat.details["tested_imputed_cell_count"] == 0


def test_differential_imputation_inference_metadata_serializes_and_reconstructs() -> (
    None
):
    result = DifferentialAnalysisWorkflow().run(_withhold_request())
    payload = result.to_payload()
    restored_caveats = result_caveats_from_payloads(payload["caveats"])
    restored = [
        caveat
        for caveat in restored_caveats
        if caveat.code == DIFFERENTIAL_IMPUTATION_WITHHOLDING_POLICY_CAVEAT_CODE
    ]

    assert len(restored) == 1
    assert restored[0].details["tested_imputed_feature_count"] == 1
    policy_payload = payload["policy_provenance"]
    assert policy_payload["missing_values"]["tested_feature_count"] == 2
    assert policy_payload["missing_values"]["tested_imputed_feature_count"] == 1
    assert policy_payload["missing_values"]["tested_imputed_cell_count"] == 1
    assert policy_payload["missing_values"]["observed_only_fit"] is False
    assert (
        policy_payload["missing_values"]["residual_df_adjusted_for_imputation"] is False
    )
    assert policy_payload["missing_values"]["inferential_status"] == (
        DIFFERENTIAL_IMPUTATION_INFERENCE_STATUS_RETAINED_IMPUTED_VALUES
    )


def test_differential_imputation_inference_metadata_applies_to_multiple_contrasts() -> (
    None
):
    result = DifferentialAnalysisWorkflow().run(
        DifferentialAnalysisRequest(
            dataset=_imputed_dataset(),
            design=_design(),
            contrasts=(
                *_contrasts(),
                Contrast(
                    name="A_vs_B",
                    numerator_condition="A",
                    denominator_condition="B",
                ),
            ),
            config=DifferentialAnalysisConfig(
                imputed_value_policy="withhold_imputed_features",
                imputed_value_max_fraction=0.20,
            ),
        )
    )

    first = result.table_for("B_vs_A")
    second = result.table_for("A_vs_B")

    assert first["inferential_status"].tolist() == second["inferential_status"].tolist()
    caveat = _caveat_by_code(
        result,
        DIFFERENTIAL_IMPUTATION_WITHHOLDING_POLICY_CAVEAT_CODE,
    )
    assert caveat.details["tested_feature_count"] == 2
    assert caveat.details["tested_imputed_feature_count"] == 1


def test_differential_imputation_policy_is_recorded_in_provenance() -> None:
    result = DifferentialAnalysisWorkflow().run(_withhold_request())

    assert result.policy_provenance is not None
    missing_values = result.policy_provenance.missing_values
    assert missing_values.imputed_value_policy == "withhold_imputed_features"
    assert missing_values.imputed_value_max_fraction == pytest.approx(0.20)
    assert missing_values.imputation_metadata_required is True
    assert missing_values.adjusted_p_value_scope == (
        "adjustment_over_tested_features_only_per_contrast"
    )
    assert missing_values.tested_feature_count == 2
    assert missing_values.withheld_feature_count == 2
    assert missing_values.tested_imputed_feature_count == 1
    assert missing_values.tested_imputed_cell_count == 1
    assert missing_values.observed_only_fit is False
    assert missing_values.residual_df_adjusted_for_imputation is False
    assert missing_values.inferential_status == (
        DIFFERENTIAL_IMPUTATION_INFERENCE_STATUS_RETAINED_IMPUTED_VALUES
    )
    assert missing_values.adjusted_p_value_denominator_feature_count == 2


def test_differential_imputation_policy_output_is_reproducible() -> None:
    first = DifferentialAnalysisWorkflow().run(_withhold_request()).table_for("B_vs_A")
    second = DifferentialAnalysisWorkflow().run(_withhold_request()).table_for("B_vs_A")

    pd.testing.assert_frame_equal(first, second)
