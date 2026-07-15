from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from phospy.api import (
    Contrast,
    DifferentialAnalysisConfig,
    DifferentialAnalysisRequest,
    DifferentialAnalysisWorkflow,
    ExperimentalDesign,
    MultipleTestingConfig,
    Organism,
    SampleDesignRecord,
)
from phospy.errors import WorkflowBoundaryError
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.science.differential.executor import (
    DifferentialAnalysisExecutor as DifferentialComputationExecutor,
)
from phospy.science.differential.models import (
    DIFFERENTIAL_RESULT_STATUS_COLUMN,
    DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN,
    DIFFERENTIAL_RESULT_STATUS_TESTED,
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_ALL_CONSTANT,
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_OTHER,
)
from phospy.science.statistics.multiple_testing import adjust_p_values
from phospy.workflows.differential.executor import (
    DifferentialAnalysisExecutor as DifferentialWorkflowExecutor,
)
from tests.support.intensity_scale_states import (
    supported_log2_intensity_scale_state,
    supported_log2_processing_state,
)
from tests.support.site_keys import protein_site_key_index, site_key_context_columns

_GENES = ["MAPK14", "AKT1", "GSK3B"]
_SITES = ["Y182", "T308", "S9"]


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
    return _dataset_from_matrix(
        matrix=_matrix(),
        genes=tuple(_GENES),
        sites=tuple(_SITES),
        site_index=site_index,
    )


def _dataset_from_matrix(
    *,
    matrix: pd.DataFrame,
    genes: tuple[str, ...],
    sites: tuple[str, ...],
    site_index: pd.Index,
) -> AnalysisReadyPhosphoDataset:
    site_metadata = pd.DataFrame(
        {
            "site_key": site_index.astype(str).tolist(),
            "display_id": [
                f"{gene};{site};" for gene, site in zip(genes, sites, strict=True)
            ],
            **site_key_context_columns(site_index),
            "gene_symbol": list(genes),
            "site": list(sites),
            "site_sequence": [("A" * 15) + site[0] + ("A" * 15) for site in sites],
            "protein_id": list(genes),
        },
        index=site_index.copy(),
    )
    return AnalysisReadyPhosphoDataset(
        phospho=matrix,
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


def test_differential_executor_receives_only_eligible_rows() -> None:
    class _ComputationExecutorSpy:
        def __init__(self) -> None:
            self.received_feature_ids: tuple[str, ...] = ()
            self._real_executor = DifferentialComputationExecutor()

        def run(self, request):
            self.received_feature_ids = tuple(
                str(feature_id) for feature_id in request.matrix.index.tolist()
            )
            return self._real_executor.run(request)

    computation_executor = _ComputationExecutorSpy()

    DifferentialAnalysisWorkflow(
        executor=DifferentialWorkflowExecutor(
            computation_executor=computation_executor,  # type: ignore[arg-type]
        )
    ).run(_request())

    assert computation_executor.received_feature_ids == tuple(
        str(feature_id) for feature_id in _site_index()[1:].tolist()
    )


def test_differential_result_preserves_original_row_alignment() -> None:
    result = DifferentialAnalysisWorkflow().run(_request())
    table = result.table_for("B_vs_A")
    expected_index = _site_index().tolist()

    assert table.index.tolist() == expected_index
    assert table.loc[:, "site_key"].tolist() == expected_index
    assert result.residual_variance_series().index.tolist() == expected_index
    assert result.feature_eligibility is not None
    assert result.feature_eligibility.index.tolist() == expected_index


def test_differential_adjusted_p_values_exclude_withheld_rows() -> None:
    base_request = _request()
    result = DifferentialAnalysisWorkflow().run(
        DifferentialAnalysisRequest(
            dataset=base_request.dataset,
            design=base_request.design,
            contrasts=base_request.contrasts,
            config=DifferentialAnalysisConfig(
                multiple_testing=MultipleTestingConfig(method="bonferroni")
            ),
        )
    )
    table = result.table_for("B_vs_A")
    tested = (
        table[DIFFERENTIAL_RESULT_STATUS_COLUMN] == DIFFERENTIAL_RESULT_STATUS_TESTED
    )

    expected_adjusted = adjust_p_values(
        table.loc[tested, "P.Value"].to_numpy(dtype=float),
        method="bonferroni",
    )
    np.testing.assert_allclose(
        table.loc[tested, "adj.P.Val"].to_numpy(dtype=float),
        expected_adjusted,
        rtol=1e-12,
        atol=1e-12,
    )
    assert table.loc[~tested, "adj.P.Val"].isna().all()


def test_differential_withheld_rows_have_null_statistics() -> None:
    result = DifferentialAnalysisWorkflow().run(_request())
    table = result.table_for("B_vs_A")
    withheld = (
        table[DIFFERENTIAL_RESULT_STATUS_COLUMN] != DIFFERENTIAL_RESULT_STATUS_TESTED
    )

    assert withheld.any()
    assert (
        table.loc[withheld, ["logFC", "t", "P.Value", "adj.P.Val"]].isna().all().all()
    )
    assert (
        table.loc[withheld, DIFFERENTIAL_RESULT_STATUS_COLUMN].astype(str).ne("").all()
    )
    assert (
        table.loc[withheld, DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN]
        .astype(str)
        .str.strip()
        .ne("")
        .all()
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


def test_differential_provenance_counts_tested_and_failed_model_fit_sites() -> None:
    genes = ("MAPK14", "AKT1", "GSK3B", "PRKACA")
    sites = ("Y182", "T308", "S9", "S339")
    site_index = protein_site_key_index(
        protein_identifiers=genes,
        sites=sites,
    )
    matrix = pd.DataFrame(
        {
            "A_1": [5.0, 1.0, 2.0, 3.0],
            "A_2": [5.0, 1.0, 2.2, 3.1],
            "B_1": [5.0, 2.0, 1.8, 4.2],
            "B_2": [5.0, 2.0, 2.0, 4.0],
        },
        index=site_index,
    )
    base_request = _request()

    result = DifferentialAnalysisWorkflow().run(
        DifferentialAnalysisRequest(
            dataset=_dataset_from_matrix(
                matrix=matrix,
                genes=genes,
                sites=sites,
                site_index=site_index,
            ),
            design=base_request.design,
            contrasts=base_request.contrasts,
            config=DifferentialAnalysisConfig(
                multiple_testing=MultipleTestingConfig(method="bonferroni")
            ),
        )
    )

    table = result.table_for("B_vs_A")
    failed_site = site_index[1]
    tested = (
        table[DIFFERENTIAL_RESULT_STATUS_COLUMN] == DIFFERENTIAL_RESULT_STATUS_TESTED
    )

    assert table[DIFFERENTIAL_RESULT_STATUS_COLUMN].tolist() == [
        DIFFERENTIAL_RESULT_STATUS_WITHHELD_ALL_CONSTANT,
        DIFFERENTIAL_RESULT_STATUS_WITHHELD_OTHER,
        DIFFERENTIAL_RESULT_STATUS_TESTED,
        DIFFERENTIAL_RESULT_STATUS_TESTED,
    ]
    assert (
        table.at[failed_site, DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN]
        == "Feature model fit failed before multiple-testing correction; residual "
        "variance was zero or non-finite."
    )
    assert table.loc[failed_site, ["logFC", "t", "P.Value", "adj.P.Val"]].isna().all()

    expected_adjusted = adjust_p_values(
        table.loc[tested, "P.Value"].to_numpy(dtype=float),
        method="bonferroni",
    )
    np.testing.assert_allclose(
        table.loc[tested, "adj.P.Val"].to_numpy(dtype=float),
        expected_adjusted,
        rtol=1e-12,
        atol=1e-12,
    )
    assert table.loc[~tested, "adj.P.Val"].isna().all()

    assert result.workflow_provenance is not None
    metrics = result.workflow_provenance["row_attrition_metrics"]
    assert metrics == {
        "input_sites": 4,
        "sites_retained_for_model_fitting": 3,
        "sites_excluded_before_testing": 1,
        "sites_with_failed_model_fit": 1,
        "sites_included_in_multiple_testing_family": 2,
    }
    row_attrition = result.workflow_provenance["row_attrition"]
    assert row_attrition["input_rows"] == 4
    assert row_attrition["final_rows"] == 2
    assert [record["reason"] for record in row_attrition["records"]] == [
        "sites_excluded_before_testing",
        "failed_model_fit",
    ]
    assert [record["removed_rows"] for record in row_attrition["records"]] == [1, 1]

    payload = result.to_payload()
    payload_provenance = payload["workflow_provenance"]
    assert payload_provenance["row_attrition"]["records"] == row_attrition["records"]


@pytest.mark.filterwarnings("ignore:overflow encountered in square:RuntimeWarning")
def test_differential_marks_non_finite_residual_variance_rows_as_failed_model_fit() -> (
    None
):
    genes = ("MAPK14", "AKT1", "GSK3B", "PRKACA", "RPS6")
    sites = ("Y182", "T308", "S9", "S339", "S235")
    site_index = protein_site_key_index(
        protein_identifiers=genes,
        sites=sites,
    )
    matrix = pd.DataFrame(
        {
            "A_1": [5.0, 1.0, 1e308, 2.0, 3.0],
            "A_2": [5.0, 1.0, -1e308, 2.2, 3.1],
            "B_1": [5.0, 2.0, 1e308, 1.8, 4.2],
            "B_2": [5.0, 2.0, -1e308, 2.0, 4.0],
        },
        index=site_index,
    )
    base_request = _request()

    result = DifferentialAnalysisWorkflow().run(
        DifferentialAnalysisRequest(
            dataset=_dataset_from_matrix(
                matrix=matrix,
                genes=genes,
                sites=sites,
                site_index=site_index,
            ),
            design=base_request.design,
            contrasts=base_request.contrasts,
            config=DifferentialAnalysisConfig(
                multiple_testing=MultipleTestingConfig(method="bonferroni")
            ),
        )
    )

    table = result.table_for("B_vs_A")
    tested = (
        table[DIFFERENTIAL_RESULT_STATUS_COLUMN] == DIFFERENTIAL_RESULT_STATUS_TESTED
    )

    assert table[DIFFERENTIAL_RESULT_STATUS_COLUMN].tolist() == [
        DIFFERENTIAL_RESULT_STATUS_WITHHELD_ALL_CONSTANT,
        DIFFERENTIAL_RESULT_STATUS_WITHHELD_OTHER,
        DIFFERENTIAL_RESULT_STATUS_WITHHELD_OTHER,
        DIFFERENTIAL_RESULT_STATUS_TESTED,
        DIFFERENTIAL_RESULT_STATUS_TESTED,
    ]
    assert (
        table.loc[site_index[2], DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN]
        == "Feature model fit failed before multiple-testing correction; residual "
        "variance was zero or non-finite."
    )
    assert table.loc[site_index[2], ["logFC", "t", "P.Value", "adj.P.Val"]].isna().all()
    expected_adjusted = adjust_p_values(
        table.loc[tested, "P.Value"].to_numpy(dtype=float),
        method="bonferroni",
    )
    np.testing.assert_allclose(
        table.loc[tested, "adj.P.Val"].to_numpy(dtype=float),
        expected_adjusted,
        rtol=1e-12,
        atol=1e-12,
    )

    assert result.workflow_provenance is not None
    assert result.workflow_provenance["row_attrition_metrics"] == {
        "input_sites": 5,
        "sites_retained_for_model_fitting": 4,
        "sites_excluded_before_testing": 1,
        "sites_with_failed_model_fit": 2,
        "sites_included_in_multiple_testing_family": 2,
    }


def test_differential_all_rows_filtered_before_model_fit_fails_before_executor() -> (
    None
):
    class _ExecutorSpy:
        calls = 0

        def run(self, request):
            self.calls += 1
            raise AssertionError("executor should not run when no rows are testable")

    matrix = pd.DataFrame(
        {
            "A_1": [5.0, 10.0, 20.0],
            "A_2": [5.0, 10.0, 20.0],
            "B_1": [5.0, 10.0, 20.0],
            "B_2": [5.0, 10.0, 20.0],
        },
        index=_site_index(),
    )
    executor = _ExecutorSpy()
    with pytest.raises(
        WorkflowBoundaryError,
        match="differential.interpreter.feature_eligibility",
    ) as exc_info:
        DifferentialAnalysisWorkflow(executor=executor).run(  # type: ignore[arg-type]
            DifferentialAnalysisRequest(
                dataset=_dataset_from_matrix(
                    matrix=matrix,
                    genes=tuple(_GENES),
                    sites=tuple(_SITES),
                    site_index=_site_index(),
                ),
                design=_request().design,
                contrasts=_request().contrasts,
            )
        )

    assert executor.calls == 0
    assert exc_info.value.details["status_counts"] == {
        DIFFERENTIAL_RESULT_STATUS_WITHHELD_ALL_CONSTANT: 3
    }
