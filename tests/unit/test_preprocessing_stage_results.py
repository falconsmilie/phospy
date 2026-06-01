from __future__ import annotations

from collections.abc import Mapping

import pandas as pd
import pytest

from phospy.errors.build import DatasetBuildError
from phospy.errors.input import PhosPyInputError
from phospy.science.datasets.builders.contracts import InterpretedDatasetBuildRequest
from phospy.science.datasets.builders.executor import DatasetBuildExecutor
from phospy.science.datasets.builders.preprocessing import (
    DatasetPreprocessor,
    build_dataset_processing_state,
)
from phospy.science.datasets.preprocessing.models import (
    PreprocessingPlan,
    PreprocessingReportRow,
    PreprocessingStageResult,
    PreprocessingState,
)
from phospy.science.datasets.preprocessing.pipeline import PreprocessingPipeline
from phospy.science.datasets.preprocessing.report_schema import PreprocessingRowAuditRow
from phospy.science.datasets.preprocessing.stage_registry import (
    PreprocessingStageMetadata,
)
from phospy.science.datasets.preprocessing.stages.comparisons import ComparisonsStage
from phospy.science.datasets.preprocessing.stages.intensity_transform import (
    IntensityTransformStage,
)
from phospy.science.datasets.preprocessing.stages.missing_data import MissingDataStage
from phospy.science.datasets.preprocessing.stages.normalisation import (
    NormalisationStage,
)
from phospy.science.datasets.preprocessing.stages.site_matrix import SiteMatrixStage
from phospy.science.datasets.preprocessing.stages.total_protein_correction import (
    TotalProteinCorrectionStage,
)
from phospy.science.datasets.processing_state import MissingDataDiagnosticsV1
from phospy.science.transformations.models import (
    IntensityScaleKind,
    QuantitativeMeaning,
)
from tests.support.intensity_scale_states import supported_linear_intensity_scale_state
from tests.support.site_keys import site_key_index_from_display_ids

_ANALYSIS_DISPLAY_IDS = ["MAPK14;Y182;", "AKT1;T308;"]
_ANALYSIS_GENE_SYMBOLS = ["MAPK14", "AKT1"]
_ANALYSIS_SITES = ["Y182", "T308"]


def _analysis_site_index() -> pd.Index:
    return site_key_index_from_display_ids(_ANALYSIS_DISPLAY_IDS)


def _analysis_site_metadata(index: pd.Index) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "site_key": index.astype(str).tolist(),
            "display_id": _ANALYSIS_DISPLAY_IDS,
            "gene_symbol": _ANALYSIS_GENE_SYMBOLS,
            "site": _ANALYSIS_SITES,
            "site_sequence": [
                ("A" * 15) + site[0] + ("A" * 15) for site in _ANALYSIS_SITES
            ],
            "localisation_confidence": [0.95, 0.9],
        },
        index=index.copy(),
    )


def _phospho() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample_a": [4.0, 8.0],
            "sample_b": [6.0, 12.0],
        },
        index=pd.Index(["row_a", "row_b"], name="source_row"),
    )


def _site_metadata(index: pd.Index) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "AKT1"],
            "site": ["Y182", "T308"],
            "site_sequence": ["SEQ_A", "SEQ_R"],
            "localisation_confidence": [0.95, 0.9],
        },
        index=index.copy(),
    )


def _sample_metadata(columns: pd.Index) -> pd.DataFrame:
    return pd.DataFrame(
        {"comparison_group": ["group1", "group2"]},
        index=columns.copy(),
    )


def _total(columns: pd.Index) -> pd.DataFrame:
    return pd.DataFrame(
        {
            str(columns[0]): [2.0, 4.0],
            str(columns[1]): [3.0, 5.0],
        },
        index=pd.Index(["MAPK14", "AKT1"], name="protein_id"),
    )


def _custom_stage_metadata(stage_key: str) -> PreprocessingStageMetadata:
    return PreprocessingStageMetadata(
        stage_key=stage_key,
        display_label=stage_key,
        provenance_stage=stage_key,
        operation_name=lambda _plan: stage_key,
        serialize_parameters=lambda _plan: {},
        consumed_input_tables=("dataset.phospho",),
        produced_output_tables=("dataset.phospho",),
    )


def _assert_stage_result_contract(result: PreprocessingStageResult) -> None:
    assert isinstance(result, PreprocessingStageResult)
    assert result.state is not None
    assert isinstance(result.diagnostics, Mapping)
    assert isinstance(tuple(result.report_rows), tuple)


def test_missing_data_stage_returns_stage_result() -> None:
    phospho = _phospho()
    phospho.loc["row_b", "sample_b"] = float("nan")
    state = PreprocessingState(
        phospho=phospho,
        site_metadata=_site_metadata(phospho.index),
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan(
            missing_data_policy="impute_row_median",
            missing_data_min_observed_values=1,
            stage_order=("missing_data",),
        ),
    )

    result = MissingDataStage().run(state)

    _assert_stage_result_contract(result)


def test_missing_data_stage_emits_typed_row_audit_report_rows() -> None:
    phospho = _phospho()
    phospho.loc["row_b", "sample_b"] = float("nan")
    state = PreprocessingState(
        phospho=phospho,
        site_metadata=_site_metadata(phospho.index),
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan(
            missing_data_policy="impute_row_median",
            missing_data_min_observed_values=1,
            stage_order=("missing_data",),
        ),
    )

    result = MissingDataStage().run(state)

    assert result.report_rows
    assert all(row.table == "row_audit" for row in result.report_rows)
    assert all(
        isinstance(row.values, PreprocessingRowAuditRow) for row in result.report_rows
    )


def test_missing_data_stage_forbid_policy_accepts_complete_matrix() -> None:
    phospho = _phospho()
    state = PreprocessingState(
        phospho=phospho,
        site_metadata=_site_metadata(phospho.index),
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan(
            missing_data_policy="forbid",
            stage_order=("missing_data",),
        ),
    )

    result = MissingDataStage().run(state)

    assert result.state is state
    assert result.diagnostics["dropped_row_count"] == 0
    assert result.diagnostics["imputed_cell_count"] == 0
    diagnostics = result.diagnostics["diagnostics"]
    assert diagnostics["input_missing_cell_count"] == 0
    assert diagnostics["output_missing_cell_count"] == 0
    assert diagnostics["imputed_cell_count"] == 0
    assert diagnostics["affected_row_count"] == 0
    assert diagnostics["affected_column_count"] == 0


def test_missing_data_stage_forbid_policy_rejects_single_missing_value() -> None:
    phospho = _phospho()
    phospho.loc["row_a", "sample_a"] = float("nan")
    state = PreprocessingState(
        phospho=phospho,
        site_metadata=_site_metadata(phospho.index),
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan(
            missing_data_policy="forbid",
            stage_order=("missing_data",),
        ),
    )

    with pytest.raises(PhosPyInputError) as exc_info:
        MissingDataStage().run(state)

    diagnostics = exc_info.value.diagnostics
    assert isinstance(diagnostics, MissingDataDiagnosticsV1)
    assert diagnostics.missing_data_policy == "forbid"
    assert diagnostics.imputation_method_id == "forbid"
    assert diagnostics.input_missing_cell_count == 1
    assert diagnostics.output_missing_cell_count == 1
    assert diagnostics.imputed_cell_count == 0
    assert diagnostics.affected_row_ids == ("row_a",)
    assert diagnostics.affected_column_ids == ("sample_a",)
    assert diagnostics.dropped_row_ids == ()
    assert diagnostics.random_seed is None
    assert diagnostics.method_parameters == {}
    assert diagnostics.matrix_scale_requirement is None
    assert diagnostics.stage_order == ("missing_data",)
    assert isinstance(diagnostics.missingness_mask_hash, str)
    assert diagnostics.missingness_mask_hash
    assert diagnostics.left_censored_assumption is False
    assert diagnostics.rows_not_imputable == ("row_a",)

    message = str(exc_info.value)
    assert "stage 'missing_data'" in message
    assert "missing_data.policy='forbid'" in message
    assert "found 1 missing values across 1 rows and 1 columns" in message
    assert "'row_a'" in message
    assert "'sample_a'" in message
    assert "choose missing_data.policy='impute_row_median'" in message


def test_missing_data_stage_forbid_policy_rejects_multiple_missing_values() -> None:
    phospho = _phospho()
    phospho.loc["row_a", "sample_a"] = float("nan")
    phospho.loc["row_a", "sample_b"] = float("nan")
    phospho.loc["row_b", "sample_b"] = float("nan")
    state = PreprocessingState(
        phospho=phospho,
        site_metadata=_site_metadata(phospho.index),
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan(
            missing_data_policy="forbid",
            stage_order=("missing_data",),
        ),
    )

    with pytest.raises(
        PhosPyInputError,
        match="found 3 missing values across 2 rows and 2 columns",
    ):
        MissingDataStage().run(state)


def test_missing_data_stage_report_rows_appear_in_final_report() -> None:
    site_index = _analysis_site_index()
    phospho = pd.DataFrame(
        {
            "sample_a": [4.0, 8.0],
            "sample_b": [6.0, float("nan")],
        },
        index=site_index,
    )
    site_metadata = _analysis_site_metadata(phospho.index)
    built = DatasetBuildExecutor().run(
        InterpretedDatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            sample_metadata=None,
            total=None,
            organism=None,
            declared_input_intensity_scale_kind=IntensityScaleKind.LINEAR,
            preprocessing_plan=PreprocessingPlan(
                missing_data_policy="impute_row_median",
                missing_data_min_observed_values=1,
                stage_order=("missing_data",),
            ),
        )
    )

    assert built.preprocessing_report is not None
    row_audit = built.preprocessing_report.row_audit
    assert not row_audit.empty
    assert "missing_data" in set(row_audit.loc[:, "stage"].astype(str))


def test_missing_data_stage_operations_report_imputation_summary_note() -> None:
    site_index = _analysis_site_index()
    phospho = pd.DataFrame(
        {
            "sample_a": [4.0, 8.0],
            "sample_b": [6.0, float("nan")],
        },
        index=site_index,
    )
    site_metadata = _analysis_site_metadata(phospho.index)

    built = DatasetBuildExecutor().run(
        InterpretedDatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            sample_metadata=None,
            total=None,
            organism=None,
            declared_input_intensity_scale_kind=IntensityScaleKind.LINEAR,
            preprocessing_plan=PreprocessingPlan(
                missing_data_policy="impute_row_median",
                missing_data_min_observed_values=1,
                stage_order=("missing_data",),
            ),
        )
    )

    assert built.preprocessing_report is not None
    operations = built.preprocessing_report.operations
    missing_data_row = operations.loc[
        operations.loc[:, "stage"].astype(str) == "missing_data"
    ].iloc[0]
    note = str(missing_data_row["notes"])
    assert "policy='impute_row_median'" in note
    assert "imputed_cells=1" in note
    assert "imputed_rows=1" in note
    assert "output_missing_cells=0" in note


def test_final_dataset_has_complete_matrix_after_missing_data_imputation() -> None:
    site_index = _analysis_site_index()
    phospho = pd.DataFrame(
        {
            "sample_a": [4.0, 8.0],
            "sample_b": [6.0, float("nan")],
        },
        index=site_index,
    )
    site_metadata = _analysis_site_metadata(phospho.index)

    built = DatasetBuildExecutor().run(
        InterpretedDatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            sample_metadata=None,
            total=None,
            organism=None,
            declared_input_intensity_scale_kind=IntensityScaleKind.LINEAR,
            preprocessing_plan=PreprocessingPlan(
                missing_data_policy="impute_row_median",
                missing_data_min_observed_values=1,
                stage_order=("missing_data",),
            ),
        )
    )

    assert int(built.phospho.isna().to_numpy().sum()) == 0
    assert built.processing_state.missing_data.complete_matrix is True


def test_intensity_transform_stage_returns_stage_result() -> None:
    phospho = _phospho()
    state = PreprocessingState(
        phospho=phospho,
        site_metadata=_site_metadata(phospho.index),
        sample_metadata=None,
        total=_total(phospho.columns),
        plan=PreprocessingPlan(
            intensity_transform_policy="log2",
            intensity_transform_pseudocount=1.0,
            stage_order=("intensity_transform",),
        ),
    )

    result = IntensityTransformStage().run(state)

    _assert_stage_result_contract(result)


def test_normalisation_stage_returns_stage_result() -> None:
    phospho = _phospho()
    state = PreprocessingState(
        phospho=phospho,
        site_metadata=_site_metadata(phospho.index),
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan(
            normalisation_policy="median_center",
            stage_order=("normalisation",),
        ),
    )

    result = NormalisationStage().run(state)

    _assert_stage_result_contract(result)


def test_comparisons_stage_returns_stage_result() -> None:
    phospho = _phospho()
    state = PreprocessingState(
        phospho=phospho,
        site_metadata=_site_metadata(phospho.index),
        sample_metadata=_sample_metadata(phospho.columns),
        total=None,
        plan=PreprocessingPlan(
            comparison_building_policy="sample_metadata_pairs",
            stage_order=("comparisons",),
        ),
    )

    result = ComparisonsStage().run(state)

    _assert_stage_result_contract(result)


def test_total_protein_stage_returns_stage_result() -> None:
    phospho = _phospho()
    state = PreprocessingState(
        phospho=phospho,
        site_metadata=_site_metadata(phospho.index),
        sample_metadata=None,
        total=_total(phospho.columns),
        plan=PreprocessingPlan(
            intensity_transform_policy="log2",
            total_protein_correction_policy="subtract_log_total",
            stage_order=("total_protein_correction",),
        ),
    )

    result = TotalProteinCorrectionStage().run(state)

    _assert_stage_result_contract(result)


def test_site_matrix_stage_returns_stage_result() -> None:
    phospho = _phospho()
    state = PreprocessingState(
        phospho=phospho,
        site_metadata=_site_metadata(phospho.index),
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan(
            site_matrix_policy="build_from_metadata",
            stage_order=("site_matrix",),
        ),
    )

    result = SiteMatrixStage().run(state)

    _assert_stage_result_contract(result)


def test_pipeline_uses_stage_owned_diagnostics_and_report_rows() -> None:
    class FakeStage:
        stage_key = "fake_stage"

        def run(self, state: PreprocessingState) -> PreprocessingStageResult:
            return PreprocessingStageResult(
                state=state,
                diagnostics={"custom_metric": 123},
                report_rows=(
                    PreprocessingReportRow(
                        table="row_audit",
                        values=PreprocessingRowAuditRow(
                            stage="fake_stage",
                            action="retained",
                            reason="test emission",
                            source_row_id="row_a",
                            site_id="row_a",
                            retained=True,
                            retained_row_id="row_a",
                            source_rows=("row_a",),
                            retained_row="row_a",
                            parameter_snapshot={"source": "fake_stage"},
                        ),
                    ),
                ),
            )

    phospho = _phospho()
    state = PreprocessingState(
        phospho=phospho,
        site_metadata=_site_metadata(phospho.index),
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan(stage_order=("fake_stage",)),
    )

    final_state, trace = PreprocessingPipeline(
        stage_registry=(FakeStage(),),
        stage_metadata_registry=(_custom_stage_metadata("fake_stage"),),
    ).run_with_trace(state)

    assert len(trace) == 1
    assert trace[0].stage == "fake_stage"
    assert trace[0].diagnostics["custom_metric"] == 123
    assert len(final_state.report_rows) == 1
    assert final_state.report_rows[0].table == "row_audit"
    assert isinstance(final_state.report_rows[0].values, PreprocessingRowAuditRow)


def test_pipeline_rejects_unsupported_stage_report_rows() -> None:
    class FakeStage:
        stage_key = "fake_stage"

        def run(self, state: PreprocessingState) -> PreprocessingStageResult:
            return PreprocessingStageResult(
                state=state,
                report_rows=(
                    PreprocessingReportRow(
                        table="custom_table",
                        values=PreprocessingRowAuditRow(
                            stage="fake_stage",
                            action="retained",
                            reason="test emission",
                            source_row_id="row_a",
                            site_id="row_a",
                            retained=True,
                            retained_row_id="row_a",
                            source_rows=("row_a",),
                            retained_row="row_a",
                            parameter_snapshot={"source": "fake_stage"},
                        ),
                    ),
                ),
            )

    phospho = _phospho()
    state = PreprocessingState(
        phospho=phospho,
        site_metadata=_site_metadata(phospho.index),
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan(stage_order=("fake_stage",)),
    )

    with pytest.raises(
        DatasetBuildError,
        match="unsupported table",
    ):
        PreprocessingPipeline(
            stage_registry=(FakeStage(),),
            stage_metadata_registry=(_custom_stage_metadata("fake_stage"),),
        ).run_with_trace(state)


def test_minimal_custom_stage_emits_supported_report_row_into_final_report() -> None:
    class FakeStage:
        stage_key = "fake_stage"

        def run(self, state: PreprocessingState) -> PreprocessingStageResult:
            return PreprocessingStageResult(
                state=state,
                report_rows=(
                    PreprocessingReportRow(
                        table="row_audit",
                        values=PreprocessingRowAuditRow(
                            stage="fake_stage",
                            action="retained",
                            reason="custom stage retained row",
                            source_row_id="MAPK14;Y182;",
                            site_id="MAPK14;Y182;",
                            retained=True,
                            retained_row_id="MAPK14;Y182;",
                            source_rows=("MAPK14;Y182;",),
                            retained_row="MAPK14;Y182;",
                            parameter_snapshot={"policy": "test"},
                        ),
                    ),
                ),
                diagnostics={"notes": "stage executed"},
            )

    site_index = _analysis_site_index()
    phospho = pd.DataFrame(
        {
            "sample_a": [4.0, 8.0],
            "sample_b": [6.0, 12.0],
        },
        index=site_index,
    )
    site_metadata = _analysis_site_metadata(phospho.index)
    executor = DatasetBuildExecutor(
        preprocessor=DatasetPreprocessor(
            pipeline=PreprocessingPipeline(
                stage_registry=(FakeStage(),),
                stage_metadata_registry=(_custom_stage_metadata("fake_stage"),),
            ),
        )
    )
    built = executor.run(
        InterpretedDatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            sample_metadata=None,
            total=None,
            organism=None,
            declared_input_intensity_scale_kind=IntensityScaleKind.LINEAR,
            preprocessing_plan=PreprocessingPlan(stage_order=("fake_stage",)),
        )
    )

    assert built.preprocessing_report is not None
    report = built.preprocessing_report
    row_audit = report.row_audit
    assert not row_audit.empty
    assert set(row_audit.loc[:, "stage"].astype(str)) == {"fake_stage"}
    assert row_audit.iloc[0]["source_row_id"] == "MAPK14;Y182;"
    assert "final_dataset_construction" in set(report.row_counts.loc[:, "stage"])
    assert "final_dataset_construction" in set(report.operations.loc[:, "stage"])


def test_executor_applies_explicit_quantitative_meaning_to_dataset_and_provenance() -> (
    None
):
    site_index = _analysis_site_index()
    phospho = pd.DataFrame(
        {
            "sample_a": [0.5, -0.2],
            "sample_b": [1.1, 0.0],
        },
        index=site_index,
    )
    site_metadata = _analysis_site_metadata(phospho.index)
    built = DatasetBuildExecutor().run(
        InterpretedDatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            sample_metadata=None,
            total=None,
            organism=None,
            declared_input_intensity_scale_kind=IntensityScaleKind.LINEAR,
            preprocessing_plan=PreprocessingPlan.default(),
            quantitative_meaning=QuantitativeMeaning.CONTRAST_LOG2_FOLD_CHANGE,
        )
    )

    assert built.intensity_scale_state.quantity is not None
    assert (
        built.intensity_scale_state.quantity.value
        == QuantitativeMeaning.CONTRAST_LOG2_FOLD_CHANGE.value
    )
    assert (
        built.processing_state.total_protein_correction.quantitative_meaning
        == QuantitativeMeaning.CONTRAST_LOG2_FOLD_CHANGE.value
    )
    assert built.preprocessing_report is not None
    final_operation = built.preprocessing_report.operations.loc[
        built.preprocessing_report.operations.loc[:, "stage"]
        == "final_dataset_construction"
    ].iloc[0]
    assert final_operation["parameters"]["quantitative_meaning"] == (
        QuantitativeMeaning.CONTRAST_LOG2_FOLD_CHANGE.value
    )
    assert built.provenance is not None
    assert (
        built.provenance.workflow_parameters["quantitative_meaning"]
        == QuantitativeMeaning.CONTRAST_LOG2_FOLD_CHANGE.value
    )


def test_pipeline_trace_preserves_intensity_transform_diagnostics() -> None:
    phospho = _phospho()
    state = PreprocessingState(
        phospho=phospho,
        site_metadata=_site_metadata(phospho.index),
        sample_metadata=None,
        total=_total(phospho.columns),
        plan=PreprocessingPlan(
            intensity_transform_policy="log2",
            intensity_transform_pseudocount=1.0,
            stage_order=("intensity_transform",),
        ),
    )

    _, trace = PreprocessingPipeline().run_with_trace(state)

    diagnostics = trace[0].diagnostics
    assert diagnostics["policy"] == "log2"
    assert diagnostics["pseudocount"] == 1.0
    assert isinstance(diagnostics.get("input_phospho_hash"), str)
    assert isinstance(diagnostics.get("output_phospho_hash"), str)


def test_pipeline_trace_preserves_normalisation_diagnostics() -> None:
    phospho = _phospho()
    state = PreprocessingState(
        phospho=phospho,
        site_metadata=_site_metadata(phospho.index),
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan(
            normalisation_policy="median_center",
            stage_order=("normalisation",),
        ),
    )

    _, trace = PreprocessingPipeline().run_with_trace(state)

    diagnostics = trace[0].diagnostics
    assert diagnostics["method"] == "median_center"
    assert diagnostics["parameters"] == {
        "applied": True,
        "centering_statistic": "median",
        "axis": "columns",
        "skipna": True,
    }
    assert diagnostics["policy"] == "median_center"
    assert diagnostics["note"] == "median centering used"
    assert diagnostics["input_matrix_shape"] == {"rows": 2, "columns": 2}
    assert diagnostics["output_matrix_shape"] == {"rows": 2, "columns": 2}
    assert diagnostics["rows_dropped"] is False
    assert diagnostics["columns_dropped"] is False
    assert diagnostics["dropped_row_count"] == 0
    assert diagnostics["dropped_column_count"] == 0
    assert "sample_a" in diagnostics["per_sample_summary_before"]
    assert "sample_a" in diagnostics["per_sample_summary_after"]
    assert isinstance(diagnostics.get("input_phospho_hash"), str)
    assert isinstance(diagnostics.get("output_phospho_hash"), str)


def test_pipeline_trace_preserves_comparison_diagnostics() -> None:
    phospho = _phospho()
    state = PreprocessingState(
        phospho=phospho,
        site_metadata=_site_metadata(phospho.index),
        sample_metadata=_sample_metadata(phospho.columns),
        total=None,
        plan=PreprocessingPlan(
            comparison_building_policy="sample_metadata_pairs",
            stage_order=("comparisons",),
        ),
    )

    _, trace = PreprocessingPipeline().run_with_trace(state)

    diagnostics = trace[0].diagnostics
    assert diagnostics["policy"] == "sample_metadata_pairs"
    assert diagnostics["sample_group_column"] == "comparison_group"
    assert diagnostics["resolved_comparison_pairs"] == [("group1", "group2")]
    assert set(diagnostics["group_labels"]) == {"group1", "group2"}
    assert isinstance(diagnostics.get("output_comparison_hash"), str)


def test_missing_data_stage_row_median_emits_structured_diagnostics() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 2.0, float("nan")],
            "sample_b": [float("nan"), 3.0, float("nan")],
            "sample_c": [4.0, 5.0, float("nan")],
        },
        index=pd.Index(["row_a", "row_b", "row_c"], name="source_row"),
    )
    state = PreprocessingState(
        phospho=phospho,
        site_metadata=pd.DataFrame(
            {
                "gene_symbol": ["MAPK14", "AKT1", "GSK3B"],
                "site": ["Y182", "T308", "S9"],
                "site_sequence": ["SEQ_A", "SEQ_R", "SEQ_C"],
            },
            index=phospho.index.copy(),
        ),
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan(
            missing_data_policy="impute_row_median",
            missing_data_min_observed_values=2,
            stage_order=("missing_data",),
        ),
    )

    result = MissingDataStage().run(state)
    diagnostics = result.diagnostics["diagnostics"]
    assert float(result.state.phospho.loc["row_a", "sample_b"]) == 2.5

    assert diagnostics["diagnostics_schema_version"] == 1
    assert diagnostics["missing_data_policy"] == "impute_row_median"
    assert diagnostics["imputation_method_id"] == "row_median"
    assert diagnostics["imputation_method_family"] == "deterministic_row_statistic"
    assert diagnostics["input_missing_cell_count"] == 4
    assert diagnostics["output_missing_cell_count"] == 0
    assert diagnostics["imputed_cell_count"] == 1
    assert diagnostics["affected_row_count"] == 2
    assert diagnostics["affected_column_count"] == 3
    assert diagnostics["affected_row_ids"] == ["row_a", "row_c"]
    assert diagnostics["affected_column_ids"] == ["sample_a", "sample_b", "sample_c"]
    assert diagnostics["imputed_row_ids"] == ["row_a"]
    assert diagnostics["imputed_column_ids"] == ["sample_b"]
    assert diagnostics["imputed_row_count"] == 1
    assert diagnostics["imputed_column_count"] == 1
    assert diagnostics["dropped_row_ids"] == ["row_c"]
    assert diagnostics["dropped_row_count"] == 1
    assert diagnostics["random_seed"] is None
    assert diagnostics["matrix_scale_requirement"] is None
    assert diagnostics["stage_order"] == ["missing_data"]
    assert isinstance(diagnostics["missingness_mask_hash"], str)
    assert isinstance(diagnostics["imputation_mask_hash"], str)
    assert diagnostics["left_censored_assumption"] is False


def test_missing_data_stage_row_median_values_remain_unchanged() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, float("nan"), float("nan")],
            "sample_b": [2.0, 10.0, float("nan")],
            "sample_c": [3.0, 20.0, 9.0],
            "sample_d": [4.0, float("nan"), float("nan")],
        },
        index=pd.Index(["row_keep", "row_impute", "row_drop"], name="source_row"),
    )
    state = PreprocessingState(
        phospho=phospho,
        site_metadata=pd.DataFrame(
            {
                "gene_symbol": ["MAPK14", "AKT1", "GSK3B"],
                "site": ["Y182", "T308", "S9"],
                "site_sequence": ["SEQ_A", "SEQ_R", "SEQ_C"],
            },
            index=phospho.index.copy(),
        ),
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan(
            missing_data_policy="impute_row_median",
            missing_data_min_observed_values=2,
            stage_order=("missing_data",),
        ),
    )

    result = MissingDataStage().run(state)

    assert result.state.phospho.index.tolist() == ["row_keep", "row_impute"]
    assert float(result.state.phospho.loc["row_keep", "sample_a"]) == 1.0
    assert float(result.state.phospho.loc["row_keep", "sample_b"]) == 2.0
    assert float(result.state.phospho.loc["row_keep", "sample_c"]) == 3.0
    assert float(result.state.phospho.loc["row_keep", "sample_d"]) == 4.0
    assert float(result.state.phospho.loc["row_impute", "sample_a"]) == 15.0
    assert float(result.state.phospho.loc["row_impute", "sample_b"]) == 10.0
    assert float(result.state.phospho.loc["row_impute", "sample_c"]) == 20.0
    assert float(result.state.phospho.loc["row_impute", "sample_d"]) == 15.0

    diagnostics = result.diagnostics["diagnostics"]
    assert diagnostics["input_missing_cell_count"] == 5
    assert diagnostics["imputed_cell_count"] == 2
    assert diagnostics["affected_row_count"] == 2
    assert diagnostics["affected_column_count"] == 3
    assert diagnostics["dropped_row_ids"] == ["row_drop"]


def test_missing_data_stage_minprob_emits_distribution_diagnostics_and_drops_rows() -> (
    None
):
    phospho = pd.DataFrame(
        {
            "sample_a": [10.0, float("nan"), float("nan"), 4.0],
            "sample_b": [9.0, 8.0, float("nan"), 6.0],
            "sample_c": [11.0, 7.0, 5.0, float("nan")],
        },
        index=pd.Index(["row_keep", "row_impute_a", "row_drop", "row_impute_c"]),
    )
    state = PreprocessingState(
        phospho=phospho,
        site_metadata=pd.DataFrame(
            {
                "gene_symbol": ["MAPK14", "AKT1", "GSK3B", "PRKACA"],
                "site": ["Y182", "T308", "S9", "S339"],
                "site_sequence": ["SEQ_A", "SEQ_R", "SEQ_C", "SEQ_D"],
            },
            index=phospho.index.copy(),
        ),
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan(
            intensity_transform_policy="log2",
            missing_data_policy="impute_minprob",
            missing_data_q=0.01,
            missing_data_width=0.3,
            missing_data_seed=12345,
            missing_data_max_missing_fraction_per_row=0.5,
            stage_order=("intensity_transform", "missing_data"),
        ),
    )

    result = MissingDataStage().run(state)
    diagnostics = result.diagnostics["diagnostics"]

    assert result.state.phospho.isna().to_numpy().sum() == 0
    assert result.state.phospho.index.tolist() == [
        "row_keep",
        "row_impute_a",
        "row_impute_c",
    ]
    assert diagnostics["imputation_method_id"] == "minprob"
    assert diagnostics["imputation_method_family"] == "left_censored_random"
    assert diagnostics["left_censored_assumption"] is True
    assert diagnostics["matrix_scale_requirement"] == "log2"
    assert diagnostics["random_seed"] == 12345
    assert diagnostics["dropped_row_ids"] == ["row_drop"]
    assert diagnostics["dropped_row_count"] == 1
    assert diagnostics["imputed_row_count"] == 2
    assert diagnostics["imputed_column_count"] == 2
    assert diagnostics["dropped_rows_above_max_missing_fraction"] == ["row_drop"]
    assert diagnostics["output_missing_cell_count"] == 0
    assert isinstance(diagnostics["imputation_mask_hash"], str)
    assert (
        diagnostics["per_column_distribution_parameters"]["sample_a"]["observed_count"]
        == 2
    )
    assert (
        diagnostics["per_column_distribution_parameters"]["sample_a"]["missing_count"]
        == 1
    )


def test_missing_data_stage_minprob_is_deterministic_for_same_seed() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [10.0, float("nan"), float("nan"), 4.0],
            "sample_b": [9.0, 8.0, float("nan"), 6.0],
            "sample_c": [11.0, 7.0, 5.0, float("nan")],
        },
        index=pd.Index(["row_keep", "row_impute_a", "row_drop", "row_impute_c"]),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "AKT1", "GSK3B", "PRKACA"],
            "site": ["Y182", "T308", "S9", "S339"],
            "site_sequence": ["SEQ_A", "SEQ_R", "SEQ_C", "SEQ_D"],
            "localisation_confidence": [0.95, 0.9, 0.92, 0.94],
        },
        index=phospho.index.copy(),
    )
    plan = PreprocessingPlan(
        intensity_transform_policy="log2",
        missing_data_policy="impute_minprob",
        missing_data_q=0.01,
        missing_data_width=0.3,
        missing_data_seed=12345,
        missing_data_max_missing_fraction_per_row=0.5,
        stage_order=("intensity_transform", "missing_data"),
    )

    first = MissingDataStage().run(
        PreprocessingState(
            phospho=phospho.copy(deep=True),
            site_metadata=site_metadata.copy(deep=True),
            sample_metadata=None,
            total=None,
            plan=plan,
        )
    )
    second = MissingDataStage().run(
        PreprocessingState(
            phospho=phospho.copy(deep=True),
            site_metadata=site_metadata.copy(deep=True),
            sample_metadata=None,
            total=None,
            plan=plan,
        )
    )

    pd.testing.assert_frame_equal(first.state.phospho, second.state.phospho)


def test_missing_data_stage_minprob_changes_with_seed() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [10.0, float("nan"), float("nan"), 4.0],
            "sample_b": [9.0, 8.0, float("nan"), 6.0],
            "sample_c": [11.0, 7.0, 5.0, float("nan")],
        },
        index=pd.Index(["row_keep", "row_impute_a", "row_drop", "row_impute_c"]),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "AKT1", "GSK3B", "PRKACA"],
            "site": ["Y182", "T308", "S9", "S339"],
            "site_sequence": ["SEQ_A", "SEQ_R", "SEQ_C", "SEQ_D"],
            "localisation_confidence": [0.95, 0.9, 0.92, 0.94],
        },
        index=phospho.index.copy(),
    )

    first = MissingDataStage().run(
        PreprocessingState(
            phospho=phospho.copy(deep=True),
            site_metadata=site_metadata.copy(deep=True),
            sample_metadata=None,
            total=None,
            plan=PreprocessingPlan(
                intensity_transform_policy="log2",
                missing_data_policy="impute_minprob",
                missing_data_q=0.01,
                missing_data_width=0.3,
                missing_data_seed=1,
                missing_data_max_missing_fraction_per_row=0.5,
                stage_order=("intensity_transform", "missing_data"),
            ),
        )
    )
    second = MissingDataStage().run(
        PreprocessingState(
            phospho=phospho.copy(deep=True),
            site_metadata=site_metadata.copy(deep=True),
            sample_metadata=None,
            total=None,
            plan=PreprocessingPlan(
                intensity_transform_policy="log2",
                missing_data_policy="impute_minprob",
                missing_data_q=0.01,
                missing_data_width=0.3,
                missing_data_seed=2,
                missing_data_max_missing_fraction_per_row=0.5,
                stage_order=("intensity_transform", "missing_data"),
            ),
        )
    )

    assert float(first.state.phospho.loc["row_impute_a", "sample_a"]) != pytest.approx(
        float(second.state.phospho.loc["row_impute_a", "sample_a"])
    )


def test_missing_data_stage_minprob_rejects_incompatible_stage_order() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [10.0, float("nan"), 4.0],
            "sample_b": [9.0, 8.0, 6.0],
            "sample_c": [11.0, 7.0, float("nan")],
        },
        index=pd.Index(["row_keep", "row_impute_a", "row_impute_c"]),
    )
    state = PreprocessingState(
        phospho=phospho,
        site_metadata=pd.DataFrame(
            {
                "gene_symbol": ["MAPK14", "AKT1", "PRKACA"],
                "site": ["Y182", "T308", "S339"],
                "site_sequence": ["SEQ_A", "SEQ_R", "SEQ_D"],
            },
            index=phospho.index.copy(),
        ),
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan(
            intensity_transform_policy="log2",
            missing_data_policy="impute_minprob",
            missing_data_q=0.01,
            missing_data_width=0.3,
            missing_data_seed=12345,
            missing_data_max_missing_fraction_per_row=1.0,
            stage_order=("missing_data",),
        ),
    )

    with pytest.raises(
        PhosPyInputError,
        match="stage_order is incompatible with minprob intensity-state requirements",
    ):
        MissingDataStage().run(state)


def test_missing_data_stage_knn_imputes_drops_and_reports_diagnostics() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 1.0, 2.0, 10.0],
            "sample_b": [1.0, 2.0, 2.0, float("nan")],
            "sample_c": [float("nan"), 3.0, 3.0, float("nan")],
        },
        index=pd.Index(["row_impute", "row_ref_1", "row_ref_2", "row_drop"]),
    )
    state = PreprocessingState(
        phospho=phospho,
        site_metadata=pd.DataFrame(
            {
                "gene_symbol": ["MAPK14", "AKT1", "GSK3B", "PRKACA"],
                "site": ["Y182", "T308", "S9", "S339"],
                "site_sequence": ["SEQ_A", "SEQ_R", "SEQ_C", "SEQ_D"],
            },
            index=phospho.index.copy(),
        ),
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan(
            missing_data_policy="impute_knn",
            missing_data_k=1,
            missing_data_distance="nan_euclidean",
            missing_data_max_missing_fraction_per_row=0.5,
            stage_order=("missing_data",),
        ),
    )

    result = MissingDataStage().run(state)
    diagnostics = result.diagnostics["diagnostics"]

    assert result.state.phospho.columns.tolist() == ["sample_a", "sample_b", "sample_c"]
    assert result.state.phospho.index.tolist() == [
        "row_impute",
        "row_ref_1",
        "row_ref_2",
    ]
    assert float(result.state.phospho.loc["row_impute", "sample_c"]) == pytest.approx(
        3.0
    )
    assert int(result.state.phospho.isna().to_numpy().sum()) == 0
    assert diagnostics["imputation_method_id"] == "knn"
    assert diagnostics["imputation_method_family"] == "nearest_neighbour"
    assert diagnostics["neighbour_count"] == 1
    assert diagnostics["distance_metric"] == "nan_euclidean"
    assert diagnostics["rows_not_imputable"] == ["row_drop"]
    assert diagnostics["dropped_row_ids"] == ["row_drop"]
    assert diagnostics["dropped_row_count"] == 1
    assert diagnostics["imputed_row_count"] == 1
    assert diagnostics["imputed_column_count"] == 1
    assert diagnostics["output_missing_cell_count"] == 0
    assert diagnostics["matrix_scale_requirement"] is None
    assert diagnostics["left_censored_assumption"] is False
    assert isinstance(diagnostics["imputation_mask_hash"], str)


def test_missing_data_stage_knn_rejects_columns_without_observed_values() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 2.0],
            "sample_b": [2.0, 3.0],
            "sample_c": [float("nan"), float("nan")],
        },
        index=pd.Index(["row_a", "row_b"]),
    )
    state = PreprocessingState(
        phospho=phospho,
        site_metadata=pd.DataFrame(
            {
                "gene_symbol": ["MAPK14", "AKT1"],
                "site": ["Y182", "T308"],
                "site_sequence": ["SEQ_A", "SEQ_R"],
            },
            index=phospho.index.copy(),
        ),
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan(
            missing_data_policy="impute_knn",
            missing_data_k=1,
            missing_data_distance="nan_euclidean",
            missing_data_max_missing_fraction_per_row=1.0,
            stage_order=("missing_data",),
        ),
    )

    with pytest.raises(
        PhosPyInputError,
        match="have no observed values after row filtering",
    ):
        MissingDataStage().run(state)


def test_processing_state_missing_data_imputed_flag_requires_provenance() -> None:
    state = build_dataset_processing_state(
        plan=PreprocessingPlan(
            missing_data_policy="impute_row_median",
            missing_data_min_observed_values=1,
            stage_order=("missing_data",),
        ),
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        preprocessing_trace=None,
    )

    assert state.missing_data.policy == "impute_row_median"
    assert state.missing_data.imputed is False
