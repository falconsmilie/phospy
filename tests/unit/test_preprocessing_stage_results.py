from __future__ import annotations

from collections.abc import Mapping

import pandas as pd

from phospy.datasets.preprocessing.models import (
    PreprocessingPlan,
    PreprocessingReportRow,
    PreprocessingStageResult,
    PreprocessingState,
)
from phospy.datasets.preprocessing.pipeline import PreprocessingPipeline
from phospy.datasets.preprocessing.stages.comparisons import ComparisonsStage
from phospy.datasets.preprocessing.stages.intensity_transform import (
    IntensityTransformStage,
)
from phospy.datasets.preprocessing.stages.missing_data import MissingDataStage
from phospy.datasets.preprocessing.stages.normalisation import NormalisationStage
from phospy.datasets.preprocessing.stages.site_matrix import SiteMatrixStage
from phospy.datasets.preprocessing.stages.total_protein_correction import (
    TotalProteinCorrectionStage,
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
            "site_sequence": ["SEQ_A", "SEQ_B"],
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
                        table="custom_table",
                        values={"custom_metric": 123},
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
        stage_registry=(FakeStage(),)
    ).run_with_trace(state)

    assert len(trace) == 1
    assert trace[0].stage == "fake_stage"
    assert trace[0].diagnostics["custom_metric"] == 123
    assert len(final_state.report_rows) == 1
    assert final_state.report_rows[0].table == "custom_table"


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
    assert diagnostics["policy"] == "median_center"
    assert diagnostics["note"] == "median centering used"
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
