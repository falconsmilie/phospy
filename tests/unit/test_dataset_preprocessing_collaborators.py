from __future__ import annotations

import pandas as pd
import pandas.testing as pdt
import pytest

from phospy.api.configs import DatasetPreprocessingConfig
from phospy.datasets.builders.preprocessing import (
    DatasetPreprocessor,
    build_dataset_processing_state,
)
from phospy.datasets.preprocessing.diagnostics import ProcessingTraceDiagnostics
from phospy.datasets.preprocessing.models import (
    DATASET_PREPROCESSING_STAGE_MISSING_DATA,
    DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION,
    DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION,
    PreprocessingPlan,
    PreprocessingStageExecution,
    PreprocessingState,
)
from phospy.datasets.preprocessing.pipeline import PreprocessingPipeline
from phospy.datasets.preprocessing.provenance_adapter import (
    PreprocessingProvenanceAdapter,
)
from phospy.datasets.preprocessing.state_builder import DatasetProcessingStateBuilder
from phospy.errors.build import DatasetBuildError
from phospy.errors.input import PhosPyInputError
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
)


def _phospho() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample_a": [1.0, 2.0, 3.0],
            "sample_b": [2.0, 1.0, 4.0],
        },
        index=pd.Index(["MAPK14;Y182;", "GSK3B;S9;", "AKT1;T308;"], name="site_id"),
    )


def _site_metadata(index: pd.Index) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "GSK3B", "AKT1"],
            "site": ["Y182", "S9", "T308"],
            "site_sequence": ["SEQ_A", "SEQ_B", "SEQ_C"],
        },
        index=index.copy(),
    )


def _trace_record(
    *,
    stage: str,
    diagnostics: dict[str, object] | None = None,
    notes: str | None = None,
    dropped_row_count: int = 0,
) -> PreprocessingStageExecution:
    return PreprocessingStageExecution(
        stage=stage,
        operation=f"{stage}.operation",
        parameters={"stage": stage},
        input_shape=(3, 2),
        output_shape=(max(3 - dropped_row_count, 0), 2),
        input_hash="input_hash",
        output_hash="output_hash",
        dropped_row_count=dropped_row_count,
        notes=notes,
        diagnostics={} if diagnostics is None else dict(diagnostics),
    )


def test_diagnostics_collaborator_parses_stage_payloads_and_row_diagnostics() -> None:
    trace = (
        _trace_record(
            stage=DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION,
            diagnostics={"requires_log_scale": True, "quantitative_meaning": "unknown"},
        ),
        _trace_record(
            stage=DATASET_PREPROCESSING_STAGE_MISSING_DATA,
            diagnostics={"imputed_cell_count": "2", "output_missing_cell_count": 0},
        ),
        _trace_record(
            stage=DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION,
            diagnostics={
                "configured": True,
                "unresolved_counts_by_reason": {"missing_accession": "2"},
                "row_diagnostics": [
                    {"row_index": 0, "row_id": "MAPK14;Y182;", "action": "filled"},
                    {"row_index": -1, "row_id": "bad"},
                    {"row_index": 2},
                ],
            },
        ),
    )

    parsed = ProcessingTraceDiagnostics.from_trace(trace)

    assert parsed.total_protein_correction is not None
    assert parsed.missing_data is not None
    assert parsed.site_sequence_resolution is not None
    assert (
        ProcessingTraceDiagnostics.resolve_optional_int(
            parsed.missing_data,
            key="imputed_cell_count",
            default=0,
        )
        == 2
    )
    assert ProcessingTraceDiagnostics.resolve_optional_mapping_int(
        parsed.site_sequence_resolution,
        key="unresolved_counts_by_reason",
    ) == {"missing_accession": 2}
    row_diagnostics = ProcessingTraceDiagnostics.resolve_site_sequence_row_diagnostics(
        parsed.site_sequence_resolution
    )
    assert len(row_diagnostics) == 1
    assert row_diagnostics[0].row_id == "MAPK14;Y182;"
    assert row_diagnostics[0].action == "filled"


def test_processing_state_builder_surfaces_invalid_diagnostics_clearly() -> None:
    plan = PreprocessingPlan.default()
    builder = DatasetProcessingStateBuilder()
    trace = (
        _trace_record(
            stage=DATASET_PREPROCESSING_STAGE_MISSING_DATA,
            diagnostics={"diagnostics_schema_version": "invalid"},
        ),
    )

    with pytest.raises(PhosPyInputError, match="diagnostics_schema_version"):
        builder.build(
            plan=plan,
            intensity_scale_state=supported_linear_intensity_scale_state(
                has_total_matrix=False
            ),
            preprocessing_trace=trace,
        )


def test_processing_state_builder_rejects_unknown_quantitative_meaning_without_defaulting() -> (
    None
):
    plan = PreprocessingPlan.default()
    builder = DatasetProcessingStateBuilder()
    trace = (
        _trace_record(
            stage=DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION,
            diagnostics={"quantitative_meaning": "not_a_supported_meaning"},
        ),
    )

    with pytest.raises(DatasetBuildError, match="quantitative_meaning must be one of:"):
        builder.build(
            plan=plan,
            intensity_scale_state=supported_linear_intensity_scale_state(
                has_total_matrix=False
            ),
            preprocessing_trace=trace,
        )


def test_provenance_adapter_builds_row_counts_and_operations() -> None:
    plan = PreprocessingPlan.from_config(DatasetPreprocessingConfig())
    adapter = PreprocessingProvenanceAdapter()
    trace = (
        _trace_record(
            stage=DATASET_PREPROCESSING_STAGE_MISSING_DATA,
            diagnostics={"imputed_cell_count": 0},
            notes="stage executed",
            dropped_row_count=1,
        ),
    )

    row_counts, operations = adapter.build_tables(
        plan=plan,
        input_row_count=3,
        output_row_count=2,
        trace=trace,
    )

    assert row_counts.iloc[0]["stage"] == "preprocessing_input"
    assert row_counts.iloc[-1]["stage"] == "preprocessing_complete"
    assert int(row_counts.iloc[1]["dropped_rows"]) == 1
    assert "stage not scheduled in preprocessing plan" in set(
        operations.loc[:, "notes"].astype(str).tolist()
    )


def test_builder_integration_assembles_pipeline_outputs_through_collaborators() -> None:
    phospho = _phospho()
    site_metadata = _site_metadata(phospho.index)
    plan = PreprocessingPlan.from_config(DatasetPreprocessingConfig())
    pipeline = PreprocessingPipeline()

    expected_state, expected_trace = pipeline.run_with_trace(
        PreprocessingState(
            phospho=phospho,
            site_metadata=site_metadata,
            sample_metadata=None,
            total=None,
            plan=plan,
        )
    )
    expected_row_counts, expected_operations = (
        PreprocessingProvenanceAdapter().build_tables(
            plan=plan,
            input_row_count=3,
            output_row_count=int(len(expected_state.phospho.index)),
            trace=expected_trace,
        )
    )
    preprocessed = DatasetPreprocessor(pipeline=pipeline).run(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=None,
        total=None,
        plan=plan,
    )
    direct_state = DatasetProcessingStateBuilder().build(
        plan=plan,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        preprocessing_trace=preprocessed.preprocessing_trace,
        final_phospho=preprocessed.phospho,
        final_site_metadata=preprocessed.site_metadata,
        final_sample_metadata=preprocessed.sample_metadata,
    )
    wrapper_state = build_dataset_processing_state(
        plan=plan,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        preprocessing_trace=preprocessed.preprocessing_trace,
        final_phospho=preprocessed.phospho,
        final_site_metadata=preprocessed.site_metadata,
        final_sample_metadata=preprocessed.sample_metadata,
    )

    pdt.assert_frame_equal(preprocessed.phospho, expected_state.phospho)
    pdt.assert_frame_equal(preprocessed.preprocessing_row_counts, expected_row_counts)
    pdt.assert_frame_equal(
        preprocessed.preprocessing_operations,
        expected_operations,
    )
    assert preprocessed.preprocessing_trace == expected_trace
    assert wrapper_state == direct_state
