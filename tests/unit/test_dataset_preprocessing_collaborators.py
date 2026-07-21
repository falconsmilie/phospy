from __future__ import annotations

import pandas as pd
import pandas.testing as pdt
import pytest

from phospy.api.configs import DatasetMissingDataConfig, DatasetPreprocessingConfig
from phospy.errors.build import DatasetBuildError
from phospy.errors.input import PhosPyInputError
from phospy.science.datasets.builders.preprocessing import (
    DatasetPreprocessor,
    build_dataset_processing_state,
)
from phospy.science.datasets.builders.transformation_resolver import (
    DatasetIntensityScaleResolver,
)
from phospy.science.datasets.preprocessing.diagnostics import ProcessingTraceDiagnostics
from phospy.science.datasets.preprocessing.models import (
    DATASET_PREPROCESSING_STAGE_MISSING_DATA,
    DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION,
    DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION,
    PreprocessingPlan,
    PreprocessingStageExecution,
    PreprocessingStageResult,
    PreprocessingState,
)
from phospy.science.datasets.preprocessing.pipeline import (
    PreprocessingPipeline,
    _resolve_state_table,
)
from phospy.science.datasets.preprocessing.provenance_adapter import (
    PreprocessingProvenanceAdapter,
)
from phospy.science.datasets.preprocessing.state_builder import (
    DatasetProcessingStateBuilder,
)
from phospy.science.transformations.models import (
    IntensityScaleKind,
    IntensityScaleState,
    MatrixIntensityScaleState,
    QuantitativeMeaning,
)
from phospy.science.transformations.transformers import IdentityTransformer
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
            "site_sequence": ["SEQ_A", "SEQ_R", "SEQ_C"],
            "localisation_confidence": [0.95, 0.9, 0.92],
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


def _supported_log2_intensity_scale_state(*, has_total_matrix: bool):
    phospho = pd.DataFrame({"sample_a": [1.0]}, index=["GENEA;S1;"])
    total = (
        None
        if not has_total_matrix
        else pd.DataFrame({"sample_a": [1.0]}, index=["GENEA"])
    )
    declared_log2 = IntensityScaleState(
        phospho=MatrixIntensityScaleState.log2(established_by="trusted.fixture"),
        total=(
            MatrixIntensityScaleState.log2(established_by="trusted.fixture")
            if has_total_matrix
            else None
        ),
    )
    return (
        DatasetIntensityScaleResolver(transformer=IdentityTransformer())
        .run(
            phospho=phospho,
            total=total,
            expected_scale_kind=IntensityScaleKind.LOG2,
            declared_input_scale_state=declared_log2,
        )
        .intensity_scale_state
    )


def test_diagnostics_collaborator_parses_stage_payloads_and_row_diagnostics() -> None:
    trace = (
        _trace_record(
            stage=DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION,
            diagnostics={"requires_log_scale": True, "quantitative_meaning": "unknown"},
        ),
        _trace_record(
            stage=DATASET_PREPROCESSING_STAGE_MISSING_DATA,
            diagnostics={"imputed_cell_count": 2, "output_missing_cell_count": 0},
        ),
        _trace_record(
            stage=DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION,
            diagnostics={
                "configured": True,
                "unresolved_counts_by_reason": {"missing_accession": 2},
                "row_diagnostics": [
                    {
                        "row_index": 0,
                        "row_id": "MAPK14;Y182;",
                        "status": "resolved",
                        "action": "filled",
                    },
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
            stage=DATASET_PREPROCESSING_STAGE_MISSING_DATA,
            key="imputed_cell_count",
            default=0,
        )
        == 2
    )
    assert ProcessingTraceDiagnostics.resolve_optional_mapping_int(
        parsed.site_sequence_resolution,
        stage=DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION,
        key="unresolved_counts_by_reason",
    ) == {"missing_accession": 2}
    row_diagnostics = ProcessingTraceDiagnostics.resolve_site_sequence_row_diagnostics(
        parsed.site_sequence_resolution
    )
    assert len(row_diagnostics) == 1
    assert row_diagnostics[0].row_id == "MAPK14;Y182;"
    assert row_diagnostics[0].action == "filled"


@pytest.mark.parametrize("invalid_value", [True, False, 1.5, "1"])
def test_integer_count_field_rejects_non_int_values(invalid_value: object) -> None:
    with pytest.raises(
        DatasetBuildError,
        match=(
            f"stage={DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION!r}.*"
            "resolved_site_count"
        ),
    ):
        ProcessingTraceDiagnostics.resolve_optional_int(
            {"resolved_site_count": invalid_value},
            stage=DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION,
            key="resolved_site_count",
            default=0,
        )


def test_site_sequence_row_diagnostics_reject_missing_required_row_id() -> None:
    with pytest.raises(
        DatasetBuildError,
        match=(
            f"stage={DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION!r}.*"
            "row_diagnostics\\[0\\]\\.row_id"
        ),
    ):
        ProcessingTraceDiagnostics.resolve_site_sequence_row_diagnostics(
            {"row_diagnostics": [{"row_index": 0}]}
        )


def test_site_sequence_diagnostics_reject_unknown_fields() -> None:
    with pytest.raises(
        DatasetBuildError,
        match=(
            f"stage {DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION!r}: "
            "unexpected_count"
        ),
    ):
        ProcessingTraceDiagnostics.validate_site_sequence_resolution_payload(
            {"configured": True, "unexpected_count": 1}
        )


def test_integer_count_field_accepts_valid_int() -> None:
    resolved = ProcessingTraceDiagnostics.resolve_optional_int(
        {"resolved_site_count": 4},
        stage=DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION,
        key="resolved_site_count",
        default=0,
    )
    assert resolved == 4


def test_optional_field_defaults_remain_for_missing_site_sequence_count() -> None:
    resolved = ProcessingTraceDiagnostics.resolve_optional_int(
        {},
        stage=DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION,
        key="resolved_site_count",
        default=0,
    )
    assert resolved == 0


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
            diagnostics={
                "diagnostics_schema_version": 1,
                "quantitative_meaning": "not_a_supported_meaning",
            },
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


def test_processing_state_builder_emits_default_total_correction_diagnostics_when_stage_absent() -> (
    None
):
    plan = PreprocessingPlan.default()
    builder = DatasetProcessingStateBuilder()

    state = builder.build(
        plan=plan,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        preprocessing_trace=None,
    )

    diagnostics = state.total_protein_correction.diagnostics
    assert diagnostics is not None
    assert diagnostics.get("diagnostics_schema_version") == 1
    assert diagnostics.get("policy") == "none"
    assert diagnostics.get("requested_policy") == "none"
    assert diagnostics.get("resolved_policy") == "none"
    assert diagnostics.get("quantitative_meaning") == "phosphosite_abundance"


@pytest.mark.parametrize(
    ("plan", "state", "explicit_quantitative_meaning", "expected"),
    [
        pytest.param(
            PreprocessingPlan.default(),
            supported_linear_intensity_scale_state(has_total_matrix=False),
            None,
            QuantitativeMeaning.PHOSPHOSITE_ABUNDANCE.value,
            id="raw-abundance-default",
        ),
        pytest.param(
            PreprocessingPlan.default(),
            _supported_log2_intensity_scale_state(has_total_matrix=False),
            None,
            QuantitativeMeaning.PHOSPHOSITE_LOG_ABUNDANCE.value,
            id="log-abundance-default",
        ),
        pytest.param(
            PreprocessingPlan.default(),
            _supported_log2_intensity_scale_state(has_total_matrix=False),
            QuantitativeMeaning.PHOSPHOSITE_LOG_ABUNDANCE,
            QuantitativeMeaning.PHOSPHOSITE_LOG_ABUNDANCE.value,
            id="log-abundance-explicit",
        ),
    ],
)
def test_processing_state_builder_quantitative_meaning_matrix(
    plan: PreprocessingPlan,
    state,
    explicit_quantitative_meaning: QuantitativeMeaning | None,
    expected: str,
) -> None:
    built = DatasetProcessingStateBuilder().build(
        plan=plan,
        intensity_scale_state=state,
        explicit_quantitative_meaning=explicit_quantitative_meaning,
        preprocessing_trace=None,
    )
    assert built.intensity_scale.quantity is not None
    assert built.intensity_scale.quantity.value == expected
    assert built.total_protein_correction.quantitative_meaning == expected
    assert built.total_protein_correction.diagnostics is not None
    assert built.total_protein_correction.diagnostics.get("quantitative_meaning") == (
        expected
    )


def test_processing_state_builder_rejects_total_correction_meaning_without_trace() -> (
    None
):
    with pytest.raises(DatasetBuildError, match="requires a total_protein_correction"):
        DatasetProcessingStateBuilder().build(
            plan=PreprocessingPlan(
                intensity_transform_policy="log2",
                total_protein_correction_policy="subtract_log_total",
            ),
            intensity_scale_state=_supported_log2_intensity_scale_state(
                has_total_matrix=True
            ),
            preprocessing_trace=None,
        )


def test_processing_state_builder_rejects_operation_derived_explicit_meaning() -> None:
    with pytest.raises(DatasetBuildError, match="may only declare direct input"):
        DatasetProcessingStateBuilder().build(
            plan=PreprocessingPlan.default(),
            intensity_scale_state=_supported_log2_intensity_scale_state(
                has_total_matrix=False
            ),
            explicit_quantitative_meaning=(
                QuantitativeMeaning.CONTRAST_LOG2_FOLD_CHANGE
            ),
            preprocessing_trace=None,
        )


def test_pipeline_rejects_malformed_stage_diagnostics_before_trace_is_recorded() -> (
    None
):
    class MalformedMissingDataStage:
        stage_key = DATASET_PREPROCESSING_STAGE_MISSING_DATA

        def run(self, state: PreprocessingState) -> PreprocessingStageResult:
            return PreprocessingStageResult(
                state=state,
                diagnostics={
                    "dropped_row_ids": (),
                    "dropped_row_count": 0,
                    "imputed_cell_count": "1",
                    "imputed_row_ids": (),
                    "notes": "stage executed",
                    "diagnostics": {"random_seed": 7},
                },
            )

    pipeline = PreprocessingPipeline(stage_registry=(MalformedMissingDataStage(),))
    state = PreprocessingState(
        phospho=_phospho(),
        site_metadata=_site_metadata(_phospho().index),
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan(
            localisation_mode="ignore",
            stage_order=(DATASET_PREPROCESSING_STAGE_MISSING_DATA,),
        ),
    )

    with pytest.raises(
        DatasetBuildError,
        match=f"stage={DATASET_PREPROCESSING_STAGE_MISSING_DATA!r}.*imputed_cell_count",
    ):
        pipeline.run_with_trace(state)


def test_pipeline_table_resolution_rejects_unknown_table_key() -> None:
    state = PreprocessingState(
        phospho=_phospho(),
        site_metadata=_site_metadata(_phospho().index),
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan.default(),
    )

    with pytest.raises(DatasetBuildError, match="unknown table key"):
        _resolve_state_table(
            state=state,
            table_name="dataset.sampl_metadata",
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
    missing_data_rows = row_counts.loc[row_counts.loc[:, "stage"] == "missing_data", :]
    assert missing_data_rows.shape[0] == 1
    assert int(missing_data_rows.iloc[0]["dropped_rows"]) == 1
    assert "stage not scheduled in preprocessing plan" in set(
        operations.loc[:, "notes"].astype(str).tolist()
    )


def test_provenance_adapter_includes_execution_summary_for_imputation_stage() -> None:
    plan = PreprocessingPlan.from_config(
        DatasetPreprocessingConfig(
            missing_data=DatasetMissingDataConfig(
                policy="impute_row_median",
                min_observed_values=1,
            )
        )
    )
    adapter = PreprocessingProvenanceAdapter()
    trace = (
        _trace_record(
            stage=DATASET_PREPROCESSING_STAGE_MISSING_DATA,
            diagnostics={"imputed_cell_count": 2, "output_missing_cell_count": 0},
            notes="stage executed",
            dropped_row_count=1,
        ),
    )

    _, operations = adapter.build_tables(
        plan=plan,
        input_row_count=3,
        output_row_count=2,
        trace=trace,
    )
    missing_data_row = operations.loc[
        operations.loc[:, "stage"] == "missing_data"
    ].iloc[0]
    parameters = missing_data_row["parameters"]
    assert isinstance(parameters, dict)
    summary = parameters.get("execution_summary")
    assert isinstance(summary, dict)
    assert summary["dropped_rows"] == 1
    assert summary["imputed_cell_count"] == 2
    assert summary["imputation_scope"] == "per_row"
    assert summary["determinism_kind"] == "deterministic"
    assert summary["diagnostic_summary"]["output_missing_cell_count"] == 0


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


def test_disabled_site_sequence_stage_does_not_emit_stage_diagnostics() -> None:
    phospho = _phospho()
    state = PreprocessingState(
        phospho=phospho,
        site_metadata=_site_metadata(phospho.index),
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan.default(),
    )

    _, trace = PreprocessingPipeline().run_with_trace(state)
    parsed = ProcessingTraceDiagnostics.from_trace(trace)

    assert parsed.site_sequence_resolution is None

    processing_state = DatasetProcessingStateBuilder().build(
        plan=state.plan,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        preprocessing_trace=trace,
        final_phospho=state.phospho,
        final_site_metadata=state.site_metadata,
        final_sample_metadata=state.sample_metadata,
    )
    assert processing_state.site_sequence_resolution.configured is False
    assert processing_state.site_sequence_resolution.mode is None
