from __future__ import annotations

import pandas as pd
import pytest

from phospy.api.configs import (
    DATASET_LOCALISATION_MODE_ALLOW_MISSING_WITH_WAIVER,
    DatasetComparisonBuildingConfig,
    DatasetIntensityTransformConfig,
    DatasetLocalisationConfig,
    DatasetNormalisationConfig,
    DatasetPreprocessingConfig,
    DatasetSiteMatrixConfig,
)
from phospy.errors.build import DatasetBuildError
from phospy.science.datasets.builders.preprocessing import DatasetPreprocessor
from phospy.science.datasets.preprocessing.models import (
    DATASET_PREPROCESSING_STAGE_BATCH_CORRECTION,
    DATASET_PREPROCESSING_STAGE_COMPARISONS,
    DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM,
    DATASET_PREPROCESSING_STAGE_LOCALISATION,
    DATASET_PREPROCESSING_STAGE_MISSING_DATA,
    DATASET_PREPROCESSING_STAGE_NORMALISATION,
    DATASET_PREPROCESSING_STAGE_SITE_MATRIX,
    DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION,
    DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION,
    PreprocessingPlan,
    PreprocessingStageResult,
    PreprocessingState,
)
from phospy.science.datasets.preprocessing.pipeline import PreprocessingPipeline
from phospy.science.datasets.preprocessing.stage_registry import (
    PreprocessingStageMetadata,
    get_preprocessing_stage_metadata,
    list_registered_preprocessing_stages,
    resolve_builder_provenance_stage_order,
    resolve_registered_preprocessing_stages,
)
from tests.support.site_keys import site_key_index_from_display_ids

_DISPLAY_IDS = ["MAPK14;Y182;", "AKT1;T308;"]


def _phospho() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample_a": [1.0, 2.0],
            "sample_b": [2.0, 4.0],
        },
        index=pd.Index(_DISPLAY_IDS, name="site_id"),
    )


def _site_metadata(index: pd.Index) -> pd.DataFrame:
    site_keys = site_key_index_from_display_ids(
        _DISPLAY_IDS,
        protein_namespace="gene_symbol",
    )
    return pd.DataFrame(
        {
            "site_key": site_keys.astype(str).tolist(),
            "display_id": _DISPLAY_IDS,
            "gene_symbol": ["MAPK14", "AKT1"],
            "site": ["Y182", "T308"],
            "site_sequence": ["SEQ_A", "SEQ_R"],
            "localisation_confidence": [0.95, 0.91],
        },
        index=index.copy(),
    )


def _sample_metadata(columns: pd.Index) -> pd.DataFrame:
    return pd.DataFrame(
        {"comparison_group": ["group_a", "group_b"]},
        index=columns.copy(),
    )


def _plan_with_multiple_stages() -> PreprocessingPlan:
    return PreprocessingPlan.from_config(
        DatasetPreprocessingConfig(
            intensity_transform=DatasetIntensityTransformConfig(
                policy="log2",
                pseudocount=1.0,
            ),
            normalisation=DatasetNormalisationConfig(policy="median_center"),
            site_matrix=DatasetSiteMatrixConfig(policy="build_from_metadata"),
            comparisons=DatasetComparisonBuildingConfig(policy="sample_metadata_pairs"),
            localisation=DatasetLocalisationConfig(
                mode=DATASET_LOCALISATION_MODE_ALLOW_MISSING_WITH_WAIVER,
                waiver_reason="test waiver",
            ),
        )
    )


def test_every_preprocessing_stage_has_registry_metadata() -> None:
    expected = {
        DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION,
        DATASET_PREPROCESSING_STAGE_LOCALISATION,
        DATASET_PREPROCESSING_STAGE_MISSING_DATA,
        DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM,
        DATASET_PREPROCESSING_STAGE_BATCH_CORRECTION,
        DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION,
        DATASET_PREPROCESSING_STAGE_SITE_MATRIX,
        DATASET_PREPROCESSING_STAGE_NORMALISATION,
        DATASET_PREPROCESSING_STAGE_COMPARISONS,
    }
    observed = {item.stage_key for item in list_registered_preprocessing_stages()}
    assert observed == expected


def test_registered_stage_keys_are_unique() -> None:
    registered = list_registered_preprocessing_stages()
    stage_keys = [metadata.stage_key for metadata in registered]
    assert len(stage_keys) == len(set(stage_keys))


def test_registered_stage_order_is_stable() -> None:
    observed = tuple(
        metadata.stage_key for metadata in list_registered_preprocessing_stages()
    )
    assert observed == (
        DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION,
        DATASET_PREPROCESSING_STAGE_LOCALISATION,
        DATASET_PREPROCESSING_STAGE_MISSING_DATA,
        DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM,
        DATASET_PREPROCESSING_STAGE_BATCH_CORRECTION,
        DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION,
        DATASET_PREPROCESSING_STAGE_SITE_MATRIX,
        DATASET_PREPROCESSING_STAGE_NORMALISATION,
        DATASET_PREPROCESSING_STAGE_COMPARISONS,
    )


def test_every_registered_stage_has_required_metadata_contract_fields() -> None:
    plan = _plan_with_multiple_stages()
    for metadata in list_registered_preprocessing_stages():
        assert metadata.stage_key
        assert metadata.display_label
        assert metadata.provenance_stage_key
        assert callable(metadata.operation_name)
        assert callable(metadata.serialize_parameters)
        assert isinstance(metadata.serialize_parameters(plan), dict)
        assert isinstance(metadata.diagnostics_metadata, dict)
        assert metadata.diagnostics_metadata


def test_duplicate_override_stage_keys_fail_registry_resolution() -> None:
    duplicate_key = DATASET_PREPROCESSING_STAGE_MISSING_DATA
    duplicate_entries = (
        PreprocessingStageMetadata(
            stage_key=duplicate_key,
            display_label="missing_data_override_a",
            operation_name=lambda _plan: "forbid",
            serialize_parameters=lambda _plan: {},
            consumed_input_tables=("dataset.phospho",),
            produced_output_tables=("dataset.phospho",),
            diagnostics_metadata={"known_diagnostics_fields": ("policy",)},
        ),
        PreprocessingStageMetadata(
            stage_key=duplicate_key,
            display_label="missing_data_override_b",
            operation_name=lambda _plan: "forbid",
            serialize_parameters=lambda _plan: {},
            consumed_input_tables=("dataset.phospho",),
            produced_output_tables=("dataset.phospho",),
            diagnostics_metadata={"known_diagnostics_fields": ("policy",)},
        ),
    )

    with pytest.raises(DatasetBuildError, match="duplicate stage key"):
        resolve_registered_preprocessing_stages(duplicate_entries)


def test_stage_metadata_rejects_unknown_consumed_table_key() -> None:
    with pytest.raises(DatasetBuildError, match="unknown table key"):
        PreprocessingStageMetadata(
            stage_key="custom_stage",
            display_label="custom_stage",
            operation_name=lambda _plan: "custom_stage",
            serialize_parameters=lambda _plan: {},
            consumed_input_tables=("dataset.sampl_metadata",),
            produced_output_tables=("dataset.phospho",),
            diagnostics_metadata={"known_diagnostics_fields": ("policy",)},
        )


def test_stage_metadata_rejects_unknown_produced_table_key() -> None:
    with pytest.raises(DatasetBuildError, match="unknown table key"):
        PreprocessingStageMetadata(
            stage_key="custom_stage",
            display_label="custom_stage",
            operation_name=lambda _plan: "custom_stage",
            serialize_parameters=lambda _plan: {},
            consumed_input_tables=("dataset.phospho",),
            produced_output_tables=("report.row_audt",),
            diagnostics_metadata={"known_diagnostics_fields": ("policy",)},
        )


def test_known_optional_missing_table_is_skipped_in_trace_fingerprints() -> None:
    class _OptionalSampleMetadataStage:
        stage_key = "optional_sample_metadata_stage"

        def run(self, state: PreprocessingState) -> PreprocessingStageResult:
            return PreprocessingStageResult(state=state)

    phospho = _phospho()
    state = PreprocessingState(
        phospho=phospho,
        site_metadata=_site_metadata(phospho.index),
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan(stage_order=("optional_sample_metadata_stage",)),
    )
    metadata = PreprocessingStageMetadata(
        stage_key="optional_sample_metadata_stage",
        display_label="optional_sample_metadata_stage",
        operation_name=lambda _plan: "optional_sample_metadata_stage",
        serialize_parameters=lambda _plan: {},
        consumed_input_tables=("dataset.sample_metadata",),
        produced_output_tables=("dataset.phospho",),
        diagnostics_metadata={"known_diagnostics_fields": ("policy",)},
    )

    _, trace = PreprocessingPipeline(
        stage_registry=(_OptionalSampleMetadataStage(),),
        stage_metadata_registry=(metadata,),
    ).run_with_trace(state)

    assert len(trace) == 1
    assert trace[0].consumed_input_tables == ()
    assert tuple(item.name for item in trace[0].produced_output_tables) == (
        "dataset.phospho",
    )


def test_pipeline_trace_parameters_and_operation_come_from_registry() -> None:
    phospho = _phospho()
    plan = _plan_with_multiple_stages()
    preprocessed = DatasetPreprocessor().run(
        phospho=phospho,
        site_metadata=_site_metadata(phospho.index),
        sample_metadata=_sample_metadata(phospho.columns),
        total=None,
        plan=plan,
    )
    trace = preprocessed.preprocessing_trace or ()
    assert trace

    for entry in trace:
        metadata = get_preprocessing_stage_metadata(entry.stage)
        assert entry.operation == metadata.operation_name(plan)
        assert dict(entry.parameters) == metadata.serialize_parameters(plan)


def test_pipeline_and_builder_use_same_stage_labels_and_operations() -> None:
    phospho = _phospho()
    plan = _plan_with_multiple_stages()
    preprocessed = DatasetPreprocessor().run(
        phospho=phospho,
        site_metadata=_site_metadata(phospho.index),
        sample_metadata=_sample_metadata(phospho.columns),
        total=None,
        plan=plan,
    )
    trace = preprocessed.preprocessing_trace or ()
    operations = preprocessed.preprocessing_operations
    assert operations is not None

    for entry in trace:
        metadata = get_preprocessing_stage_metadata(entry.stage)
        matching = operations.loc[operations.loc[:, "stage"] == metadata.display_label]
        assert matching.shape[0] == 1
        assert str(matching.iloc[0]["operation"]) == entry.operation


def test_provenance_operations_are_resolved_from_registry_metadata() -> None:
    phospho = _phospho()
    plan = _plan_with_multiple_stages()
    preprocessed = DatasetPreprocessor().run(
        phospho=phospho,
        site_metadata=_site_metadata(phospho.index),
        sample_metadata=_sample_metadata(phospho.columns),
        total=None,
        plan=plan,
    )
    operations = preprocessed.preprocessing_operations
    assert operations is not None

    expected_metadata = {
        metadata.display_label: metadata
        for metadata in resolve_builder_provenance_stage_order(plan)
    }
    for record in operations.to_dict(orient="records"):
        label = str(record["stage"])
        metadata = expected_metadata.get(label)
        assert metadata is not None
        assert str(record["operation"]) == metadata.operation_name(plan)


def test_unknown_stage_metadata_fails_with_clear_error() -> None:
    class _UnknownStage:
        stage_key = "unknown_stage"

        def run(self, state: PreprocessingState) -> PreprocessingStageResult:
            return PreprocessingStageResult(state=state)

    phospho = _phospho()
    state = PreprocessingState(
        phospho=phospho,
        site_metadata=_site_metadata(phospho.index),
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan(stage_order=("unknown_stage",)),
    )

    with pytest.raises(
        DatasetBuildError,
        match="metadata is not registered for stage 'unknown_stage'",
    ):
        PreprocessingPipeline(stage_registry=(_UnknownStage(),)).run_with_trace(state)


def test_registered_stage_factories_expose_run_method() -> None:
    for metadata in list_registered_preprocessing_stages():
        if metadata.stage_factory is None:
            continue
        stage = metadata.stage_factory()
        run_method = getattr(stage, "run", None)
        assert callable(run_method)


def test_custom_stage_registration_is_stage_owned() -> None:
    class FakeStage:
        stage_key = "fake_stage"

        def run(self, state: PreprocessingState) -> PreprocessingStageResult:
            return PreprocessingStageResult(
                state=state,
                diagnostics={
                    "notes": "stage executed",
                    "diagnostics": {"policy": "fake"},
                },
            )

    fake_contract = PreprocessingStageMetadata(
        stage_key="fake_stage",
        display_label="fake_stage",
        provenance_stage="fake_stage",
        operation_name=lambda _plan: "fake",
        serialize_parameters=lambda _plan: {"mode": "test"},
        consumed_input_tables=("dataset.phospho",),
        produced_output_tables=("dataset.phospho",),
        stage_factory=FakeStage,
        diagnostics_metadata={"known_diagnostics_fields": ("policy",)},
    )
    state = PreprocessingState(
        phospho=_phospho(),
        site_metadata=_site_metadata(_phospho().index),
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan(stage_order=("fake_stage",)),
    )

    _, trace = PreprocessingPipeline(
        stage_contract_registry=(fake_contract,),
    ).run_with_trace(state)

    assert len(trace) == 1
    assert trace[0].stage == "fake_stage"
    assert trace[0].operation == "fake"
    assert trace[0].parameters == {"mode": "test"}
    assert trace[0].diagnostics["policy"] == "fake"
