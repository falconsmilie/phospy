from __future__ import annotations

import pandas as pd
import pytest

from phospy.api.configs import (
    DatasetComparisonBuildingConfig,
    DatasetIntensityTransformConfig,
    DatasetNormalisationConfig,
    DatasetPreprocessingConfig,
    DatasetSiteMatrixConfig,
)
from phospy.datasets.builders.preprocessing import DatasetPreprocessor
from phospy.datasets.preprocessing.models import (
    DATASET_PREPROCESSING_STAGE_COMPARISONS,
    DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM,
    DATASET_PREPROCESSING_STAGE_MISSING_DATA,
    DATASET_PREPROCESSING_STAGE_NORMALISATION,
    DATASET_PREPROCESSING_STAGE_SITE_MATRIX,
    DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION,
    DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION,
    PreprocessingPlan,
    PreprocessingStageResult,
    PreprocessingState,
)
from phospy.datasets.preprocessing.pipeline import PreprocessingPipeline
from phospy.datasets.preprocessing.stage_registry import (
    get_preprocessing_stage_metadata,
    list_registered_preprocessing_stages,
)
from phospy.errors.build import DatasetBuildError


def _phospho() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample_a": [1.0, 2.0],
            "sample_b": [2.0, 4.0],
        },
        index=pd.Index(["MAPK14;Y182;", "AKT1;T308;"], name="site_id"),
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
        )
    )


def test_every_preprocessing_stage_has_registry_metadata() -> None:
    expected = {
        DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION,
        DATASET_PREPROCESSING_STAGE_MISSING_DATA,
        DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM,
        DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION,
        DATASET_PREPROCESSING_STAGE_SITE_MATRIX,
        DATASET_PREPROCESSING_STAGE_NORMALISATION,
        DATASET_PREPROCESSING_STAGE_COMPARISONS,
    }
    observed = {item.stage_key for item in list_registered_preprocessing_stages()}
    assert observed == expected


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
