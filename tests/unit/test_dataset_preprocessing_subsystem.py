from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd
import pandas.testing as pdt
import pytest

from phospy.api.configs import (
    DATASET_LOCALISATION_MODE_ALLOW_MISSING_WITH_WAIVER,
    DatasetComparisonBuildingConfig,
    DatasetIntensityTransformConfig,
    DatasetLocalisationConfig,
    DatasetMissingDataConfig,
    DatasetNormalisationConfig,
    DatasetPreprocessingConfig,
    DatasetRuvReadinessConfig,
    DatasetSiteMatrixConfig,
    DatasetSiteSequenceResolutionConfig,
    DatasetTotalProteinCorrectionConfig,
)
from phospy.api.requests import DatasetBuildRequest
from phospy.errors.input import PhosPyInputError
from phospy.science.datasets.builders.contracts import (
    InterpretedDatasetBuildRequest,
    PreprocessedDatasetBuildTables,
)
from phospy.science.datasets.builders.executor import DatasetBuildExecutor
from phospy.science.datasets.builders.interpreter import DatasetBuildRequestInterpreter
from phospy.science.datasets.builders.preprocessing import (
    DatasetPreprocessor,
    build_dataset_processing_state,
)
from phospy.science.datasets.builders.transformation_resolver import (
    ResolvedIntensityScale,
)
from phospy.science.datasets.preprocessing.models import (
    DATASET_PREPROCESSING_STAGE_LOCALISATION,
    PREPROCESSING_STAGE_ORDER_RATIONALE_MINPROB_MISSING_DATA,
    PREPROCESSING_STAGE_ORDER_RATIONALE_NON_MINPROB_MISSING_DATA,
    PreprocessingPlan,
    PreprocessingStageResult,
    PreprocessingState,
)
from phospy.science.datasets.preprocessing.pipeline import PreprocessingPipeline
from phospy.science.datasets.preprocessing.stage_registry import (
    PreprocessingStageMetadata,
)
from phospy.science.references.models import Organism
from phospy.science.sites.site_keys import (
    build_protein_scoped_site_key,
    encode_site_key,
)
from phospy.science.transformations.models import (
    IntensityScaleKind,
    QuantitativeMeaning,
)
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
)

ROOT = Path(__file__).resolve().parents[2]


def _phospho() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample_a": [1.0, 2.0, 3.0],
            "sample_b": [2.0, 1.0, 4.0],
            "sample_c": [3.0, 1.5, 5.0],
        },
        index=["MAPK14;Y182;", "GSK3B;S9;", "AKT1;T308;"],
    )


def _site_metadata(index: pd.Index | None = None) -> pd.DataFrame:
    data = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "GSK3B", "AKT1"],
            "site": ["Y182", "S9", "T308"],
            "site_sequence": [
                "LDFGLARHTDDEMTGYVATRWYRAPEIMLNW",
                "AAAAAAAAAAAAAAASAAAAAAAAAAAAAAA",
                "RPHFPQFSYSASGTA",
            ],
            "localisation_confidence": [0.95, 0.9, 0.92],
        },
        index=["MAPK14;Y182;", "GSK3B;S9;", "AKT1;T308;"],
    )
    return data if index is None else data.loc[index]


def _sample_metadata(columns: pd.Index) -> pd.DataFrame:
    return pd.DataFrame(
        {"batch": [1, 1, 2]},
        index=columns,
    )


def _comparison_sample_metadata(columns: pd.Index) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "comparison_group": ["group1", "group1", "group4", "group4"],
        },
        index=columns,
    )


def _ruv_site_metadata(index: pd.Index) -> pd.DataFrame:
    metadata = _site_metadata(index).copy(deep=True)
    metadata.loc[:, "is_control_feature"] = [
        position == 0 for position in range(len(metadata.index))
    ]
    return metadata


def _ruv_sample_metadata(columns: pd.Index) -> pd.DataFrame:
    sample_count = len(columns)
    midpoint = max(sample_count // 2, 1)
    return pd.DataFrame(
        {
            "replicate_group": [
                "group_a" if position < midpoint else "group_b"
                for position in range(sample_count)
            ],
            "batch": [
                f"batch_{(position % 2) + 1}" for position in range(sample_count)
            ],
        },
        index=columns.copy(),
    )


def _build_processing_state_from_preprocessor(
    *,
    phospho: pd.DataFrame,
    site_metadata: pd.DataFrame,
    sample_metadata: pd.DataFrame | None,
    config: DatasetPreprocessingConfig,
):
    if config.localisation.mode == "require_threshold":
        config = replace(
            config,
            localisation=DatasetLocalisationConfig(
                mode=DATASET_LOCALISATION_MODE_ALLOW_MISSING_WITH_WAIVER,
                waiver_reason="test helper waiver",
            ),
        )
    plan = PreprocessingPlan.from_config(config)
    preprocessed = DatasetPreprocessor().run(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=sample_metadata,
        total=None,
        plan=plan,
    )
    state = build_dataset_processing_state(
        plan=plan,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        preprocessing_trace=preprocessed.preprocessing_trace,
        final_phospho=preprocessed.phospho,
        final_site_metadata=preprocessed.site_metadata,
        final_sample_metadata=preprocessed.sample_metadata,
    )
    return preprocessed, state


def _internal_site_matrix_config(
    *,
    policy: str = "build_from_metadata",
    duplicate_site_policy: str = "max_mean_signal",
    missing_data_policy: str = "drop_any_missing",
    minimum_observed_values: int | None = None,
) -> DatasetSiteMatrixConfig:
    """Construct a site-matrix config bypassing public constructor validation.

    Subsystem tests use this to exercise internal preprocessing-lane behaviors
    that are intentionally unsupported on the public config surface.
    """

    config = object.__new__(DatasetSiteMatrixConfig)
    object.__setattr__(config, "policy", policy)
    object.__setattr__(config, "duplicate_site_policy", duplicate_site_policy)
    object.__setattr__(config, "missing_data_policy", missing_data_policy)
    object.__setattr__(config, "minimum_observed_values", minimum_observed_values)
    return config


def _plan_without_missing_stage(
    config: DatasetPreprocessingConfig,
) -> PreprocessingPlan:
    if config.localisation.mode == "require_threshold":
        config = replace(
            config,
            localisation=DatasetLocalisationConfig(
                mode=DATASET_LOCALISATION_MODE_ALLOW_MISSING_WITH_WAIVER,
                waiver_reason="test helper waiver",
            ),
        )
    plan = PreprocessingPlan.from_config(config)
    return replace(
        plan,
        stage_order=tuple(
            stage
            for stage in plan.stage_order
            if stage not in {"missing_data", DATASET_PREPROCESSING_STAGE_LOCALISATION}
        ),
    )


def _total(columns: pd.Index) -> pd.DataFrame:
    return pd.DataFrame(
        {
            str(columns[0]): [10.0, 12.0],
            str(columns[1]): [10.5, 12.5],
            str(columns[2]): [11.0, 13.0],
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


def test_preprocessing_pipeline_applies_plan_order() -> None:
    calls: list[str] = []

    class StageA:
        stage_key = "stage_a"

        def run(self, state: PreprocessingState) -> PreprocessingStageResult:
            calls.append(self.stage_key)
            return PreprocessingStageResult(state=state)

    class StageB:
        stage_key = "stage_b"

        def run(self, state: PreprocessingState) -> PreprocessingStageResult:
            calls.append(self.stage_key)
            return PreprocessingStageResult(state=state)

    state = PreprocessingState(
        phospho=_phospho(),
        site_metadata=_site_metadata(),
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan(
            missing_data_policy="forbid",
            missing_data_min_observed_values=None,
            total_protein_correction_policy="none",
            site_matrix_policy="as_input",
            comparison_building_policy="none",
            stage_order=("stage_b", "stage_a"),
        ),
    )

    pipeline = PreprocessingPipeline(
        stage_registry=(StageA(), StageB()),
        stage_metadata_registry=(
            _custom_stage_metadata("stage_a"),
            _custom_stage_metadata("stage_b"),
        ),
    )
    observed = pipeline.run(state)

    assert observed is state
    assert calls == ["stage_b", "stage_a"]


def test_preprocessing_pipeline_passes_stage_state_forward() -> None:
    observed_first_value: list[float] = []

    class AddOneStage:
        stage_key = "add_one"

        def run(self, state: PreprocessingState) -> PreprocessingStageResult:
            return PreprocessingStageResult(
                state=replace(state, phospho=state.phospho + 1.0)
            )

    class InspectStage:
        stage_key = "inspect"

        def run(self, state: PreprocessingState) -> PreprocessingStageResult:
            observed_first_value.append(float(state.phospho.iloc[0, 0]))
            return PreprocessingStageResult(
                state=replace(
                    state, site_metadata=state.site_metadata.assign(seen=True)
                )
            )

    phospho = _phospho()
    sample_metadata = _sample_metadata(phospho.columns)
    total = _total(phospho.columns)
    state = PreprocessingState(
        phospho=phospho,
        site_metadata=_site_metadata(),
        sample_metadata=sample_metadata,
        total=total,
        plan=PreprocessingPlan(
            missing_data_policy="forbid",
            missing_data_min_observed_values=None,
            total_protein_correction_policy="none",
            site_matrix_policy="as_input",
            comparison_building_policy="none",
            stage_order=("add_one", "inspect"),
        ),
    )
    pipeline = PreprocessingPipeline(
        stage_registry=(AddOneStage(), InspectStage()),
        stage_metadata_registry=(
            _custom_stage_metadata("add_one"),
            _custom_stage_metadata("inspect"),
        ),
    )

    observed = pipeline.run(state)

    assert observed_first_value == [2.0]
    assert "seen" in observed.site_metadata.columns
    assert observed.sample_metadata is sample_metadata
    assert observed.total is total


def test_dataset_preprocessor_regression_forbid_policy_is_passthrough() -> None:
    phospho = _phospho()
    site_metadata = _site_metadata()
    sample_metadata = _sample_metadata(phospho.columns)
    total = _total(phospho.columns)

    preprocessed = DatasetPreprocessor().run(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=sample_metadata,
        total=total,
        plan=PreprocessingPlan.from_config(DatasetPreprocessingConfig()),
    )

    assert preprocessed.phospho is phospho
    assert preprocessed.site_metadata is site_metadata
    assert preprocessed.sample_metadata is sample_metadata
    assert preprocessed.total is total
    assert preprocessed.row_audit is not None
    assert preprocessed.row_audit.empty


def test_dataset_preprocessor_regression_impute_row_median_policy() -> None:
    phospho = _phospho()
    phospho.loc["MAPK14;Y182;", "sample_a"] = float("nan")
    phospho.loc["GSK3B;S9;", :] = float("nan")
    site_metadata = _site_metadata()
    sample_metadata = _sample_metadata(phospho.columns)
    total = _total(phospho.columns)

    preprocessed = DatasetPreprocessor().run(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=sample_metadata,
        total=total,
        plan=PreprocessingPlan.from_config(
            DatasetPreprocessingConfig(
                missing_data=DatasetMissingDataConfig(
                    policy="impute_row_median",
                    min_observed_values=2,
                )
            )
        ),
    )

    expected_phospho = pd.DataFrame(
        {
            "sample_a": [2.5, 3.0],
            "sample_b": [2.0, 4.0],
            "sample_c": [3.0, 5.0],
        },
        index=["MAPK14;Y182;", "AKT1;T308;"],
    )
    expected_site_metadata = _site_metadata(pd.Index(["MAPK14;Y182;", "AKT1;T308;"]))

    pdt.assert_frame_equal(preprocessed.phospho, expected_phospho)
    pdt.assert_frame_equal(preprocessed.site_metadata, expected_site_metadata)
    assert preprocessed.sample_metadata is sample_metadata
    assert preprocessed.total is total
    dropped = preprocessed.row_audit.loc[
        (preprocessed.row_audit.loc[:, "stage"] == "missing_data")
        & (preprocessed.row_audit.loc[:, "action"] == "dropped")
    ]
    assert set(dropped.loc[:, "source_row_id"].astype(str)) == {"GSK3B;S9;"}
    imputed = preprocessed.row_audit.loc[
        (preprocessed.row_audit.loc[:, "stage"] == "missing_data")
        & (preprocessed.row_audit.loc[:, "action"] == "imputed")
    ]
    assert set(imputed.loc[:, "source_row_id"].astype(str)) == {"MAPK14;Y182;"}
    assert "imputed_columns" in imputed.iloc[0]["parameter_snapshot"]


def test_build_dataset_processing_state_row_median_diagnostics_include_row_medians_used() -> (
    None
):
    phospho = _phospho()
    phospho.loc["MAPK14;Y182;", "sample_a"] = float("nan")
    phospho.loc["GSK3B;S9;", :] = float("nan")
    preprocessed, state = _build_processing_state_from_preprocessor(
        phospho=phospho,
        site_metadata=_site_metadata(),
        sample_metadata=_sample_metadata(phospho.columns),
        config=DatasetPreprocessingConfig(
            missing_data=DatasetMissingDataConfig(
                policy="impute_row_median",
                min_observed_values=2,
            )
        ),
    )

    diagnostics = state.missing_data.diagnostics
    assert diagnostics is not None
    diagnostics_payload = diagnostics.to_payload()
    assert diagnostics_payload["imputation_method_id"] == "row_median"
    assert diagnostics_payload["row_medians_used"] == {"MAPK14;Y182;": 2.5}

    row_median = float(diagnostics_payload["row_medians_used"]["MAPK14;Y182;"])
    assert state.missing_data.complete_matrix is True
    assert state.missing_data.imputed is True
    assert state.missing_data.diagnostics is not None

    imputed_value = float(
        state.missing_data.diagnostics["row_medians_used"]["MAPK14;Y182;"]
    )
    assert imputed_value == pytest.approx(row_median)
    assert preprocessed.phospho.loc["MAPK14;Y182;", "sample_a"] == pytest.approx(
        row_median
    )


@pytest.mark.parametrize(
    "missing_data_config",
    [
        DatasetMissingDataConfig(
            policy="impute_minprob",
            q=0.01,
            width=0.3,
            seed=12345,
            max_missing_fraction_per_row=0.5,
        ),
        DatasetMissingDataConfig(
            policy="impute_knn",
            k=1,
            distance="nan_euclidean",
            max_missing_fraction_per_row=0.5,
        ),
    ],
)
def test_build_dataset_processing_state_non_row_median_diagnostics_have_empty_row_medians_used(
    missing_data_config: DatasetMissingDataConfig,
) -> None:
    phospho = _phospho()
    phospho.loc["MAPK14;Y182;", "sample_a"] = float("nan")
    phospho.loc["GSK3B;S9;", :] = float("nan")
    _, state = _build_processing_state_from_preprocessor(
        phospho=phospho,
        site_metadata=_site_metadata(),
        sample_metadata=_sample_metadata(phospho.columns),
        config=DatasetPreprocessingConfig(
            intensity_transform=DatasetIntensityTransformConfig(
                policy="log2",
                pseudocount=1.0,
            )
            if missing_data_config.policy == "impute_minprob"
            else DatasetIntensityTransformConfig(policy="identity"),
            missing_data=missing_data_config,
        ),
    )

    diagnostics = state.missing_data.diagnostics
    assert diagnostics is not None
    assert diagnostics.to_payload()["row_medians_used"] == {}


def test_dataset_missing_data_config_rejects_minprob_without_seed() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="missing_data.seed must be an int",
    ):
        DatasetMissingDataConfig(
            policy="impute_minprob",
            q=0.01,
            width=0.3,
            max_missing_fraction_per_row=0.5,
        )


def test_dataset_missing_data_config_rejects_minprob_invalid_q() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="missing_data.q must satisfy 0 < q < 0.5",
    ):
        DatasetMissingDataConfig(
            policy="impute_minprob",
            q=0.5,
            width=0.3,
            seed=12345,
            max_missing_fraction_per_row=0.5,
        )


def test_dataset_missing_data_config_rejects_minprob_invalid_width() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="missing_data.width must satisfy 0 < width <= 1.0",
    ):
        DatasetMissingDataConfig(
            policy="impute_minprob",
            q=0.01,
            width=0.0,
            seed=12345,
            max_missing_fraction_per_row=0.5,
        )


def test_dataset_missing_data_config_rejects_knn_without_k() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="missing_data.k must be an int",
    ):
        DatasetMissingDataConfig(
            policy="impute_knn",
            distance="nan_euclidean",
            max_missing_fraction_per_row=0.5,
        )


def test_dataset_missing_data_config_rejects_knn_unsupported_distance() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="missing_data.distance must be 'nan_euclidean'",
    ):
        DatasetMissingDataConfig(
            policy="impute_knn",
            k=2,
            distance="euclidean",
            max_missing_fraction_per_row=0.5,
        )


def test_ruv_readiness_disabled_reports_not_configured() -> None:
    _, state = _build_processing_state_from_preprocessor(
        phospho=_phospho(),
        site_metadata=_site_metadata(),
        sample_metadata=None,
        config=DatasetPreprocessingConfig(),
    )

    assert state.ruv_readiness.enabled is False
    assert state.ruv_readiness.ready is False
    assert "not configured" in set(state.ruv_readiness.reasons)


def test_ruv_readiness_enabled_reports_missing_control_feature_column() -> None:
    _, state = _build_processing_state_from_preprocessor(
        phospho=_phospho(),
        site_metadata=_site_metadata(),
        sample_metadata=_ruv_sample_metadata(_phospho().columns),
        config=DatasetPreprocessingConfig(
            ruv_readiness=DatasetRuvReadinessConfig(enabled=True)
        ),
    )

    assert state.ruv_readiness.enabled is True
    assert state.ruv_readiness.ready is False
    assert "control feature column missing" in set(state.ruv_readiness.reasons)


def test_ruv_readiness_enabled_reports_missing_sample_metadata() -> None:
    _, state = _build_processing_state_from_preprocessor(
        phospho=_phospho(),
        site_metadata=_ruv_site_metadata(_phospho().index),
        sample_metadata=None,
        config=DatasetPreprocessingConfig(
            ruv_readiness=DatasetRuvReadinessConfig(enabled=True)
        ),
    )

    assert state.ruv_readiness.enabled is True
    assert state.ruv_readiness.ready is False
    assert "sample metadata unavailable" in set(state.ruv_readiness.reasons)


def test_ruv_readiness_enabled_reports_missing_replicate_group_column() -> None:
    phospho = _phospho()
    sample_metadata = pd.DataFrame(
        {"batch": ["batch_1", "batch_1", "batch_2"]},
        index=phospho.columns.copy(),
    )
    _, state = _build_processing_state_from_preprocessor(
        phospho=phospho,
        site_metadata=_ruv_site_metadata(phospho.index),
        sample_metadata=sample_metadata,
        config=DatasetPreprocessingConfig(
            ruv_readiness=DatasetRuvReadinessConfig(enabled=True)
        ),
    )

    assert state.ruv_readiness.enabled is True
    assert state.ruv_readiness.ready is False
    assert "replicate group column missing" in set(state.ruv_readiness.reasons)


def test_ruv_readiness_enabled_reports_ready_when_requirements_are_met() -> None:
    phospho = _phospho()
    _, state = _build_processing_state_from_preprocessor(
        phospho=phospho,
        site_metadata=_ruv_site_metadata(phospho.index),
        sample_metadata=_ruv_sample_metadata(phospho.columns),
        config=DatasetPreprocessingConfig(
            ruv_readiness=DatasetRuvReadinessConfig(enabled=True)
        ),
    )

    assert state.ruv_readiness.enabled is True
    assert state.ruv_readiness.ready is True
    assert state.ruv_readiness.reasons == ()
    assert state.ruv_readiness.matrix_complete is True
    assert state.ruv_readiness.control_feature_count >= 1
    assert state.ruv_readiness.replicate_group_count >= 2
    assert (state.ruv_readiness.batch_count or 0) >= 1


@pytest.mark.parametrize(
    ("missing_data_config", "expected_method_id"),
    [
        (
            DatasetMissingDataConfig(
                policy="impute_row_median",
                min_observed_values=2,
            ),
            "row_median",
        ),
        (
            DatasetMissingDataConfig(
                policy="impute_minprob",
                q=0.01,
                width=0.3,
                seed=12345,
                max_missing_fraction_per_row=0.5,
            ),
            "minprob",
        ),
        (
            DatasetMissingDataConfig(
                policy="impute_knn",
                k=1,
                distance="nan_euclidean",
                max_missing_fraction_per_row=0.5,
            ),
            "knn",
        ),
    ],
)
def test_ruv_readiness_reports_imputation_method_id_from_missing_data_diagnostics(
    missing_data_config: DatasetMissingDataConfig,
    expected_method_id: str,
) -> None:
    phospho = _phospho()
    phospho.loc["MAPK14;Y182;", "sample_a"] = float("nan")
    phospho.loc["GSK3B;S9;", :] = float("nan")
    config = DatasetPreprocessingConfig(
        intensity_transform=DatasetIntensityTransformConfig(
            policy="log2",
            pseudocount=1.0,
        )
        if missing_data_config.policy == "impute_minprob"
        else DatasetIntensityTransformConfig(policy="identity"),
        missing_data=missing_data_config,
        ruv_readiness=DatasetRuvReadinessConfig(enabled=True),
    )
    _, state = _build_processing_state_from_preprocessor(
        phospho=phospho,
        site_metadata=_ruv_site_metadata(phospho.index),
        sample_metadata=_ruv_sample_metadata(phospho.columns),
        config=config,
    )

    assert state.ruv_readiness.imputation_method_id == expected_method_id
    assert state.ruv_readiness.missingness_mask_preserved is True


def test_dataset_preprocessor_knn_imputes_expected_values_and_preserves_labels() -> (
    None
):
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 1.0, 2.0, 10.0],
            "sample_b": [1.0, 2.0, 2.0, float("nan")],
            "sample_c": [float("nan"), 3.0, 3.0, float("nan")],
        },
        index=pd.Index(["row_impute", "row_ref_1", "row_ref_2", "row_drop"]),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "AKT1", "GSK3B", "PRKACA"],
            "site": ["Y182", "T308", "S9", "S339"],
            "site_sequence": ["SEQ_A", "SEQ_R", "SEQ_C", "SEQ_D"],
            "localisation_confidence": [0.95, 0.96, 0.9, 0.88],
        },
        index=phospho.index.copy(),
    )
    preprocessed = DatasetPreprocessor().run(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan.from_config(
            DatasetPreprocessingConfig(
                missing_data=DatasetMissingDataConfig(
                    policy="impute_knn",
                    k=1,
                    distance="nan_euclidean",
                    max_missing_fraction_per_row=0.5,
                )
            )
        ),
    )

    expected = pd.DataFrame(
        {
            "sample_a": [1.0, 1.0, 2.0],
            "sample_b": [1.0, 2.0, 2.0],
            "sample_c": [3.0, 3.0, 3.0],
        },
        index=pd.Index(["row_impute", "row_ref_1", "row_ref_2"]),
    )
    pdt.assert_frame_equal(preprocessed.phospho, expected)
    assert preprocessed.phospho.index.tolist() == [
        "row_impute",
        "row_ref_1",
        "row_ref_2",
    ]
    assert preprocessed.phospho.columns.tolist() == ["sample_a", "sample_b", "sample_c"]
    assert int(preprocessed.phospho.isna().to_numpy().sum()) == 0
    assert preprocessed.row_audit is not None
    dropped = preprocessed.row_audit.loc[
        (preprocessed.row_audit.loc[:, "stage"] == "missing_data")
        & (preprocessed.row_audit.loc[:, "action"] == "dropped")
    ]
    assert set(dropped.loc[:, "source_row_id"].astype(str)) == {"row_drop"}


def test_preprocessing_plan_rejects_minprob_without_log2_transform() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="missing_data.policy='impute_minprob' requires",
    ):
        PreprocessingPlan.from_config(
            DatasetPreprocessingConfig(
                intensity_transform=DatasetIntensityTransformConfig(policy="identity"),
                missing_data=DatasetMissingDataConfig(
                    policy="impute_minprob",
                    q=0.01,
                    width=0.3,
                    seed=12345,
                    max_missing_fraction_per_row=0.5,
                ),
            )
        )


def test_preprocessing_plan_orders_minprob_after_intensity_transform() -> None:
    plan = PreprocessingPlan.from_config(
        DatasetPreprocessingConfig(
            intensity_transform=DatasetIntensityTransformConfig(
                policy="log2",
                pseudocount=1.0,
            ),
            missing_data=DatasetMissingDataConfig(
                policy="impute_minprob",
                q=0.01,
                width=0.3,
                seed=12345,
                max_missing_fraction_per_row=0.5,
            ),
            total_protein_correction=DatasetTotalProteinCorrectionConfig(
                policy="subtract_log_total"
            ),
            site_matrix=DatasetSiteMatrixConfig(policy="build_from_metadata"),
        )
    )

    assert plan.stage_order == (
        "localisation_confidence",
        "intensity_transform",
        "missing_data",
        "total_protein_correction",
        "site_matrix",
    )
    assert plan.stage_order_resolution[1].stage == "intensity_transform"
    assert plan.stage_order_resolution[2].stage == "missing_data"
    assert (
        plan.stage_order_resolution[2].rationale
        == PREPROCESSING_STAGE_ORDER_RATIONALE_MINPROB_MISSING_DATA
    )


@pytest.mark.parametrize(
    "missing_data",
    (
        DatasetMissingDataConfig(policy="impute_row_median", min_observed_values=1),
        DatasetMissingDataConfig(
            policy="impute_knn",
            k=1,
            distance="nan_euclidean",
            max_missing_fraction_per_row=0.5,
        ),
    ),
)
def test_preprocessing_plan_orders_non_minprob_before_log2_transform(
    missing_data: DatasetMissingDataConfig,
) -> None:
    plan = PreprocessingPlan.from_config(
        DatasetPreprocessingConfig(
            intensity_transform=DatasetIntensityTransformConfig(
                policy="log2",
                pseudocount=1.0,
            ),
            missing_data=missing_data,
        )
    )

    assert plan.stage_order[:3] == (
        "localisation_confidence",
        "missing_data",
        "intensity_transform",
    )
    assert plan.stage_order_resolution[1].stage == "missing_data"
    assert (
        plan.stage_order_resolution[1].rationale
        == PREPROCESSING_STAGE_ORDER_RATIONALE_NON_MINPROB_MISSING_DATA
    )


def test_preprocessing_plan_orders_fasta_resolution_before_minprob_with_log2(
    tmp_path: Path,
) -> None:
    fasta_path = tmp_path / "proteins.fasta"
    fasta_path.write_text(">P1\nAAAASAAAA\n", encoding="utf-8")
    plan = PreprocessingPlan.from_config(
        DatasetPreprocessingConfig(
            site_sequence_resolution=DatasetSiteSequenceResolutionConfig(
                fasta_path=str(fasta_path),
                flank_size=2,
            ),
            intensity_transform=DatasetIntensityTransformConfig(
                policy="log2",
                pseudocount=1.0,
            ),
            missing_data=DatasetMissingDataConfig(
                policy="impute_minprob",
                q=0.01,
                width=0.3,
                seed=12345,
                max_missing_fraction_per_row=0.5,
            ),
        )
    )

    assert plan.stage_order == (
        "site_sequence_resolution",
        "localisation_confidence",
        "intensity_transform",
        "missing_data",
    )


def test_preprocessing_plan_rejects_minprob_without_log2_when_fasta_enabled() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="missing_data.policy='impute_minprob' requires",
    ):
        PreprocessingPlan.from_config(
            DatasetPreprocessingConfig(
                site_sequence_resolution=DatasetSiteSequenceResolutionConfig(
                    fasta_path="local.fasta",
                ),
                intensity_transform=DatasetIntensityTransformConfig(policy="identity"),
                missing_data=DatasetMissingDataConfig(
                    policy="impute_minprob",
                    q=0.01,
                    width=0.3,
                    seed=12345,
                    max_missing_fraction_per_row=0.5,
                ),
            )
        )


def test_dataset_preprocessor_minprob_is_deterministic_for_same_seed() -> None:
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
            "localisation_confidence": [0.95, 0.96, 0.9, 0.88],
        },
        index=phospho.index.copy(),
    )
    config = DatasetPreprocessingConfig(
        intensity_transform=DatasetIntensityTransformConfig(
            policy="log2", pseudocount=1.0
        ),
        missing_data=DatasetMissingDataConfig(
            policy="impute_minprob",
            q=0.01,
            width=0.3,
            seed=12345,
            max_missing_fraction_per_row=0.5,
        ),
    )
    plan = PreprocessingPlan.from_config(config)

    first = DatasetPreprocessor().run(
        phospho=phospho.copy(deep=True),
        site_metadata=site_metadata.copy(deep=True),
        sample_metadata=None,
        total=None,
        plan=plan,
    )
    second = DatasetPreprocessor().run(
        phospho=phospho.copy(deep=True),
        site_metadata=site_metadata.copy(deep=True),
        sample_metadata=None,
        total=None,
        plan=plan,
    )

    pdt.assert_frame_equal(first.phospho, second.phospho)
    assert first.phospho.isna().to_numpy().sum() == 0
    assert first.phospho.index.tolist() == ["row_keep", "row_impute_a", "row_impute_c"]


def test_dataset_preprocessor_minprob_changes_with_seed() -> None:
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
            "localisation_confidence": [0.95, 0.96, 0.9, 0.88],
        },
        index=phospho.index.copy(),
    )

    first_plan = PreprocessingPlan.from_config(
        DatasetPreprocessingConfig(
            intensity_transform=DatasetIntensityTransformConfig(
                policy="log2",
                pseudocount=1.0,
            ),
            missing_data=DatasetMissingDataConfig(
                policy="impute_minprob",
                q=0.01,
                width=0.3,
                seed=1,
                max_missing_fraction_per_row=0.5,
            ),
        )
    )
    second_plan = PreprocessingPlan.from_config(
        DatasetPreprocessingConfig(
            intensity_transform=DatasetIntensityTransformConfig(
                policy="log2",
                pseudocount=1.0,
            ),
            missing_data=DatasetMissingDataConfig(
                policy="impute_minprob",
                q=0.01,
                width=0.3,
                seed=2,
                max_missing_fraction_per_row=0.5,
            ),
        )
    )

    first = DatasetPreprocessor().run(
        phospho=phospho.copy(deep=True),
        site_metadata=site_metadata.copy(deep=True),
        sample_metadata=None,
        total=None,
        plan=first_plan,
    )
    second = DatasetPreprocessor().run(
        phospho=phospho.copy(deep=True),
        site_metadata=site_metadata.copy(deep=True),
        sample_metadata=None,
        total=None,
        plan=second_plan,
    )

    assert first.phospho.index.tolist() == second.phospho.index.tolist()
    assert float(first.phospho.loc["row_impute_a", "sample_a"]) != pytest.approx(
        float(second.phospho.loc["row_impute_a", "sample_a"])
    )


def test_preprocessing_plan_orders_intensity_transform_before_total_correction() -> (
    None
):
    plan = PreprocessingPlan.from_config(
        DatasetPreprocessingConfig(
            intensity_transform=DatasetIntensityTransformConfig(
                policy="log2",
                pseudocount=1.0,
            ),
            total_protein_correction=DatasetTotalProteinCorrectionConfig(
                policy="subtract_log_total"
            ),
            site_matrix=DatasetSiteMatrixConfig(policy="build_from_metadata"),
        )
    )
    assert plan.stage_order == (
        "localisation_confidence",
        "missing_data",
        "intensity_transform",
        "total_protein_correction",
        "site_matrix",
    )


def test_preprocessing_plan_orders_comparisons_after_upstream_stages() -> None:
    plan = PreprocessingPlan.from_config(
        DatasetPreprocessingConfig(
            intensity_transform=DatasetIntensityTransformConfig(
                policy="log2",
                pseudocount=1.0,
            ),
            total_protein_correction=DatasetTotalProteinCorrectionConfig(
                policy="subtract_log_total"
            ),
            site_matrix=DatasetSiteMatrixConfig(policy="build_from_metadata"),
            comparisons=DatasetComparisonBuildingConfig(policy="sample_metadata_pairs"),
        )
    )
    assert plan.stage_order == (
        "localisation_confidence",
        "missing_data",
        "intensity_transform",
        "total_protein_correction",
        "site_matrix",
        "comparisons",
    )


def test_dataset_preprocessor_applies_subtract_log_total_policy() -> None:
    phospho = _phospho()
    site_metadata = _site_metadata()
    total = pd.DataFrame(
        {
            "sample_a": [0.5, 2.0, 1.0],
            "sample_b": [1.0, 1.0, 1.0],
            "sample_c": [1.5, 0.5, 1.0],
        },
        index=pd.Index(["MAPK14", "GSK3B", "AKT1"], name="protein_id"),
    )

    preprocessed = DatasetPreprocessor().run(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=None,
        total=total,
        plan=PreprocessingPlan.from_config(
            DatasetPreprocessingConfig(
                intensity_transform=DatasetIntensityTransformConfig(
                    policy="log2",
                    pseudocount=1.0,
                ),
                total_protein_correction=DatasetTotalProteinCorrectionConfig(
                    policy="subtract_log_total"
                ),
            )
        ),
    )

    total_by_site = total.reindex(site_metadata.loc[:, "gene_symbol"].tolist())
    total_by_site.index = phospho.index.copy()
    expected = np.log2(phospho + 1.0) - np.log2(total_by_site + 1.0)
    pdt.assert_frame_equal(preprocessed.phospho, expected)
    assert preprocessed.site_metadata is site_metadata
    assert preprocessed.total is not None
    pdt.assert_frame_equal(preprocessed.total, np.log2(total + 1.0))


def test_dataset_preprocessor_total_correction_uses_log_ratio_formula() -> None:
    phospho = pd.DataFrame(
        {"sample_a": [15.0]},
        index=pd.Index(["MAPK14;Y182;"], name="site_id"),
    )
    total = pd.DataFrame(
        {"sample_a": [3.0]},
        index=pd.Index(["MAPK14"], name="protein_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14"],
            "site": ["Y182"],
            "site_sequence": ["SEQ_A"],
            "localisation_confidence": [0.95],
        },
        index=phospho.index.copy(),
    )

    preprocessed = DatasetPreprocessor().run(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=None,
        total=total,
        plan=PreprocessingPlan.from_config(
            DatasetPreprocessingConfig(
                intensity_transform=DatasetIntensityTransformConfig(
                    policy="log2",
                    pseudocount=1.0,
                ),
                total_protein_correction=DatasetTotalProteinCorrectionConfig(
                    policy="subtract_log_total"
                ),
            )
        ),
    )

    observed = float(preprocessed.phospho.iloc[0, 0])
    assert observed == pytest.approx(2.0)
    assert observed != pytest.approx(np.log2(15.0 - 3.0 + 1.0))


def test_dataset_preprocessor_builds_site_matrix_from_metadata_policy() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [2.0, 8.0, 5.0, 1.0],
            "sample_b": [2.5, 8.5, 5.5, float("nan")],
        },
        index=pd.Index(
            ["row_a", "row_b", "row_c", "row_d"],
            name="input_row",
        ),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "MAPK14", "AKT1", "GSK3B"],
            "site": ["Y182", "Y182", "T308", "S9"],
            "site_sequence": ["SEQ_A", "SEQ_R", "", "SEQ_D"],
            "source_uid": ["UID_A", "UID_B", "UID_C", "UID_D"],
            "localisation_confidence": [0.95, 0.92, 0.91, 0.9],
        },
        index=phospho.index.copy(),
    )

    preprocessed = DatasetPreprocessor().run(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=None,
        total=None,
        plan=_plan_without_missing_stage(
            DatasetPreprocessingConfig(
                site_matrix=DatasetSiteMatrixConfig(policy="build_from_metadata")
            )
        ),
    )

    expected_phospho = pd.DataFrame(
        {
            "sample_a": [8.0],
            "sample_b": [8.5],
        },
        index=pd.Index(["MAPK14;Y182;"], name="input_row"),
    )
    expected_site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14"],
            "site": ["Y182"],
            "site_sequence": ["SEQ_R"],
            "source_uid": ["UID_B"],
            "localisation_confidence": [0.92],
        },
        index=expected_phospho.index.copy(),
    )

    pdt.assert_frame_equal(preprocessed.phospho, expected_phospho)
    pdt.assert_frame_equal(preprocessed.site_metadata, expected_site_metadata)


def test_dataset_preprocessor_site_matrix_retain_missing_policy_keeps_partial_rows() -> (
    None
):
    phospho = pd.DataFrame(
        {
            "sample_a": [2.0, 8.0],
            "sample_b": [float("nan"), 8.5],
        },
        index=pd.Index(["row_a", "row_b"], name="input_row"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "AKT1"],
            "site": ["Y182", "T308"],
            "site_sequence": ["SEQ_A", "SEQ_R"],
            "localisation_confidence": [0.95, 0.96],
        },
        index=phospho.index.copy(),
    )

    preprocessed = DatasetPreprocessor().run(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=None,
        total=None,
        plan=_plan_without_missing_stage(
            DatasetPreprocessingConfig(
                site_matrix=_internal_site_matrix_config(
                    policy="build_from_metadata",
                    missing_data_policy="retain_missing",
                )
            )
        ),
    )

    assert preprocessed.phospho.index.tolist() == ["AKT1;T308;", "MAPK14;Y182;"]
    assert pd.isna(preprocessed.phospho.loc["MAPK14;Y182;", "sample_b"])
    dropped = preprocessed.row_audit.loc[
        (preprocessed.row_audit.loc[:, "stage"] == "site_matrix")
        & (preprocessed.row_audit.loc[:, "action"] == "dropped")
    ]
    assert dropped.empty


def test_dataset_preprocessor_site_matrix_supports_min_observed_and_duplicate_aggregate_mean() -> (
    None
):
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 3.0, float("nan")],
            "sample_b": [2.0, 4.0, float("nan")],
            "sample_c": [float("nan"), float("nan"), 9.0],
        },
        index=pd.Index(["row_a", "row_b", "row_c"], name="input_row"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "MAPK14", "AKT1"],
            "site": ["Y182", "Y182", "T308"],
            "site_sequence": ["SEQ_A", "SEQ_R", "SEQ_C"],
            "uid": ["A", "B", "C"],
            "localisation_confidence": [0.95, 0.92, 0.97],
        },
        index=phospho.index.copy(),
    )

    preprocessed = DatasetPreprocessor().run(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=None,
        total=None,
        plan=_plan_without_missing_stage(
            DatasetPreprocessingConfig(
                site_matrix=_internal_site_matrix_config(
                    policy="build_from_metadata",
                    missing_data_policy="require_min_observed_values",
                    minimum_observed_values=1,
                    duplicate_site_policy="aggregate_mean",
                )
            )
        ),
    )

    assert preprocessed.phospho.index.tolist() == ["AKT1;T308;", "MAPK14;Y182;"]
    assert preprocessed.phospho.loc["MAPK14;Y182;", "sample_a"] == pytest.approx(2.0)
    assert preprocessed.phospho.loc["MAPK14;Y182;", "sample_b"] == pytest.approx(3.0)
    assert pd.isna(preprocessed.phospho.loc["MAPK14;Y182;", "sample_c"])
    assert pd.isna(preprocessed.site_metadata.loc["MAPK14;Y182;", "uid"])


def test_dataset_preprocessor_site_matrix_duplicate_first_policy_keeps_first_row() -> (
    None
):
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 3.0],
            "sample_b": [2.0, 4.0],
        },
        index=pd.Index(["row_a", "row_b"], name="input_row"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "MAPK14"],
            "site": ["Y182", "Y182"],
            "site_sequence": ["SEQ_A", "SEQ_R"],
            "uid": ["A", "B"],
            "localisation_confidence": [0.95, 0.92],
        },
        index=phospho.index.copy(),
    )

    preprocessed = DatasetPreprocessor().run(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan.from_config(
            DatasetPreprocessingConfig(
                site_matrix=_internal_site_matrix_config(
                    policy="build_from_metadata",
                    duplicate_site_policy="first",
                    missing_data_policy="retain_missing",
                )
            )
        ),
    )

    assert preprocessed.phospho.index.tolist() == ["MAPK14;Y182;"]
    assert preprocessed.phospho.loc["MAPK14;Y182;", "sample_a"] == pytest.approx(1.0)
    assert preprocessed.site_metadata.loc["MAPK14;Y182;", "uid"] == "A"


def test_dataset_preprocessor_site_matrix_duplicate_aggregate_median_policy() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 30.0, 3.0],
            "sample_b": [2.0, 40.0, 4.0],
        },
        index=pd.Index(["row_a", "row_b", "row_c"], name="input_row"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "MAPK14", "AKT1"],
            "site": ["Y182", "Y182", "T308"],
            "site_sequence": ["SEQ_A", "SEQ_R", "SEQ_C"],
            "uid": ["A", "B", "C"],
            "localisation_confidence": [0.95, 0.92, 0.97],
        },
        index=phospho.index.copy(),
    )

    preprocessed = DatasetPreprocessor().run(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan.from_config(
            DatasetPreprocessingConfig(
                site_matrix=_internal_site_matrix_config(
                    policy="build_from_metadata",
                    duplicate_site_policy="aggregate_median",
                    missing_data_policy="retain_missing",
                )
            )
        ),
    )

    assert preprocessed.phospho.index.tolist() == ["AKT1;T308;", "MAPK14;Y182;"]
    assert preprocessed.phospho.loc["MAPK14;Y182;", "sample_a"] == pytest.approx(15.5)
    assert preprocessed.phospho.loc["MAPK14;Y182;", "sample_b"] == pytest.approx(21.0)
    assert pd.isna(preprocessed.site_metadata.loc["MAPK14;Y182;", "uid"])


def test_dataset_preprocessor_rejects_site_matrix_min_observed_above_sample_count() -> (
    None
):
    phospho = pd.DataFrame(
        {"sample_a": [1.0], "sample_b": [2.0]},
        index=pd.Index(["row_a"], name="input_row"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14"],
            "site": ["Y182"],
            "site_sequence": ["SEQ_A"],
            "localisation_confidence": [0.95],
        },
        index=phospho.index.copy(),
    )

    with pytest.raises(
        PhosPyInputError,
        match="minimum_observed_values cannot exceed phospho sample count",
    ):
        DatasetPreprocessor().run(
            phospho=phospho,
            site_metadata=site_metadata,
            sample_metadata=None,
            total=None,
            plan=PreprocessingPlan.from_config(
                DatasetPreprocessingConfig(
                    site_matrix=_internal_site_matrix_config(
                        policy="build_from_metadata",
                        missing_data_policy="require_min_observed_values",
                        minimum_observed_values=3,
                    )
                )
            ),
        )


def test_dataset_preprocessor_rejects_site_matrix_duplicate_rows_in_error_mode() -> (
    None
):
    phospho = pd.DataFrame(
        {"sample_a": [1.0, 3.0], "sample_b": [2.0, 4.0]},
        index=pd.Index(["row_a", "row_b"], name="input_row"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "MAPK14"],
            "site": ["Y182", "Y182"],
            "site_sequence": ["SEQ_A", "SEQ_R"],
            "localisation_confidence": [0.95, 0.92],
        },
        index=phospho.index.copy(),
    )

    with pytest.raises(
        PhosPyInputError,
        match="duplicate_site_policy='error'",
    ):
        DatasetPreprocessor().run(
            phospho=phospho,
            site_metadata=site_metadata,
            sample_metadata=None,
            total=None,
            plan=PreprocessingPlan.from_config(
                DatasetPreprocessingConfig(
                    site_matrix=DatasetSiteMatrixConfig(
                        policy="build_from_metadata",
                        duplicate_site_policy="error",
                    )
                )
            ),
        )


def test_dataset_preprocessor_rejects_site_matrix_build_without_site_sequence() -> None:
    phospho = _phospho()
    site_metadata = _site_metadata().drop(columns=["site_sequence"])

    with pytest.raises(
        PhosPyInputError,
        match=(
            "site-matrix construction produced no retained rows after filtering; "
            "input_rows=3, dropped_missing_sequence=3"
        ),
    ):
        DatasetPreprocessor().run(
            phospho=phospho,
            site_metadata=site_metadata,
            sample_metadata=None,
            total=None,
            plan=PreprocessingPlan.from_config(
                DatasetPreprocessingConfig(
                    site_matrix=DatasetSiteMatrixConfig(policy="build_from_metadata")
                )
            ),
        )


def test_dataset_preprocessor_rejects_correction_when_total_columns_mismatch() -> None:
    phospho = _phospho()
    site_metadata = _site_metadata()
    total = pd.DataFrame(
        {
            "sample_a": [1.0, 2.0, 3.0],
            "sample_b": [1.0, 2.0, 3.0],
        },
        index=pd.Index(["MAPK14", "GSK3B", "AKT1"], name="protein_id"),
    )

    with pytest.raises(
        PhosPyInputError,
        match="requires total columns to exactly match phospho columns",
    ):
        DatasetPreprocessor().run(
            phospho=phospho,
            site_metadata=site_metadata,
            sample_metadata=None,
            total=total,
            plan=PreprocessingPlan.from_config(
                DatasetPreprocessingConfig(
                    intensity_transform=DatasetIntensityTransformConfig(
                        policy="log2",
                        pseudocount=1.0,
                    ),
                    total_protein_correction=DatasetTotalProteinCorrectionConfig(
                        policy="subtract_log_total"
                    ),
                )
            ),
        )


def test_dataset_preprocessor_rejects_correction_when_proteins_are_unmatched() -> None:
    phospho = _phospho()
    site_metadata = _site_metadata()
    total = _total(phospho.columns)

    with pytest.raises(
        PhosPyInputError,
        match="requires complete phosphosite-to-total mapping under the configured identity policy",
    ):
        DatasetPreprocessor().run(
            phospho=phospho,
            site_metadata=site_metadata,
            sample_metadata=None,
            total=total,
            plan=PreprocessingPlan.from_config(
                DatasetPreprocessingConfig(
                    intensity_transform=DatasetIntensityTransformConfig(
                        policy="log2",
                        pseudocount=1.0,
                    ),
                    total_protein_correction=DatasetTotalProteinCorrectionConfig(
                        policy="subtract_log_total"
                    ),
                )
            ),
        )


def test_dataset_preprocessor_subtract_log_total_matches_historical_baseline_fixture() -> (
    None
):
    phospho_fixture = pd.read_csv(
        ROOT
        / "tests"
        / "fixtures"
        / "rewrite_parity"
        / "protein_correction"
        / "reference_input_phospho.csv"
    ).set_index("site_id")
    site_metadata_fixture = pd.read_csv(
        ROOT
        / "tests"
        / "fixtures"
        / "rewrite_parity"
        / "protein_correction"
        / "reference_input_site_metadata.csv"
    ).set_index("site_id")
    total_fixture = pd.read_csv(
        ROOT
        / "tests"
        / "fixtures"
        / "rewrite_parity"
        / "protein_correction"
        / "reference_input_total.csv"
    ).set_index("protein_id")
    corrected_fixture = pd.read_csv(
        ROOT
        / "tests"
        / "fixtures"
        / "rewrite_parity"
        / "protein_correction"
        / "reference_corrected_matrix.csv"
    ).set_index("site_id")

    phospho_columns = tuple(phospho_fixture.columns.astype(str))
    phospho = phospho_fixture.astype(float).copy()
    phospho.columns = pd.Index(phospho_columns)
    phospho.index = pd.Index(phospho_fixture.index.astype(str), name="site_id")
    site_metadata = site_metadata_fixture.copy()
    site_metadata.index = pd.Index(
        site_metadata_fixture.index.astype(str),
        name="site_id",
    )
    total = total_fixture.astype(float).copy()
    total.columns = pd.Index(phospho_columns)
    total.index = pd.Index(total_fixture.index.astype(str), name="protein_id")

    expected = corrected_fixture.loc[:, list(phospho_columns)].astype(float).copy()
    expected.index = pd.Index(corrected_fixture.index.astype(str), name="site_id")

    preprocessed = DatasetPreprocessor().run(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=None,
        total=total,
        plan=_plan_without_missing_stage(
            DatasetPreprocessingConfig(
                intensity_transform=DatasetIntensityTransformConfig(
                    policy="log2",
                    pseudocount=1.0,
                ),
                total_protein_correction=DatasetTotalProteinCorrectionConfig(
                    policy="subtract_log_total"
                ),
            )
        ),
    )

    pdt.assert_frame_equal(preprocessed.phospho, expected)


def test_dataset_preprocessor_site_matrix_build_matches_historical_baseline_fixture() -> (
    None
):
    corrected_fixture = pd.read_csv(
        ROOT
        / "tests"
        / "fixtures"
        / "rewrite_parity"
        / "site_matrix"
        / "reference_phospho_corrected.csv"
    )
    expected_matrix_fixture = pd.read_csv(
        ROOT
        / "tests"
        / "fixtures"
        / "rewrite_parity"
        / "site_matrix"
        / "reference_expected_matrix.csv"
    )
    expected_input_fixture = pd.read_csv(
        ROOT
        / "tests"
        / "fixtures"
        / "rewrite_parity"
        / "site_matrix"
        / "reference_expected_input.csv"
    )

    corrected_cols = tuple(f"phospho_corrected_{position}" for position in range(1, 7))
    phospho = corrected_fixture.loc[:, list(corrected_cols)].astype(float).copy()
    phospho.index = pd.Index(
        corrected_fixture.loc[:, "uid"].astype(str),
        name="source_uid",
    )

    site_tokens = (
        corrected_fixture.loc[:, "gene_p_site"]
        .astype(str)
        .str.split("_", n=1, expand=True)
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": corrected_fixture.loc[:, "gene_names"].astype(str).tolist(),
            "site": site_tokens.loc[:, 1].astype(str).tolist(),
            "site_sequence": corrected_fixture.loc[:, "centralized_sequence"]
            .astype(str)
            .tolist(),
            "localisation_confidence": [0.95] * len(phospho.index),
        },
        index=phospho.index.copy(),
    )

    preprocessed = DatasetPreprocessor().run(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=None,
        total=None,
        plan=_plan_without_missing_stage(
            DatasetPreprocessingConfig(
                site_matrix=DatasetSiteMatrixConfig(policy="build_from_metadata")
            )
        ),
    )

    expected_phospho = (
        expected_matrix_fixture.set_index(expected_matrix_fixture.columns[0])
        .loc[:, list(corrected_cols)]
        .astype(float)
    )
    expected_phospho.index = pd.Index(
        expected_phospho.index.astype(str), name="source_uid"
    )
    expected_site_metadata = pd.DataFrame(
        {
            "gene_symbol": expected_input_fixture.loc[:, "gene_names"]
            .astype(str)
            .tolist(),
            "site": expected_input_fixture.loc[:, "p_site"].astype(str).tolist(),
            "site_sequence": expected_input_fixture.loc[:, "centralized_sequence"]
            .astype(str)
            .tolist(),
            "localisation_confidence": [0.95] * int(expected_input_fixture.shape[0]),
        },
        index=pd.Index(
            expected_input_fixture.loc[:, "site_id"].astype(str),
            name="source_uid",
        ),
    )

    pdt.assert_frame_equal(preprocessed.phospho, expected_phospho)
    pdt.assert_frame_equal(
        preprocessed.site_metadata.loc[
            :, ["gene_symbol", "site", "site_sequence", "localisation_confidence"]
        ],
        expected_site_metadata.loc[
            :, ["gene_symbol", "site", "site_sequence", "localisation_confidence"]
        ],
    )


def test_dataset_preprocessor_builds_inferred_comparisons_from_sample_groups() -> None:
    phospho = pd.DataFrame(
        {
            "sample_1": [8.0, 2.0],
            "sample_2": [8.0, 4.0],
            "sample_3": [5.0, 1.0],
            "sample_4": [5.0, 1.0],
        },
        index=pd.Index(["MAPK14;Y182;", "AKT1;T308;"], name="site_id"),
    )
    site_metadata = _site_metadata(phospho.index)
    sample_metadata = _comparison_sample_metadata(phospho.columns)

    preprocessed = DatasetPreprocessor().run(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=sample_metadata,
        total=None,
        plan=PreprocessingPlan.from_config(
            DatasetPreprocessingConfig(
                comparisons=DatasetComparisonBuildingConfig(
                    policy="sample_metadata_pairs"
                )
            )
        ),
    )

    assert preprocessed.comparisons is not None
    expected = pd.DataFrame(
        {"p_group1_group4": [3.0, 2.0]},
        index=phospho.index.copy(),
    )
    pdt.assert_frame_equal(preprocessed.comparisons, expected)


def test_dataset_preprocessor_builds_explicit_comparison_pairs() -> None:
    phospho = pd.DataFrame(
        {
            "sample_1": [8.0],
            "sample_2": [8.0],
            "sample_3": [5.0],
            "sample_4": [5.0],
        },
        index=pd.Index(["MAPK14;Y182;"], name="site_id"),
    )
    site_metadata = _site_metadata(phospho.index)
    sample_metadata = _comparison_sample_metadata(phospho.columns)

    preprocessed = DatasetPreprocessor().run(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=sample_metadata,
        total=None,
        plan=PreprocessingPlan.from_config(
            DatasetPreprocessingConfig(
                comparisons=DatasetComparisonBuildingConfig(
                    policy="sample_metadata_pairs",
                    pairs=(("group4", "group1"),),
                )
            )
        ),
    )

    assert preprocessed.comparisons is not None
    expected = pd.DataFrame(
        {"p_group4_group1": [-3.0]},
        index=phospho.index.copy(),
    )
    pdt.assert_frame_equal(preprocessed.comparisons, expected)


def test_dataset_preprocessor_rejects_comparison_building_without_sample_metadata() -> (
    None
):
    with pytest.raises(
        PhosPyInputError,
        match="policy='sample_metadata_pairs' requires sample_metadata input data",
    ):
        DatasetPreprocessor().run(
            phospho=_phospho().iloc[:, :2].copy(deep=True),
            site_metadata=_site_metadata().iloc[:3].copy(deep=True),
            sample_metadata=None,
            total=None,
            plan=PreprocessingPlan.from_config(
                DatasetPreprocessingConfig(
                    comparisons=DatasetComparisonBuildingConfig(
                        policy="sample_metadata_pairs"
                    )
                )
            ),
        )


def test_dataset_preprocessor_rejects_comparison_building_without_group_column() -> (
    None
):
    phospho = _phospho().iloc[:, :2].copy(deep=True)
    site_metadata = _site_metadata()
    sample_metadata = pd.DataFrame({"batch": [1, 2]}, index=phospho.columns)

    with pytest.raises(
        PhosPyInputError,
        match="requires sample_metadata column 'comparison_group'",
    ):
        DatasetPreprocessor().run(
            phospho=phospho,
            site_metadata=site_metadata,
            sample_metadata=sample_metadata,
            total=None,
            plan=PreprocessingPlan.from_config(
                DatasetPreprocessingConfig(
                    comparisons=DatasetComparisonBuildingConfig(
                        policy="sample_metadata_pairs"
                    )
                )
            ),
        )


def test_dataset_preprocessor_comparison_building_matches_reference_pairwise_expectation() -> (
    None
):
    phospho = pd.DataFrame(
        {
            "sample_a": [7.0],
            "sample_b": [4.0],
        },
        index=pd.Index(["PRKACA;S339;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["PRKACA"],
            "site": ["S339"],
            "site_sequence": ["AAAAAA"],
            "localisation_confidence": [0.95],
        },
        index=phospho.index.copy(),
    )
    sample_metadata = pd.DataFrame(
        {"comparison_group": ["sample_a", "sample_b"]},
        index=phospho.columns.copy(),
    )
    expected_fixture = pd.read_csv(
        ROOT
        / "tests"
        / "fixtures"
        / "rewrite_parity"
        / "comparison_building"
        / "reference_pairwise_expected.csv"
    ).set_index("site_id")

    preprocessed = DatasetPreprocessor().run(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=sample_metadata,
        total=None,
        plan=PreprocessingPlan.from_config(
            DatasetPreprocessingConfig(
                comparisons=DatasetComparisonBuildingConfig(
                    policy="sample_metadata_pairs",
                    pairs=(("sample_a", "sample_b"),),
                )
            )
        ),
    )

    assert preprocessed.comparisons is not None
    expected = expected_fixture.astype(float).copy()
    expected.index = phospho.index.copy()
    pdt.assert_frame_equal(preprocessed.comparisons, expected)
    assert preprocessed.comparison_group_stats is not None
    assert preprocessed.comparison_pair_stats is not None
    assert {"n", "mean", "sd", "sem"}.issubset(
        set(preprocessed.comparison_group_stats.columns)
    )
    assert {
        "left_n",
        "right_n",
        "left_mean",
        "right_mean",
        "left_sd",
        "right_sd",
        "left_sem",
        "right_sem",
        "effect_size",
    }.issubset(set(preprocessed.comparison_pair_stats.columns))
    pair_stats = preprocessed.comparison_pair_stats.copy(deep=True)
    assert pair_stats.shape[0] == preprocessed.comparisons.shape[0]
    paired_with_matrix = pair_stats.merge(
        preprocessed.comparisons.reset_index()
        .rename(columns={"index": "site_id"})
        .rename(columns={"p_sample_a_sample_b": "expected_effect_size"}),
        how="inner",
        on="site_id",
    )
    assert paired_with_matrix.shape[0] == preprocessed.comparisons.shape[0]
    assert (
        paired_with_matrix.loc[:, "effect_size"]
        == paired_with_matrix.loc[:, "expected_effect_size"]
    ).all()


def test_dataset_preprocessor_comparison_stats_follow_pandas_single_sample_convention() -> (
    None
):
    phospho = pd.DataFrame(
        {
            "sample_1": [8.0, 2.0],
            "sample_2": [10.0, 4.0],
            "sample_3": [5.0, 1.0],
        },
        index=pd.Index(["MAPK14;Y182;", "AKT1;T308;"], name="site_id"),
    )
    site_metadata = _site_metadata(phospho.index)
    sample_metadata = pd.DataFrame(
        {"comparison_group": ["group1", "group1", "group2"]},
        index=phospho.columns.copy(),
    )

    preprocessed = DatasetPreprocessor().run(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=sample_metadata,
        total=None,
        plan=PreprocessingPlan.from_config(
            DatasetPreprocessingConfig(
                comparisons=DatasetComparisonBuildingConfig(
                    policy="sample_metadata_pairs"
                )
            )
        ),
    )

    assert preprocessed.comparisons is not None
    expected = pd.DataFrame(
        {"p_group1_group2": [4.0, 2.0]},
        index=phospho.index.copy(),
    )
    pdt.assert_frame_equal(preprocessed.comparisons, expected)
    assert preprocessed.comparison_group_stats is not None
    assert preprocessed.comparison_pair_stats is not None

    group_stats = preprocessed.comparison_group_stats
    group2_rows = group_stats.loc[group_stats["group"] == "group2"].sort_values(
        "site_id"
    )
    assert group2_rows["n"].tolist() == [1, 1]
    assert group2_rows["mean"].tolist() == [1.0, 5.0]
    assert group2_rows["sd"].isna().all()
    assert group2_rows["sem"].isna().all()

    pair_stats = preprocessed.comparison_pair_stats.sort_values("site_id")
    assert pair_stats["comparison"].tolist() == ["p_group1_group2", "p_group1_group2"]
    assert pair_stats["left_n"].tolist() == [2, 2]
    assert pair_stats["right_n"].tolist() == [1, 1]
    assert pair_stats["left_mean"].tolist() == [3.0, 9.0]
    assert pair_stats["right_mean"].tolist() == [1.0, 5.0]
    assert pair_stats["right_sd"].isna().all()
    assert pair_stats["right_sem"].isna().all()
    assert pair_stats["effect_size"].tolist() == [2.0, 4.0]


def test_executor_delegates_preprocessing_to_internal_subsystem() -> None:
    calls: list[str] = []
    phospho = _phospho().iloc[:2].copy(deep=True)
    encoded_site_keys = [
        encode_site_key(
            build_protein_scoped_site_key(
                organism="rat",
                protein_namespace="protein_id",
                protein_identifier=protein_id,
                residue=residue,
                position=position,
                isoform_id=None,
                field_name="test.site_key",
                error_type=ValueError,
            )
        )
        for protein_id, residue, position in (("P1", "Y", 182), ("P2", "S", 9))
    ]
    phospho.index = pd.Index(encoded_site_keys, name="site_key")
    site_metadata = _site_metadata().iloc[:2].copy(deep=True)
    site_metadata.index = phospho.index.copy()
    site_metadata.loc[:, "display_id"] = phospho.index.tolist()
    site_metadata.loc[:, "site_key"] = phospho.index.tolist()
    sample_metadata = _sample_metadata(phospho.columns)
    total = _total(phospho.columns)

    interpreted = InterpretedDatasetBuildRequest(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=sample_metadata,
        total=total,
        organism=Organism.RAT,
        declared_input_intensity_scale_kind=IntensityScaleKind.LINEAR,
        preprocessing_plan=PreprocessingPlan.default(),
    )

    preprocessed_tables = PreprocessedDatasetBuildTables(
        phospho=phospho + 10.0,
        site_metadata=site_metadata.assign(processed=True),
        sample_metadata=sample_metadata.assign(processed=True),
        total=total + 5.0,
    )

    class PreprocessorSpy:
        def run(
            self,
            *,
            phospho: pd.DataFrame,
            site_metadata: pd.DataFrame,
            sample_metadata: pd.DataFrame | None,
            total: pd.DataFrame | None,
            plan: PreprocessingPlan,
        ) -> PreprocessedDatasetBuildTables:
            calls.append("preprocessor")
            assert phospho is interpreted.phospho
            assert site_metadata is interpreted.site_metadata
            assert sample_metadata is interpreted.sample_metadata
            assert total is interpreted.total
            assert plan is interpreted.preprocessing_plan
            return preprocessed_tables

    class ResolverSpy:
        def run(
            self,
            *,
            phospho: pd.DataFrame,
            total: pd.DataFrame | None,
            expected_scale_kind: object | None = None,
        ) -> ResolvedIntensityScale:
            calls.append("resolver")
            assert phospho is preprocessed_tables.phospho
            assert total is preprocessed_tables.total
            assert expected_scale_kind is not None
            return ResolvedIntensityScale(
                phospho=phospho,
                total=total,
                intensity_scale_state=supported_linear_intensity_scale_state(
                    has_total_matrix=True
                ),
            )

    built = DatasetBuildExecutor(
        preprocessor=PreprocessorSpy(),
        intensity_scale_resolver=ResolverSpy(),
    ).run(interpreted)

    pdt.assert_frame_equal(built.phospho, preprocessed_tables.phospho)
    pdt.assert_frame_equal(built.site_metadata, preprocessed_tables.site_metadata)
    assert built._sample_metadata is preprocessed_tables.sample_metadata
    pdt.assert_frame_equal(built.total, preprocessed_tables.total)
    assert calls == ["preprocessor", "resolver"]


def test_dataset_interpreter_does_not_apply_preprocessing_science() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, float("nan")],
            "sample_b": [2.0, float("nan")],
        },
        index=["MAPK14;Y182;", "GSK3B;S9;"],
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "GSK3B"],
            "site": ["Y182", "S9"],
            "site_sequence": [
                "LDFGLARHTDDEMTGYVATRWYRAPEIMLNW",
                "AAAAAAAAAAAAAAASAAAAAAAAAAAAAAA",
            ],
            "localisation_confidence": [0.95, 0.9],
        },
        index=phospho.index,
    )
    request = DatasetBuildRequest(
        phospho=phospho,
        site_metadata=site_metadata,
        preprocessing_config=DatasetPreprocessingConfig(
            missing_data=DatasetMissingDataConfig(
                policy="impute_row_median",
                min_observed_values=2,
            ),
        ),
    )

    interpreted = DatasetBuildRequestInterpreter().run(request)

    assert interpreted.phospho.isna().to_numpy().sum() == 2
    assert interpreted.phospho.index.tolist() == ["MAPK14;Y182;", "GSK3B;S9;"]


def test_dataset_builder_request_quantitative_meaning_propagates_to_provenance() -> (
    None
):
    phospho = pd.DataFrame(
        {"sample_a": [0.25, -0.5], "sample_b": [1.0, 0.0]},
        index=pd.Index(["MAPK14;Y182;", "AKT1;T308;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "AKT1"],
            "site": ["Y182", "T308"],
            "site_sequence": ["SEQ_A", "SEQ_R"],
            "localisation_confidence": [0.95, 0.96],
        },
        index=phospho.index.copy(),
    )

    built = DatasetBuildExecutor().run(
        DatasetBuildRequestInterpreter().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=site_metadata,
                input_intensity_scale="log2",
                quantitative_meaning=QuantitativeMeaning.CONTRAST_LOG2_FOLD_CHANGE.value,
            )
        )
    )

    expected = QuantitativeMeaning.CONTRAST_LOG2_FOLD_CHANGE.value
    assert built.intensity_scale_state.quantity is not None
    assert built.intensity_scale_state.quantity.value == expected
    assert built.provenance is not None
    assert built.provenance.workflow_parameters["quantitative_meaning"] == expected


def test_dataset_interpreter_defers_reference_site_sequence_fill_when_fasta_is_configured() -> (
    None
):
    calls: list[dict[str, object]] = []

    class SiteSequenceDeriverSpy:
        def run(
            self,
            site_metadata: pd.DataFrame,
            *,
            organism: Organism | None,
            allow_partial: bool = False,
            derive_missing_from_reference: bool = True,
        ) -> pd.DataFrame:
            calls.append(
                {
                    "organism": organism,
                    "allow_partial": allow_partial,
                    "derive_missing_from_reference": derive_missing_from_reference,
                }
            )
            return site_metadata

    request = DatasetBuildRequest(
        phospho=pd.DataFrame(
            {"sample_a": [1.0]},
            index=pd.Index(["MAPK14;Y182;"], name="site_id"),
        ),
        site_metadata=pd.DataFrame(
            {
                "gene_symbol": ["MAPK14"],
                "site": ["Y182"],
                "protein_accession": ["P1"],
                "site_sequence": [pd.NA],
            },
            index=pd.Index(["MAPK14;Y182;"], name="site_id"),
        ),
        organism=Organism.RAT,
        input_intensity_scale="linear",
        preprocessing_config=DatasetPreprocessingConfig(
            site_sequence_resolution=DatasetSiteSequenceResolutionConfig(
                fasta_path="local.fasta",
                flank_size=2,
            ),
        ),
    )

    interpreted = DatasetBuildRequestInterpreter(
        site_sequence_deriver=SiteSequenceDeriverSpy()
    ).run(request)

    assert len(calls) == 1
    assert calls[0]["organism"] == Organism.RAT
    assert calls[0]["allow_partial"] is True
    assert calls[0]["derive_missing_from_reference"] is False
    assert pd.isna(interpreted.site_metadata.loc["MAPK14;Y182;", "site_sequence"])


def test_dataset_interpreter_keeps_reference_site_sequence_fill_enabled_without_fasta() -> (
    None
):
    calls: list[dict[str, object]] = []

    class SiteSequenceDeriverSpy:
        def run(
            self,
            site_metadata: pd.DataFrame,
            *,
            organism: Organism | None,
            allow_partial: bool = False,
            derive_missing_from_reference: bool = True,
        ) -> pd.DataFrame:
            calls.append(
                {
                    "organism": organism,
                    "allow_partial": allow_partial,
                    "derive_missing_from_reference": derive_missing_from_reference,
                }
            )
            return site_metadata

    request = DatasetBuildRequest(
        phospho=pd.DataFrame(
            {"sample_a": [1.0]},
            index=pd.Index(["MAPK14;Y182;"], name="site_id"),
        ),
        site_metadata=pd.DataFrame(
            {
                "gene_symbol": ["MAPK14"],
                "site": ["Y182"],
                "protein_accession": ["P1"],
            },
            index=pd.Index(["MAPK14;Y182;"], name="site_id"),
        ),
        organism=Organism.RAT,
        input_intensity_scale="linear",
    )

    DatasetBuildRequestInterpreter(site_sequence_deriver=SiteSequenceDeriverSpy()).run(
        request
    )

    assert len(calls) == 1
    assert calls[0]["organism"] == Organism.RAT
    assert calls[0]["allow_partial"] is False
    assert calls[0]["derive_missing_from_reference"] is True


def test_preprocessing_plan_defaults_keep_identity_transform_and_no_normalisation() -> (
    None
):
    plan = PreprocessingPlan.from_config(DatasetPreprocessingConfig())
    assert plan.intensity_transform_policy == "identity"
    assert plan.intensity_transform_pseudocount == pytest.approx(1.0)
    assert plan.normalisation_policy == "none"
    assert "intensity_transform" not in plan.stage_order
    assert "normalisation" not in plan.stage_order


def test_default_preprocessing_report_records_noop_normalisation_method() -> None:
    phospho = _phospho()
    site_metadata = _site_metadata()

    preprocessed = DatasetPreprocessor().run(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan.from_config(DatasetPreprocessingConfig()),
    )

    normalisation_operation = preprocessed.preprocessing_operations.loc[
        preprocessed.preprocessing_operations.loc[:, "stage"] == "normalisation"
    ]
    assert normalisation_operation.shape[0] == 1
    assert normalisation_operation.iloc[0]["operation"] == "none"
    assert normalisation_operation.iloc[0]["parameters"] == {"applied": False}
    assert normalisation_operation.iloc[0]["notes"] == (
        "stage not scheduled in preprocessing plan"
    )


def test_dataset_preprocessor_applies_log2_intensity_transform_policy() -> None:
    phospho = _phospho()
    site_metadata = _site_metadata()

    preprocessed = DatasetPreprocessor().run(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan.from_config(
            DatasetPreprocessingConfig(
                intensity_transform=DatasetIntensityTransformConfig(
                    policy="log2",
                    pseudocount=1.0,
                )
            )
        ),
    )

    expected = pd.DataFrame(
        {
            "sample_a": [1.0, 1.584962500721156, 2.0],
            "sample_b": [1.584962500721156, 1.0, 2.321928094887362],
            "sample_c": [2.0, 1.3219280948873624, 2.584962500721156],
        },
        index=phospho.index.copy(),
    )
    pdt.assert_frame_equal(preprocessed.phospho, expected)
    assert preprocessed.phospho.index.equals(phospho.index)
    assert preprocessed.phospho.columns.equals(phospho.columns)
    assert preprocessed.site_metadata is site_metadata


def test_dataset_preprocessor_rejects_log2_transform_when_values_are_invalid() -> None:
    phospho = _phospho()
    phospho.loc["MAPK14;Y182;", "sample_a"] = -1.0

    with pytest.raises(
        PhosPyInputError,
        match="requires all non-missing values plus pseudocount to be greater than 0",
    ):
        DatasetPreprocessor().run(
            phospho=phospho,
            site_metadata=_site_metadata(),
            sample_metadata=None,
            total=None,
            plan=PreprocessingPlan.from_config(
                DatasetPreprocessingConfig(
                    intensity_transform=DatasetIntensityTransformConfig(
                        policy="log2",
                        pseudocount=1.0,
                    )
                )
            ),
        )


def test_dataset_preprocessor_applies_median_center_normalisation_policy() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 2.0, 3.0],
            "sample_b": [4.0, 2.0, 0.0],
        },
        index=_phospho().index.copy(),
    )

    preprocessed = DatasetPreprocessor().run(
        phospho=phospho,
        site_metadata=_site_metadata(phospho.index),
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan.from_config(
            DatasetPreprocessingConfig(
                normalisation=DatasetNormalisationConfig(policy="median_center")
            )
        ),
    )

    expected = pd.DataFrame(
        {
            "sample_a": [-1.0, 0.0, 1.0],
            "sample_b": [2.0, 0.0, -2.0],
        },
        index=phospho.index.copy(),
    )
    pdt.assert_frame_equal(preprocessed.phospho, expected)
    assert preprocessed.phospho.index.equals(phospho.index)
    assert preprocessed.phospho.columns.equals(phospho.columns)


def test_dataset_preprocessor_quantile_normalisation_equalises_column_distributions() -> (
    None
):
    phospho = pd.DataFrame(
        {
            "sample_a": [5.0, 2.0, 3.0, 4.0],
            "sample_b": [4.0, 1.0, 2.5, 2.0],
            "sample_c": [8.0, 6.0, 7.0, 9.0],
        },
        index=pd.Index(["A;S1;", "B;S1;", "C;S1;", "D;S1;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["A", "B", "C", "D"],
            "site": ["S1", "S1", "S1", "S1"],
            "site_sequence": ["SEQ_A", "SEQ_R", "SEQ_C", "SEQ_D"],
            "localisation_confidence": [0.95, 0.94, 0.93, 0.92],
        },
        index=phospho.index.copy(),
    )

    preprocessed = DatasetPreprocessor().run(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan.from_config(
            DatasetPreprocessingConfig(
                normalisation=DatasetNormalisationConfig(policy="quantile")
            )
        ),
    )

    sorted_a = sorted(preprocessed.phospho.loc[:, "sample_a"].tolist())
    sorted_b = sorted(preprocessed.phospho.loc[:, "sample_b"].tolist())
    sorted_c = sorted(preprocessed.phospho.loc[:, "sample_c"].tolist())
    assert sorted_a == pytest.approx(sorted_b)
    assert sorted_b == pytest.approx(sorted_c)
    assert preprocessed.phospho.index.equals(phospho.index)
    assert preprocessed.phospho.columns.equals(phospho.columns)
    assert preprocessed.site_metadata is site_metadata


def test_dataset_preprocessor_quantile_normalisation_is_deterministic_with_ties_and_missing() -> (
    None
):
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 1.0, 3.0, float("nan")],
            "sample_b": [4.0, 2.0, 2.0, 1.0],
            "sample_c": [7.0, 8.0, 9.0, 10.0],
        },
        index=pd.Index(["A;S1;", "B;S1;", "C;S1;", "D;S1;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["A", "B", "C", "D"],
            "site": ["S1", "S1", "S1", "S1"],
            "site_sequence": ["SEQ_A", "SEQ_R", "SEQ_C", "SEQ_D"],
            "localisation_confidence": [0.95, 0.94, 0.93, 0.92],
        },
        index=phospho.index.copy(),
    )
    config = DatasetPreprocessingConfig(
        normalisation=DatasetNormalisationConfig(policy="quantile")
    )

    first = DatasetPreprocessor().run(
        phospho=phospho.copy(deep=True),
        site_metadata=site_metadata.copy(deep=True),
        sample_metadata=None,
        total=None,
        plan=_plan_without_missing_stage(config),
    )
    second = DatasetPreprocessor().run(
        phospho=phospho.copy(deep=True),
        site_metadata=site_metadata.copy(deep=True),
        sample_metadata=None,
        total=None,
        plan=_plan_without_missing_stage(config),
    )

    pdt.assert_frame_equal(first.phospho, second.phospho)
    assert first.phospho.loc["A;S1;", "sample_a"] == pytest.approx(
        first.phospho.loc["B;S1;", "sample_a"]
    )
    assert pd.isna(first.phospho.loc["D;S1;", "sample_a"])


def test_dataset_preprocessor_site_matrix_min_observed_filter_can_keep_all_rows() -> (
    None
):
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 10.0],
            "sample_b": [2.0, float("nan")],
            "sample_c": [float("nan"), 12.0],
        },
        index=pd.Index(["row_a", "row_b"], name="input_row"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "AKT1"],
            "site": ["Y182", "T308"],
            "site_sequence": ["SEQ_A", "SEQ_R"],
            "localisation_confidence": [0.95, 0.96],
        },
        index=phospho.index.copy(),
    )

    preprocessed = DatasetPreprocessor().run(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=None,
        total=None,
        plan=_plan_without_missing_stage(
            DatasetPreprocessingConfig(
                site_matrix=_internal_site_matrix_config(
                    policy="build_from_metadata",
                    missing_data_policy="require_min_observed_values",
                    minimum_observed_values=2,
                    duplicate_site_policy="aggregate_mean",
                )
            )
        ),
    )

    assert preprocessed.phospho.index.tolist() == ["AKT1;T308;", "MAPK14;Y182;"]
    dropped = preprocessed.row_audit.loc[
        (preprocessed.row_audit.loc[:, "stage"] == "site_matrix")
        & (preprocessed.row_audit.loc[:, "action"] == "dropped")
    ]
    assert dropped.empty


def test_dataset_preprocessor_site_matrix_min_observed_filter_can_remove_all_rows() -> (
    None
):
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, float("nan")],
            "sample_b": [float("nan"), 5.0],
            "sample_c": [float("nan"), float("nan")],
        },
        index=pd.Index(["row_a", "row_b"], name="input_row"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "AKT1"],
            "site": ["Y182", "T308"],
            "site_sequence": ["SEQ_A", "SEQ_R"],
            "localisation_confidence": [0.95, 0.96],
        },
        index=phospho.index.copy(),
    )

    with pytest.raises(
        PhosPyInputError,
        match="site-matrix construction produced no retained rows after filtering",
    ):
        DatasetPreprocessor().run(
            phospho=phospho,
            site_metadata=site_metadata,
            sample_metadata=None,
            total=None,
            plan=_plan_without_missing_stage(
                DatasetPreprocessingConfig(
                    site_matrix=_internal_site_matrix_config(
                        policy="build_from_metadata",
                        missing_data_policy="require_min_observed_values",
                        minimum_observed_values=2,
                        duplicate_site_policy="aggregate_mean",
                    )
                )
            ),
        )


def test_dataset_preprocessor_trace_tracks_stagewise_transform_progression() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, float("nan"), 4.0],
            "sample_b": [2.0, 3.0, 5.0],
            "sample_c": [3.0, 5.0, 6.0],
        },
        index=pd.Index(["SITE_A", "SITE_B", "SITE_C"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "AKT1", "GSK3B"],
            "site": ["Y182", "T308", "S9"],
            "site_sequence": ["SEQ_A", "SEQ_R", "SEQ_C"],
            "localisation_confidence": [0.95, 0.96, 0.9],
        },
        index=phospho.index.copy(),
    )

    preprocessed = DatasetPreprocessor().run(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan.from_config(
            DatasetPreprocessingConfig(
                missing_data=DatasetMissingDataConfig(
                    policy="impute_row_median",
                    min_observed_values=2,
                ),
                intensity_transform=DatasetIntensityTransformConfig(
                    policy="log2",
                    pseudocount=1.0,
                ),
                normalisation=DatasetNormalisationConfig(policy="median_center"),
            )
        ),
    )

    stage_names = [stage.stage for stage in preprocessed.preprocessing_trace]
    assert stage_names == [
        "localisation_confidence",
        "missing_data",
        "intensity_transform",
        "normalisation",
    ]
    for stage in preprocessed.preprocessing_trace:
        assert isinstance(stage.phospho_input_hash, str)
        assert isinstance(stage.phospho_output_hash, str)
        assert len(stage.phospho_input_hash) == 64
        assert len(stage.phospho_output_hash) == 64
        assert stage.phospho_input_hash != stage.phospho_output_hash


@pytest.mark.parametrize(
    ("intensity_transform", "normalisation", "expected"),
    [
        (
            DatasetIntensityTransformConfig(policy="log2", pseudocount=1.0),
            DatasetNormalisationConfig(policy="none"),
            "intensity_transform.policy='log2' requires numeric phospho columns",
        ),
        (
            DatasetIntensityTransformConfig(policy="identity", pseudocount=1.0),
            DatasetNormalisationConfig(policy="median_center"),
            "normalisation.policy='median_center' requires numeric phospho columns",
        ),
    ],
)
def test_dataset_preprocessor_rejects_non_numeric_phospho_columns_for_new_methods(
    intensity_transform: DatasetIntensityTransformConfig,
    normalisation: DatasetNormalisationConfig,
    expected: str,
) -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 2.0],
            "metadata_note": ["a", "b"],
        },
        index=pd.Index(["MAPK14;Y182;", "AKT1;T308;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "AKT1"],
            "site": ["Y182", "T308"],
            "site_sequence": ["SEQ_A", "SEQ_R"],
            "localisation_confidence": [0.95, 0.96],
        },
        index=phospho.index.copy(),
    )

    with pytest.raises(PhosPyInputError, match=expected):
        DatasetPreprocessor().run(
            phospho=phospho,
            site_metadata=site_metadata,
            sample_metadata=None,
            total=None,
            plan=PreprocessingPlan.from_config(
                DatasetPreprocessingConfig(
                    intensity_transform=intensity_transform,
                    normalisation=normalisation,
                )
            ),
        )
