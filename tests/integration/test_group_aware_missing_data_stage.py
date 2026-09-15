from __future__ import annotations

import json
from copy import deepcopy

import numpy as np
import pandas as pd
import pandas.testing as pdt
import pytest

from phospy import AnalysisReadyDatasetBuilder, KinaseWorkflow
from phospy.advanced import (
    DatasetMissingDataConfig,
    KinasePredictionConfig,
    KinaseScoringConfig,
    ReferenceContextCompatibilityPolicy,
)
from phospy.api import (
    DatasetBuildRequest,
    DatasetPreprocessingConfig,
    KinaseWorkflowRequest,
    Organism,
    ReferenceBundle,
)
from phospy.errors.input import PhosPyInputError
from phospy.errors.validation import DatasetValidationError
from phospy.io.bundles.kinase import (
    KinaseWorkflowConfigSnapshot,
    load_kinase_workflow_bundle,
    save_kinase_workflow_bundle,
)
from phospy.provenance.models import DeterminismKind
from phospy.science.datasets.builders.preprocessing import (
    DatasetPreprocessor,
    build_dataset_processing_state,
)
from phospy.science.datasets.preprocessing.models import (
    PreprocessingPlan,
    PreprocessingState,
)
from phospy.science.datasets.preprocessing.stages.missing_data import (
    knn as knn_module,
)
from phospy.science.datasets.preprocessing.stages.missing_data import (
    stage as missing_data_stage_module,
)
from phospy.science.datasets.preprocessing.stages.missing_data.group_aware_routing import (
    route_group_aware_missingness,
)
from phospy.science.datasets.preprocessing.stages.missing_data.models import (
    GroupAwarePolicyOutcome,
    KnnPolicyOutcome,
    MinProbPolicyOutcome,
)
from phospy.science.datasets.preprocessing.stages.missing_data.stage import (
    MissingDataStage,
    execute_group_aware_imputation,
)
from phospy.science.datasets.processing_state import MissingDataDiagnosticsV2
from phospy.science.transformations.models import IntensityScaleKind
from tests.support.analysis_ready_dataset_factories import (
    trusted_analysis_ready_dataset_from_tables,
)
from tests.support.intensity_scale_states import supported_log2_intensity_scale_state
from tests.support.site_keys import protein_site_key_index


def _inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    display_ids = ["GENE1;S1;", "GENE2;S2;", "GENE3;S3;", "GENE4;S4;"]
    index = protein_site_key_index(
        protein_identifiers=["P1", "P2", "P3", "P4"],
        sites=["S1", "S2", "S3", "S4"],
        protein_namespace="protein_id",
        organism="rat",
    )
    columns = pd.Index(["A1", "A2", "A3", "B1", "B2", "B3"], name="sample")
    phospho = pd.DataFrame(
        [
            [1.0, 2.0, np.nan, 5.0, 6.0, 7.0],
            [np.nan, np.nan, np.nan, 8.0, 9.0, 10.0],
            [2.0, 3.0, 4.0, 6.0, 7.0, 8.0],
            [1.0, np.nan, np.nan, 4.0, 5.0, 6.0],
        ],
        index=index,
        columns=columns,
    )
    site_metadata = pd.DataFrame(
        {
            "site_key": index.astype(str).tolist(),
            "display_id": display_ids,
            "gene_symbol": ["GENE1", "GENE2", "GENE3", "GENE4"],
            "site": ["S1", "S2", "S3", "S4"],
            "site_sequence": [
                "AAAAAAAAAAAAAAASAAAAAAAAAAAAAAA",
                "CCCCCCCCCCCCCCCSCCCCCCCCCCCCCCC",
                "DDDDDDDDDDDDDDDSDDDDDDDDDDDDDDD",
                "EEEEEEEEEEEEEEESEEEEEEEEEEEEEEE",
            ],
            "protein_id": ["P1", "P2", "P3", "P4"],
            "localisation_confidence": [0.9, 0.9, 0.9, 0.9],
        },
        index=index.copy(),
    )
    sample_metadata = pd.DataFrame(
        {"condition": ["A", "A", "A", "B", "B", "B"]},
        index=columns.copy(),
    )
    return phospho, site_metadata, sample_metadata


def _group_aware_config(seed: int) -> DatasetPreprocessingConfig:
    return DatasetPreprocessingConfig(
        missing_data=DatasetMissingDataConfig(
            policy="impute_group_aware",
            group_column="condition",
            min_partial_observed_fraction=0.5,
            min_reference_observed_fraction=0.75,
            q=0.01,
            width=0.3,
            seed=seed,
            k=2,
            distance="nan_euclidean",
        )
    )


def _build(seed: int, *, sample_metadata: pd.DataFrame | None = None):
    phospho, site_metadata, default_sample_metadata = _inputs()
    return AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            sample_metadata=(
                default_sample_metadata if sample_metadata is None else sample_metadata
            ),
            organism=Organism.RAT,
            preprocessing_config=_group_aware_config(seed),
            input_intensity_scale="log2",
            allow_suspicious_declared_input_intensity_scale=True,
        )
    )


def _missing_data_stage(dataset):
    assert dataset.provenance is not None
    return next(
        stage
        for stage in dataset.provenance.preprocessing_stages
        if stage.stage == "missing_data"
    )


def test_group_aware_stage_produces_complete_mixed_mechanism_matrix(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    phospho, _, _ = _inputs()
    mechanism_inputs: list[pd.DataFrame] = []
    real_knn = missing_data_stage_module.impute_knn_targets
    real_minprob = missing_data_stage_module.impute_minprob_targets

    def _capture_knn_input(phospho: pd.DataFrame, **kwargs: object):
        mechanism_inputs.append(phospho.copy(deep=True))
        return real_knn(phospho, **kwargs)

    def _capture_minprob_input(phospho: pd.DataFrame, **kwargs: object):
        mechanism_inputs.append(phospho.copy(deep=True))
        return real_minprob(phospho, **kwargs)

    monkeypatch.setattr(
        missing_data_stage_module, "impute_knn_targets", _capture_knn_input
    )
    monkeypatch.setattr(
        missing_data_stage_module,
        "impute_minprob_targets",
        _capture_minprob_input,
    )
    dataset = _build(123)
    knn_row, minprob_row, _, dropped_row = phospho.index
    original_retained = phospho.drop(index=dropped_row)

    assert dropped_row not in dataset.phospho.index
    assert len(mechanism_inputs) == 2
    expected_knn_target_mask_hash = (
        "b34fa2ba8391e3915e7308786d3c89aafb4227e4f45f851f028feb0f2d1e52d7"
    )
    expected_minprob_target_mask_hash = (
        "18efa8fe46ab30e8b284162af38a84f51c220d6a9bd3842b588650abbd1a0a56"
    )
    expected_knn_imputation_mask_hash = (
        "4c9086793071a4613c5283f7831aa707c8c74bab7c8cfe95f18e86251373f5f1"
    )
    expected_minprob_imputation_mask_hash = (
        "5d7770b23d948dbb6d6e136f98fa658badb29bdf6ddc4551a9986f28320f9faa"
    )
    expected_imputation_mask_hash = (
        "ed61ed60a299677b4a23af967fe76999aa14c390a8d7a9b855ff89eb4f2f32e8"
    )
    pdt.assert_frame_equal(mechanism_inputs[0], original_retained)
    pdt.assert_frame_equal(mechanism_inputs[1], original_retained)
    assert not dataset.phospho.isna().to_numpy().any()
    assert dataset.phospho.at[knn_row, "A3"] == 4.0
    assert dataset.phospho.loc[minprob_row, ["A1", "A2", "A3"]].notna().all()

    retained_input = phospho.loc[dataset.phospho.index]
    observed = retained_input.notna()
    np.testing.assert_array_equal(
        dataset.phospho.to_numpy()[observed.to_numpy()],
        retained_input.to_numpy()[observed.to_numpy()],
    )
    observation_mask = dataset.imputation_observed_mask_dataframe()
    assert observation_mask is not None
    pdt.assert_frame_equal(
        observation_mask,
        retained_input.notna(),
        check_names=False,
    )

    assert dataset.preprocessing_report is not None
    row_audit = dataset.preprocessing_report.row_audit
    assert set(row_audit["source_row_id"]) == {
        knn_row,
        minprob_row,
        dropped_row,
    }
    assert set(row_audit["action"]) == {"dropped", "imputed"}

    stage = _missing_data_stage(dataset)
    assert stage.determinism is DeterminismKind.SEEDED_STOCHASTIC
    assert stage.random_seed == 123
    assert {item.name for item in stage.consumed_input_tables} == {
        "dataset.phospho",
        "dataset.site_metadata",
        "dataset.sample_metadata",
    }
    diagnostics = stage.diagnostics or {}
    parameters = diagnostics.get("method_parameters", {})
    assert parameters["knn_target_cell_count"] == 1
    assert parameters["minprob_target_cell_count"] == 3
    assert parameters["mechanism_input"] == "original_retained_matrix"
    typed_diagnostics = dataset.processing_state.missing_data.diagnostics
    assert isinstance(typed_diagnostics, MissingDataDiagnosticsV2)
    assert typed_diagnostics.group_aware.observed_group_sizes == {"A": 3, "B": 3}
    assert typed_diagnostics.group_aware.knn_imputed_cell_count == 1
    assert typed_diagnostics.group_aware.minprob_imputed_cell_count == 3
    assert typed_diagnostics.group_aware.minprob_left_censored_assumption is True
    assert [
        record.to_payload() for record in typed_diagnostics.group_aware.routed_rows
    ] == [
        {
            "row_id": str(knn_row),
            "knn_imputed_columns": ["A3"],
            "minprob_imputed_columns": [],
            "knn_group_labels": ["A"],
            "minprob_group_labels": [],
        },
        {
            "row_id": str(minprob_row),
            "knn_imputed_columns": [],
            "minprob_imputed_columns": ["A1", "A2", "A3"],
            "knn_group_labels": [],
            "minprob_group_labels": ["A"],
        },
    ]
    assert (
        typed_diagnostics.group_aware.knn_target_mask_hash
        == expected_knn_target_mask_hash
    )
    assert (
        typed_diagnostics.group_aware.minprob_target_mask_hash
        == expected_minprob_target_mask_hash
    )
    assert (
        typed_diagnostics.group_aware.knn_imputation_mask_hash
        == expected_knn_imputation_mask_hash
    )
    assert (
        typed_diagnostics.group_aware.minprob_imputation_mask_hash
        == expected_minprob_imputation_mask_hash
    )
    assert typed_diagnostics.imputation_mask_hash == expected_imputation_mask_hash

    operations = dataset.preprocessing_report.operations
    missing_operation = operations.loc[operations["stage"] == "missing_data"].iloc[0]
    operation_parameters = missing_operation["parameters"]
    assert operation_parameters == {
        "missing_data_policy": "impute_group_aware",
        "missing_data_min_observed_values": None,
        "missing_data_q": 0.01,
        "missing_data_width": 0.3,
        "missing_data_seed": 123,
        "missing_data_k": 2,
        "missing_data_distance": "nan_euclidean",
        "missing_data_max_missing_fraction_per_row": None,
        "missing_data_no_overlap_policy": "error",
        "missing_data_group_column": "condition",
        "missing_data_min_partial_observed_fraction": 0.5,
        "missing_data_min_reference_observed_fraction": 0.75,
        "missing_data_input_scale": "log2",
        "missing_data_input_scale_source": "method_required",
        "missing_data_imputation_operation_order": "no_intensity_transform",
        "execution_summary": operation_parameters["execution_summary"],
    }
    assert (
        operation_parameters["execution_summary"]["imputation_scope"] == "group_aware"
    )


def test_group_aware_stage_preserves_canonical_identity_through_audit() -> None:
    columns = pd.Index([" a1 ", "a2 ", " b1", "b2"], name="sample")
    index = protein_site_key_index(
        protein_identifiers=["P1", "P2"],
        sites=["S1", "S2"],
        protein_namespace="protein_id",
        organism="rat",
    )
    phospho = pd.DataFrame(
        [[1.0, np.nan, np.nan, np.nan], [2.0, 3.0, 4.0, 5.0]],
        index=index,
        columns=columns,
    )
    site_metadata = pd.DataFrame({"site": ["S1", "S2"]}, index=index.copy())
    sample_metadata = pd.DataFrame(
        {"condition": ["B", "A", "B", "A"]},
        index=pd.Index(["b2", "a2", "b1", "a1"], name="sample"),
    )
    state = PreprocessingState(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=sample_metadata,
        total=None,
        plan=PreprocessingPlan(
            missing_data_policy="impute_group_aware",
            missing_data_group_column="condition",
            missing_data_min_partial_observed_fraction=0.5,
            missing_data_min_reference_observed_fraction=0.5,
            missing_data_q=0.01,
            missing_data_width=0.3,
            missing_data_seed=123,
            missing_data_k=1,
            missing_data_distance="nan_euclidean",
            stage_order=("missing_data",),
        ),
    )
    result = MissingDataStage().run(state)

    assert phospho.columns.equals(columns)
    assert result.state.phospho.columns.equals(columns)
    assert not result.state.phospho.isna().to_numpy().any()
    diagnostics = MissingDataDiagnosticsV2.from_mapping(
        result.diagnostics["diagnostics"], field_name="diagnostics"
    )
    assert diagnostics.affected_column_ids == ("a2 ", " b1", "b2")
    assert diagnostics.imputed_column_ids == ("a2 ", " b1", "b2")
    assert diagnostics.per_column_distribution_parameters is not None
    assert set(diagnostics.per_column_distribution_parameters) == set(columns)
    assert [record.to_payload() for record in diagnostics.group_aware.routed_rows] == [
        {
            "row_id": str(index[0]),
            "knn_imputed_columns": ["a2 "],
            "minprob_imputed_columns": [" b1", "b2"],
            "knn_group_labels": ["A"],
            "minprob_group_labels": ["B"],
        }
    ]
    method_parameters = diagnostics.method_parameters
    assert method_parameters["resolved_group_samples"] == {
        "A": [" a1 ", "a2 "],
        "B": [" b1", "b2"],
    }
    assert result.state.row_audit is not None
    audit = result.state.row_audit.iloc[0]["parameter_snapshot"]
    assert audit["knn_imputed_columns"] == ("a2 ",)
    assert audit["minprob_imputed_columns"] == (" b1", "b2")
    assert audit["knn_affected_groups"] == ["A"]
    assert audit["minprob_affected_groups"] == ["B"]
    assert audit["mechanisms"] == [
        "partial_observation_knn",
        "asymmetric_absence_minprob",
    ]
    assert audit["resolved_group_samples"] == {
        "A": [" a1 ", "a2 "],
        "B": [" b1", "b2"],
    }


def test_group_aware_original_labels_survive_final_dataset_binding() -> None:
    phospho, site_metadata, sample_metadata = _inputs()
    columns = pd.Index(
        [" A1 ", "A2 ", " A3", "B1 ", " B2", "B3 "],
        name="sample",
    )
    phospho.columns = columns
    sample_metadata.index = columns.copy()
    site_metadata = site_metadata.assign(
        organism="rat",
        protein_namespace="protein_id",
        protein_identifier=["P1", "P2", "P3", "P4"],
    )
    plan = PreprocessingPlan(
        missing_data_policy="impute_group_aware",
        missing_data_group_column="condition",
        missing_data_min_partial_observed_fraction=0.5,
        missing_data_min_reference_observed_fraction=0.75,
        missing_data_q=0.01,
        missing_data_width=0.3,
        missing_data_seed=123,
        missing_data_k=2,
        missing_data_distance="nan_euclidean",
        stage_order=("missing_data",),
    )

    preprocessed = DatasetPreprocessor().run(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=sample_metadata,
        total=None,
        plan=plan,
        initial_quantitative_scale_kind=IntensityScaleKind.LOG2,
    )
    processing_state = build_dataset_processing_state(
        plan=plan,
        intensity_scale_state=supported_log2_intensity_scale_state(
            has_total_matrix=False
        ),
        preprocessing_trace=preprocessed.preprocessing_trace,
        final_phospho=preprocessed.phospho,
        final_site_metadata=preprocessed.site_metadata,
        final_sample_metadata=preprocessed.sample_metadata,
    )

    assert phospho.columns.equals(columns)
    assert preprocessed.phospho.columns.equals(columns)
    assert not preprocessed.phospho.isna().to_numpy().any()
    diagnostics = processing_state.missing_data.diagnostics
    assert isinstance(diagnostics, MissingDataDiagnosticsV2)
    assert diagnostics.method_parameters["resolved_group_samples"] == {
        "A": [" A1 ", "A2 ", " A3"],
        "B": ["B1 ", " B2", "B3 "],
    }
    assert [record.to_payload() for record in diagnostics.group_aware.routed_rows] == [
        {
            "row_id": str(phospho.index[0]),
            "knn_imputed_columns": [" A3"],
            "minprob_imputed_columns": [],
            "knn_group_labels": ["A"],
            "minprob_group_labels": [],
        },
        {
            "row_id": str(phospho.index[1]),
            "knn_imputed_columns": [],
            "minprob_imputed_columns": [" A1 ", "A2 ", " A3"],
            "knn_group_labels": [],
            "minprob_group_labels": ["A"],
        },
    ]
    audit_by_row = {
        row["source_row_id"]: row["parameter_snapshot"]
        for _, row in preprocessed.row_audit.iterrows()
        if row["action"] == "imputed"
    }
    assert audit_by_row[phospho.index[0]]["knn_imputed_columns"] == (" A3",)
    assert audit_by_row[phospho.index[0]]["knn_affected_groups"] == ["A"]
    assert audit_by_row[phospho.index[1]]["minprob_imputed_columns"] == (
        " A1 ",
        "A2 ",
        " A3",
    )
    assert audit_by_row[phospho.index[1]]["minprob_affected_groups"] == ["A"]

    dataset = trusted_analysis_ready_dataset_from_tables(
        phospho=preprocessed.phospho,
        site_metadata=preprocessed.site_metadata,
        sample_metadata=preprocessed.sample_metadata,
        intensity_scale_state=processing_state.intensity_scale,
        processing_state=processing_state,
        imputation_observation_mask=preprocessed.imputation_observation_mask,
        organism=Organism.RAT,
    )

    assert dataset.phospho.columns.equals(columns)
    assert dataset.sample_metadata.index.equals(columns)
    assert dataset.processing_state == processing_state


def test_group_aware_numerical_execution_returns_dedicated_typed_outcome() -> None:
    phospho, site_metadata, sample_metadata = _inputs()
    state = PreprocessingState(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=sample_metadata,
        total=None,
        plan=PreprocessingPlan(
            missing_data_policy="impute_group_aware",
            missing_data_group_column="condition",
            missing_data_min_partial_observed_fraction=0.5,
            missing_data_min_reference_observed_fraction=0.75,
            missing_data_q=0.01,
            missing_data_width=0.3,
            missing_data_seed=123,
            missing_data_k=2,
            missing_data_distance="nan_euclidean",
            stage_order=("missing_data",),
        ),
    )
    routing = route_group_aware_missingness(
        phospho=phospho,
        sample_metadata=sample_metadata,
        group_column="condition",
        min_partial_observed_fraction=0.5,
        min_reference_observed_fraction=0.75,
    )

    outcome = execute_group_aware_imputation(
        state=state,
        routing=routing,
        q=0.01,
        width=0.3,
        seed=123,
        k=2,
        distance="nan_euclidean",
        no_overlap_policy="error",
    )

    assert isinstance(outcome, GroupAwarePolicyOutcome)
    assert not isinstance(outcome, (KnnPolicyOutcome, MinProbPolicyOutcome))
    assert outcome.knn_imputed_cell_count == 1
    assert outcome.minprob_imputed_cell_count == 3
    assert not outcome.phospho.isna().to_numpy().any()


def test_group_aware_row_audit_records_both_mechanisms_for_one_row(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    phospho = pd.DataFrame(
        [
            [1.0, np.nan, np.nan, np.nan, 5.0, 6.0],
            [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
        ],
        index=pd.Index(["mixed", "donor"], name="site_key"),
        columns=pd.Index(["a1", "a2", "b1", "b2", "c1", "c2"], name="sample"),
    )
    original = phospho.copy(deep=True)
    mechanism_inputs: list[pd.DataFrame] = []
    real_knn = missing_data_stage_module.impute_knn_targets
    real_minprob = missing_data_stage_module.impute_minprob_targets

    def _capture_knn_input(phospho: pd.DataFrame, **kwargs: object):
        mechanism_inputs.append(phospho.copy(deep=True))
        return real_knn(phospho, **kwargs)

    def _capture_minprob_input(phospho: pd.DataFrame, **kwargs: object):
        mechanism_inputs.append(phospho.copy(deep=True))
        return real_minprob(phospho, **kwargs)

    monkeypatch.setattr(
        missing_data_stage_module, "impute_knn_targets", _capture_knn_input
    )
    monkeypatch.setattr(
        missing_data_stage_module,
        "impute_minprob_targets",
        _capture_minprob_input,
    )
    state = PreprocessingState(
        phospho=phospho,
        site_metadata=pd.DataFrame({"site": ["S1", "S2"]}, index=phospho.index.copy()),
        sample_metadata=pd.DataFrame(
            {"condition": ["A", "A", "B", "B", "C", "C"]},
            index=phospho.columns.copy(),
        ),
        total=None,
        plan=PreprocessingPlan(
            missing_data_policy="impute_group_aware",
            missing_data_group_column="condition",
            missing_data_min_partial_observed_fraction=0.5,
            missing_data_min_reference_observed_fraction=0.75,
            missing_data_q=0.01,
            missing_data_width=0.3,
            missing_data_seed=42,
            missing_data_k=1,
            missing_data_distance="nan_euclidean",
            stage_order=("missing_data",),
        ),
    )

    result = MissingDataStage().run(state)
    pdt.assert_frame_equal(phospho, original)
    assert len(mechanism_inputs) == 2
    pdt.assert_frame_equal(mechanism_inputs[0], original)
    pdt.assert_frame_equal(mechanism_inputs[1], original)
    assert not result.state.phospho.isna().to_numpy().any()
    observed = original.notna().to_numpy()
    np.testing.assert_array_equal(
        result.state.phospho.to_numpy()[observed],
        original.to_numpy()[observed],
    )
    assert result.state.row_audit is not None
    audit = result.state.row_audit.iloc[0]
    assert audit["reason"] == "group-aware row retained/imputed"
    assert audit["parameter_snapshot"]["mechanisms"] == [
        "partial_observation_knn",
        "asymmetric_absence_minprob",
    ]
    assert audit["parameter_snapshot"]["knn_imputed_cell_count_for_row"] == 1
    assert audit["parameter_snapshot"]["minprob_imputed_cell_count_for_row"] == 2
    assert audit["parameter_snapshot"]["knn_affected_groups"] == ["A"]
    assert audit["parameter_snapshot"]["minprob_affected_groups"] == ["B"]
    assert audit["parameter_snapshot"]["resolved_group_samples"] == {
        "A": ["a1", "a2"],
        "B": ["b1", "b2"],
        "C": ["c1", "c2"],
    }


def test_group_aware_row_audit_distinguishes_unsupported_routes() -> None:
    phospho = pd.DataFrame(
        [
            [1.0, np.nan, np.nan, 4.0],
            [np.nan, np.nan, 3.0, np.nan],
            [1.0, 2.0, 3.0, 4.0],
        ],
        index=pd.Index(["partial", "absence", "donor"], name="site_key"),
        columns=pd.Index(["a1", "a2", "b1", "b2"], name="sample"),
    )
    state = PreprocessingState(
        phospho=phospho,
        site_metadata=pd.DataFrame(
            {"site": ["S1", "S2", "S3"]}, index=phospho.index.copy()
        ),
        sample_metadata=pd.DataFrame(
            {"condition": ["A", "A", "B", "B"]},
            index=phospho.columns.copy(),
        ),
        total=None,
        plan=PreprocessingPlan(
            missing_data_policy="impute_group_aware",
            missing_data_group_column="condition",
            missing_data_min_partial_observed_fraction=0.75,
            missing_data_min_reference_observed_fraction=0.75,
            missing_data_q=0.01,
            missing_data_width=0.3,
            missing_data_seed=42,
            missing_data_k=1,
            missing_data_distance="nan_euclidean",
            stage_order=("missing_data",),
        ),
    )

    result = MissingDataStage().run(state)
    assert result.state.row_audit is not None
    reasons = dict(
        zip(
            result.state.row_audit["source_row_id"],
            result.state.row_audit["reason"],
            strict=True,
        )
    )
    assert reasons == {
        "partial": "unsupported partial coverage",
        "absence": ("unsupported partial coverage and unsupported fully-missing group"),
    }
    snapshots = dict(
        zip(
            result.state.row_audit["source_row_id"],
            result.state.row_audit["parameter_snapshot"],
            strict=True,
        )
    )
    assert snapshots["partial"]["observed_finite_count_by_group"] == {"A": 1, "B": 1}
    assert snapshots["partial"]["observed_fraction_by_group"] == {
        "A": 0.5,
        "B": 0.5,
    }
    assert snapshots["absence"]["observed_finite_count_by_group"] == {"A": 0, "B": 1}
    assert snapshots["absence"]["observed_fraction_by_group"] == {
        "A": 0.0,
        "B": 0.5,
    }
    typed = MissingDataDiagnosticsV2.from_mapping(
        result.diagnostics["diagnostics"], field_name="diagnostics"
    )
    assert typed.group_aware.unsupported_partial_row_ids == ("partial", "absence")
    assert typed.group_aware.unsupported_absence_row_ids == ("absence",)


def test_group_aware_seed_reproducibility_changes_only_minprob_cells() -> None:
    phospho, _, _ = _inputs()
    knn_row, minprob_row = phospho.index[:2]
    first = _build(123)
    repeated = _build(123)
    changed_seed = _build(456)

    pdt.assert_frame_equal(first.phospho, repeated.phospho)
    first_mask = first.imputation_observed_mask_dataframe()
    repeated_mask = repeated.imputation_observed_mask_dataframe()
    changed_mask = changed_seed.imputation_observed_mask_dataframe()
    pdt.assert_frame_equal(
        first_mask,
        repeated_mask,
    )
    pdt.assert_frame_equal(
        first_mask,
        changed_mask,
    )
    assert first.phospho.at[knn_row, "A3"] == changed_seed.phospho.at[knn_row, "A3"]
    assert not np.array_equal(
        first.phospho.loc[minprob_row, ["A1", "A2"]].to_numpy(),
        changed_seed.phospho.loc[minprob_row, ["A1", "A2"]].to_numpy(),
    )
    first_parameters = (_missing_data_stage(first).diagnostics or {})[
        "method_parameters"
    ]
    changed_parameters = (_missing_data_stage(changed_seed).diagnostics or {})[
        "method_parameters"
    ]
    assert (
        first_parameters["knn_target_mask_hash"]
        == changed_parameters["knn_target_mask_hash"]
    )
    assert (
        first_parameters["minprob_target_mask_hash"]
        == changed_parameters["minprob_target_mask_hash"]
    )


def test_group_aware_minprob_only_rows_do_not_inflate_knn_work_estimate(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    phospho = pd.DataFrame(
        [
            [10.0, np.nan, 10.2, 11.0, 11.1, 10.9],
            [np.nan, np.nan, np.nan, 11.0, 11.1, 10.9],
            [np.nan, np.nan, np.nan, 12.0, 12.1, 11.9],
            [np.nan, np.nan, np.nan, 13.0, 13.1, 12.9],
            [9.9, 10.1, 10.2, 10.8, 11.0, 10.9],
        ],
        index=pd.Index(
            ["knn_target", "minprob_1", "minprob_2", "minprob_3", "donor"],
            name="site_key",
        ),
        columns=pd.Index(["a1", "a2", "a3", "b1", "b2", "b3"], name="sample"),
    )
    state = PreprocessingState(
        phospho=phospho,
        site_metadata=pd.DataFrame(
            {"site": ["S1", "S2", "S3", "S4", "S5"]},
            index=phospho.index.copy(),
        ),
        sample_metadata=pd.DataFrame(
            {"condition": ["A", "A", "A", "B", "B", "B"]},
            index=phospho.columns.copy(),
        ),
        total=None,
        plan=PreprocessingPlan(
            missing_data_policy="impute_group_aware",
            missing_data_group_column="condition",
            missing_data_min_partial_observed_fraction=0.5,
            missing_data_min_reference_observed_fraction=0.75,
            missing_data_q=0.01,
            missing_data_width=0.3,
            missing_data_seed=42,
            missing_data_k=1,
            missing_data_distance="nan_euclidean",
            stage_order=("missing_data",),
        ),
    )
    one_knn_target_row_work = phospho.shape[0] * phospho.shape[1]
    monkeypatch.setattr(
        knn_module,
        "KNN_MAX_DISTANCE_FEATURE_OPERATIONS",
        one_knn_target_row_work,
    )

    result = MissingDataStage().run(state)
    diagnostics = MissingDataDiagnosticsV2.from_mapping(
        result.diagnostics["diagnostics"], field_name="diagnostics"
    )

    assert diagnostics.group_aware.knn_routed_cell_count == 1
    assert diagnostics.group_aware.minprob_routed_cell_count == 9
    assert not result.state.phospho.isna().to_numpy().any()


def test_group_aware_processing_state_survives_supported_bundle_round_trip(
    tmp_path,
) -> None:
    dataset = _build(123)
    display_ids = dataset.site_metadata["display_id"].astype(str).tolist()
    reference_sequences = dataset.site_metadata.set_index("display_id").loc[
        :, ["site_sequence"]
    ]
    references = ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["K1", "K1", "K2", "K2"],
                "substrate_site": [
                    display_ids[0],
                    display_ids[1],
                    display_ids[1],
                    display_ids[2],
                ],
            }
        ),
        site_sequences=reference_sequences,
    )
    request = KinaseWorkflowRequest(
        dataset=dataset,
        references=references,
        scoring_config=KinaseScoringConfig(
            reliability_profile="custom",
            min_substrates=2,
            reference_context_compatibility_policy=(
                ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT
            ),
        ),
        prediction_config=KinasePredictionConfig(
            top_k=2,
            deterministic_max_selected_kinases=2,
            adaptive_ensemble_runs=2,
        ),
        activity_config=None,
    )
    result = KinaseWorkflow().run(request)
    bundle_root = tmp_path / "group_aware_kinase_bundle"
    save_kinase_workflow_bundle(
        result,
        bundle_root,
        config_snapshot=KinaseWorkflowConfigSnapshot.from_request(request),
        output_format="csv",
    )

    loaded = load_kinase_workflow_bundle(bundle_root).result.dataset
    diagnostics = loaded.processing_state.missing_data.diagnostics
    original_diagnostics = dataset.processing_state.missing_data.diagnostics
    assert loaded.processing_state.missing_data.policy.value == "impute_group_aware"
    assert isinstance(diagnostics, MissingDataDiagnosticsV2)
    assert isinstance(original_diagnostics, MissingDataDiagnosticsV2)
    assert diagnostics.to_payload() == original_diagnostics.to_payload()
    assert (
        diagnostics.group_aware.knn_target_mask_hash
        == original_diagnostics.group_aware.knn_target_mask_hash
    )
    assert (
        diagnostics.group_aware.minprob_target_mask_hash
        == original_diagnostics.group_aware.minprob_target_mask_hash
    )
    assert (
        diagnostics.group_aware.knn_imputation_mask_hash
        == original_diagnostics.group_aware.knn_imputation_mask_hash
    )
    assert (
        diagnostics.group_aware.minprob_imputation_mask_hash
        == original_diagnostics.group_aware.minprob_imputation_mask_hash
    )
    assert diagnostics.imputation_mask_hash == original_diagnostics.imputation_mask_hash
    assert diagnostics.group_aware.rejected_rows == (
        original_diagnostics.group_aware.rejected_rows
    )
    assert diagnostics.group_aware.routed_rows == (
        original_diagnostics.group_aware.routed_rows
    )

    manifest_path = bundle_root / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    valid_manifest = deepcopy(manifest)

    tampered_retained_count = deepcopy(valid_manifest)
    tampered_group_aware = tampered_retained_count["dataset"]["metadata"][
        "processing_state"
    ]["missing_data"]["diagnostics"]["group_aware"]
    tampered_group_aware["retained_row_count"] += 1
    manifest_path.write_text(
        json.dumps(tampered_retained_count, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    with pytest.raises(DatasetValidationError, match="retained_row_count"):
        load_kinase_workflow_bundle(bundle_root)

    tampered_hash = deepcopy(valid_manifest)
    tampered_diagnostics = tampered_hash["dataset"]["metadata"]["processing_state"][
        "missing_data"
    ]["diagnostics"]
    tampered_diagnostics["group_aware"]["knn_target_mask_hash"] = "tampered"
    tampered_diagnostics["method_parameters"]["knn_target_mask_hash"] = "tampered"
    manifest_path.write_text(
        json.dumps(tampered_hash, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    with pytest.raises(DatasetValidationError, match="knn_target_mask_hash"):
        load_kinase_workflow_bundle(bundle_root)

    tampered_groups = deepcopy(valid_manifest)
    tampered_method_parameters = tampered_groups["dataset"]["metadata"][
        "processing_state"
    ]["missing_data"]["diagnostics"]["method_parameters"]
    tampered_method_parameters["resolved_group_samples"] = {
        "A": ["B1", "B2", "B3"],
        "B": ["A1", "A2", "A3"],
    }
    tampered_routed_rows = tampered_groups["dataset"]["metadata"]["processing_state"][
        "missing_data"
    ]["diagnostics"]["group_aware"]["routed_rows"]
    tampered_routed_rows[0]["knn_group_labels"] = ["B"]
    tampered_routed_rows[1]["minprob_group_labels"] = ["B"]
    manifest_path.write_text(
        json.dumps(tampered_groups, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    with pytest.raises(DatasetValidationError, match="resolved_group_samples"):
        load_kinase_workflow_bundle(bundle_root)

    manifest = deepcopy(valid_manifest)
    manifest["dataset"]["metadata"]["processing_state"]["missing_data"]["policy"] = (
        "impute_knn"
    )
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    with pytest.raises(
        PhosPyInputError,
        match=r"missing_data\.policy must match.*diagnostics\.missing_data_policy",
    ):
        load_kinase_workflow_bundle(bundle_root)


@pytest.mark.parametrize("metadata_case", ["missing", "misaligned"])
def test_group_aware_stage_requires_aligned_sample_metadata(metadata_case: str) -> None:
    phospho, site_metadata, sample_metadata = _inputs()
    resolved_metadata = None
    if metadata_case == "misaligned":
        resolved_metadata = sample_metadata.drop(index="B3")

    with pytest.raises(PhosPyInputError, match="sample_metadata"):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=site_metadata,
                sample_metadata=resolved_metadata,
                organism=Organism.RAT,
                preprocessing_config=_group_aware_config(123),
                input_intensity_scale="log2",
                allow_suspicious_declared_input_intensity_scale=True,
            )
        )


def test_standalone_missing_data_policy_does_not_consume_sample_metadata() -> None:
    phospho, site_metadata, _ = _inputs()
    dataset = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            sample_metadata=None,
            organism=Organism.RAT,
            preprocessing_config=DatasetPreprocessingConfig(
                missing_data=DatasetMissingDataConfig(
                    policy="impute_row_median",
                    min_observed_values=1,
                    input_scale="log2",
                )
            ),
            input_intensity_scale="log2",
            allow_suspicious_declared_input_intensity_scale=True,
        )
    )

    stage = _missing_data_stage(dataset)
    assert "dataset.sample_metadata" not in {
        item.name for item in stage.consumed_input_tables
    }


def test_group_aware_knn_no_overlap_error_identifies_cell_context() -> None:
    phospho, site_metadata, sample_metadata = _inputs()
    phospho = phospho.iloc[:2].copy(deep=True)
    phospho.iloc[0] = [1.0, 2.0, np.nan, np.nan, np.nan, np.nan]
    phospho.iloc[1] = [np.nan, np.nan, 4.0, 5.0, 6.0, 7.0]
    site_metadata = site_metadata.loc[phospho.index]
    config = DatasetPreprocessingConfig(
        missing_data=DatasetMissingDataConfig(
            policy="impute_group_aware",
            group_column="condition",
            min_partial_observed_fraction=0.3,
            min_reference_observed_fraction=0.3,
            q=0.01,
            width=0.3,
            seed=123,
            k=1,
            distance="nan_euclidean",
        )
    )

    with pytest.raises(
        PhosPyInputError,
        match=r"no eligible donor.*affected row=.*affected column='A3'",
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=site_metadata,
                sample_metadata=sample_metadata,
                organism=Organism.RAT,
                preprocessing_config=config,
                input_intensity_scale="log2",
                allow_suspicious_declared_input_intensity_scale=True,
            )
        )
