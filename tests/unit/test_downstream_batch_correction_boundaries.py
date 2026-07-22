from __future__ import annotations

import ast
import inspect
from dataclasses import fields, replace
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

from phospy import AnalysisReadyPhosphoDataset
from phospy.api import (
    BatchCovariate,
    Contrast,
    DifferentialAnalysisConfig,
    DifferentialAnalysisRequest,
    EnrichmentConfig,
    EnrichmentWorkflowRequest,
    KinaseActivityConfig,
    KinasePredictionConfig,
    KinaseScoringConfig,
    KinaseWorkflowRequest,
    Organism,
    PtmSetCollection,
    ReferenceBundle,
    SampleDesignRecord,
    SignalomeConfig,
    SignalomeWorkflowRequest,
)
from phospy.api.requests import ExperimentalDesign
from phospy.api.results import (
    KinasePredictionResult,
    KinaseScoringResult,
    KinaseWorkflowResult,
)
from phospy.errors import PhosPyInputError
from phospy.science.datasets.models import DatasetPreprocessingReport
from phospy.science.datasets.preprocessing.batch_correction import (
    BatchCorrectionDiagnostics,
    BatchCorrectionPolicy,
    BatchCorrectionReport,
    LinearResidualizeBatchCorrectionEngine,
)
from phospy.science.datasets.processing_state import RuvReadinessState
from phospy.validation.datasets.preprocessing import (
    reject_external_corrected_output_after_downstream_preprocessing,
)
from phospy.workflows.batch_correction import BatchCorrectionWorkflow
from phospy.workflows.differential.interpreter import DifferentialAnalysisInterpreter
from phospy.workflows.differential.validator import DifferentialAnalysisValidator
from tests.support.analysis_ready_dataset_factories import (
    trusted_analysis_ready_dataset_from_tables,
)
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
    supported_log2_intensity_scale_state,
    supported_log2_processing_state,
)
from tests.support.site_keys import protein_site_key_index, site_key_context_columns

_DOWNSTREAM_WORKFLOW_DIRS = (
    Path("src/phospy/workflows/differential"),
    Path("src/phospy/workflows/kinase"),
    Path("src/phospy/workflows/enrichment"),
    Path("src/phospy/workflows/signalome"),
)
_FORBIDDEN_CORRECTION_SURFACE_TOKENS = (
    "batch_correction",
    "control_site",
    "negative_control",
    "n_unwanted",
    "ruv",
    "sps",
    "unwanted_factor",
)


def _corrected_dataset(
    *,
    log2_scale: bool,
    ruv_readiness: RuvReadinessState | None = None,
) -> AnalysisReadyPhosphoDataset:
    genes = ["MAPK14", "GSK3B"]
    sites = ["Y182", "S9"]
    display_ids = ["MAPK14;Y182;", "GSK3B;S9;"]
    site_index = protein_site_key_index(protein_identifiers=genes, sites=sites)
    phospho = pd.DataFrame(
        {
            "A_1": [10.0, 20.0],
            "A_2": [10.2, 20.2],
            "B_1": [11.0, 19.8],
            "B_2": [11.1, 20.0],
        },
        index=site_index,
    )
    site_metadata = pd.DataFrame(
        {
            "site_key": site_index.astype(str).tolist(),
            "display_id": display_ids,
            **site_key_context_columns(site_index),
            "gene_symbol": genes,
            "site": sites,
            "site_sequence": [
                ("A" * 15) + str(site).strip().upper()[0] + ("A" * 15) for site in sites
            ],
            "protein_id": genes,
            "localisation_confidence": [0.95, 0.91],
        },
        index=site_index.copy(),
    )
    sample_metadata = pd.DataFrame(
        {
            "batch": ["run_1", "run_2", "run_1", "run_2"],
            "condition": ["A", "A", "B", "B"],
            "replicate_group": ["r1", "r2", "r1", "r2"],
        },
        index=phospho.columns.copy(),
    )
    preprocessing_report = DatasetPreprocessingReport.from_rows(
        batch_correction=_applied_batch_correction_report(phospho.shape)
    )
    intensity_scale_state = (
        supported_log2_intensity_scale_state(has_total_matrix=False)
        if log2_scale
        else supported_linear_intensity_scale_state(has_total_matrix=False)
    )
    processing_state = (
        supported_log2_processing_state(has_total_matrix=False)
        if log2_scale
        else supported_linear_processing_state(has_total_matrix=False)
    )
    if ruv_readiness is not None:
        processing_state = replace(processing_state, ruv_readiness=ruv_readiness)
    return trusted_analysis_ready_dataset_from_tables(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=sample_metadata,
        organism=Organism.RAT,
        intensity_scale_state=intensity_scale_state,
        processing_state=processing_state,
        preprocessing_report=preprocessing_report,
    )


def _applied_batch_correction_report(shape: tuple[int, int]) -> BatchCorrectionReport:
    return BatchCorrectionReport(
        status="applied",
        policy=BatchCorrectionPolicy(
            method="linear_residualize_batch",
            batch_column="batch",
            condition_column="condition",
        ),
        diagnostics=BatchCorrectionDiagnostics(
            number_of_batches=2,
            batch_levels=("run_1", "run_2"),
            condition_levels=("A", "B"),
            confounding_check_status="passed",
            matrix_shape_before=shape,
            matrix_shape_after=shape,
        ),
    )


def _differential_design(
    *, include_batch_covariate: bool = False
) -> ExperimentalDesign:
    fixed_effects = (BatchCovariate(),) if include_batch_covariate else ()
    return ExperimentalDesign(
        samples=(
            SampleDesignRecord(
                sample_id="A_1",
                condition="A",
                biological_replicate_id="A_r1",
                batch="run_1",
            ),
            SampleDesignRecord(
                sample_id="A_2",
                condition="A",
                biological_replicate_id="A_r2",
                batch="run_2",
            ),
            SampleDesignRecord(
                sample_id="B_1",
                condition="B",
                biological_replicate_id="B_r1",
                batch="run_1",
            ),
            SampleDesignRecord(
                sample_id="B_2",
                condition="B",
                biological_replicate_id="B_r2",
                batch="run_2",
            ),
        ),
        fixed_effects=fixed_effects,
    )


def _contrast() -> tuple[Contrast, ...]:
    return (
        Contrast(
            name="B_vs_A",
            numerator_condition="B",
            denominator_condition="A",
        ),
    )


def _reference_bundle(dataset: AnalysisReadyPhosphoDataset) -> ReferenceBundle:
    metadata = dataset.site_metadata
    display_ids = metadata.loc[:, "display_id"].astype(str).tolist()
    return ReferenceBundle(
        organism=dataset.organism,
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["MAP2K6", "MAP2K6"],
                "substrate_site": display_ids,
            }
        ),
        site_sequences=pd.DataFrame(
            {"site_sequence": metadata.loc[:, "site_sequence"].astype(str).tolist()},
            index=pd.Index(display_ids, name="site_id"),
        ),
    )


def _kinase_result(dataset: AnalysisReadyPhosphoDataset) -> KinaseWorkflowResult:
    site_index = dataset.phospho.index.copy()
    score_matrix = pd.DataFrame(
        {"MAP2K6": [1.0, 0.5]},
        index=site_index.copy(),
    )
    prediction_matrix = pd.DataFrame(
        {"MAP2K6": [0.9, 0.2]},
        index=site_index.copy(),
    )
    return KinaseWorkflowResult(
        dataset=dataset,
        references=_reference_bundle(dataset),
        scoring_result=KinaseScoringResult(
            profile_scores=score_matrix,
            rank_weighted_fusion_scores=score_matrix.copy(deep=True),
        ),
        prediction_result=KinasePredictionResult(pred_mat=prediction_matrix),
        activity_result=None,
    )


def _ready_ruv_state() -> RuvReadinessState:
    return RuvReadinessState(
        enabled=True,
        ready=True,
        reasons=(),
        control_feature_column="is_control_feature",
        replicate_group_column="replicate_group",
        batch_column="batch",
        control_feature_count=1,
        replicate_group_count=2,
        batch_count=2,
        requires_complete_matrix=True,
        matrix_complete=True,
        imputation_method_id=None,
        missingness_mask_preserved=False,
    )


def test_downstream_requests_accept_corrected_analysis_ready_boundaries() -> None:
    differential_dataset = _corrected_dataset(log2_scale=True)
    kinase_dataset = _corrected_dataset(log2_scale=False)
    signalome_dataset = _corrected_dataset(log2_scale=False)
    enrichment_dataset = _corrected_dataset(log2_scale=False)

    differential_request = DifferentialAnalysisRequest(
        dataset=differential_dataset,
        design=_differential_design(),
        contrasts=_contrast(),
    )
    kinase_request = KinaseWorkflowRequest(
        dataset=kinase_dataset,
        references=_reference_bundle(kinase_dataset),
    )
    signalome_request = SignalomeWorkflowRequest(
        kinase_result=_kinase_result(signalome_dataset)
    )
    enrichment_request = EnrichmentWorkflowRequest(
        identifier_column="site_key",
        identifier_kind="site_key",
        set_collection=PtmSetCollection(
            sets={"corrected_input_sites": tuple(enrichment_dataset.phospho.index)},
            identifier_kind="site_key",
        ),
        input_table=enrichment_dataset.site_metadata.loc[:, ["site_key"]],
        background_universe=tuple(enrichment_dataset.phospho.index.astype(str)),
    )

    assert differential_request.dataset is differential_dataset
    assert kinase_request.dataset is kinase_dataset
    assert signalome_request.kinase_result.dataset is signalome_dataset
    assert enrichment_request.input_table is not None


@pytest.mark.parametrize(
    "stage",
    [
        "differential_analysis_preparation",
        "kinase_activity_preparation",
        "enrichment_preparation",
        "signalome_preparation",
    ],
)
def test_external_corrected_output_boundary_rejects_downstream_preparation_aliases(
    stage: str,
) -> None:
    with pytest.raises(
        PhosPyInputError,
        match=(
            "external corrected output cannot be integrated after downstream "
            "preprocessing stages.*only matrix-changing preprocessing input.*"
            "SpsRuvBatchCorrectionConfig"
        ),
    ):
        reject_external_corrected_output_after_downstream_preprocessing(
            ("missing_data", "batch_correction", stage)
        )


@pytest.mark.parametrize(
    "request_factory",
    [
        lambda: DifferentialAnalysisRequest(
            dataset=_corrected_dataset(log2_scale=True),
            design=_differential_design(),
            contrasts=_contrast(),
        ),
        lambda: KinaseWorkflowRequest(
            dataset=_corrected_dataset(log2_scale=False),
            references=_reference_bundle(_corrected_dataset(log2_scale=False)),
        ),
        lambda: SignalomeWorkflowRequest(
            kinase_result=_kinase_result(_corrected_dataset(log2_scale=False))
        ),
        lambda: EnrichmentWorkflowRequest(
            identifier_column="site_key",
            identifier_kind="site_key",
            set_collection=PtmSetCollection(
                sets={"sites": ("site_a",)},
                identifier_kind="site_key",
            ),
            selected_identifiers=("site_a",),
            background_universe=("site_a",),
        ),
    ],
)
def test_downstream_requests_do_not_accept_direct_ruv_or_sps_config_channels(
    request_factory: Any,
) -> None:
    request = request_factory()
    request_type = type(request)

    with pytest.raises(TypeError, match="unexpected keyword argument"):
        request_type(
            **{field.name: getattr(request, field.name) for field in fields(request)},
            ruv_config=object(),
        )
    with pytest.raises(TypeError, match="unexpected keyword argument"):
        request_type(
            **{field.name: getattr(request, field.name) for field in fields(request)},
            sps_correction_config=object(),
        )


@pytest.mark.parametrize(
    "contract_type",
    [
        DifferentialAnalysisRequest,
        DifferentialAnalysisConfig,
        KinaseWorkflowRequest,
        KinaseScoringConfig,
        KinasePredictionConfig,
        KinaseActivityConfig,
        SignalomeWorkflowRequest,
        SignalomeConfig,
        EnrichmentWorkflowRequest,
        EnrichmentConfig,
    ],
)
def test_downstream_public_contracts_do_not_expose_correction_parameters(
    contract_type: type[object],
) -> None:
    field_names = {field.name for field in fields(contract_type)}
    signature_names = set(inspect.signature(contract_type).parameters)
    exposed_names = field_names | signature_names

    forbidden = {
        name
        for name in exposed_names
        if any(token in name.lower() for token in _FORBIDDEN_CORRECTION_SURFACE_TOKENS)
    }
    assert forbidden == set()


def test_downstream_workflows_do_not_import_batch_correction_executors() -> None:
    forbidden_imports: list[str] = []
    for root in _DOWNSTREAM_WORKFLOW_DIRS:
        for path in root.rglob("*.py"):
            tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
            for node in ast.walk(tree):
                if isinstance(node, ast.ImportFrom):
                    module = node.module or ""
                    imported_names = {alias.name for alias in node.names}
                    if module in {
                        "phospy.workflows.batch_correction",
                        "phospy.workflows.batch_correction.workflow",
                        "phospy.science.batch_correction",
                        "phospy.science.batch_correction.executor",
                    }:
                        forbidden_imports.append(f"{path}:{module}")
                    if (
                        module
                        == "phospy.science.datasets.preprocessing.batch_correction"
                        and imported_names
                        & {
                            "BatchCorrectionEngine",
                            "LinearResidualizeBatchCorrectionEngine",
                        }
                    ):
                        forbidden_imports.append(f"{path}:{module}")
                if isinstance(node, ast.Import):
                    for alias in node.names:
                        if alias.name in {
                            "phospy.workflows.batch_correction",
                            "phospy.science.batch_correction",
                        }:
                            forbidden_imports.append(f"{path}:{alias.name}")

    assert forbidden_imports == []


def test_differential_batch_covariate_remains_model_term_not_preprocessing_correction(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def fail_if_batch_correction_runs(*args: object, **kwargs: object) -> object:
        raise AssertionError("downstream differential must not run batch correction")

    monkeypatch.setattr(
        LinearResidualizeBatchCorrectionEngine,
        "run",
        fail_if_batch_correction_runs,
    )
    monkeypatch.setattr(BatchCorrectionWorkflow, "run", fail_if_batch_correction_runs)
    request = DifferentialAnalysisRequest(
        dataset=_corrected_dataset(log2_scale=True),
        design=_differential_design(include_batch_covariate=True),
        contrasts=_contrast(),
    )

    interpreted = DifferentialAnalysisInterpreter().run(
        DifferentialAnalysisValidator().run(request)
    )

    assert interpreted.execution_design is not None
    batch_covariates = [
        covariate
        for covariate in interpreted.execution_design.covariate_columns
        if covariate.name == "batch"
    ]
    assert len(batch_covariates) == 1
    assert batch_covariates[0].kind == "batch"
    assert batch_covariates[0].columns == ("batch[run_2]",)
    assert "batch" in interpreted.execution_design.formula
    assert (
        interpreted.dataset_preprocessing_report is request.dataset.preprocessing_report
    )
    assert interpreted.computation_request.design.to_dataframe().columns.tolist() == [
        "A",
        "B",
        "batch[run_2]",
    ]


def test_ruv_readiness_is_report_only_and_does_not_modify_dataset_matrix() -> None:
    dataset = _corrected_dataset(log2_scale=True, ruv_readiness=_ready_ruv_state())
    original = dataset.phospho.copy(deep=True)

    request = DifferentialAnalysisRequest(
        dataset=dataset,
        design=_differential_design(include_batch_covariate=True),
        contrasts=_contrast(),
    )
    validated = DifferentialAnalysisValidator().run(request)

    pd.testing.assert_frame_equal(dataset.phospho, original)
    pd.testing.assert_frame_equal(validated.dataset.phospho, original)
    assert validated.dataset.processing_state.ruv_readiness.ready is True
    assert validated.dataset.preprocessing_report is not None
    assert validated.dataset.preprocessing_report.batch_correction is not None


def test_linear_residualize_batch_contract_remains_distinct_from_sps_ruv() -> None:
    report = _applied_batch_correction_report((2, 4))
    payload = report.to_payload()
    payload_text = repr(payload).lower()

    assert report.method == "linear_residualize_batch"
    assert report.design_preservation_policy == "preserve_condition_effects"
    assert report.confounding_check_status == "passed"
    assert "ruv" not in payload_text
    assert "sps" not in payload_text
