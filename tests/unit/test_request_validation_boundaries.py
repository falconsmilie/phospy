from __future__ import annotations

from dataclasses import is_dataclass

import pandas as pd
import pytest

from phospy import (
    AnalysisReadyDatasetBuilder,
    AnalysisReadyPhosphoDataset,
    DifferentialAnalysisWorkflow,
    KinaseWorkflow,
    SignalomeWorkflow,
)
from phospy.api import (
    Contrast,
    DatasetBuildRequest,
    DifferentialAnalysisRequest,
    ExperimentalDesign,
    KinaseWorkflowRequest,
    Organism,
    ReferenceBundle,
    SampleDesignRecord,
    SignalomeWorkflowRequest,
)
from phospy.api.results import (
    KinasePredictionResult,
    KinaseScoringResult,
    KinaseWorkflowResult,
)
from phospy.errors import PhosPyInputError, WorkflowValidationError
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
    supported_log2_intensity_scale_state,
    supported_log2_processing_state,
)
from tests.support.site_keys import protein_site_key_index, site_key_context_columns


class _UnexpectedStage:
    def __init__(self, name: str) -> None:
        self.name = name
        self.called = False

    def run(self, request: object) -> object:
        self.called = True
        raise AssertionError(f"{self.name} must not run after request validation fails")


def _assert_not_called(*stages: _UnexpectedStage) -> None:
    assert all(not stage.called for stage in stages)


def _analysis_ready_dataset(
    *,
    log2_scale: bool,
    include_protein_id: bool = True,
) -> AnalysisReadyPhosphoDataset:
    genes = ["MAPK14", "GSK3B"]
    sites = ["Y182", "S9"]
    display_ids = ["MAPK14;Y182;", "GSK3B;S9;"]
    site_index = protein_site_key_index(protein_identifiers=genes, sites=sites)
    phospho = pd.DataFrame(
        {
            "A_1": [1.0, 2.0],
            "A_2": [1.1, 2.1],
            "B_1": [2.1, 2.0],
            "B_2": [2.0, 2.2],
        },
        index=site_index,
    )
    site_metadata_payload: dict[str, object] = {
        "site_key": site_index.astype(str).tolist(),
        "display_id": display_ids,
        **site_key_context_columns(site_index),
        "gene_symbol": genes,
        "site": sites,
        "site_sequence": [
            ("A" * 15) + str(site).strip().upper()[0] + ("A" * 15) for site in sites
        ],
    }
    if include_protein_id:
        site_metadata_payload["protein_id"] = genes
    site_metadata = pd.DataFrame(site_metadata_payload, index=site_index.copy())

    if log2_scale:
        intensity_scale_state = supported_log2_intensity_scale_state(
            has_total_matrix=False
        )
        processing_state = supported_log2_processing_state(has_total_matrix=False)
    else:
        intensity_scale_state = supported_linear_intensity_scale_state(
            has_total_matrix=False
        )
        processing_state = supported_linear_processing_state(has_total_matrix=False)

    return AnalysisReadyPhosphoDataset(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        intensity_scale_state=intensity_scale_state,
        processing_state=processing_state,
    )


def _differential_design() -> ExperimentalDesign:
    return ExperimentalDesign(
        samples=(
            SampleDesignRecord(
                sample_id="A_1",
                condition="A",
                biological_replicate_id="A_r1",
            ),
            SampleDesignRecord(
                sample_id="A_2",
                condition="A",
                biological_replicate_id="A_r2",
            ),
            SampleDesignRecord(
                sample_id="B_1",
                condition="B",
                biological_replicate_id="B_r1",
            ),
            SampleDesignRecord(
                sample_id="B_2",
                condition="B",
                biological_replicate_id="B_r2",
            ),
        )
    )


def _reference_bundle(dataset: AnalysisReadyPhosphoDataset) -> ReferenceBundle:
    display_ids = dataset.site_metadata.loc[:, "display_id"].astype(str).tolist()
    sequences = dataset.site_metadata.loc[:, "site_sequence"].astype(str).tolist()
    return ReferenceBundle(
        organism=dataset.organism,
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["MAP2K6"],
                "substrate_site": [display_ids[0]],
            }
        ),
        site_sequences=pd.DataFrame(
            {"site_sequence": [sequences[0]]},
            index=pd.Index([display_ids[0]], name="site_id"),
        ),
    )


def _kinase_result_missing_signalome_protein_id() -> KinaseWorkflowResult:
    dataset = _analysis_ready_dataset(log2_scale=False, include_protein_id=False)
    site_index = dataset.phospho.index.copy()
    prediction_matrix = pd.DataFrame({"MAP2K6": [0.9, 0.2]}, index=site_index)
    score_matrix = pd.DataFrame({"MAP2K6": [1.5, 0.5]}, index=site_index.copy())
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


@pytest.mark.parametrize(
    "request_type",
    [
        DatasetBuildRequest,
        DifferentialAnalysisRequest,
        KinaseWorkflowRequest,
        SignalomeWorkflowRequest,
    ],
)
def test_request_dataclasses_document_inert_payload_boundary(
    request_type: type[object],
) -> None:
    assert is_dataclass(request_type)
    assert "__post_init__" not in request_type.__dict__
    assert "Construction stores the payload only" in (request_type.__doc__ or "")


def test_request_constructors_accept_payloads_without_scientific_validation() -> None:
    DatasetBuildRequest(
        phospho=object(),  # type: ignore[arg-type]
        site_metadata=object(),  # type: ignore[arg-type]
        site_resolution_mode="not-a-mode",
    )
    DifferentialAnalysisRequest(
        dataset=object(),  # type: ignore[arg-type]
        design=object(),  # type: ignore[arg-type]
        contrasts=object(),  # type: ignore[arg-type]
        config=object(),  # type: ignore[arg-type]
    )
    KinaseWorkflowRequest(
        dataset=object(),  # type: ignore[arg-type]
        references=object(),  # type: ignore[arg-type]
        reference_display_ambiguity_policy="not-a-policy",  # type: ignore[arg-type]
    )
    SignalomeWorkflowRequest(
        kinase_result=object(),  # type: ignore[arg-type]
        config=object(),  # type: ignore[arg-type]
    )


def test_invalid_dataset_build_request_fails_at_builder_validation_stage() -> None:
    interpreter = _UnexpectedStage("dataset builder interpreter")
    executor = _UnexpectedStage("dataset builder executor")
    request = DatasetBuildRequest(site_resolution_mode="not-a-mode")

    with pytest.raises(PhosPyInputError, match="site_resolution_mode"):
        AnalysisReadyDatasetBuilder(
            interpreter=interpreter,  # type: ignore[arg-type]
            executor=executor,  # type: ignore[arg-type]
        ).run(request)

    _assert_not_called(interpreter, executor)


def test_invalid_differential_request_fails_at_workflow_validation_stage() -> None:
    interpreter = _UnexpectedStage("differential interpreter")
    executor = _UnexpectedStage("differential executor")
    request = DifferentialAnalysisRequest(
        dataset=_analysis_ready_dataset(log2_scale=True),
        design=_differential_design(),
        contrasts=(
            Contrast(
                name="B_vs_C",
                numerator_condition="B",
                denominator_condition="C",
            ),
        ),
    )

    with pytest.raises(WorkflowValidationError, match="denominator condition"):
        DifferentialAnalysisWorkflow(
            interpreter=interpreter,  # type: ignore[arg-type]
            executor=executor,  # type: ignore[arg-type]
        ).run(request)

    _assert_not_called(interpreter, executor)


def test_invalid_kinase_request_fails_at_workflow_validation_stage() -> None:
    interpreter = _UnexpectedStage("kinase interpreter")
    executor = _UnexpectedStage("kinase executor")
    request = KinaseWorkflowRequest(
        dataset=_analysis_ready_dataset(log2_scale=False),
        reference_display_ambiguity_policy="not-a-policy",  # type: ignore[arg-type]
    )

    with pytest.raises(
        WorkflowValidationError,
        match="reference_display_ambiguity_policy",
    ):
        KinaseWorkflow(
            interpreter=interpreter,  # type: ignore[arg-type]
            executor=executor,  # type: ignore[arg-type]
        ).run(request)

    _assert_not_called(interpreter, executor)


def test_invalid_signalome_request_fails_at_workflow_validation_stage() -> None:
    interpreter = _UnexpectedStage("signalome interpreter")
    executor = _UnexpectedStage("signalome executor")
    request = SignalomeWorkflowRequest(
        kinase_result=_kinase_result_missing_signalome_protein_id()
    )

    with pytest.raises(WorkflowValidationError) as exc_info:
        SignalomeWorkflow(
            interpreter=interpreter,  # type: ignore[arg-type]
            executor=executor,  # type: ignore[arg-type]
        ).run(request)

    message = str(exc_info.value)
    assert "protein_id" in message
    assert "signalome workflow request kinase_result.dataset.site_metadata" in message
    _assert_not_called(interpreter, executor)
