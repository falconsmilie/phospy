from __future__ import annotations

import pandas as pd
import pytest

from phospy import AnalysisReadyDatasetBuilder, KinaseWorkflow, SignalomeWorkflow
from phospy.api import (
    DatasetBuildRequest,
    KinaseWorkflowRequest,
    Organism,
    ReferenceBundle,
    ReferencePreset,
    SignalomeWorkflowRequest,
)
from phospy.api.results import (
    KinasePredictionResult,
    KinaseScoringResult,
    KinaseWorkflowResult,
)
from phospy.errors import (
    PhosPyInputError,
    ReferenceValidationError,
    SignalomeScaleError,
    WorkflowValidationError,
)
from tests.support.signalome_config import build_signalome_config


def _dataset_request() -> DatasetBuildRequest:
    return DatasetBuildRequest(
        phospho=pd.DataFrame(
            {"sample_a": [1.0]},
            index=pd.Index(["MAPK14;Y182;"], name="site_id"),
        ),
        site_metadata=pd.DataFrame(
            {
                "gene_symbol": ["MAPK14"],
                "site": ["Y182"],
                "site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"],
            },
            index=pd.Index(["MAPK14;Y182;"], name="site_id"),
        ),
        organism=Organism.RAT,
    )


def _dataset():
    request = _dataset_request()
    return AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=request.phospho,
            site_metadata=request.site_metadata,
            organism=request.organism,
        )
    )


def _kinase_request() -> KinaseWorkflowRequest:
    return KinaseWorkflowRequest(
        dataset=_dataset(),
        references=ReferencePreset.AUTO,
    )


def _kinase_result() -> KinaseWorkflowResult:
    return KinaseWorkflowResult(
        dataset=_dataset(),
        references=ReferenceBundle(
            organism=Organism.RAT,
            kinase_substrate_map=pd.DataFrame(
                {"kinase": ["MAP2K6"], "substrate_site": ["MAPK14;Y182;"]}
            ),
            site_sequences=pd.DataFrame(
                {"site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"]},
                index=pd.Index(["MAPK14;Y182;"], name="site_id"),
            ),
        ),
        scoring_result=KinaseScoringResult(
            profile_scores=pd.DataFrame(
                {"MAP2K6": [1.0]},
                index=pd.Index(["MAPK14;Y182;"], name="site_id"),
            )
        ),
        prediction_result=KinasePredictionResult(
            pred_mat=pd.DataFrame(
                {"MAP2K6": [0.8]},
                index=pd.Index(["MAPK14;Y182;"], name="site_id"),
            )
        ),
    )


def _signalome_request() -> SignalomeWorkflowRequest:
    return SignalomeWorkflowRequest(
        kinase_result=_kinase_result(),
        config=build_signalome_config(),
    )


class _PassThroughValidator:
    def run(self, request: object) -> object:
        return request


class _PassThroughInterpreter:
    def run(self, request: object) -> object:
        return request


class _UnexpectedExecutor:
    def run(self, resolved: object) -> object:
        raise RuntimeError("unexpected internal defect")


def test_builder_propagates_invalid_input_table_failures() -> None:
    class InvalidInputValidator:
        def run(self, request: object) -> object:
            raise PhosPyInputError("dataset.phospho has an invalid input table shape")

    builder = AnalysisReadyDatasetBuilder(
        validator=InvalidInputValidator(),
        interpreter=_PassThroughInterpreter(),
        executor=_UnexpectedExecutor(),
    )
    with pytest.raises(PhosPyInputError, match="invalid input table shape"):
        builder.run(_dataset_request())


def test_kinase_workflow_propagates_invalid_user_config_failures() -> None:
    class InvalidConfigValidator:
        def run(self, request: object) -> object:
            raise WorkflowValidationError(
                "invalid user config: scoring_config.min_substrates must be >= 2"
            )

    workflow = KinaseWorkflow(
        validator=InvalidConfigValidator(),
        interpreter=_PassThroughInterpreter(),
        executor=_UnexpectedExecutor(),
    )
    with pytest.raises(WorkflowValidationError, match="invalid user config"):
        workflow.run(_kinase_request())


def test_kinase_workflow_propagates_invalid_reference_schema_failures() -> None:
    class InvalidReferenceInterpreter:
        def run(self, request: object) -> object:
            raise ReferenceValidationError(
                "reference schema validation failed: missing substrate_site column"
            )

    workflow = KinaseWorkflow(
        validator=_PassThroughValidator(),
        interpreter=InvalidReferenceInterpreter(),
        executor=_UnexpectedExecutor(),
    )
    with pytest.raises(ReferenceValidationError, match="reference schema validation"):
        workflow.run(_kinase_request())


def test_signalome_workflow_propagates_resource_limit_guard_failures() -> None:
    class ResourceLimitExecutor:
        def run(self, resolved: object) -> object:
            raise SignalomeScaleError(
                "signalome scale guard blocked execution: site count exceeds max_full_candidate_scoring_sites"
            )

    workflow = SignalomeWorkflow(
        validator=_PassThroughValidator(),
        interpreter=_PassThroughInterpreter(),
        executor=ResourceLimitExecutor(),
    )
    with pytest.raises(SignalomeScaleError, match="scale guard blocked execution"):
        workflow.run(_signalome_request())


def test_builder_does_not_wrap_unexpected_internal_exceptions() -> None:
    builder = AnalysisReadyDatasetBuilder(
        validator=_PassThroughValidator(),
        interpreter=_PassThroughInterpreter(),
        executor=_UnexpectedExecutor(),
    )
    with pytest.raises(RuntimeError, match="unexpected internal defect"):
        builder.run(_dataset_request())
