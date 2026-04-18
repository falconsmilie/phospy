from __future__ import annotations

import inspect

import pandas as pd
import pytest

from phospy.api.configs import KinasePredictionConfig, KinaseScoringConfig
from phospy.api.requests import (
    DatasetBuildRequest,
    KinaseWorkflowRequest,
    SignalomeWorkflowRequest,
)
from phospy.api.results import (
    KinaseWorkflowResult,
    SignalomeWorkflowResult,
)
from phospy.datasets.models import AnalysisReadyPhosphoDataset
from phospy.errors import ReferenceCompatibilityError
from phospy.references.models import Organism, ReferencePreset
from phospy.transformations.models import TransformationState
from phospy.validation.datasets.analysis_ready import AnalysisReadyDatasetValidator
from phospy.workflows.kinase.public import KinaseWorkflow
from phospy.workflows.kinase.validator import KinaseWorkflowValidator


def _dataset() -> AnalysisReadyPhosphoDataset:
    index = pd.Index(["MAPK14;Y182;"], name="site_id")
    return AnalysisReadyPhosphoDataset(
        phospho=pd.DataFrame({"sample_a": [1.0]}, index=index),
        site_metadata=pd.DataFrame(
            {
                "gene_symbol": ["MAPK14"],
                "site": ["Y182"],
                "site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"],
            },
            index=index,
        ),
        organism=Organism.RAT,
        transformation_state=TransformationState.raw(has_total_matrix=False),
    )


def test_request_dtos_do_not_own_runtime_validation_logic() -> None:
    assert "__post_init__" not in DatasetBuildRequest.__dict__
    assert "__post_init__" not in KinaseWorkflowRequest.__dict__
    assert "__post_init__" not in SignalomeWorkflowRequest.__dict__


def test_result_ownership_keeps_kinase_result_as_typed_container() -> None:
    assert "__post_init__" not in KinaseWorkflowResult.__dict__


def test_result_ownership_keeps_signalome_result_contract_check_narrow() -> None:
    source = inspect.getsource(SignalomeWorkflowResult.__post_init__)
    assert "signalome_result.expanded_signalome" in source
    assert "must be KinaseWorkflowResult" not in source
    assert "must be SignalomeAssignments" not in source
    assert "must be SignalomeModules" not in source
    assert "must be KinaseNetwork" not in source


def test_reference_compatibility_is_owned_by_reference_resolution() -> None:
    request = KinaseWorkflowRequest(
        dataset=_dataset(),
        references=ReferencePreset.HUMAN,
        scoring_config=KinaseScoringConfig(min_substrates=2),
        prediction_config=KinasePredictionConfig(top_k=5, ensemble_size=5),
        activity_config=None,
    )
    assert KinaseWorkflowValidator().run(request) is request
    with pytest.raises(ReferenceCompatibilityError):
        KinaseWorkflow().run(request)


def test_dataset_validation_composition_is_outside_validation_subdomains() -> None:
    dataset_validator_source = inspect.getsource(AnalysisReadyDatasetValidator)
    dataset_post_init_source = inspect.getsource(
        AnalysisReadyPhosphoDataset.__post_init__
    )
    assert "TransformationStateValidator" not in dataset_validator_source
    assert (
        "transformation_state"
        not in inspect.signature(AnalysisReadyDatasetValidator.run).parameters
    )
    assert "_DATASET_VALIDATOR.run(" in dataset_post_init_source
    assert "_TRANSFORMATION_STATE_VALIDATOR.run(" in dataset_post_init_source
