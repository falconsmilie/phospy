from __future__ import annotations

import pandas as pd
import pytest

from phospy.api.configs import (
    KinasePredictionConfig,
    KinaseScoringConfig,
    LocalisationRequirement,
)
from phospy.api.requests import KinaseWorkflowRequest, SignalomeWorkflowRequest
from phospy.api.results import KinaseWorkflowResult
from phospy.errors.validation import WorkflowValidationError
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.science.prediction.models import KinasePredictionResult, KinaseScoringResult
from phospy.science.references.models import Organism, ReferenceBundle
from phospy.workflows.kinase.validator import KinaseWorkflowValidator
from phospy.workflows.signalome.validator import SignalomeWorkflowValidator
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
)
from tests.support.signalome_config import build_signalome_config
from tests.support.site_keys import site_key_index_from_display_ids

pytestmark = pytest.mark.integration


def _dataset(
    *, localisation_probability: object | None = None
) -> AnalysisReadyPhosphoDataset:
    display_ids = ["MAPK14;Y182;"]
    site_index = site_key_index_from_display_ids(
        display_ids,
        protein_namespace="gene_symbol",
    )
    phospho = pd.DataFrame(
        {"sample_a": [1.0], "sample_b": [1.2]},
        index=site_index,
    )
    site_metadata = pd.DataFrame(
        {
            "site_key": site_index.tolist(),
            "display_id": display_ids,
            "gene_symbol": ["MAPK14"],
            "site": ["Y182"],
            "site_sequence": [("A" * 15) + "Y" + ("A" * 15)],
            "protein_id": ["P53778"],
        },
        index=site_index.copy(),
    )
    if localisation_probability is not None:
        site_metadata.loc[:, "localisation_probability"] = [localisation_probability]
    return AnalysisReadyPhosphoDataset(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
    )


def _references() -> ReferenceBundle:
    return ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {"kinase": ["MAP2K6"], "substrate_site": ["MAPK14;Y182;"]}
        ),
        site_sequences=pd.DataFrame(
            {"site_sequence": [("A" * 15) + "Y" + ("A" * 15)]},
            index=pd.Index(["MAPK14;Y182;"], name="site_id"),
        ),
    )


def _kinase_result(dataset: AnalysisReadyPhosphoDataset) -> KinaseWorkflowResult:
    prediction_matrix = pd.DataFrame(
        {"MAP2K6": [0.9]},
        index=dataset.phospho.index.copy(),
    )
    score_matrix = pd.DataFrame(
        {"MAP2K6": [1.5]},
        index=dataset.phospho.index.copy(),
    )
    return KinaseWorkflowResult(
        dataset=dataset,
        references=_references(),
        scoring_result=KinaseScoringResult(
            profile_scores=score_matrix,
            rank_weighted_fusion_scores=score_matrix,
        ),
        prediction_result=KinasePredictionResult(pred_mat=prediction_matrix),
        activity_result=None,
    )


def test_analysis_ready_dataset_with_unknown_localisation_allows_default_kinase_validation() -> (
    None
):
    request = KinaseWorkflowRequest(
        dataset=_dataset(),
        references=_references(),
        scoring_config=KinaseScoringConfig(min_substrates=2),
        prediction_config=KinasePredictionConfig(
            top_k=5,
            deterministic_max_selected_kinases=5,
            adaptive_ensemble_runs=5,
        ),
        activity_config=None,
    )
    validated = KinaseWorkflowValidator().run(request)
    assert validated is request


def test_analysis_ready_dataset_with_unknown_localisation_rejects_localisation_strict_kinase_validation() -> (
    None
):
    request = KinaseWorkflowRequest(
        dataset=_dataset(),
        references=_references(),
        scoring_config=KinaseScoringConfig(
            min_substrates=2,
            localisation_requirement=LocalisationRequirement(require_present=True),
        ),
        prediction_config=KinasePredictionConfig(
            top_k=5,
            deterministic_max_selected_kinases=5,
            adaptive_ensemble_runs=5,
        ),
        activity_config=None,
    )
    with pytest.raises(WorkflowValidationError, match="missing required column"):
        KinaseWorkflowValidator().run(request)


def test_analysis_ready_dataset_with_unknown_localisation_rejects_localisation_strict_signalome_validation() -> (
    None
):
    request = SignalomeWorkflowRequest(
        kinase_result=_kinase_result(_dataset()),
        config=build_signalome_config(
            substrate_support_cutoff=0.5,
            localisation_requirement=LocalisationRequirement(require_present=True),
        ),
    )
    with pytest.raises(WorkflowValidationError, match="missing required column"):
        SignalomeWorkflowValidator().run(request)


def test_analysis_ready_dataset_with_below_threshold_localisation_rejects_threshold_policy() -> (
    None
):
    request = KinaseWorkflowRequest(
        dataset=_dataset(localisation_probability=0.6),
        references=_references(),
        scoring_config=KinaseScoringConfig(
            min_substrates=2,
            localisation_requirement=LocalisationRequirement(minimum_probability=0.75),
        ),
        prediction_config=KinasePredictionConfig(
            top_k=5,
            deterministic_max_selected_kinases=5,
            adaptive_ensemble_runs=5,
        ),
        activity_config=None,
    )
    with pytest.raises(WorkflowValidationError, match="must be >= 0.750"):
        KinaseWorkflowValidator().run(request)
