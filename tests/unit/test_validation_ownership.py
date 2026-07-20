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
from phospy.errors import ReferenceCompatibilityError
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.science.references.models import Organism, ReferenceBundle, ReferencePreset
from phospy.science.references.resolution import ReferenceResolver
from phospy.science.sites.site_keys import (
    build_protein_scoped_site_key,
    encode_site_key,
)
from phospy.validation.datasets.analysis_ready import (
    AnalysisReadyDatasetModelBoundaryValidator,
)
from phospy.validation.ownership import VALIDATION_RULE_OWNERS
from phospy.validation.references.bundle import ReferenceBundleValidator
from phospy.validation.references.compatibility import ReferenceCompatibilityValidator
from phospy.workflows.kinase.public import KinaseWorkflow
from phospy.workflows.kinase.validator import KinaseWorkflowValidator
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
)
from tests.support.site_keys import site_key_context_columns


def _site_key() -> str:
    return encode_site_key(
        build_protein_scoped_site_key(
            organism="rat",
            protein_namespace="protein_id",
            protein_identifier="P28482",
            residue="Y",
            position=182,
            field_name="tests.unit.test_validation_ownership.site_key",
            error_type=ValueError,
        )
    )


def _dataset() -> AnalysisReadyPhosphoDataset:
    index = pd.Index([_site_key()], name="site_key")
    return AnalysisReadyPhosphoDataset(
        phospho=pd.DataFrame({"sample_a": [1.0]}, index=index),
        site_metadata=pd.DataFrame(
            {
                "site_key": index.tolist(),
                "display_id": ["MAPK14;Y182;"],
                **site_key_context_columns(index),
                "gene_symbol": ["MAPK14"],
                "site": ["Y182"],
                "site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"],
            },
            index=index,
        ),
        organism=Organism.RAT,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
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


def test_reference_compatibility_is_enforced_at_workflow_runtime_boundary() -> None:
    request = KinaseWorkflowRequest(
        dataset=_dataset(),
        references=ReferencePreset.HUMAN,
        scoring_config=KinaseScoringConfig(min_substrates=2),
        prediction_config=KinasePredictionConfig(
            top_k=5,
            deterministic_max_selected_kinases=5,
            adaptive_ensemble_runs=5,
        ),
        activity_config=None,
    )
    assert KinaseWorkflowValidator().run(request) is request
    with pytest.raises(ReferenceCompatibilityError):
        KinaseWorkflow().run(request)


def test_reference_resolver_delegates_compatibility_and_does_not_duplicate_checks() -> (
    None
):
    resolver_source = inspect.getsource(ReferenceResolver.run)
    assert "self._compatibility_validator.run(" in resolver_source
    assert "resolve_preset_organism" in resolver_source
    assert "run_bundle_organism" not in resolver_source
    assert "ReferenceBundleValidator" not in resolver_source
    assert "requested reference preset must match" not in resolver_source
    assert "ReferencePreset.AUTO requires dataset.organism" not in resolver_source


def test_reference_compatibility_validator_is_single_owner_for_compatibility_rules() -> (
    None
):
    source = inspect.getsource(ReferenceCompatibilityValidator)
    assert "requested reference preset must match" in source
    assert "ReferencePreset.AUTO requires dataset.organism" in source
    assert "references.organism must match dataset.organism" in source


def test_reference_bundle_contract_validation_has_single_owner() -> None:
    bundle_source = inspect.getsource(ReferenceBundle._init_reference_bundle)
    validator_source = inspect.getsource(ReferenceBundleValidator.run)
    resolver_source = inspect.getsource(ReferenceResolver.run)
    assert "ReferenceBundleValidator().run(" in bundle_source
    assert "KinaseSubstrateReference(" not in bundle_source
    assert "SiteSequenceReference(" not in bundle_source
    assert "KinaseSubstrateReference(" in validator_source
    assert "SiteSequenceReference(" in validator_source
    assert "ReferenceBundleValidator" not in resolver_source


def test_dataset_validation_composition_is_outside_validation_subdomains() -> None:
    dataset_validator_source = inspect.getsource(
        AnalysisReadyDatasetModelBoundaryValidator
    )
    dataset_constructor_source = inspect.getsource(
        AnalysisReadyPhosphoDataset.__init__
    ) + inspect.getsource(AnalysisReadyPhosphoDataset._init_analysis_ready_tables)
    assert "AnalysisReadyPhosphoDataset(" in dataset_validator_source
    assert "PhosphoIntensityMatrix(" not in dataset_validator_source
    assert "SiteMetadataTable(" not in dataset_validator_source
    assert "TotalProteinMatrix(" not in dataset_validator_source
    assert (
        "intensity_scale_state"
        in inspect.signature(AnalysisReadyDatasetModelBoundaryValidator.run).parameters
    )
    assert "PhosphoIntensityMatrix(" in dataset_constructor_source
    assert "SiteMetadataTable(" in dataset_constructor_source
    assert "SampleMetadataTable(" in dataset_constructor_source
    assert "TotalProteinMatrix(" in dataset_constructor_source
    assert "_INTENSITY_SCALE_STATE_VALIDATOR.run(" in dataset_constructor_source


def test_major_validation_rules_have_documented_owners() -> None:
    documented = {entry.rule: entry.owner for entry in VALIDATION_RULE_OWNERS}
    assert len(documented) == len(VALIDATION_RULE_OWNERS)
    assert documented["dataset build request input source types"]
    assert (
        documented["dataset build preprocessing config policy"]
        == "DatasetPreprocessingConfigValidator.run"
    )
    assert (
        documented["kinase workflow request config policy"]
        == "KinaseWorkflowConfigValidator.run"
    )
    assert (
        documented["signalome workflow request config policy"]
        == "SignalomeConfigValidator.run"
    )
    assert (
        documented["reference input compatibility (preset/bundle vs dataset organism)"]
        == "ReferenceCompatibilityValidator.run"
    )
    assert documented["reference bundle structural contract"]
    assert (
        documented["analysis-ready dataset structural contract"]
        == "AnalysisReadyPhosphoDataset.__init__"
    )
    assert documented["dataset/intensity-scale-state coherence"]
    assert (
        documented["sequence-aware workflow centred site-sequence context"]
        == "enforce_centred_site_sequence_context"
    )
    assert documented["signalome result expanded_signalome field type/ownership"]
