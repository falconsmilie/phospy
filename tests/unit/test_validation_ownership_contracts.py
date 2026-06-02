from __future__ import annotations

import inspect

from phospy.validation.datasets import (
    protein_scoped_site_identity as site_key_validation,
)
from phospy.validation.datasets import site_metadata as site_metadata_validation
from phospy.validation.workflows import identity as workflow_identity_validation
from phospy.workflows.differential.validator import DifferentialAnalysisValidator
from phospy.workflows.kinase.validator import KinaseWorkflowValidator
from phospy.workflows.signalome.validator import SignalomeWorkflowValidator


def test_workflow_validators_compose_shared_site_identity_validator() -> None:
    shared_source = inspect.getsource(
        site_metadata_validation.enforce_site_identity_rows
    )
    site_key_index_source = inspect.getsource(
        site_key_validation.enforce_analysis_ready_site_key_index
    )
    site_key_alignment_source = inspect.getsource(
        site_key_validation.enforce_site_key_column_matches_index
    )
    workflow_shared_source = inspect.getsource(
        workflow_identity_validation.enforce_workflow_site_identity_contract
    )
    differential_source = inspect.getsource(DifferentialAnalysisValidator.run)
    kinase_source = inspect.getsource(KinaseWorkflowValidator.run)
    signalome_source = inspect.getsource(
        SignalomeWorkflowValidator._require_explicit_site_metadata_protein_identity
    )

    assert "build_phosphosite_identity(" in shared_source
    assert "require_site_key_index(" in site_key_index_source
    assert "enforce_site_key_column(" in site_key_alignment_source
    assert "enforce_analysis_ready_site_key_index(" in workflow_shared_source
    assert "enforce_site_key_column_matches_index(" in workflow_shared_source
    assert "enforce_display_id_column(" in workflow_shared_source
    assert "require_exact_index_match(" in workflow_shared_source
    assert "enforce_site_identity_rows(" not in workflow_shared_source
    assert "validate_no_conflicting_identity_collisions(" not in workflow_shared_source
    assert "enforce_workflow_site_identity_contract(" in differential_source
    assert "enforce_workflow_site_identity_contract(" in kinase_source
    assert "enforce_workflow_site_identity_contract(" in signalome_source
    assert "build_phosphosite_identity(" not in kinase_source
    assert "build_phosphosite_identity(" not in differential_source
    assert "build_phosphosite_identity(" not in signalome_source
    assert "allow_opaque_site_values=True" not in kinase_source
    assert "allow_opaque_site_values=True" not in signalome_source


def test_sequence_aware_workflow_validators_compose_shared_centred_context_validator() -> (
    None
):
    shared_source = inspect.getsource(
        site_metadata_validation.enforce_centred_site_sequence_context
    )
    workflow_shared_source = inspect.getsource(
        workflow_identity_validation.enforce_workflow_site_identity_contract
    )
    kinase_source = inspect.getsource(KinaseWorkflowValidator.run)
    signalome_source = inspect.getsource(
        SignalomeWorkflowValidator._require_explicit_site_metadata_protein_identity
    )

    assert "site_sequence_column" in shared_source
    assert "requires centred sequence context" in shared_source
    assert "enforce_centred_site_sequence_context(" in workflow_shared_source
    assert "enforce_centred_site_sequence_context(" not in kinase_source
    assert "enforce_centred_site_sequence_context(" not in signalome_source
