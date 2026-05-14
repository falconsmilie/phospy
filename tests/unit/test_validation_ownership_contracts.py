from __future__ import annotations

import inspect

from phospy.validation.datasets import site_metadata as site_metadata_validation
from phospy.workflows.kinase.validator import KinaseWorkflowValidator
from phospy.workflows.signalome.validator import SignalomeWorkflowValidator


def test_workflow_validators_compose_shared_site_identity_validator() -> None:
    shared_source = inspect.getsource(
        site_metadata_validation.enforce_site_identity_rows
    )
    kinase_source = inspect.getsource(KinaseWorkflowValidator.run)
    signalome_source = inspect.getsource(
        SignalomeWorkflowValidator._require_explicit_site_metadata_protein_identity
    )

    assert "build_phosphosite_identity(" in shared_source
    assert "enforce_site_identity_rows(" in kinase_source
    assert "enforce_site_identity_rows(" in signalome_source
    assert "build_phosphosite_identity(" not in kinase_source
    assert "build_phosphosite_identity(" not in signalome_source
    assert "allow_opaque_site_values=True" not in kinase_source
    assert "allow_opaque_site_values=True" not in signalome_source


def test_sequence_aware_workflow_validators_compose_shared_centred_context_validator() -> (
    None
):
    shared_source = inspect.getsource(
        site_metadata_validation.enforce_centred_site_sequence_context
    )
    kinase_source = inspect.getsource(KinaseWorkflowValidator.run)
    signalome_source = inspect.getsource(
        SignalomeWorkflowValidator._require_explicit_site_metadata_protein_identity
    )

    assert "site_sequence_column" in shared_source
    assert "requires centred sequence context" in shared_source
    assert "enforce_centred_site_sequence_context(" in kinase_source
    assert "enforce_centred_site_sequence_context(" in signalome_source
