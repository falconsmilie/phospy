from __future__ import annotations

import dataclasses
import inspect

from phospy.workflows.kinase import provenance
from phospy.workflows.kinase.contracts import (
    ResolvedKinaseSiteUniverses,
    ResolvedKinaseWorkflowRequest,
)
from phospy.workflows.kinase.reference_projection import (
    KinaseReferenceProjectionSummary,
)


def test_reference_projection_summary_is_typed_resolved_request_state() -> None:
    fields = {
        field.name: field.type
        for field in dataclasses.fields(ResolvedKinaseWorkflowRequest)
    }

    assert fields["reference_projection_summary"] in {
        KinaseReferenceProjectionSummary | None,
        "KinaseReferenceProjectionSummary | None",
    }


def test_site_universes_do_not_accept_source_reference_identifier_namespace() -> None:
    field_names = {
        field.name for field in dataclasses.fields(ResolvedKinaseSiteUniverses)
    }

    assert "reference_supported_membership_sites" in field_names
    assert "source_reference_substrate_identifiers" not in field_names
    assert "reference_projection_summary" not in field_names


def test_reference_projection_attrition_is_not_generic_membership_attrition() -> None:
    source = inspect.getsource(provenance._reference_attrition_records)

    assert '"attrition_type": "reference_attrition"' in source
    assert "source_reference_substrate_identifiers" in source
    assert "membership_attrition" not in source
    assert "reference_supported_membership_sites" not in source
