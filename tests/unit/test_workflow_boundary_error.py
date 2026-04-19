from __future__ import annotations

from phospy.errors.workflows import WorkflowBoundaryError


def test_workflow_boundary_error_supports_message_only_construction() -> None:
    error = WorkflowBoundaryError("workflow boundary message")

    assert str(error) == "workflow boundary message"
    assert error.seam is None
    assert error.next_action is None
    assert error.details == {}


def test_workflow_boundary_error_exposes_structured_diagnostics() -> None:
    error = WorkflowBoundaryError(
        seam="kinase.interpreter.reference_coverage",
        next_action="use compatible references",
        details={"dataset_sites": 2, "overlap_sites": 0},
        message_prefix="kinase workflow boundary validation failed",
    )

    assert error.seam == "kinase.interpreter.reference_coverage"
    assert error.next_action == "use compatible references"
    assert error.details == {"dataset_sites": 2, "overlap_sites": 0}
    assert "seam=kinase.interpreter.reference_coverage" in str(error)
    assert "dataset_sites=2" in str(error)
    assert "next_action=use compatible references" in str(error)


def test_workflow_boundary_error_copies_details_mapping() -> None:
    details = {"shared_sites": 0}
    error = WorkflowBoundaryError(
        seam="signalome.interpreter.site_alignment",
        next_action="align score and prediction matrices",
        details=details,
    )

    details["shared_sites"] = 99

    assert error.details == {"shared_sites": 0}
