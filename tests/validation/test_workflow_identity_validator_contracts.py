from __future__ import annotations

from collections.abc import Callable

import pytest

from phospy import AnalysisReadyPhosphoDataset
from phospy.errors import WorkflowValidationError
from phospy.workflows.differential.validator import DifferentialAnalysisValidator
from phospy.workflows.kinase.validator import KinaseWorkflowValidator
from phospy.workflows.signalome.validator import SignalomeWorkflowValidator
from tests.support.workflow_identity_coherence import (
    build_duplicate_display_differential_request,
    build_duplicate_display_kinase_request,
    build_duplicate_display_signalome_request,
    corrupt_dataset_to_display_index,
    drop_site_metadata_column,
    set_site_sequence_values,
)

_RequestFactory = Callable[[], object]
_DatasetGetter = Callable[[object], AnalysisReadyPhosphoDataset]
_ValidatorRunner = Callable[[object], object]


def _workflow_validator_cases() -> tuple[
    tuple[str, _ValidatorRunner, _RequestFactory, _DatasetGetter],
    ...,
]:
    return (
        (
            "differential",
            DifferentialAnalysisValidator().run,
            build_duplicate_display_differential_request,
            lambda request: request.dataset,
        ),
        (
            "kinase",
            KinaseWorkflowValidator().run,
            build_duplicate_display_kinase_request,
            lambda request: request.dataset,
        ),
        (
            "signalome",
            SignalomeWorkflowValidator().run,
            build_duplicate_display_signalome_request,
            lambda request: request.kinase_result.dataset,
        ),
    )


@pytest.mark.parametrize(
    ("workflow_name", "run_validator", "request_factory", "dataset_getter"),
    _workflow_validator_cases(),
)
@pytest.mark.parametrize(
    ("column_name", "pattern"),
    (
        ("site_key", "missing required columns: site_key"),
        ("display_id", "missing required columns: display_id"),
    ),
)
def test_workflow_validators_require_site_key_and_display_id(
    workflow_name: str,
    run_validator: _ValidatorRunner,
    request_factory: _RequestFactory,
    dataset_getter: _DatasetGetter,
    column_name: str,
    pattern: str,
) -> None:
    request = request_factory()
    dataset = dataset_getter(request)
    drop_site_metadata_column(dataset, column_name)

    with pytest.raises(WorkflowValidationError, match=pattern):
        run_validator(request)


@pytest.mark.parametrize(
    ("workflow_name", "run_validator", "request_factory", "dataset_getter"),
    _workflow_validator_cases(),
)
def test_workflow_validators_reject_display_indexed_datasets(
    workflow_name: str,
    run_validator: _ValidatorRunner,
    request_factory: _RequestFactory,
    dataset_getter: _DatasetGetter,
) -> None:
    request = request_factory()
    corrupt_dataset_to_display_index(dataset_getter(request))

    with pytest.raises(WorkflowValidationError, match="display-indexed"):
        run_validator(request)


@pytest.mark.parametrize(
    ("workflow_name", "run_validator", "request_factory", "dataset_getter"),
    (
        (
            "kinase",
            KinaseWorkflowValidator().run,
            build_duplicate_display_kinase_request,
            lambda request: request.dataset,
        ),
        (
            "signalome",
            SignalomeWorkflowValidator().run,
            build_duplicate_display_signalome_request,
            lambda request: request.kinase_result.dataset,
        ),
    ),
)
def test_kinase_and_signalome_validators_enforce_sequence_context(
    workflow_name: str,
    run_validator: _ValidatorRunner,
    request_factory: _RequestFactory,
    dataset_getter: _DatasetGetter,
) -> None:
    request = request_factory()
    set_site_sequence_values(dataset_getter(request), ["AAAAAAAA", "AAAAAAAA"])

    with pytest.raises(
        WorkflowValidationError,
        match="requires centred sequence context",
    ):
        run_validator(request)
