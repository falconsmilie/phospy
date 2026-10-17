from __future__ import annotations

from types import SimpleNamespace

import pytest

from phospy.datasets.processing_state import DatasetPreprocessingReport
from phospy.errors.workflows import PhosPyWorkflowError
from phospy.kinase.eligibility import KinaseEligibilityReportComposer
from phospy.kinase.provenance import KinaseProvenanceBuilder
from phospy.preparation_records.report_schema import PreprocessingOperationRow
from phospy.provenance.scientific_policy_models import (
    ScientificPolicyId,
    ScientificPolicyRecord,
)


def _request_with_preprocessing_report(
    report: DatasetPreprocessingReport,
) -> SimpleNamespace:
    return SimpleNamespace(dataset=SimpleNamespace(preprocessing_report=report))


def _duplicate_policy(
    *, value: str = "mean", description: str = "attached"
) -> ScientificPolicyRecord:
    return ScientificPolicyRecord(
        id=ScientificPolicyId.DUPLICATE_SITE_RESOLUTION,
        name="dataset_duplicate_site_resolution",
        version="custom-version",
        description=description,
        parameters={"duplicate_site_policy": value, "evidence": "preserve-me"},
        assumptions=("dataset-owned",),
        output_scale="dataset sites",
        quantitative_meaning="site abundance",
    )


def _request_with_policy_records(
    *policies: ScientificPolicyRecord,
    plan_value: str = "mean",
) -> SimpleNamespace:
    provenance = SimpleNamespace(
        workflow_parameters={
            "preprocessing_plan": {
                "stage_order": ["missing_data", "site_matrix"],
                "site_matrix_duplicate_site_policy": plan_value,
            }
        },
        scientific_policies=policies,
    )
    return SimpleNamespace(dataset=SimpleNamespace(provenance=provenance))


def test_localisation_counts_use_dataset_owned_last_operation_semantics() -> None:
    report = DatasetPreprocessingReport.from_rows(
        operation_rows=(
            PreprocessingOperationRow(
                step_order=1,
                stage="localisation_confidence",
                operation="filter",
                parameters={},
                input_rows=12,
                output_rows=10,
                notes="first",
            ),
            PreprocessingOperationRow(
                step_order=2,
                stage="missing_data",
                operation="forbid",
                parameters={},
                input_rows=10,
                output_rows=9,
                notes="middle",
            ),
            PreprocessingOperationRow(
                step_order=3,
                stage="localisation_confidence",
                operation="filter",
                parameters={},
                input_rows=8,
                output_rows=9,
                notes="last",
            ),
        )
    )

    counts = KinaseEligibilityReportComposer._resolve_localisation_counts(
        _request_with_preprocessing_report(report),  # type: ignore[arg-type]
        final_dataset_sites=7,
    )

    assert counts == (7, 0)


def test_localisation_counts_remain_absent_without_durable_evidence() -> None:
    report = DatasetPreprocessingReport.from_rows(
        operation_rows=(
            PreprocessingOperationRow(
                step_order=1,
                stage="missing_data",
                operation="forbid",
                parameters={},
                input_rows=5,
                output_rows=4,
                notes="unrelated",
            ),
        )
    )

    counts = KinaseEligibilityReportComposer._resolve_localisation_counts(
        _request_with_preprocessing_report(report),  # type: ignore[arg-type]
        final_dataset_sites=4,
    )

    assert counts == (None, None)


def test_duplicate_site_policy_reuses_exact_dataset_record_and_payload() -> None:
    attached = _duplicate_policy(description="non-default attached payload")

    resolved = KinaseProvenanceBuilder._resolve_duplicate_site_resolution_policy(
        request=_request_with_policy_records(attached),  # type: ignore[arg-type]
    )

    assert resolved is attached
    assert resolved.to_payload() == attached.to_payload()


@pytest.mark.parametrize(
    ("policies", "plan_value", "message"),
    [
        ((), "mean", "observed 0"),
        (
            (_duplicate_policy(), _duplicate_policy(description="duplicate")),
            "mean",
            "observed 2",
        ),
        ((_duplicate_policy(value="median"),), "mean", "disagree"),
        (
            (
                ScientificPolicyRecord(
                    id=ScientificPolicyId.CANDIDATE_SUBSTRATE_SELECTION,
                    name="unrelated",
                    version="1",
                    description="unrelated",
                    parameters={"duplicate_site_policy": "mean"},
                    assumptions=(),
                ),
            ),
            "mean",
            "observed 0",
        ),
    ],
)
def test_duplicate_site_policy_fails_closed_for_missing_or_conflicting_evidence(
    policies: tuple[ScientificPolicyRecord, ...],
    plan_value: str,
    message: str,
) -> None:
    with pytest.raises(PhosPyWorkflowError, match=message) as exc_info:
        KinaseProvenanceBuilder._resolve_duplicate_site_resolution_policy(
            request=_request_with_policy_records(  # type: ignore[arg-type]
                *policies,
                plan_value=plan_value,
            )
        )

    assert "seam=kinase.provenance.duplicate_site_resolution_policy" in str(
        exc_info.value
    )
