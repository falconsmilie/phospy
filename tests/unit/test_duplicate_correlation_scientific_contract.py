from __future__ import annotations

from dataclasses import FrozenInstanceError

import pytest

from phospy.errors import PhosPyInputError
from phospy.provenance.models import TableFingerprint
from phospy.science.differential.models.duplicate_correlation import (
    DUPLICATE_CORRELATION_BLOCK_TREATMENT_CONSENSUS_CORRELATION,
    DUPLICATE_CORRELATION_COVARIANCE_STRUCTURE_COMPOUND_SYMMETRY,
    DUPLICATE_CORRELATION_ESTIMATOR_FEATURE_WISE_REML,
    DUPLICATE_CORRELATION_ESTIMATOR_POLICY_VERSION,
    DUPLICATE_CORRELATION_GLS_FIT_STATUS_FIT,
    DUPLICATE_CORRELATION_METHOD_REML_FISHER_TRIMMED_MEAN,
    DUPLICATE_CORRELATION_TRIM_FRACTION,
    DUPLICATE_CORRELATION_WORKFLOW_PROVENANCE_VERSION,
    DuplicateCorrelationBlockStructureSummary,
    DuplicateCorrelationBoundarySummary,
    DuplicateCorrelationConsensusResult,
    DuplicateCorrelationConsensusSummary,
    DuplicateCorrelationConvergenceSummary,
    DuplicateCorrelationFailureReason,
    DuplicateCorrelationFeatureEstimate,
    DuplicateCorrelationFeatureStatus,
    DuplicateCorrelationReasonCount,
    DuplicateCorrelationWorkflowProvenance,
    duplicate_correlation_workflow_provenance_from_payload,
)


def _matrix_fingerprint() -> TableFingerprint:
    return TableFingerprint(
        name="differential.approved_model_matrix",
        rows=2,
        columns=4,
        index_name="site_key",
        column_names=("A_1", "A_2", "B_1", "B_2"),
        dtypes=("float64", "float64", "float64", "float64"),
        exact_hash_algorithm="sha256",
        exact_hash_value="abc123",
        tolerance_hash_algorithm="sha256-tolerance-v1",
        tolerance_hash_value="def456",
    )


def _design_fingerprint() -> TableFingerprint:
    return TableFingerprint(
        name="differential.non_block_fixed_effect_design",
        rows=5,
        columns=2,
        index_name="sample",
        column_names=("A", "B"),
        dtypes=("float64", "float64"),
        exact_hash_algorithm="sha256",
        exact_hash_value="design123",
        tolerance_hash_algorithm="sha256-tolerance-v1",
        tolerance_hash_value="design456",
    )


def _block_assignment_fingerprint() -> TableFingerprint:
    return TableFingerprint(
        name="differential.duplicate_correlation_block_assignment",
        rows=5,
        columns=2,
        index_name="sample",
        column_names=("sample_id", "block_id"),
        dtypes=("object", "object"),
        exact_hash_algorithm="sha256",
        exact_hash_value="block123",
        tolerance_hash_algorithm="sha256-tolerance-v1",
        tolerance_hash_value="block456",
    )


def _block_summary() -> DuplicateCorrelationBlockStructureSummary:
    return DuplicateCorrelationBlockStructureSummary(
        block_id_field_name="block_id",
        sample_count=5,
        block_count=3,
        repeated_block_count=2,
        singleton_block_count=1,
        correlated_pair_count=2,
        block_levels=("donor_1", "donor_2", "donor_3"),
        minimum_block_size=1,
        maximum_block_size=2,
    )


def _success_consensus_result() -> DuplicateCorrelationConsensusResult:
    return DuplicateCorrelationConsensusResult(
        method=DUPLICATE_CORRELATION_METHOD_REML_FISHER_TRIMMED_MEAN,
        trim_fraction=DUPLICATE_CORRELATION_TRIM_FRACTION,
        success=True,
        consensus_correlation=0.25,
        eligible_feature_count=2,
        estimated_feature_count=1,
        failed_feature_count=1,
        non_finite_feature_count=1,
        feature_estimates=(
            DuplicateCorrelationFeatureEstimate(
                feature_id="site_1",
                status=DuplicateCorrelationFeatureStatus.ESTIMATED,
                correlation=0.25,
            ),
            DuplicateCorrelationFeatureEstimate(
                feature_id="site_2",
                status=DuplicateCorrelationFeatureStatus.NON_FINITE,
                failure_reason=(
                    DuplicateCorrelationFailureReason.NUMERICAL_NON_CONVERGENCE
                ),
            ),
        ),
    )


def _success_consensus_summary() -> DuplicateCorrelationConsensusSummary:
    return _success_consensus_result().to_summary()


def _failure_reason_counts() -> tuple[DuplicateCorrelationReasonCount, ...]:
    return (
        DuplicateCorrelationReasonCount(
            reason=DuplicateCorrelationFailureReason.NUMERICAL_NON_CONVERGENCE,
            count=1,
        ),
    )


def _convergence_summary() -> DuplicateCorrelationConvergenceSummary:
    return DuplicateCorrelationConvergenceSummary(
        converged_feature_count=1,
        boundary_feature_count=0,
        failed_optimisation_feature_count=0,
        non_finite_objective_or_estimate_feature_count=1,
    )


def _boundary_summary() -> DuplicateCorrelationBoundarySummary:
    return DuplicateCorrelationBoundarySummary(
        lower_correlation_bound=-0.9999,
        upper_correlation_bound=0.9999,
        lower_boundary_feature_count=0,
        upper_boundary_feature_count=0,
        positive_definite_tolerance=1.0e-10,
        fisher_boundary_tolerance=1.0e-10,
    )


def _workflow_provenance(**overrides: object) -> DuplicateCorrelationWorkflowProvenance:
    values: dict[str, object] = {
        "model": "duplicate_correlation",
        "provenance_version": DUPLICATE_CORRELATION_WORKFLOW_PROVENANCE_VERSION,
        "requested_paired_design_policy": "duplicate_correlation",
        "normalised_paired_design_policy": "duplicate_correlation",
        "block_treatment": (
            DUPLICATE_CORRELATION_BLOCK_TREATMENT_CONSENSUS_CORRELATION
        ),
        "covariance_structure": (
            DUPLICATE_CORRELATION_COVARIANCE_STRUCTURE_COMPOUND_SYMMETRY
        ),
        "estimator": DUPLICATE_CORRELATION_ESTIMATOR_FEATURE_WISE_REML,
        "estimator_policy_version": DUPLICATE_CORRELATION_ESTIMATOR_POLICY_VERSION,
        "trim_fraction": DUPLICATE_CORRELATION_TRIM_FRACTION,
        "matrix_authority": "workflow approved differential fitting matrix",
        "analysis_matrix_fingerprint": _matrix_fingerprint(),
        "authoritative_matrix_fingerprint": _matrix_fingerprint(),
        "design_authority": "experimental-design interpreter",
        "design_fingerprint": _design_fingerprint(),
        "block_authority": "workflow validation domain",
        "block_assignment_fingerprint": _block_assignment_fingerprint(),
        "estimator_authority": "science.differential duplicate-correlation estimator",
        "gls_authority": "science.differential GLS fitter",
        "failure_authority": "workflow validation plus estimator typed failures",
        "block_structure": _block_summary(),
        "consensus": _success_consensus_summary(),
        "attempted_feature_count": 2,
        "trimmed_feature_count_each_tail": 0,
        "retained_feature_count_after_trimming": 1,
        "failure_reason_counts": _failure_reason_counts(),
        "convergence_summary": _convergence_summary(),
        "boundary_summary": _boundary_summary(),
        "sample_count": 5,
        "block_count": 3,
        "repeated_block_count": 2,
        "singleton_block_count": 1,
        "minimum_block_size": 1,
        "maximum_block_size": 2,
        "design_rank": 2,
        "gls_fit_status": DUPLICATE_CORRELATION_GLS_FIT_STATUS_FIT,
        "imputed_values_participated": True,
        "imputed_feature_count": 2,
        "imputed_cell_count": 3,
    }
    values.update(overrides)
    return DuplicateCorrelationWorkflowProvenance(**values)  # type: ignore[arg-type]


def test_duplicate_correlation_internal_models_accept_valid_contract() -> None:
    provenance = _workflow_provenance()

    payload = provenance.to_payload()

    assert provenance.requested_paired_design_policy == "duplicate_correlation"
    assert provenance.normalised_paired_design_policy == "duplicate_correlation"
    assert provenance.covariance_structure == "compound_symmetry"
    assert provenance.estimator == "feature-wise REML"
    assert provenance.consensus.consensus_correlation == 0.25
    assert payload["authoritative_matrix_fingerprint"] == (
        _matrix_fingerprint_payload()
    )
    assert payload["analysis_matrix_fingerprint"] == _matrix_fingerprint_payload()
    assert "feature_estimates" not in provenance.consensus.to_payload()
    assert "feature_estimates" not in payload["consensus"]


def test_workflow_provenance_is_frozen() -> None:
    provenance = _workflow_provenance()

    with pytest.raises(FrozenInstanceError):
        provenance.gls_fit_status = "failed"  # type: ignore[misc]


def test_workflow_provenance_payload_round_trips_without_feature_estimates() -> None:
    provenance = _workflow_provenance()
    payload = provenance.to_payload()

    restored = duplicate_correlation_workflow_provenance_from_payload(payload)

    assert restored == provenance
    assert restored.to_payload() == payload
    assert "feature_estimates" not in payload["consensus"]


@pytest.mark.parametrize(
    "kwargs",
    (
        {"sample_count": 0},
        {"sample_count": -1},
        {"block_count": 0},
        {"block_count": 4},
        {"repeated_block_count": 2, "singleton_block_count": 2},
        {
            "repeated_block_count": 0,
            "singleton_block_count": 3,
            "correlated_pair_count": 1,
        },
        {"repeated_block_count": 2, "correlated_pair_count": 1},
        {"sample_count": 4, "repeated_block_count": 2, "singleton_block_count": 1},
    ),
)
def test_block_structure_rejects_impossible_counts(kwargs: dict[str, int]) -> None:
    values = {
        "block_id_field_name": "block_id",
        "sample_count": 5,
        "block_count": 3,
        "repeated_block_count": 1,
        "singleton_block_count": 2,
        "correlated_pair_count": 1,
        "block_levels": ("a", "b", "c"),
    }
    values.update(kwargs)

    with pytest.raises(PhosPyInputError):
        DuplicateCorrelationBlockStructureSummary(**values)  # type: ignore[arg-type]


@pytest.mark.parametrize(
    "block_levels",
    (
        ("a", "", "c"),
        ("a", "b", "a"),
    ),
)
def test_block_structure_rejects_invalid_block_levels(
    block_levels: tuple[str, ...],
) -> None:
    with pytest.raises(PhosPyInputError):
        DuplicateCorrelationBlockStructureSummary(
            block_id_field_name="block_id",
            sample_count=5,
            block_count=3,
            repeated_block_count=1,
            singleton_block_count=2,
            correlated_pair_count=1,
            block_levels=block_levels,
        )


@pytest.mark.parametrize("trim_fraction", (-0.1, 0.2, 0.5, float("nan")))
def test_consensus_rejects_non_canonical_trim_values(trim_fraction: float) -> None:
    with pytest.raises(PhosPyInputError, match="trim_fraction"):
        DuplicateCorrelationConsensusResult(
            method=DUPLICATE_CORRELATION_METHOD_REML_FISHER_TRIMMED_MEAN,
            trim_fraction=trim_fraction,
            success=True,
            consensus_correlation=0.1,
            eligible_feature_count=1,
            estimated_feature_count=1,
            failed_feature_count=0,
            non_finite_feature_count=0,
        )


@pytest.mark.parametrize("correlation", (-1.0, 1.0, float("inf")))
def test_feature_estimate_rejects_invalid_correlations(correlation: float) -> None:
    with pytest.raises(PhosPyInputError, match="correlation"):
        DuplicateCorrelationFeatureEstimate(
            feature_id="site_1",
            status=DuplicateCorrelationFeatureStatus.ESTIMATED,
            correlation=correlation,
        )


@pytest.mark.parametrize("method", ("", "feature_reml_untrimmed_mean"))
def test_consensus_rejects_invalid_method_identifier(method: str) -> None:
    with pytest.raises(PhosPyInputError, match="method"):
        DuplicateCorrelationConsensusResult(
            method=method,
            trim_fraction=DUPLICATE_CORRELATION_TRIM_FRACTION,
            success=True,
            consensus_correlation=0.1,
            eligible_feature_count=1,
            estimated_feature_count=1,
            failed_feature_count=0,
            non_finite_feature_count=0,
        )


def test_consensus_rejects_invalid_failure_reason() -> None:
    with pytest.raises(PhosPyInputError, match="failure_reason"):
        DuplicateCorrelationConsensusResult(
            method=DUPLICATE_CORRELATION_METHOD_REML_FISHER_TRIMMED_MEAN,
            trim_fraction=DUPLICATE_CORRELATION_TRIM_FRACTION,
            success=False,
            consensus_correlation=None,
            eligible_feature_count=1,
            estimated_feature_count=0,
            failed_feature_count=1,
            non_finite_feature_count=0,
            failure_reason="not_a_supported_reason",  # type: ignore[arg-type]
        )


def test_consensus_rejects_contradictory_success_failure_state() -> None:
    with pytest.raises(PhosPyInputError, match="failure reason"):
        DuplicateCorrelationConsensusResult(
            method=DUPLICATE_CORRELATION_METHOD_REML_FISHER_TRIMMED_MEAN,
            trim_fraction=DUPLICATE_CORRELATION_TRIM_FRACTION,
            success=True,
            consensus_correlation=0.1,
            eligible_feature_count=1,
            estimated_feature_count=1,
            failed_feature_count=0,
            non_finite_feature_count=0,
            failure_reason=DuplicateCorrelationFailureReason.NUMERICAL_NON_CONVERGENCE,
        )

    with pytest.raises(PhosPyInputError, match="consensus correlation"):
        DuplicateCorrelationConsensusResult(
            method=DUPLICATE_CORRELATION_METHOD_REML_FISHER_TRIMMED_MEAN,
            trim_fraction=DUPLICATE_CORRELATION_TRIM_FRACTION,
            success=False,
            consensus_correlation=0.1,
            eligible_feature_count=1,
            estimated_feature_count=0,
            failed_feature_count=1,
            non_finite_feature_count=0,
            failure_reason=DuplicateCorrelationFailureReason.NUMERICAL_NON_CONVERGENCE,
        )


def test_consensus_rejects_inconsistent_feature_counts() -> None:
    with pytest.raises(PhosPyInputError, match="must equal eligible_feature_count"):
        DuplicateCorrelationConsensusResult(
            method=DUPLICATE_CORRELATION_METHOD_REML_FISHER_TRIMMED_MEAN,
            trim_fraction=DUPLICATE_CORRELATION_TRIM_FRACTION,
            success=True,
            consensus_correlation=0.1,
            eligible_feature_count=3,
            estimated_feature_count=1,
            failed_feature_count=1,
            non_finite_feature_count=0,
        )


def test_retained_feature_estimates_must_match_summary_counts() -> None:
    with pytest.raises(PhosPyInputError, match="feature_estimates length"):
        DuplicateCorrelationConsensusResult(
            method=DUPLICATE_CORRELATION_METHOD_REML_FISHER_TRIMMED_MEAN,
            trim_fraction=DUPLICATE_CORRELATION_TRIM_FRACTION,
            success=True,
            consensus_correlation=0.1,
            eligible_feature_count=2,
            estimated_feature_count=1,
            failed_feature_count=1,
            non_finite_feature_count=0,
            feature_estimates=(
                DuplicateCorrelationFeatureEstimate(
                    feature_id="site_1",
                    status=DuplicateCorrelationFeatureStatus.ESTIMATED,
                    correlation=0.1,
                ),
            ),
        )


def test_workflow_provenance_rejects_feature_retaining_consensus_result() -> None:
    with pytest.raises(PhosPyInputError, match="DuplicateCorrelationConsensusSummary"):
        _workflow_provenance(
            consensus=_success_consensus_result(),  # type: ignore[arg-type]
        )


def test_workflow_provenance_requires_authoritative_matrix_fingerprint() -> None:
    with pytest.raises(PhosPyInputError, match="authoritative_matrix_fingerprint"):
        _workflow_provenance(
            authoritative_matrix_fingerprint="not-a-fingerprint",  # type: ignore[arg-type]
        )


def test_workflow_provenance_rejects_contradictory_imputation_counts() -> None:
    with pytest.raises(PhosPyInputError, match="imputed counts"):
        _workflow_provenance(
            imputed_values_participated=False,
            imputed_feature_count=0,
            imputed_cell_count=1,
        )


def test_workflow_provenance_reconstruction_rejects_tampered_policy() -> None:
    payload = _workflow_provenance().to_payload()
    payload["normalised_paired_design_policy"] = "fixed_block"

    with pytest.raises(PhosPyInputError, match="normalised_paired_design_policy"):
        duplicate_correlation_workflow_provenance_from_payload(payload)


def test_workflow_provenance_reconstruction_rejects_tampered_fingerprint() -> None:
    payload = _workflow_provenance().to_payload()
    fingerprint_payload = dict(_matrix_fingerprint_payload())
    fingerprint_payload["exact_hash_value"] = "tampered"
    payload["authoritative_matrix_fingerprint"] = fingerprint_payload

    with pytest.raises(PhosPyInputError, match="analysis_matrix_fingerprint"):
        duplicate_correlation_workflow_provenance_from_payload(payload)


def test_workflow_provenance_reconstruction_rejects_failure_count_mismatch() -> None:
    payload = _workflow_provenance().to_payload()
    payload["failure_reason_counts"] = []

    with pytest.raises(PhosPyInputError, match="failure_reason_counts"):
        duplicate_correlation_workflow_provenance_from_payload(payload)


def test_workflow_provenance_reconstruction_rejects_gls_failure_status() -> None:
    payload = _workflow_provenance().to_payload()
    payload["gls_fit_status"] = "failed"

    with pytest.raises(PhosPyInputError, match="gls_fit_status"):
        duplicate_correlation_workflow_provenance_from_payload(payload)


def test_feature_level_estimator_failure_is_not_final_row_attrition() -> None:
    consensus = _success_consensus_result()

    assert consensus.eligible_feature_count == 2
    assert consensus.failed_feature_count == 1
    assert all(
        estimate.status != DuplicateCorrelationFeatureStatus.ESTIMATED
        for estimate in consensus.feature_estimates
        if estimate.failure_reason is not None
    )


def _matrix_fingerprint_payload() -> dict[str, object]:
    return {
        "name": "differential.approved_model_matrix",
        "rows": 2,
        "columns": 4,
        "index_name": "site_key",
        "column_names": ["A_1", "A_2", "B_1", "B_2"],
        "dtypes": ["float64", "float64", "float64", "float64"],
        "exact_hash_algorithm": "sha256",
        "exact_hash_value": "abc123",
        "tolerance_hash_algorithm": "sha256-tolerance-v1",
        "tolerance_hash_value": "def456",
        "index_structure": None,
        "column_index_structure": None,
    }
