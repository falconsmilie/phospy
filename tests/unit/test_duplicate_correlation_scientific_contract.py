from __future__ import annotations

import pytest

from phospy.errors import PhosPyInputError
from phospy.provenance.models import TableFingerprint
from phospy.science.differential.models.duplicate_correlation import (
    DUPLICATE_CORRELATION_METHOD_REML_FISHER_TRIMMED_MEAN,
    DUPLICATE_CORRELATION_TRIM_FRACTION,
    DuplicateCorrelationBlockStructureSummary,
    DuplicateCorrelationConsensusResult,
    DuplicateCorrelationConsensusSummary,
    DuplicateCorrelationFailureReason,
    DuplicateCorrelationFeatureEstimate,
    DuplicateCorrelationFeatureStatus,
    DuplicateCorrelationWorkflowProvenance,
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


def _block_summary() -> DuplicateCorrelationBlockStructureSummary:
    return DuplicateCorrelationBlockStructureSummary(
        block_id_field_name="block_id",
        sample_count=5,
        block_count=3,
        repeated_block_count=2,
        singleton_block_count=1,
        correlated_pair_count=2,
        block_levels=("donor_1", "donor_2", "donor_3"),
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


def test_duplicate_correlation_internal_models_accept_valid_contract() -> None:
    provenance = DuplicateCorrelationWorkflowProvenance(
        model="duplicate_correlation",
        matrix_authority="workflow approved differential fitting matrix",
        authoritative_matrix_fingerprint=_matrix_fingerprint(),
        design_authority="experimental-design interpreter",
        block_authority="workflow validation domain",
        estimator_authority="science.differential duplicate-correlation estimator",
        gls_authority="science.differential GLS fitter",
        failure_authority="workflow validation plus estimator typed failures",
        block_structure=_block_summary(),
        consensus=_success_consensus_summary(),
        imputed_values_participated=True,
        imputed_feature_count=2,
        imputed_cell_count=3,
    )

    payload = provenance.to_payload()

    assert provenance.consensus.consensus_correlation == 0.25
    assert payload["authoritative_matrix_fingerprint"] == (
        _matrix_fingerprint_payload()
    )
    assert "feature_estimates" not in provenance.consensus.to_payload()
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
        DuplicateCorrelationWorkflowProvenance(
            model="duplicate_correlation",
            matrix_authority="workflow approved differential fitting matrix",
            authoritative_matrix_fingerprint=_matrix_fingerprint(),
            design_authority="experimental-design interpreter",
            block_authority="workflow validation domain",
            estimator_authority="science.differential duplicate-correlation estimator",
            gls_authority="science.differential GLS fitter",
            failure_authority="workflow validation plus estimator typed failures",
            block_structure=_block_summary(),
            consensus=_success_consensus_result(),  # type: ignore[arg-type]
            imputed_values_participated=False,
        )


def test_workflow_provenance_requires_authoritative_matrix_fingerprint() -> None:
    with pytest.raises(PhosPyInputError, match="authoritative_matrix_fingerprint"):
        DuplicateCorrelationWorkflowProvenance(
            model="duplicate_correlation",
            matrix_authority="workflow approved differential fitting matrix",
            authoritative_matrix_fingerprint="not-a-fingerprint",  # type: ignore[arg-type]
            design_authority="experimental-design interpreter",
            block_authority="workflow validation domain",
            estimator_authority="science.differential duplicate-correlation estimator",
            gls_authority="science.differential GLS fitter",
            failure_authority="workflow validation plus estimator typed failures",
            block_structure=_block_summary(),
            consensus=_success_consensus_summary(),
            imputed_values_participated=False,
        )


def test_workflow_provenance_rejects_contradictory_imputation_counts() -> None:
    with pytest.raises(PhosPyInputError, match="imputed counts"):
        DuplicateCorrelationWorkflowProvenance(
            model="duplicate_correlation",
            matrix_authority="workflow approved differential fitting matrix",
            authoritative_matrix_fingerprint=_matrix_fingerprint(),
            design_authority="experimental-design interpreter",
            block_authority="workflow validation domain",
            estimator_authority="science.differential duplicate-correlation estimator",
            gls_authority="science.differential GLS fitter",
            failure_authority="workflow validation plus estimator typed failures",
            block_structure=_block_summary(),
            consensus=_success_consensus_summary(),
            imputed_values_participated=False,
            imputed_cell_count=1,
        )


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
