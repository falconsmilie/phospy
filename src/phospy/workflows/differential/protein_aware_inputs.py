"""Private protein-aware differential input resolution.

This module stops at execution-ready request construction. It does not activate
the public protein-aware workflow branch and does not fit the model.
"""

from __future__ import annotations

import math
from collections import Counter
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any, NoReturn, cast

import numpy as np
import numpy.typing as npt
import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.errors.validation import DatasetValidationError
from phospy.errors.workflows import WorkflowBoundaryError
from phospy.science.configs.differential import (
    PAIRED_DESIGN_POLICY_DUPLICATE_CORRELATION,
    SUPPORTED_DIFFERENTIAL_PROTEIN_AWARE_MODEL_METHODS,
    DifferentialProteinAwareModelMethod,
)
from phospy.science.configs.preprocessing.total_protein import (
    DATASET_PROTEIN_AWARE_PREPARATION_POLICY_PREPARE_MODEL_INPUTS,
)
from phospy.science.datasets.internal_view import (
    DatasetInternalView,
    ProteinAwarePreparationInternalView,
)
from phospy.science.datasets.preprocessing.policy_models import (
    TotalProteinCorrectionPolicy,
)
from phospy.science.datasets.preprocessing.protein_aware_alignment import (
    ProteinAwarePreparationEligibility,
)
from phospy.science.datasets.preprocessing.protein_aware_models import (
    PROTEIN_AWARE_PREPARATION_SCHEMA_VERSION,
)
from phospy.science.differential.linear_model import (
    DIFFERENTIAL_LINEAR_MODEL_MAX_CONDITION_NUMBER,
    DifferentialDesignDecompositionError,
    decompose_differential_design,
)
from phospy.science.differential.models import (
    DIFFERENTIAL_RESULT_REASON_PROTEIN_AUGMENTED_DESIGN_ILL_CONDITIONED,
    DIFFERENTIAL_RESULT_REASON_PROTEIN_AUGMENTED_DESIGN_NON_POSITIVE_RESIDUAL_DOF,
    DIFFERENTIAL_RESULT_REASON_PROTEIN_AUGMENTED_DESIGN_RANK_DEFICIENT,
    DIFFERENTIAL_RESULT_REASON_PROTEIN_CONTRAST_NON_ESTIMABLE,
    DIFFERENTIAL_RESULT_REASON_PROTEIN_COVARIATE_NON_FINITE,
    DIFFERENTIAL_RESULT_REASON_PROTEIN_COVARIATE_ZERO_VARIANCE,
    DIFFERENTIAL_RESULT_REASON_PROTEIN_PREPARATION_EXCLUDED,
    DIFFERENTIAL_RESULT_REASON_PROTEIN_PREPARATION_FALLBACK,
    DIFFERENTIAL_RESULT_STATUS_COLUMN,
    DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN,
    DIFFERENTIAL_RESULT_STATUS_TESTED,
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_AUGMENTED_DESIGN_INVALID,
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_CONTRAST_NON_ESTIMABLE,
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_COVARIATE_INVALID,
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_PREPARATION_INELIGIBLE,
    ContrastMatrix,
    DesignMatrix,
    EmpiricalBayesConfig,
)
from phospy.science.differential.models.protein_aware import (
    PROTEIN_AWARE_DIFFERENTIAL_MATCHED_PAIR_COLUMNS,
    ProteinAwareDifferentialComputationRequest,
)
from phospy.science.differential.protein_covariate_adjusted import (
    PROTEIN_AWARE_CENTERING_POLICY,
    PROTEIN_AWARE_COVARIATE_COEFFICIENT_NAME,
)
from phospy.science.statistics.multiple_testing import MultipleTestingCorrection
from phospy.science.transformations.models import IntensityScaleKind
from phospy.workflows._pandas_typing import (
    dataframe_column,
    dataframe_copy,
    dataframe_loc,
    index_snapshot,
)
from phospy.workflows.differential.eligibility import (
    DifferentialPreFitEligibilityResolver,
    differential_status_counts,
    filter_matrix_for_feature_ids,
)
from phospy.workflows.differential.models import (
    DifferentialFeatureEligibilityInputs,
    ProteinAwareDifferentialResolvedInputs,
    ValidatedDifferentialAnalysisRequest,
)

_FloatArray = npt.NDArray[np.float64]
_PROTEIN_VARIANCE_RELATIVE_TOLERANCE = float(np.finfo(np.float64).eps ** 0.75)


@dataclass(frozen=True, slots=True)
class _CandidateSite:
    site_key: str
    protein_identifier: str
    total_protein_row_key: str


@dataclass(frozen=True, slots=True)
class _ProteinCovariateSummary:
    raw_mean: float
    raw_standard_deviation: float
    centered_variance: float
    centered_vector: _FloatArray


@dataclass(frozen=True, slots=True)
class _ScaledSvdDiagnostics:
    rank: int | None
    residual_degrees_of_freedom: float | None
    condition_number: float | None


@dataclass(frozen=True, slots=True)
class _ProteinGroupResolution:
    total_protein_row_key: str
    status: str
    reason: str
    failure_message: str
    raw_mean: float | None
    raw_standard_deviation: float | None
    centered_variance: float | None
    rank: int | None
    residual_degrees_of_freedom: float | None
    condition_number: float | None
    max_condition_number: float

    @property
    def is_testable(self) -> bool:
        return self.status == DIFFERENTIAL_RESULT_STATUS_TESTED


class ProteinAwareDifferentialInputResolver:
    """Resolve validated ordinary differential inputs into protein-aware inputs."""

    def __init__(
        self,
        *,
        pre_fit_eligibility_resolver: DifferentialPreFitEligibilityResolver
        | None = None,
    ) -> None:
        self._pre_fit_eligibility_resolver = (
            pre_fit_eligibility_resolver or DifferentialPreFitEligibilityResolver()
        )

    def run(
        self,
        request: ValidatedDifferentialAnalysisRequest,
    ) -> ProteinAwareDifferentialResolvedInputs:
        if not isinstance(request, ValidatedDifferentialAnalysisRequest):
            _raise_boundary(
                seam="differential.protein_aware_inputs.validated_request_type",
                next_action=(
                    "validate the public differential request before resolving "
                    "protein-aware execution inputs"
                ),
            )

        method_id = _resolve_selected_method(request)
        dataset_view = request.dataset_view or DatasetInternalView(request.dataset)
        sidecar = _require_supported_sidecar(dataset_view)
        total = dataset_view.total
        if total is None:
            _raise_boundary(
                seam="differential.protein_aware_inputs.total_matrix_missing",
                next_action=(
                    "provide an analysis-ready dataset with an owned total-protein "
                    "matrix before selecting protein-aware differential analysis"
                ),
                details={"method": method_id},
            )

        _validate_sidecar_binding(dataset_view)
        _validate_log2_states(request)
        _validate_sidecar_log2_evidence(sidecar)
        _validate_no_prior_total_protein_subtraction(request)
        _validate_no_duplicate_correlation_policy(request)
        _validate_no_actual_technical_replicate_aggregation(request)

        matrix = dataframe_loc(
            dataset_view.phospho,
            slice(None),
            list(request.analysis_sample_ids),
        )
        full_site_ids = tuple(str(value) for value in matrix.index.tolist())
        sample_order = tuple(str(value) for value in request.analysis_sample_ids)
        _validate_sample_alignment(
            request=request,
            total=total,
            protein_covariates=sidecar.protein_covariate_matrix,
            sample_order=sample_order,
        )
        _validate_reserved_protein_covariate_name(request)

        ordinary_eligibility = self._pre_fit_eligibility_resolver.run(
            dataset_view=dataset_view,
            matrix=matrix,
            analysis_sample_ids=request.analysis_sample_ids,
            design=request.design,
            contrasts=request.contrasts,
            policy=request.config.imputed_value_policy,
            max_fraction=request.config.imputed_value_max_fraction,
            minimum_condition_replicates=request.config.minimum_condition_replicates,
        ).feature_eligibility_inputs

        protein_covariates = _selected_protein_covariates(
            sidecar.protein_covariate_matrix,
            sample_order=sample_order,
        )
        site_eligibility = sidecar.site_eligibility
        matched_pairs = sidecar.matched_pairs
        candidate_sites = _candidate_sites(
            full_site_ids=full_site_ids,
            ordinary_eligibility=ordinary_eligibility,
            site_eligibility=site_eligibility,
            matched_pairs=matched_pairs,
        )
        group_resolutions = _resolve_protein_group_eligibility(
            candidate_sites=candidate_sites,
            protein_covariates=protein_covariates,
            base_design=request.design_matrix,
            base_contrasts=request.contrast_matrix,
        )
        feature_eligibility = _build_feature_eligibility_inputs(
            full_site_ids=full_site_ids,
            ordinary_eligibility=ordinary_eligibility,
            site_eligibility=site_eligibility,
            candidate_sites=candidate_sites,
            group_resolutions=group_resolutions,
            method_id=str(method_id),
            matrix_index=matrix.index,
        )
        tested_site_ids = feature_eligibility.testable_feature_ids
        status_counts = _status_count_items(feature_eligibility.result_status)
        reason_counts = _reason_count_items(feature_eligibility.feature_metadata)
        eligibility_counts = _eligibility_count_items(
            full_site_ids=full_site_ids,
            ordinary_eligibility=ordinary_eligibility,
            candidate_sites=candidate_sites,
            tested_site_ids=tested_site_ids,
        )
        if not tested_site_ids:
            _raise_boundary(
                seam="differential.protein_aware_inputs.all_sites_withheld",
                next_action=(
                    "provide at least one phosphosite with ordinary eligibility, "
                    "protein-aware preparation eligibility, usable protein "
                    "covariates, and an admissible augmented design"
                ),
                details={
                    "eligibility_counts": dict(eligibility_counts),
                    "status_counts": dict(status_counts),
                    "reason_counts": dict(reason_counts),
                },
            )

        tested_matched_pairs = _matched_pairs_for_sites(
            candidate_sites=candidate_sites,
            site_ids=tested_site_ids,
        )
        tested_row_keys = _unique_row_keys(tested_matched_pairs)
        resolved_covariates = dataframe_loc(
            protein_covariates,
            list(tested_row_keys),
            list(sample_order),
        )
        computation_request = _build_computation_request(
            matrix=matrix,
            tested_site_ids=tested_site_ids,
            base_design=request.design_matrix,
            base_contrasts=request.contrast_matrix,
            matched_pairs=tested_matched_pairs,
            resolved_protein_covariates=resolved_covariates,
            sample_order=sample_order,
            empirical_bayes=request.config.empirical_bayes,
            multiple_testing_method=request.config.multiple_testing.method,
            method_id=method_id,
        )
        return ProteinAwareDifferentialResolvedInputs(
            computation_request=computation_request,
            feature_eligibility_inputs=feature_eligibility,
            matched_pairs=tested_matched_pairs,
            candidate_matched_pairs=_matched_pairs_for_sites(
                candidate_sites=candidate_sites,
                site_ids=tuple(candidate.site_key for candidate in candidate_sites),
            ),
            resolved_protein_covariates=resolved_covariates,
            site_eligibility_metadata=feature_eligibility.feature_metadata,
            full_site_ids=full_site_ids,
            tested_site_ids=tested_site_ids,
            sample_order=sample_order,
            base_design=request.design_matrix,
            base_contrasts=request.contrast_matrix,
            method_id=method_id,
            preparation_policy=str(sidecar.report.preparation_policy),
            protein_mapping_policy=str(sidecar.report.protein_mapping_policy),
            eligibility_counts=eligibility_counts,
            status_counts=status_counts,
            reason_counts=reason_counts,
        )


def _resolve_selected_method(
    request: ValidatedDifferentialAnalysisRequest,
) -> DifferentialProteinAwareModelMethod:
    model_config = request.config.protein_aware_model
    if model_config is None:
        _raise_boundary(
            seam="differential.protein_aware_inputs.model_not_selected",
            next_action=(
                "call the protein-aware input resolver only for a validated request "
                "with DifferentialAnalysisConfig.protein_aware_model set"
            ),
        )
    method_id = model_config.method
    if method_id not in SUPPORTED_DIFFERENTIAL_PROTEIN_AWARE_MODEL_METHODS:
        supported = ", ".join(
            repr(value) for value in SUPPORTED_DIFFERENTIAL_PROTEIN_AWARE_MODEL_METHODS
        )
        _raise_boundary(
            seam="differential.protein_aware_inputs.method",
            next_action="select a supported protein-aware differential method",
            details={"method": str(method_id), "supported_methods": supported},
        )
    return method_id


def _require_supported_sidecar(
    dataset_view: DatasetInternalView,
) -> ProteinAwarePreparationInternalView:
    sidecar = dataset_view.protein_aware_preparation
    if sidecar is None:
        _raise_boundary(
            seam="differential.protein_aware_inputs.sidecar_missing",
            next_action=(
                "build the analysis-ready dataset with "
                "DatasetProteinAwarePreparationConfig(policy='prepare_model_inputs') "
                "before selecting protein-aware differential analysis"
            ),
        )
    report = sidecar.report
    if report.schema_version != PROTEIN_AWARE_PREPARATION_SCHEMA_VERSION:
        _raise_boundary(
            seam="differential.protein_aware_inputs.sidecar_schema",
            next_action=(
                "rebuild the dataset-owned protein-aware preparation sidecar with "
                "the supported schema version"
            ),
            details={
                "schema_version": str(report.schema_version),
                "supported_schema_version": PROTEIN_AWARE_PREPARATION_SCHEMA_VERSION,
            },
        )
    if (
        str(report.preparation_policy)
        != DATASET_PROTEIN_AWARE_PREPARATION_POLICY_PREPARE_MODEL_INPUTS
    ):
        _raise_boundary(
            seam="differential.protein_aware_inputs.sidecar_policy",
            next_action=(
                "rebuild the dataset-owned protein-aware preparation sidecar with "
                "policy='prepare_model_inputs'"
            ),
            details={"preparation_policy": str(report.preparation_policy)},
        )
    return sidecar


def _validate_sidecar_binding(dataset_view: DatasetInternalView) -> None:
    try:
        dataset_view.validate_protein_aware_preparation_binding()
    except DatasetValidationError as exc:
        _raise_boundary(
            seam="differential.protein_aware_inputs.sidecar_binding",
            next_action=(
                "rebuild protein-aware preparation from the exact analysis-ready "
                "dataset that will be used for differential analysis"
            ),
            details={"binding_error": str(exc)},
        )


def _validate_log2_states(request: ValidatedDifferentialAnalysisRequest) -> None:
    scale_state = request.dataset.intensity_scale_state
    phospho_state = scale_state.phospho
    if (
        not scale_state.is_established
        or phospho_state.kind is not IntensityScaleKind.LOG2
        or not phospho_state.transformed
    ):
        _raise_boundary(
            seam="differential.protein_aware_inputs.phospho_scale",
            next_action=(
                "provide established log2 phosphosite intensities before selecting "
                "protein-aware differential analysis"
            ),
            details={
                "phospho_scale": str(phospho_state.kind.value),
                "established": bool(scale_state.is_established),
            },
        )
    total_state = scale_state.total
    if (
        total_state is None
        or total_state.kind is not IntensityScaleKind.LOG2
        or not total_state.transformed
    ):
        _raise_boundary(
            seam="differential.protein_aware_inputs.total_scale",
            next_action=(
                "provide established log2 total-protein intensities before selecting "
                "protein-aware differential analysis"
            ),
            details={
                "total_scale": None if total_state is None else total_state.kind.value,
                "established": bool(scale_state.is_established),
            },
        )


def _validate_sidecar_log2_evidence(
    sidecar: ProteinAwarePreparationInternalView,
) -> None:
    diagnostics = sidecar.report.transformation_state
    if diagnostics is None:
        _raise_boundary(
            seam="differential.protein_aware_inputs.prepared_covariate_scale",
            next_action=(
                "rebuild protein-aware preparation with transformation diagnostics "
                "proving established log2 total-protein covariates"
            ),
            details={"transformation_state": None},
        )
    if not bool(diagnostics.compatible):
        _raise_boundary(
            seam="differential.protein_aware_inputs.prepared_covariate_scale",
            next_action=(
                "rebuild protein-aware preparation from established log2 phosphosite "
                "and total-protein matrices"
            ),
            details={
                "transformation_state_compatible": False,
                "phospho_transformation_state": _state_payload_for_details(
                    diagnostics.phospho_transformation_state
                ),
                "total_protein_transformation_state": _state_payload_for_details(
                    diagnostics.total_protein_transformation_state
                ),
            },
        )
    _validate_sidecar_matrix_log2_evidence(
        diagnostics.phospho_transformation_state,
        seam="differential.protein_aware_inputs.prepared_phospho_scale",
        details_key="phospho_transformation_state",
        next_action=(
            "rebuild protein-aware preparation with established log2 phosphosite "
            "transformation evidence"
        ),
    )
    _validate_sidecar_matrix_log2_evidence(
        diagnostics.total_protein_transformation_state,
        seam="differential.protein_aware_inputs.prepared_covariate_scale",
        details_key="total_protein_transformation_state",
        next_action=(
            "rebuild protein-aware preparation with established log2 total-protein "
            "covariate transformation evidence"
        ),
    )


def _validate_sidecar_matrix_log2_evidence(
    state: Mapping[str, object] | None,
    *,
    seam: str,
    details_key: str,
    next_action: str,
) -> None:
    established_by = None if state is None else state.get("established_by")
    established_by_text = established_by if isinstance(established_by, str) else ""
    if (
        state is not None
        and str(state.get("kind")) == IntensityScaleKind.LOG2.value
        and state.get("transformed") is True
        and established_by_text.strip()
    ):
        return
    _raise_boundary(
        seam=seam,
        next_action=next_action,
        details={details_key: _state_payload_for_details(state)},
    )


def _state_payload_for_details(
    state: Mapping[str, object] | None,
) -> dict[str, object] | None:
    return None if state is None else dict(state)


def _validate_no_prior_total_protein_subtraction(
    request: ValidatedDifferentialAnalysisRequest,
) -> None:
    correction = request.dataset.processing_state.total_protein_correction
    if bool(correction.applied) or str(correction.policy) != str(
        TotalProteinCorrectionPolicy.NONE
    ):
        _raise_boundary(
            seam="differential.protein_aware_inputs.prior_total_protein_subtraction",
            next_action=(
                "use uncorrected phosphosite log-abundances and total-protein "
                "log-abundances for the protein-covariate-adjusted estimator"
            ),
            details={
                "total_protein_correction_policy": str(correction.policy),
                "total_protein_correction_applied": bool(correction.applied),
            },
        )


def _validate_no_duplicate_correlation_policy(
    request: ValidatedDifferentialAnalysisRequest,
) -> None:
    if str(request.config.paired_design_policy) == str(
        PAIRED_DESIGN_POLICY_DUPLICATE_CORRELATION
    ):
        _raise_boundary(
            seam="differential.protein_aware_inputs.duplicate_correlation",
            next_action=(
                "use paired_design_policy='fixed_block' or paired_design_policy='reject' "
                "with protein-aware differential analysis"
            ),
            details={"paired_design_policy": str(request.config.paired_design_policy)},
        )


def _validate_no_actual_technical_replicate_aggregation(
    request: ValidatedDifferentialAnalysisRequest,
) -> None:
    plan = request.technical_replicate_aggregation_plan
    if plan is not None and plan.requires_aggregation:
        _raise_boundary(
            seam="differential.protein_aware_inputs.technical_replicate_aggregation",
            next_action=(
                "aggregate technical replicates before dataset construction or use a "
                "validated request whose technical-replicate plan is a no-op"
            ),
            details={
                "technical_replicate_group_count": int(len(plan.groups)),
                "aggregate_phospho": bool(plan.aggregate_phospho),
                "aggregate_total_protein": bool(plan.aggregate_total_protein),
            },
        )


def _validate_sample_alignment(
    *,
    request: ValidatedDifferentialAnalysisRequest,
    total: pd.DataFrame,
    protein_covariates: pd.DataFrame,
    sample_order: tuple[str, ...],
) -> None:
    duplicates = _duplicate_values(sample_order)
    if duplicates:
        _raise_boundary(
            seam="differential.protein_aware_inputs.sample_alignment",
            next_action="provide unique analysis sample identifiers",
            details={"duplicate_analysis_sample_ids": duplicates[:5]},
        )
    design_samples = tuple(str(value) for value in request.design_matrix.frame.index)
    if design_samples != sample_order:
        _raise_boundary(
            seam="differential.protein_aware_inputs.sample_alignment",
            next_action=(
                "preserve validator-resolved sample order between the design matrix "
                "and protein-aware input resolution"
            ),
            details={
                "analysis_sample_ids": list(sample_order),
                "design_sample_ids": list(design_samples),
            },
        )
    total_sample_set = {str(value) for value in total.columns}
    protein_covariate_sample_set = {str(value) for value in protein_covariates.columns}
    missing_total_samples = [
        sample_id for sample_id in sample_order if sample_id not in total_sample_set
    ]
    missing_covariate_samples = [
        sample_id
        for sample_id in sample_order
        if sample_id not in protein_covariate_sample_set
    ]
    if missing_total_samples or missing_covariate_samples:
        _raise_boundary(
            seam="differential.protein_aware_inputs.sample_alignment",
            next_action=(
                "rebuild protein-aware preparation so total-protein covariates cover "
                "the validated analysis samples"
            ),
            details={
                "missing_total_samples": missing_total_samples[:5],
                "missing_protein_covariate_samples": missing_covariate_samples[:5],
            },
        )


def _validate_reserved_protein_covariate_name(
    request: ValidatedDifferentialAnalysisRequest,
) -> None:
    if (
        PROTEIN_AWARE_COVARIATE_COEFFICIENT_NAME
        not in request.design_matrix.frame.columns
    ):
        return
    _raise_boundary(
        seam="differential.protein_aware_inputs.reserved_covariate_name",
        next_action=(
            "rename the user-supplied design coefficient before selecting the "
            "protein-aware estimator"
        ),
        details={"reserved_name": PROTEIN_AWARE_COVARIATE_COEFFICIENT_NAME},
    )


def _selected_protein_covariates(
    protein_covariates: pd.DataFrame,
    *,
    sample_order: tuple[str, ...],
) -> pd.DataFrame:
    selected = dataframe_loc(
        protein_covariates,
        slice(None),
        list(sample_order),
    )
    return pd.DataFrame(
        selected.to_numpy(dtype=float),
        index=pd.Index(
            [str(value) for value in selected.index.tolist()],
            name=selected.index.name,
        ),
        columns=pd.Index(sample_order, name=selected.columns.name),
    )


def _candidate_sites(
    *,
    full_site_ids: tuple[str, ...],
    ordinary_eligibility: DifferentialFeatureEligibilityInputs,
    site_eligibility: pd.DataFrame,
    matched_pairs: pd.DataFrame,
) -> tuple[_CandidateSite, ...]:
    ordinary_status = {
        str(site_key): str(status)
        for site_key, status in zip(
            ordinary_eligibility.result_status.index,
            ordinary_eligibility.result_status.to_numpy(dtype=object),
            strict=True,
        )
    }
    site_eligibility_by_site = _frame_by_site(site_eligibility)
    matched_pairs_by_site = _frame_by_site(matched_pairs)
    candidates: list[_CandidateSite] = []
    missing_eligibility_sites = [
        site_key
        for site_key in full_site_ids
        if site_key not in site_eligibility_by_site
    ]
    if missing_eligibility_sites:
        _raise_boundary(
            seam="differential.protein_aware_inputs.site_eligibility_alignment",
            next_action=(
                "rebuild protein-aware preparation so every phosphosite row has "
                "sidecar eligibility"
            ),
            details={"missing_site_keys": missing_eligibility_sites[:5]},
        )
    for site_key in full_site_ids:
        if ordinary_status[site_key] != DIFFERENTIAL_RESULT_STATUS_TESTED:
            continue
        eligibility = str(site_eligibility_by_site[site_key]["eligibility"])
        if (
            eligibility
            != ProteinAwarePreparationEligibility.ELIGIBLE_FOR_PROTEIN_AWARE_PREPARATION.value
        ):
            continue
        pair = matched_pairs_by_site.get(site_key)
        if pair is None:
            _raise_boundary(
                seam="differential.protein_aware_inputs.sidecar_binding",
                next_action=(
                    "rebuild protein-aware preparation so every eligible site has "
                    "one matched total-protein row"
                ),
                details={"eligible_site_without_matched_pair": site_key},
            )
        candidates.append(
            _CandidateSite(
                site_key=site_key,
                protein_identifier=str(pair["protein_identifier"]),
                total_protein_row_key=str(pair["total_protein_row_key"]),
            )
        )
    return tuple(candidates)


def _frame_by_site(frame: pd.DataFrame) -> dict[str, dict[str, object]]:
    if "site_key" not in frame.columns:
        _raise_boundary(
            seam="differential.protein_aware_inputs.sidecar_binding",
            next_action="rebuild protein-aware preparation with a site_key column",
        )
    rows: dict[str, dict[str, object]] = {}
    for record in frame.to_dict(orient="records"):
        site_key = str(record["site_key"])
        if site_key in rows:
            _raise_boundary(
                seam="differential.protein_aware_inputs.sidecar_binding",
                next_action=(
                    "rebuild protein-aware preparation so site_key values are unique"
                ),
                details={"duplicate_site_key": site_key},
            )
        rows[site_key] = cast(dict[str, object], record)
    return rows


def _resolve_protein_group_eligibility(
    *,
    candidate_sites: tuple[_CandidateSite, ...],
    protein_covariates: pd.DataFrame,
    base_design: DesignMatrix,
    base_contrasts: ContrastMatrix,
) -> dict[str, _ProteinGroupResolution]:
    resolutions: dict[str, _ProteinGroupResolution] = {}
    for row_key in dict.fromkeys(
        candidate.total_protein_row_key for candidate in candidate_sites
    ):
        protein_vector = (
            dataframe_loc(
                protein_covariates,
                [row_key],
                slice(None),
            )
            .iloc[0]
            .to_numpy(dtype=float)
        )
        summary = _summarize_protein_covariate(cast(_FloatArray, protein_vector))
        resolutions[row_key] = _evaluate_protein_group(
            row_key=row_key,
            summary=summary,
            base_design=base_design,
            base_contrasts=base_contrasts,
        )
    return resolutions


def _evaluate_protein_group(
    *,
    row_key: str,
    summary: _ProteinCovariateSummary,
    base_design: DesignMatrix,
    base_contrasts: ContrastMatrix,
) -> _ProteinGroupResolution:
    sample_count = int(base_design.frame.shape[0])
    coefficient_count = int(base_design.frame.shape[1]) + 1
    if not np.isfinite(summary.centered_vector).all():
        return _group_failure(
            row_key=row_key,
            status=DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_COVARIATE_INVALID,
            reason=DIFFERENTIAL_RESULT_REASON_PROTEIN_COVARIATE_NON_FINITE,
            message="protein covariate contains non-finite values",
            sample_count=sample_count,
            coefficient_count=coefficient_count,
            summary=summary,
        )
    centered_sum_squares = float(summary.centered_vector @ summary.centered_vector)
    scale = max(float(np.linalg.norm(summary.centered_vector)), 1.0)
    tolerance = (
        _PROTEIN_VARIANCE_RELATIVE_TOLERANCE
        * float(max(sample_count, coefficient_count))
        * scale
    )
    if (
        not math.isfinite(centered_sum_squares)
        or centered_sum_squares <= tolerance
        or summary.centered_variance <= 0.0
    ):
        return _group_failure(
            row_key=row_key,
            status=DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_COVARIATE_INVALID,
            reason=DIFFERENTIAL_RESULT_REASON_PROTEIN_COVARIATE_ZERO_VARIANCE,
            message="protein covariate has zero or numerically unusable variance",
            sample_count=sample_count,
            coefficient_count=coefficient_count,
            summary=summary,
        )
    augmented_design = _augmented_design(
        base_design=base_design.frame,
        centered_protein=summary.centered_vector,
    )
    augmented_contrasts = _augmented_contrasts(base_contrasts.frame)
    try:
        decomposition = decompose_differential_design(augmented_design.to_numpy(float))
    except DifferentialDesignDecompositionError as exc:
        diagnostics = _scaled_svd_diagnostics(augmented_design.to_numpy(float))
        return _group_failure(
            row_key=row_key,
            status=DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_AUGMENTED_DESIGN_INVALID,
            reason=_decomposition_failure_reason(str(exc)),
            message=str(exc),
            sample_count=int(augmented_design.shape[0]),
            coefficient_count=int(augmented_design.shape[1]),
            summary=summary,
            rank=diagnostics.rank,
            residual_degrees_of_freedom=diagnostics.residual_degrees_of_freedom,
            condition_number=diagnostics.condition_number,
        )
    invalid_contrast_positions = decomposition.invalid_contrast_positions(
        augmented_contrasts.to_numpy(dtype=float)
    )
    if invalid_contrast_positions:
        contrast_names = tuple(
            str(augmented_contrasts.columns[position])
            for position in invalid_contrast_positions
        )
        return _group_failure(
            row_key=row_key,
            status=DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_CONTRAST_NON_ESTIMABLE,
            reason=DIFFERENTIAL_RESULT_REASON_PROTEIN_CONTRAST_NON_ESTIMABLE,
            message=(
                "protein-aware augmented contrast is non-estimable or zero-length; "
                f"contrasts={contrast_names}"
            ),
            sample_count=decomposition.sample_count,
            coefficient_count=decomposition.coefficient_count,
            summary=summary,
            rank=decomposition.rank,
            residual_degrees_of_freedom=decomposition.residual_degrees_of_freedom,
            condition_number=decomposition.condition_number,
            max_condition_number=decomposition.max_condition_number,
        )
    return _ProteinGroupResolution(
        total_protein_row_key=row_key,
        status=DIFFERENTIAL_RESULT_STATUS_TESTED,
        reason="",
        failure_message="",
        raw_mean=_finite_or_none(summary.raw_mean),
        raw_standard_deviation=_finite_or_none(summary.raw_standard_deviation),
        centered_variance=_finite_or_none(summary.centered_variance),
        rank=decomposition.rank,
        residual_degrees_of_freedom=decomposition.residual_degrees_of_freedom,
        condition_number=decomposition.condition_number,
        max_condition_number=decomposition.max_condition_number,
    )


def _summarize_protein_covariate(
    protein_vector: _FloatArray,
) -> _ProteinCovariateSummary:
    vector = np.asarray(protein_vector, dtype=np.float64)
    if not np.isfinite(vector).all():
        return _ProteinCovariateSummary(
            raw_mean=math.nan,
            raw_standard_deviation=math.nan,
            centered_variance=math.nan,
            centered_vector=np.full(vector.shape, np.nan, dtype=np.float64),
        )
    raw_mean = float(np.mean(vector))
    raw_standard_deviation = float(np.std(vector, ddof=1)) if vector.size > 1 else 0.0
    centered = cast(_FloatArray, np.asarray(vector - raw_mean, dtype=np.float64))
    centered_variance = float(np.var(centered, ddof=1)) if centered.size > 1 else 0.0
    return _ProteinCovariateSummary(
        raw_mean=raw_mean,
        raw_standard_deviation=raw_standard_deviation,
        centered_variance=centered_variance,
        centered_vector=centered,
    )


def _group_failure(
    *,
    row_key: str,
    status: str,
    reason: str,
    message: str,
    sample_count: int,
    coefficient_count: int,
    summary: _ProteinCovariateSummary,
    rank: int | None = None,
    residual_degrees_of_freedom: float | None = None,
    condition_number: float | None = None,
    max_condition_number: float = DIFFERENTIAL_LINEAR_MODEL_MAX_CONDITION_NUMBER,
) -> _ProteinGroupResolution:
    del sample_count, coefficient_count
    return _ProteinGroupResolution(
        total_protein_row_key=row_key,
        status=status,
        reason=reason,
        failure_message=message,
        raw_mean=_finite_or_none(summary.raw_mean),
        raw_standard_deviation=_finite_or_none(summary.raw_standard_deviation),
        centered_variance=_finite_or_none(summary.centered_variance),
        rank=rank,
        residual_degrees_of_freedom=residual_degrees_of_freedom,
        condition_number=condition_number,
        max_condition_number=max_condition_number,
    )


def _augmented_design(
    *,
    base_design: pd.DataFrame,
    centered_protein: _FloatArray,
) -> pd.DataFrame:
    frame = dataframe_copy(base_design)
    frame.loc[:, PROTEIN_AWARE_COVARIATE_COEFFICIENT_NAME] = centered_protein
    return frame


def _augmented_contrasts(base_contrasts: pd.DataFrame) -> pd.DataFrame:
    protein_row = pd.DataFrame(
        np.zeros((1, int(base_contrasts.shape[1])), dtype=float),
        index=pd.Index(
            [PROTEIN_AWARE_COVARIATE_COEFFICIENT_NAME],
            name=base_contrasts.index.name,
        ),
        columns=base_contrasts.columns.copy(),
    )
    return pd.concat([dataframe_copy(base_contrasts), protein_row], axis=0)


def _scaled_svd_diagnostics(design_values: _FloatArray) -> _ScaledSvdDiagnostics:
    values = np.asarray(design_values, dtype=np.float64)
    sample_count = int(values.shape[0])
    coefficient_count = int(values.shape[1])
    try:
        column_scales = np.asarray(np.linalg.norm(values, axis=0), dtype=float)
        usable = np.isfinite(column_scales) & (column_scales > 0.0)
        if not np.any(usable):
            return _ScaledSvdDiagnostics(
                rank=0,
                residual_degrees_of_freedom=float(sample_count),
                condition_number=math.inf,
            )
        scaled = values[:, usable] / column_scales[usable][np.newaxis, :]
        singular_values = np.asarray(
            np.linalg.svd(scaled, compute_uv=False),
            dtype=np.float64,
        )
    except np.linalg.LinAlgError:
        return _ScaledSvdDiagnostics(
            rank=None,
            residual_degrees_of_freedom=None,
            condition_number=None,
        )
    if singular_values.size == 0 or not np.isfinite(singular_values).all():
        return _ScaledSvdDiagnostics(
            rank=None,
            residual_degrees_of_freedom=None,
            condition_number=None,
        )
    largest = float(singular_values[0])
    tolerance = (
        np.finfo(np.float64).eps * float(max(sample_count, coefficient_count)) * largest
    )
    rank = int(np.count_nonzero(singular_values > tolerance))
    smallest = float(singular_values[-1])
    condition_number = math.inf if smallest == 0.0 else float(largest / smallest)
    return _ScaledSvdDiagnostics(
        rank=rank,
        residual_degrees_of_freedom=float(sample_count - rank),
        condition_number=condition_number,
    )


def _decomposition_failure_reason(message: str) -> str:
    if "residual degrees of freedom must be positive" in message:
        return DIFFERENTIAL_RESULT_REASON_PROTEIN_AUGMENTED_DESIGN_NON_POSITIVE_RESIDUAL_DOF
    if "rank deficient" in message:
        return DIFFERENTIAL_RESULT_REASON_PROTEIN_AUGMENTED_DESIGN_RANK_DEFICIENT
    if "too ill-conditioned" in message:
        return DIFFERENTIAL_RESULT_REASON_PROTEIN_AUGMENTED_DESIGN_ILL_CONDITIONED
    return DIFFERENTIAL_RESULT_REASON_PROTEIN_AUGMENTED_DESIGN_ILL_CONDITIONED


def _build_feature_eligibility_inputs(
    *,
    full_site_ids: tuple[str, ...],
    ordinary_eligibility: DifferentialFeatureEligibilityInputs,
    site_eligibility: pd.DataFrame,
    candidate_sites: tuple[_CandidateSite, ...],
    group_resolutions: dict[str, _ProteinGroupResolution],
    method_id: str,
    matrix_index: pd.Index,
) -> DifferentialFeatureEligibilityInputs:
    ordinary_metadata = ordinary_eligibility.feature_metadata
    ordinary_status = ordinary_eligibility.result_status.astype(str)
    ordinary_reasons = dataframe_column(
        ordinary_metadata,
        DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN,
    ).astype(str)
    site_eligibility_by_site = _frame_by_site(site_eligibility)
    candidates_by_site = {
        candidate.site_key: candidate for candidate in candidate_sites
    }

    statuses: list[str] = []
    reasons: list[str] = []
    metadata_columns: dict[str, list[Any]] = {
        "protein_aware_method_id": [],
        "protein_aware_centering_policy": [],
        "protein_aware_candidate": [],
        "protein_aware_tested": [],
        "protein_aware_preparation_eligibility": [],
        "protein_aware_preparation_reasons": [],
        "protein_aware_mapping_status": [],
        "protein_identifier": [],
        "total_protein_row_key": [],
        "protein_covariate_raw_mean": [],
        "protein_covariate_raw_standard_deviation": [],
        "protein_covariate_centered_variance": [],
        "protein_augmented_design_rank": [],
        "protein_augmented_design_residual_degrees_of_freedom": [],
        "protein_augmented_design_condition_number": [],
        "protein_augmented_design_max_condition_number": [],
        "protein_aware_failure_message": [],
    }

    for site_key in full_site_ids:
        ordinary_site_status = str(ordinary_status.loc[site_key])
        final_status = ordinary_site_status
        final_reason = str(ordinary_reasons.loc[site_key])
        prep_row = site_eligibility_by_site[site_key]
        prep_eligibility = str(prep_row["eligibility"])
        prep_reasons = _normalize_reasons(prep_row.get("reasons", ()))
        mapping_status = str(prep_row["mapping_status"])
        protein_identifier = _optional_string(prep_row.get("protein_identifier"))
        total_row_key = _optional_string(prep_row.get("total_protein_row_key"))
        failure_message = ""
        raw_mean: float | None = None
        raw_std: float | None = None
        centered_variance: float | None = None
        rank: int | None = None
        residual_dof: float | None = None
        condition_number: float | None = None
        max_condition_number = DIFFERENTIAL_LINEAR_MODEL_MAX_CONDITION_NUMBER
        candidate = candidates_by_site.get(site_key)
        protein_tested = False

        if ordinary_site_status == DIFFERENTIAL_RESULT_STATUS_TESTED:
            if (
                prep_eligibility
                != ProteinAwarePreparationEligibility.ELIGIBLE_FOR_PROTEIN_AWARE_PREPARATION.value
            ):
                final_status = (
                    DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_PREPARATION_INELIGIBLE
                )
                final_reason = _preparation_ineligible_reason(prep_eligibility)
                failure_message = (
                    "protein-aware preparation reported "
                    f"{prep_eligibility!r} for this phosphosite"
                )
            elif candidate is not None:
                protein_identifier = candidate.protein_identifier
                total_row_key = candidate.total_protein_row_key
                resolution = group_resolutions[candidate.total_protein_row_key]
                raw_mean = resolution.raw_mean
                raw_std = resolution.raw_standard_deviation
                centered_variance = resolution.centered_variance
                rank = resolution.rank
                residual_dof = resolution.residual_degrees_of_freedom
                condition_number = resolution.condition_number
                max_condition_number = resolution.max_condition_number
                if resolution.is_testable:
                    protein_tested = True
                else:
                    final_status = resolution.status
                    final_reason = resolution.reason
                    failure_message = resolution.failure_message

        statuses.append(final_status)
        reasons.append(final_reason)
        metadata_columns["protein_aware_method_id"].append(method_id)
        metadata_columns["protein_aware_centering_policy"].append(
            PROTEIN_AWARE_CENTERING_POLICY
        )
        metadata_columns["protein_aware_candidate"].append(candidate is not None)
        metadata_columns["protein_aware_tested"].append(protein_tested)
        metadata_columns["protein_aware_preparation_eligibility"].append(
            prep_eligibility
        )
        metadata_columns["protein_aware_preparation_reasons"].append(prep_reasons)
        metadata_columns["protein_aware_mapping_status"].append(mapping_status)
        metadata_columns["protein_identifier"].append(protein_identifier)
        metadata_columns["total_protein_row_key"].append(total_row_key)
        metadata_columns["protein_covariate_raw_mean"].append(raw_mean)
        metadata_columns["protein_covariate_raw_standard_deviation"].append(raw_std)
        metadata_columns["protein_covariate_centered_variance"].append(
            centered_variance
        )
        metadata_columns["protein_augmented_design_rank"].append(rank)
        metadata_columns["protein_augmented_design_residual_degrees_of_freedom"].append(
            residual_dof
        )
        metadata_columns["protein_augmented_design_condition_number"].append(
            condition_number
        )
        metadata_columns["protein_augmented_design_max_condition_number"].append(
            max_condition_number
        )
        metadata_columns["protein_aware_failure_message"].append(failure_message)

    metadata_columns[DIFFERENTIAL_RESULT_STATUS_COLUMN] = list(statuses)
    metadata_columns[DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN] = list(reasons)
    feature_metadata = pd.concat(
        [
            dataframe_copy(ordinary_metadata).drop(
                columns=[
                    DIFFERENTIAL_RESULT_STATUS_COLUMN,
                    DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN,
                ],
                errors="ignore",
            ),
            pd.DataFrame(
                metadata_columns,
                index=index_snapshot(ordinary_metadata.index),
            ),
        ],
        axis=1,
    )
    result_status = pd.Series(
        statuses,
        index=index_snapshot(matrix_index),
        name=DIFFERENTIAL_RESULT_STATUS_COLUMN,
        dtype=str,
    )
    testable_site_ids = tuple(
        site_key
        for site_key, status in zip(full_site_ids, statuses, strict=True)
        if status == DIFFERENTIAL_RESULT_STATUS_TESTED
    )
    return DifferentialFeatureEligibilityInputs(
        feature_metadata=feature_metadata,
        result_status=result_status,
        testable_feature_ids=testable_site_ids,
        attach_to_result_tables=True,
    )


def _preparation_ineligible_reason(eligibility: str) -> str:
    if eligibility == ProteinAwarePreparationEligibility.FALLBACK_TO_PHOSPHO_ONLY.value:
        return DIFFERENTIAL_RESULT_REASON_PROTEIN_PREPARATION_FALLBACK
    return DIFFERENTIAL_RESULT_REASON_PROTEIN_PREPARATION_EXCLUDED


def _matched_pairs_for_sites(
    *,
    candidate_sites: tuple[_CandidateSite, ...],
    site_ids: tuple[str, ...],
) -> pd.DataFrame:
    candidate_by_site = {candidate.site_key: candidate for candidate in candidate_sites}
    rows = []
    for site_id in site_ids:
        candidate = candidate_by_site[site_id]
        rows.append(
            {
                "site_key": candidate.site_key,
                "protein_identifier": candidate.protein_identifier,
                "total_protein_row_key": candidate.total_protein_row_key,
            }
        )
    return pd.DataFrame(
        rows,
        columns=list(PROTEIN_AWARE_DIFFERENTIAL_MATCHED_PAIR_COLUMNS),
    )


def _unique_row_keys(matched_pairs: pd.DataFrame) -> tuple[str, ...]:
    return tuple(
        dict.fromkeys(
            str(value)
            for value in dataframe_column(matched_pairs, "total_protein_row_key")
        )
    )


def _build_computation_request(
    *,
    matrix: pd.DataFrame,
    tested_site_ids: tuple[str, ...],
    base_design: DesignMatrix,
    base_contrasts: ContrastMatrix,
    matched_pairs: pd.DataFrame,
    resolved_protein_covariates: pd.DataFrame,
    sample_order: tuple[str, ...],
    empirical_bayes: EmpiricalBayesConfig,
    multiple_testing_method: MultipleTestingCorrection,
    method_id: DifferentialProteinAwareModelMethod,
) -> ProteinAwareDifferentialComputationRequest:
    try:
        return ProteinAwareDifferentialComputationRequest(
            phosphosite_matrix=filter_matrix_for_feature_ids(
                matrix=matrix,
                feature_ids=tested_site_ids,
            ),
            base_design=base_design,
            base_contrasts=base_contrasts,
            sample_order=sample_order,
            matched_pairs=matched_pairs,
            resolved_protein_covariates=resolved_protein_covariates,
            empirical_bayes=empirical_bayes,
            multiple_testing_method=multiple_testing_method,
            method_id=method_id,
        )
    except PhosPyInputError as exc:
        _raise_boundary(
            seam="differential.protein_aware_inputs.computation_request",
            next_action=(
                "resolve protein-aware feature eligibility and aligned covariates "
                "before constructing the science computation request"
            ),
            details={"input_error": str(exc)},
        )


def _status_count_items(result_status: pd.Series) -> tuple[tuple[str, int], ...]:
    return tuple(differential_status_counts(result_status).items())


def _reason_count_items(feature_metadata: pd.DataFrame) -> tuple[tuple[str, int], ...]:
    reasons = (
        dataframe_column(feature_metadata, DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN)
        .astype(str)
        .tolist()
    )
    counts: Counter[str] = Counter()
    for reason in reasons:
        if reason:
            counts[reason] += 1
    return tuple((reason, int(count)) for reason, count in counts.items())


def _eligibility_count_items(
    *,
    full_site_ids: tuple[str, ...],
    ordinary_eligibility: DifferentialFeatureEligibilityInputs,
    candidate_sites: tuple[_CandidateSite, ...],
    tested_site_ids: tuple[str, ...],
) -> tuple[tuple[str, int], ...]:
    return (
        ("total_site_count", len(full_site_ids)),
        (
            "ordinary_testable_site_count",
            len(ordinary_eligibility.testable_feature_ids),
        ),
        ("protein_preparation_candidate_site_count", len(candidate_sites)),
        ("protein_aware_tested_site_count", len(tested_site_ids)),
    )


def _normalize_reasons(value: object) -> tuple[str, ...]:
    if value is None:
        return ()
    if isinstance(value, str):
        text = value.strip()
        return () if not text else (text,)
    if isinstance(value, tuple | list):
        return tuple(str(item) for item in value)
    return (str(value),)


def _optional_string(value: object) -> str | None:
    if value is None or bool(pd.isna(cast(Any, value))):
        return None
    text = str(value)
    return text if text else None


def _duplicate_values(values: tuple[str, ...]) -> list[str]:
    counts = Counter(values)
    return [value for value in dict.fromkeys(values) if counts[value] > 1]


def _finite_or_none(value: float) -> float | None:
    return value if math.isfinite(value) else None


def _raise_boundary(
    *,
    seam: str,
    next_action: str,
    details: dict[str, object] | None = None,
) -> NoReturn:
    raise WorkflowBoundaryError(
        seam=seam,
        next_action=next_action,
        details=details,
        message_prefix="differential workflow boundary validation failed",
    )


__all__ = ["ProteinAwareDifferentialInputResolver"]
