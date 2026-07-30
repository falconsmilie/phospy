"""Module-selection stability diagnostics for signalome clustering."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from hashlib import sha256

import numpy as np

from phospy.science.signalomes.clustering.candidate_scoring import (
    _ProfileDegeneracySummary,
)
from phospy.science.signalomes.clustering.candidate_selection import (
    _ModuleSelectionComputation,
    filter_cluster_candidates,
    select_best_candidate_count,
)
from phospy.science.signalomes.clustering.policies import (
    SIGNALOME_MODULE_SELECTION_STABILITY_ASSIGNMENT_METRIC,
    SIGNALOME_MODULE_SELECTION_STABILITY_DEFAULT_MAX_SITES,
    SIGNALOME_MODULE_SELECTION_STABILITY_DEFAULT_PERTURBATIONS,
    SIGNALOME_MODULE_SELECTION_STABILITY_PERTURBATION_SCALE,
    SIGNALOME_MODULE_SELECTION_STABILITY_SEED_POLICY_CALLER_FIXED,
    SIGNALOME_MODULE_SELECTION_STABILITY_SEED_POLICY_INPUT_DERIVED,
    SIGNALOME_MODULE_SELECTION_STABILITY_SEED_POLICY_NOT_APPLICABLE,
    SIGNALOME_MODULE_SELECTION_STABILITY_THRESHOLD_DELTA,
)
from phospy.science.signalomes.models import (
    SIGNALOME_MODULE_SELECTION_STABILITY_METHOD,
    SIGNALOME_MODULE_SELECTION_STABILITY_STATUS_NOT_COMPUTABLE,
    SIGNALOME_MODULE_SELECTION_STABILITY_STATUS_STABLE,
    SIGNALOME_MODULE_SELECTION_STABILITY_STATUS_UNSTABLE,
    SIGNALOME_MODULE_SELECTION_STABILITY_VERSION,
    SignalomeClusterCandidateScore,
    SignalomeModuleSelectionAssignmentSimilaritySummary,
    SignalomeModuleSelectionStabilityReport,
    SignalomeModuleSelectionStabilityStatus,
    SignalomeModuleSelectionThresholdSensitivity,
    SignalomeModuleSelectionThresholdSensitivityRecord,
)

_BASE_LIMITATIONS: tuple[str, ...] = (
    "This report is descriptive sensitivity analysis for automatic module-count "
    "selection, not a statistical significance test.",
    "Selected-count frequencies are counts across seeded perturbations, not "
    "probabilities or calibrated confidence values.",
    "Assignment similarity is a descriptive partition-agreement score, not a "
    "confidence probability.",
    "Seeded score perturbations do not replace biological replication, an "
    "independent validation cohort, or a pre-specified inferential design.",
)


@dataclass(frozen=True, slots=True)
class _ResolvedStabilitySeed:
    value: int
    policy: str


def build_not_computable_module_selection_stability_report(
    *,
    reason: str,
    selected_module_count: int,
    input_site_count: int,
    input_dimension_count: int,
    perturbation_count: int = 0,
    limitations: tuple[str, ...] = (),
) -> SignalomeModuleSelectionStabilityReport:
    """Build a typed not-computable stability report without fabricated scores."""

    return SignalomeModuleSelectionStabilityReport(
        evaluation_method=SIGNALOME_MODULE_SELECTION_STABILITY_METHOD,
        evaluation_version=SIGNALOME_MODULE_SELECTION_STABILITY_VERSION,
        seed_policy=SIGNALOME_MODULE_SELECTION_STABILITY_SEED_POLICY_NOT_APPLICABLE,
        random_seed=None,
        perturbation_count=int(perturbation_count),
        selected_count_frequency={},
        assignment_similarity_metric=SIGNALOME_MODULE_SELECTION_STABILITY_ASSIGNMENT_METRIC,
        assignment_similarity=SignalomeModuleSelectionAssignmentSimilaritySummary(
            metric=SIGNALOME_MODULE_SELECTION_STABILITY_ASSIGNMENT_METRIC,
            evaluated_perturbations=0,
            minimum=None,
            median=None,
            mean=None,
            maximum=None,
        ),
        threshold_sensitivity=SignalomeModuleSelectionThresholdSensitivity(
            method="primary_fallback_threshold_grid",
            version=SIGNALOME_MODULE_SELECTION_STABILITY_VERSION,
            records=(),
            selected_count_frequency={},
            disagrees_with_selected_count=False,
        ),
        status=SIGNALOME_MODULE_SELECTION_STABILITY_STATUS_NOT_COMPUTABLE,
        limitations=_BASE_LIMITATIONS + tuple(str(value) for value in limitations),
        not_computable_reason=str(reason),
        base_selected_module_count=int(selected_module_count),
        input_site_count=int(input_site_count),
        input_dimension_count=int(input_dimension_count),
    )


def evaluate_module_selection_stability(
    *,
    scoring_values: np.ndarray,
    base_selection: _ModuleSelectionComputation,
    profile_degeneracy: _ProfileDegeneracySummary,
    primary_threshold: float,
    fallback_threshold: float,
    max_clusters: int,
    perturbation_count: int,
    random_seed: int | None,
    max_stability_sites: int,
    select_core: Callable[[np.ndarray], _ModuleSelectionComputation],
) -> SignalomeModuleSelectionStabilityReport:
    """Evaluate seeded perturbation and threshold sensitivity for auto-selection."""

    values = np.asarray(scoring_values, dtype=float)
    selected_module_count = int(base_selection.diagnostics.selected_module_count)
    input_site_count = int(values.shape[0])
    input_dimension_count = int(values.shape[1]) if values.ndim == 2 else 0
    if values.ndim != 2:
        return build_not_computable_module_selection_stability_report(
            reason="scoring_values must be a two-dimensional matrix",
            selected_module_count=selected_module_count,
            input_site_count=input_site_count,
            input_dimension_count=input_dimension_count,
            perturbation_count=perturbation_count,
        )

    not_computable_reason = _not_computable_reason(
        values=values,
        profile_degeneracy=profile_degeneracy,
        perturbation_count=perturbation_count,
        max_stability_sites=max_stability_sites,
    )
    if not_computable_reason is not None:
        return build_not_computable_module_selection_stability_report(
            reason=not_computable_reason,
            selected_module_count=selected_module_count,
            input_site_count=input_site_count,
            input_dimension_count=input_dimension_count,
            perturbation_count=perturbation_count,
        )
    if not base_selection.diagnostics.candidate_scores:
        return build_not_computable_module_selection_stability_report(
            reason=(
                "candidate module counts were not evaluated, so stability cannot "
                "be distinguished from deterministic fallback behavior"
            ),
            selected_module_count=selected_module_count,
            input_site_count=input_site_count,
            input_dimension_count=input_dimension_count,
            perturbation_count=perturbation_count,
        )

    base_labels = _labels_for_selection(
        base_selection,
        selected_module_count=selected_module_count,
        input_site_count=input_site_count,
    )
    if base_labels is None:
        return build_not_computable_module_selection_stability_report(
            reason=(
                "base automatic selection did not retain labels for the selected "
                "module count"
            ),
            selected_module_count=selected_module_count,
            input_site_count=input_site_count,
            input_dimension_count=input_dimension_count,
            perturbation_count=perturbation_count,
        )

    seed = _resolve_stability_seed(
        scoring_values=values,
        random_seed=random_seed,
    )
    generator = np.random.default_rng(seed.value)
    selected_count_frequency: dict[int, int] = {}
    assignment_similarities: list[float] = []
    perturbation_scales = _perturbation_scales(values)
    for _position in range(int(perturbation_count)):
        perturbed_values = values + generator.normal(
            loc=0.0,
            scale=perturbation_scales,
            size=values.shape,
        )
        perturbed_selection = select_core(perturbed_values)
        perturbed_count = int(perturbed_selection.diagnostics.selected_module_count)
        selected_count_frequency[perturbed_count] = (
            selected_count_frequency.get(perturbed_count, 0) + 1
        )
        perturbed_labels = _labels_for_selection(
            perturbed_selection,
            selected_module_count=perturbed_count,
            input_site_count=input_site_count,
        )
        if perturbed_labels is None:
            continue
        assignment_similarities.append(
            _pairwise_coassignment_agreement(base_labels, perturbed_labels)
        )

    threshold_sensitivity = _threshold_sensitivity(
        candidate_scores=base_selection.diagnostics.candidate_scores,
        base_selected_module_count=selected_module_count,
        primary_threshold=primary_threshold,
        fallback_threshold=fallback_threshold,
    )
    assignment_summary = _assignment_similarity_summary(assignment_similarities)
    status = _stability_status(
        selected_module_count=selected_module_count,
        perturbation_count=int(perturbation_count),
        selected_count_frequency=selected_count_frequency,
        assignment_summary=assignment_summary,
        threshold_sensitivity=threshold_sensitivity,
    )
    return SignalomeModuleSelectionStabilityReport(
        evaluation_method=SIGNALOME_MODULE_SELECTION_STABILITY_METHOD,
        evaluation_version=SIGNALOME_MODULE_SELECTION_STABILITY_VERSION,
        seed_policy=seed.policy,
        random_seed=seed.value,
        perturbation_count=int(perturbation_count),
        selected_count_frequency=selected_count_frequency,
        assignment_similarity_metric=SIGNALOME_MODULE_SELECTION_STABILITY_ASSIGNMENT_METRIC,
        assignment_similarity=assignment_summary,
        threshold_sensitivity=threshold_sensitivity,
        status=status,
        limitations=_limitations_for_status(
            status=status,
            threshold_sensitivity=threshold_sensitivity,
        ),
        not_computable_reason=None,
        base_selected_module_count=selected_module_count,
        input_site_count=input_site_count,
        input_dimension_count=input_dimension_count,
    )


def _not_computable_reason(
    *,
    values: np.ndarray,
    profile_degeneracy: _ProfileDegeneracySummary,
    perturbation_count: int,
    max_stability_sites: int,
) -> str | None:
    site_count = int(values.shape[0])
    dimension_count = int(values.shape[1])
    if int(perturbation_count) < 1:
        return "perturbation_count must be at least one"
    if site_count < 3:
        return "fewer than three phosphosite profiles are available"
    if dimension_count < 2:
        return "fewer than two score dimensions are available"
    non_degenerate_count = site_count - int(profile_degeneracy.excluded_count)
    if non_degenerate_count < 3:
        return (
            "fewer than three non-degenerate phosphosite profiles are available "
            "after correlation-degeneracy filtering"
        )
    if site_count > int(max_stability_sites):
        return (
            "input exceeds the automatic stability perturbation guard; "
            f"site_count={site_count}; max_stability_sites={int(max_stability_sites)}"
        )
    return None


def _resolve_stability_seed(
    *,
    scoring_values: np.ndarray,
    random_seed: int | None,
) -> _ResolvedStabilitySeed:
    if random_seed is not None:
        seed_value = int(random_seed)
        if seed_value < 0 or seed_value >= 2**32:
            raise ValueError("module-selection stability random_seed must fit uint32")
        return _ResolvedStabilitySeed(
            value=seed_value,
            policy=SIGNALOME_MODULE_SELECTION_STABILITY_SEED_POLICY_CALLER_FIXED,
        )

    values = np.ascontiguousarray(np.asarray(scoring_values, dtype="<f8"))
    digest = sha256()
    digest.update(SIGNALOME_MODULE_SELECTION_STABILITY_METHOD.encode("utf-8"))
    digest.update(SIGNALOME_MODULE_SELECTION_STABILITY_VERSION.encode("utf-8"))
    digest.update(str(values.shape).encode("utf-8"))
    digest.update(values.view(np.uint8).tobytes())
    seed_value = int.from_bytes(digest.digest()[:8], byteorder="little") % (2**32)
    return _ResolvedStabilitySeed(
        value=seed_value,
        policy=SIGNALOME_MODULE_SELECTION_STABILITY_SEED_POLICY_INPUT_DERIVED,
    )


def _perturbation_scales(values: np.ndarray) -> np.ndarray:
    column_scales = np.std(values, axis=0)
    finite_positive = column_scales[np.isfinite(column_scales) & (column_scales > 0.0)]
    fallback = 1.0 if finite_positive.size == 0 else float(np.median(finite_positive))
    column_scales = np.where(
        np.isfinite(column_scales) & (column_scales > 0.0),
        column_scales,
        fallback,
    )
    return (SIGNALOME_MODULE_SELECTION_STABILITY_PERTURBATION_SCALE * column_scales)[
        np.newaxis,
        :,
    ]


def _labels_for_selection(
    selection: _ModuleSelectionComputation,
    *,
    selected_module_count: int,
    input_site_count: int,
) -> np.ndarray | None:
    if int(selected_module_count) <= 1:
        return np.zeros(int(input_site_count), dtype=int)
    labels = selection.candidate_labels.get(int(selected_module_count))
    if labels is None:
        return None
    resolved = np.asarray(labels, dtype=int)
    if resolved.shape != (int(input_site_count),):
        return None
    return resolved


def _pairwise_coassignment_agreement(
    left_labels: np.ndarray,
    right_labels: np.ndarray,
) -> float:
    left = np.asarray(left_labels, dtype=int)
    right = np.asarray(right_labels, dtype=int)
    if left.shape != right.shape:
        raise ValueError("label vectors must have identical shape")
    if left.size <= 1:
        return 1.0
    pair_mask = np.triu(np.ones((left.size, left.size), dtype=bool), k=1)
    left_pairs = left[:, None] == left[None, :]
    right_pairs = right[:, None] == right[None, :]
    return float(np.mean(left_pairs[pair_mask] == right_pairs[pair_mask]))


def _assignment_similarity_summary(
    similarities: list[float],
) -> SignalomeModuleSelectionAssignmentSimilaritySummary:
    if not similarities:
        return SignalomeModuleSelectionAssignmentSimilaritySummary(
            metric=SIGNALOME_MODULE_SELECTION_STABILITY_ASSIGNMENT_METRIC,
            evaluated_perturbations=0,
            minimum=None,
            median=None,
            mean=None,
            maximum=None,
        )
    values = np.asarray(similarities, dtype=float)
    return SignalomeModuleSelectionAssignmentSimilaritySummary(
        metric=SIGNALOME_MODULE_SELECTION_STABILITY_ASSIGNMENT_METRIC,
        evaluated_perturbations=int(values.size),
        minimum=float(np.min(values)),
        median=float(np.median(values)),
        mean=float(np.mean(values)),
        maximum=float(np.max(values)),
    )


def _threshold_sensitivity(
    *,
    candidate_scores: dict[int, SignalomeClusterCandidateScore],
    base_selected_module_count: int,
    primary_threshold: float,
    fallback_threshold: float,
) -> SignalomeModuleSelectionThresholdSensitivity:
    primary_thresholds = _threshold_values(primary_threshold)
    fallback_thresholds = _threshold_values(fallback_threshold)
    records: list[SignalomeModuleSelectionThresholdSensitivityRecord] = []
    frequencies: dict[int, int] = {}
    for primary in primary_thresholds:
        for fallback in fallback_thresholds:
            selected_count, threshold_used = _selected_count_for_thresholds(
                candidate_scores=candidate_scores,
                primary_threshold=primary,
                fallback_threshold=fallback,
            )
            records.append(
                SignalomeModuleSelectionThresholdSensitivityRecord(
                    primary_threshold=primary,
                    fallback_threshold=fallback,
                    selected_module_count=selected_count,
                    threshold_used=threshold_used,
                )
            )
            frequencies[selected_count] = frequencies.get(selected_count, 0) + 1
    return SignalomeModuleSelectionThresholdSensitivity(
        method="primary_fallback_threshold_grid",
        version=SIGNALOME_MODULE_SELECTION_STABILITY_VERSION,
        records=tuple(records),
        selected_count_frequency=frequencies,
        disagrees_with_selected_count=any(
            record.selected_module_count != int(base_selected_module_count)
            for record in records
        ),
    )


def _threshold_values(threshold: float) -> tuple[float, ...]:
    values = {
        _clamp_threshold(
            float(threshold) - SIGNALOME_MODULE_SELECTION_STABILITY_THRESHOLD_DELTA
        ),
        _clamp_threshold(float(threshold)),
        _clamp_threshold(
            float(threshold) + SIGNALOME_MODULE_SELECTION_STABILITY_THRESHOLD_DELTA
        ),
    }
    return tuple(sorted(values))


def _clamp_threshold(value: float) -> float:
    return round(min(1.0, max(0.0, float(value))), 10)


def _selected_count_for_thresholds(
    *,
    candidate_scores: dict[int, SignalomeClusterCandidateScore],
    primary_threshold: float,
    fallback_threshold: float,
) -> tuple[int, float | None]:
    primary_candidates = filter_cluster_candidates(
        candidate_scores,
        threshold=float(primary_threshold),
    )
    if primary_candidates:
        return select_best_candidate_count(primary_candidates), float(primary_threshold)
    fallback_candidates = filter_cluster_candidates(
        candidate_scores,
        threshold=float(fallback_threshold),
    )
    if fallback_candidates:
        return select_best_candidate_count(fallback_candidates), float(
            fallback_threshold
        )
    return 1, None


def _stability_status(
    *,
    selected_module_count: int,
    perturbation_count: int,
    selected_count_frequency: dict[int, int],
    assignment_summary: SignalomeModuleSelectionAssignmentSimilaritySummary,
    threshold_sensitivity: SignalomeModuleSelectionThresholdSensitivity,
) -> SignalomeModuleSelectionStabilityStatus:
    if assignment_summary.evaluated_perturbations != int(perturbation_count):
        return SIGNALOME_MODULE_SELECTION_STABILITY_STATUS_NOT_COMPUTABLE
    if selected_count_frequency.get(int(selected_module_count), 0) != int(
        perturbation_count
    ):
        return SIGNALOME_MODULE_SELECTION_STABILITY_STATUS_UNSTABLE
    if threshold_sensitivity.disagrees_with_selected_count:
        return SIGNALOME_MODULE_SELECTION_STABILITY_STATUS_UNSTABLE
    if assignment_summary.minimum is None or assignment_summary.minimum < 0.95:
        return SIGNALOME_MODULE_SELECTION_STABILITY_STATUS_UNSTABLE
    return SIGNALOME_MODULE_SELECTION_STABILITY_STATUS_STABLE


def _limitations_for_status(
    *,
    status: str,
    threshold_sensitivity: SignalomeModuleSelectionThresholdSensitivity,
) -> tuple[str, ...]:
    limitations: list[str] = list(_BASE_LIMITATIONS)
    if status == SIGNALOME_MODULE_SELECTION_STABILITY_STATUS_UNSTABLE:
        limitations.append(
            "At least one seeded perturbation or threshold-grid point changed the "
            "selected count or assignment partition; downstream module summaries "
            "should be interpreted as configuration-sensitive descriptive output."
        )
    if threshold_sensitivity.disagrees_with_selected_count:
        limitations.append(
            "Threshold sensitivity changed the selected module count for at least "
            "one primary/fallback threshold-grid point."
        )
    return tuple(limitations)


__all__ = [
    "SIGNALOME_MODULE_SELECTION_STABILITY_ASSIGNMENT_METRIC",
    "SIGNALOME_MODULE_SELECTION_STABILITY_DEFAULT_MAX_SITES",
    "SIGNALOME_MODULE_SELECTION_STABILITY_DEFAULT_PERTURBATIONS",
    "SIGNALOME_MODULE_SELECTION_STABILITY_PERTURBATION_SCALE",
    "SIGNALOME_MODULE_SELECTION_STABILITY_SEED_POLICY_CALLER_FIXED",
    "SIGNALOME_MODULE_SELECTION_STABILITY_SEED_POLICY_INPUT_DERIVED",
    "SIGNALOME_MODULE_SELECTION_STABILITY_SEED_POLICY_NOT_APPLICABLE",
    "SIGNALOME_MODULE_SELECTION_STABILITY_THRESHOLD_DELTA",
    "build_not_computable_module_selection_stability_report",
    "evaluate_module_selection_stability",
]
