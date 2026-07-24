"""Resolved-input eligibility validation for kinase workflow execution."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import NoReturn, TypedDict, cast

import pandas as pd

from phospy.contracts.configs import (
    KINASE_ATTRITION_POLICY_ON_VIOLATION_ERROR,
    KINASE_RELIABILITY_PROFILE_PRODUCTION,
)
from phospy.contracts.configs.kinase import (
    validate_kinase_production_reliability_invariants,
)
from phospy.errors.workflows import PhosPyWorkflowError, WorkflowBoundaryError
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.science.prediction.motif_scoring import (
    KINASE_LIBRARY_RESIDUE_CLASS_SER_THR,
    KINASE_LIBRARY_RESIDUE_CLASS_TYR,
)
from phospy.science.references.kinase_library import KinaseLibraryResource
from phospy.science.references.models import ReferenceBundle
from phospy.workflows.kinase.attrition_metrics import (
    KinaseAttritionMetrics,
    KinaseAttritionPolicyViolation,
    build_kinase_attrition_metrics_from_overlap,
    build_kinase_attrition_provenance_payload,
    evaluate_kinase_attrition_policy,
)
from phospy.workflows.kinase.contracts import ResolvedKinaseExecutionConfig
from phospy.workflows.kinase.scoring_mode_contracts import (
    kinase_scoring_mode_input_contract,
)


class _OverlapSummary(TypedDict):
    dataset_sites: int
    reference_sites: int
    overlap_sites: int
    overlap_site_ids: set[str]
    reference_kinases: int
    kinases_with_overlap: int
    max_quantified_sites_per_kinase: int
    per_kinase_quantified: pd.Series


@dataclass(slots=True)
class ResolvedKinaseInputs:
    """Interpreter-resolved inputs before final execution contract construction."""

    dataset: AnalysisReadyPhosphoDataset
    dataset_phospho: pd.DataFrame
    references: ReferenceBundle
    kinase_substrate_map: pd.DataFrame
    site_sequences: pd.DataFrame
    site_identity_map: pd.DataFrame
    scoring_site_index: pd.Index
    activity_phospho_matrix: pd.DataFrame
    execution_config: ResolvedKinaseExecutionConfig
    reference_site_count: int
    kinase_library_resource: KinaseLibraryResource | None = None
    site_sequence_merge_diagnostics: dict[str, object] = field(default_factory=dict)
    attrition_metrics: KinaseAttritionMetrics | None = None
    attrition_policy_violations: tuple[KinaseAttritionPolicyViolation, ...] = ()


class ResolvedKinaseEligibilityValidator:
    """Validate eligibility checks that require resolved references and config."""

    _KINASE_COLUMN = "kinase"
    _SUBSTRATE_COLUMN = "substrate_site"

    def run(self, resolved_inputs: ResolvedKinaseInputs) -> None:
        overlap_counts = self._summarize_overlap(
            dataset=resolved_inputs.dataset_phospho,
            kinase_substrate_map=resolved_inputs.kinase_substrate_map,
        )
        mode_contract = kinase_scoring_mode_input_contract(
            resolved_inputs.execution_config.scoring_mode
        )
        self._validate_reliability_profile(resolved_inputs)
        requires_profile_context = mode_contract.requires_substrate_reference_overlap
        if requires_profile_context:
            self._validate_reference_coverage(
                overlap_counts=overlap_counts,
                resolved_inputs=resolved_inputs,
            )
            self._validate_eligible_kinases(
                overlap_counts=overlap_counts,
                resolved_inputs=resolved_inputs,
            )
        attrition_metrics = build_kinase_attrition_metrics_from_overlap(
            total_dataset_sites=int(overlap_counts["dataset_sites"]),
            reference_overlap_site_ids=overlap_counts["overlap_site_ids"],
            sequence_supported_site_index=resolved_inputs.scoring_site_index,
        )
        resolved_inputs.attrition_metrics = attrition_metrics
        resolved_inputs.attrition_policy_violations = self._enforce_attrition_policy(
            attrition_metrics=attrition_metrics,
            resolved_inputs=resolved_inputs,
        )
        self._validate_scoring_site_support(
            scoring_site_index=resolved_inputs.scoring_site_index,
            dataset=resolved_inputs.dataset_phospho,
            site_sequences=resolved_inputs.site_sequences,
        )
        if requires_profile_context:
            self._validate_scored_site_support(attrition_metrics=attrition_metrics)
        self._validate_kinase_library_resource_usability(
            resource=resolved_inputs.kinase_library_resource,
            site_sequences=resolved_inputs.site_sequences.reindex(
                resolved_inputs.scoring_site_index
            ),
            resolved_inputs=resolved_inputs,
        )

    @classmethod
    def _summarize_overlap(
        cls,
        *,
        dataset: pd.DataFrame,
        kinase_substrate_map: pd.DataFrame,
    ) -> _OverlapSummary:
        dataset_sites = set(dataset.index.astype(str).tolist())
        reference_sites = set(
            kinase_substrate_map.loc[:, cls._SUBSTRATE_COLUMN].astype(str).tolist()
        )
        overlapping_sites = dataset_sites.intersection(reference_sites)
        overlapping_map = kinase_substrate_map[
            kinase_substrate_map.loc[:, cls._SUBSTRATE_COLUMN]
            .astype(str)
            .isin(overlapping_sites)
        ]
        per_kinase_quantified = cast(
            pd.Series,
            overlapping_map.groupby(cls._KINASE_COLUMN, sort=False)[
                cls._SUBSTRATE_COLUMN
            ]
            .nunique()
            .astype("int64"),
        )
        max_quantified_sites_per_kinase = (
            0
            if per_kinase_quantified.empty
            else int(per_kinase_quantified.to_numpy(dtype="int64", copy=False).max())
        )
        return {
            "dataset_sites": len(dataset_sites),
            "reference_sites": len(reference_sites),
            "overlap_sites": len(overlapping_sites),
            "overlap_site_ids": overlapping_sites,
            "reference_kinases": int(
                kinase_substrate_map.loc[:, cls._KINASE_COLUMN].nunique()
            ),
            "kinases_with_overlap": int(per_kinase_quantified.size),
            "max_quantified_sites_per_kinase": max_quantified_sites_per_kinase,
            "per_kinase_quantified": per_kinase_quantified,
        }

    def _validate_reference_coverage(
        self,
        *,
        overlap_counts: _OverlapSummary,
        resolved_inputs: ResolvedKinaseInputs,
    ) -> None:
        overlap_sites = int(overlap_counts["overlap_sites"])
        if overlap_sites > 0:
            return
        self._raise_boundary_error(
            seam="kinase.interpreter.reference_coverage",
            next_action=(
                "use references whose display IDs match dataset display_id metadata "
                "or verify dataset site_key/display_id identity mapping"
            ),
            dataset_sites=overlap_counts["dataset_sites"],
            reference_sites=int(resolved_inputs.reference_site_count),
            overlap_sites=overlap_sites,
            scoring_config_min_substrates=(
                resolved_inputs.execution_config.scoring_min_substrates
            ),
        )

    def _validate_eligible_kinases(
        self,
        *,
        overlap_counts: _OverlapSummary,
        resolved_inputs: ResolvedKinaseInputs,
    ) -> None:
        per_kinase_quantified = overlap_counts["per_kinase_quantified"]
        if not isinstance(per_kinase_quantified, pd.Series):
            raise PhosPyWorkflowError(
                "kinase resolved eligibility validator expected "
                "overlap_counts['per_kinase_quantified'] to be a pandas Series, "
                f"got {type(per_kinase_quantified).__name__}"
            )
        config = resolved_inputs.execution_config
        eligible_mask = per_kinase_quantified >= int(config.scoring_min_substrates)
        eligible_kinases = per_kinase_quantified[eligible_mask]
        if not eligible_kinases.empty:
            return
        self._raise_boundary_error(
            seam="kinase.interpreter.eligible_kinases",
            next_action=(
                "lower scoring_config.min_substrates or provide references with "
                "deeper overlap for the current dataset "
                "(scientific floor: min_substrates >= 2)"
            ),
            reference_kinases=overlap_counts["reference_kinases"],
            kinases_with_overlap=overlap_counts["kinases_with_overlap"],
            eligible_kinases=int(eligible_kinases.size),
            max_quantified_sites_per_kinase=overlap_counts[
                "max_quantified_sites_per_kinase"
            ],
            scoring_config_min_substrates=config.scoring_min_substrates,
            prediction_config_deterministic_max_selected_kinases=(
                config.prediction_deterministic_max_selected_kinases
            ),
            prediction_config_adaptive_ensemble_runs=(
                config.prediction_adaptive_ensemble_runs
            ),
            prediction_config_mode=config.prediction_mode,
        )

    def _validate_scoring_site_support(
        self,
        *,
        scoring_site_index: pd.Index,
        dataset: pd.DataFrame,
        site_sequences: pd.DataFrame,
    ) -> None:
        if not scoring_site_index.empty:
            return
        self._raise_boundary_error(
            seam="kinase.interpreter.sequence_support",
            next_action=(
                "ensure references.site_sequences display IDs can be projected to "
                "dataset site_key rows"
            ),
            dataset_sites=int(dataset.index.size),
            reference_sequence_sites=int(site_sequences.index.size),
            sequence_supported_sites=0,
        )

    def _validate_scored_site_support(
        self,
        *,
        attrition_metrics: KinaseAttritionMetrics,
    ) -> None:
        if int(attrition_metrics.scored_sites) > 0:
            return
        self._raise_boundary_error(
            seam="kinase.interpreter.scored_site_support",
            next_action=(
                "ensure at least one dataset site has both kinase-substrate "
                "reference overlap and usable site-sequence support"
            ),
            **attrition_metrics.to_payload(),
        )

    def _validate_kinase_library_resource_usability(
        self,
        *,
        resource: KinaseLibraryResource | None,
        site_sequences: pd.DataFrame,
        resolved_inputs: ResolvedKinaseInputs,
    ) -> None:
        if resource is None:
            return
        resource_residue_classes = {
            matrix.residue_class.value for matrix in resource.matrices
        }
        site_residue_classes = {
            residue_class
            for residue_class in (
                _central_residue_class(
                    value,
                    upstream=int(resource.sequence_window.upstream_residues),
                    downstream=int(resource.sequence_window.downstream_residues),
                )
                for value in site_sequences.loc[:, "site_sequence"].tolist()
            )
            if residue_class is not None
        }
        if site_residue_classes.intersection(resource_residue_classes):
            return
        self._raise_boundary_error(
            seam="kinase.interpreter.kinase_library_resource_usability",
            next_action=(
                "provide a Kinase Library resource with residue_class lanes that "
                "match at least one resolved scoring-site sequence"
            ),
            scoring_mode=resolved_inputs.execution_config.scoring_mode,
            resource_residue_classes=tuple(sorted(resource_residue_classes)),
            scoring_site_residue_classes=tuple(sorted(site_residue_classes)),
            scoring_site_count=int(site_sequences.shape[0]),
        )

    def _validate_reliability_profile(
        self,
        resolved_inputs: ResolvedKinaseInputs,
    ) -> None:
        config = resolved_inputs.execution_config
        if (
            config.effective_reliability_profile
            is not KINASE_RELIABILITY_PROFILE_PRODUCTION
        ):
            return
        try:
            validate_kinase_production_reliability_invariants(
                min_substrates=int(config.scoring_min_substrates),
                profile_self_inclusion_policy=config.profile_self_inclusion_policy,
                localisation_requirement=config.localisation_requirement,
                attrition_policy=config.attrition_policy,
                field_name="kinase.execution_config",
                error_type=WorkflowBoundaryError,
            )
        except WorkflowBoundaryError as exc:
            self._raise_boundary_error(
                seam="kinase.interpreter.production_reliability_profile",
                next_action=(
                    "use KinaseScoringConfig.production(...) with explicit "
                    "coverage thresholds, leave-one-out profile scoring, and "
                    "strict site-level localisation"
                ),
                requested_reliability_profile=(
                    None
                    if config.requested_reliability_profile is None
                    else str(config.requested_reliability_profile)
                ),
                effective_reliability_profile=str(config.effective_reliability_profile),
                scoring_config_min_substrates=int(config.scoring_min_substrates),
                profile_self_inclusion_policy=str(config.profile_self_inclusion_policy),
                localisation_requirement_policy=str(
                    config.localisation_requirement.policy
                ),
                localisation_requirement_require_present=bool(
                    config.localisation_requirement.require_present
                ),
                localisation_requirement_minimum_probability=(
                    None
                    if config.localisation_requirement.minimum_probability is None
                    else float(config.localisation_requirement.minimum_probability)
                ),
                attrition_policy_on_violation=str(config.attrition_policy.on_violation),
                minimum_reference_overlap_fraction=float(
                    config.attrition_policy.minimum_reference_overlap_fraction
                ),
                minimum_sequence_supported_fraction=float(
                    config.attrition_policy.minimum_sequence_supported_fraction
                ),
                minimum_scored_fraction=float(
                    config.attrition_policy.minimum_scored_fraction
                ),
                validation_error=str(exc),
            )

    def _enforce_attrition_policy(
        self,
        *,
        attrition_metrics: KinaseAttritionMetrics,
        resolved_inputs: ResolvedKinaseInputs,
    ) -> tuple[KinaseAttritionPolicyViolation, ...]:
        policy = resolved_inputs.execution_config.attrition_policy
        violations = evaluate_kinase_attrition_policy(
            metrics=attrition_metrics,
            policy=policy,
        )
        if not violations:
            return ()
        if policy.on_violation != KINASE_ATTRITION_POLICY_ON_VIOLATION_ERROR:
            return violations
        first_violation = violations[0]
        attrition_provenance = build_kinase_attrition_provenance_payload(
            metrics=attrition_metrics,
            policy=policy,
            violations=violations,
        )
        raise WorkflowBoundaryError(
            first_violation.message,
            seam="kinase.interpreter.attrition_policy",
            next_action=(
                "relax scoring_config.attrition_policy thresholds or provide "
                "references and site-sequence support with stronger dataset coverage"
            ),
            details={
                "violations": [violation.to_payload() for violation in violations],
                "attrition_provenance": attrition_provenance,
                "policy_outcome": attrition_provenance["policy_outcome"],
                **attrition_metrics.to_payload(),
            },
        )

    @staticmethod
    def _raise_boundary_error(
        *,
        seam: str,
        next_action: str,
        **details: object,
    ) -> NoReturn:
        raise WorkflowBoundaryError(
            seam=seam,
            next_action=next_action,
            details=details,
            message_prefix="kinase workflow boundary validation failed",
        )


def _central_residue_class(
    value: object,
    *,
    upstream: int,
    downstream: int,
) -> str | None:
    if not isinstance(value, str):
        return None
    sequence = value.strip().upper()
    if not sequence:
        return None
    window_size = int(upstream) + 1 + int(downstream)
    if len(sequence) == window_size:
        centre_index = int(upstream)
    elif len(sequence) >= window_size and len(sequence) % 2 == 1:
        centre_index = len(sequence) // 2
    else:
        return None
    if centre_index < 0 or centre_index >= len(sequence):
        return None
    residue = sequence[centre_index]
    if residue in {"S", "T"}:
        return KINASE_LIBRARY_RESIDUE_CLASS_SER_THR
    if residue == "Y":
        return KINASE_LIBRARY_RESIDUE_CLASS_TYR
    return None


__all__ = [
    "ResolvedKinaseEligibilityValidator",
    "ResolvedKinaseInputs",
]
