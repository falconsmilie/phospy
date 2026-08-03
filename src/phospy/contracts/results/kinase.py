"""Public kinase workflow result contracts."""
# pyright: reportMissingTypeStubs=false

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field

import pandas as pd

from phospy.contracts.result_caveats import ResultCaveat, validate_result_caveats
from phospy.errors.input import PhosPyInputError
from phospy.frames.ownership import export_optional_dataframe
from phospy.provenance.immutability import freeze_json_mapping, thaw_json_mapping
from phospy.provenance.models import RunProvenance
from phospy.science.activities.models import (
    ActivityMethodDiagnostics,
    KinaseActivityResult,
    KseaZScoreActivityDiagnostics,
    SsgseaSubstrateEnrichmentActivityDiagnostics,
    WeightedSubstrateActivityDiagnostics,
)
from phospy.science.datasets.models import (
    AnalysisReadyPhosphoDataset,
    DatasetPreprocessingReport,
)
from phospy.science.prediction.models import KinasePredictionResult, KinaseScoringResult
from phospy.science.references.models import ReferenceBundle
from phospy.science.tables.kinase import KinaseSubstrateContributionTable

KinaseWorkflowCaveat = ResultCaveat


@dataclass(frozen=True, slots=True)
class KinaseWorkflowPreprocessingAttritionSummary:
    """Preprocessing-owned site attrition counters composed into kinase results."""

    input_rows: int
    rows_removed_during_preprocessing: int
    rows_removed_invalid_or_missing_site_identifiers: int
    duplicate_sites_merged_or_resolved: int
    output_rows: int
    sequence_complete_sites: int | None = None


@dataclass(frozen=True, slots=True)
class KinaseWorkflowScoringAttritionSummary:
    """Kinase workflow site-eligibility counters after preprocessing."""

    rows_removed_invalid_or_missing_site_identifiers: int
    final_quantitative_sites_entering_scoring: int
    sites_with_valid_site_sequence: int
    sites_without_usable_site_sequence: int
    sites_eligible_for_motif_scoring: int
    sites_with_kinase_substrate_reference_profile_evidence: int
    sites_contributing_to_final_fused_prediction_scoring_output: int
    sites_contributing_to_activity_scoring: int | None


@dataclass(frozen=True, slots=True)
class KinaseWorkflowSiteAttritionSummary:
    """Compact, user-facing kinase site attrition summary."""

    preprocessing: KinaseWorkflowPreprocessingAttritionSummary
    scoring: KinaseWorkflowScoringAttritionSummary


_KINASE_ATTRITION_POLICY_OUTCOMES = frozenset({"passed", "warned", "failed"})


@dataclass(frozen=True, slots=True)
class KinaseWorkflowAttritionProvenance:
    """Structured kinase attrition metrics, policy, and policy outcome."""

    metrics: Mapping[str, object]
    policy: Mapping[str, object]
    policy_outcome: str
    policy_violations: tuple[Mapping[str, object], ...] = ()
    warning_messages: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        metrics = _require_mapping(
            self.metrics,
            field_name="kinase_workflow_result.attrition_provenance.metrics",
        )
        policy = _require_mapping(
            self.policy,
            field_name="kinase_workflow_result.attrition_provenance.policy",
        )
        policy_outcome = _require_policy_outcome(self.policy_outcome)
        policy_violations = tuple(
            _require_mapping(
                violation,
                field_name=(
                    "kinase_workflow_result.attrition_provenance."
                    f"policy_violations[{index}]"
                ),
            )
            for index, violation in enumerate(self.policy_violations)
        )
        warning_messages = tuple(
            _require_non_empty_text(
                message,
                field_name=(
                    "kinase_workflow_result.attrition_provenance.warning_messages[]"
                ),
            )
            for message in self.warning_messages
        )
        object.__setattr__(self, "metrics", metrics)
        object.__setattr__(self, "policy", policy)
        object.__setattr__(self, "policy_outcome", policy_outcome)
        object.__setattr__(self, "policy_violations", policy_violations)
        object.__setattr__(self, "warning_messages", warning_messages)

    def to_payload(self) -> dict[str, object]:
        return {
            "metrics": thaw_json_mapping(
                self.metrics,
                field_name="kinase_workflow_result.attrition_provenance.metrics",
            ),
            "policy": thaw_json_mapping(
                self.policy,
                field_name="kinase_workflow_result.attrition_provenance.policy",
            ),
            "policy_outcome": self.policy_outcome,
            "policy_violations": [
                thaw_json_mapping(
                    violation,
                    field_name=(
                        "kinase_workflow_result.attrition_provenance."
                        f"policy_violations[{index}]"
                    ),
                )
                for index, violation in enumerate(self.policy_violations)
            ],
            "warning_messages": list(self.warning_messages),
        }


@dataclass(frozen=True, slots=True)
class KinaseEligibilityReport:
    """Compact, user-facing kinase workflow eligibility counters.

    `eligible_kinases` counts kinases whose projected, sequence-supported
    quantified substrates meet `KinaseScoringConfig.min_substrates`.
    `excluded_kinases_below_min_substrates` counts overlapping kinases with too
    few usable substrates. The report is count-based; it does not add per-kinase
    weak-support flags to scoring result tables.
    """

    total_dataset_sites: int
    sequence_complete_sites: int
    localisation_eligible_sites: int | None
    reference_overlap_sites: int
    excluded_no_reference_match: int
    excluded_low_localisation: int | None
    eligible_kinases: int
    excluded_kinases_below_min_substrates: int


@dataclass(frozen=True, slots=True, init=False)
class KinaseWorkflowResult:
    """Top-level public kinase workflow result.

    This is a workflow-owned container, not a direct user-construction
    validator. Workflow execution is the supported construction path for
    scientific coherence across `dataset`, `references`, scoring, prediction,
    optional activity, eligibility, attrition, and provenance. The nested stage
    result objects own their public table schemas; this container keeps the
    workflow-assembled objects together without re-running workflow validation.
    """

    dataset: AnalysisReadyPhosphoDataset
    references: ReferenceBundle
    scoring_result: KinaseScoringResult
    prediction_result: KinasePredictionResult
    eligibility_report: KinaseEligibilityReport | None = None
    site_attrition_summary: KinaseWorkflowSiteAttritionSummary | None = None
    attrition_provenance: KinaseWorkflowAttritionProvenance | None = None
    activity_result: KinaseActivityResult | None = None
    provenance: RunProvenance | None = None
    caveats: tuple[ResultCaveat, ...] = ()
    _substrate_contributions: pd.DataFrame | None = field(
        init=False,
        repr=False,
        compare=False,
    )

    def __init__(
        self,
        dataset: AnalysisReadyPhosphoDataset,
        references: ReferenceBundle,
        scoring_result: KinaseScoringResult,
        prediction_result: KinasePredictionResult,
        eligibility_report: KinaseEligibilityReport | None = None,
        site_attrition_summary: KinaseWorkflowSiteAttritionSummary | None = None,
        attrition_provenance: KinaseWorkflowAttritionProvenance | None = None,
        activity_result: KinaseActivityResult | None = None,
        provenance: RunProvenance | None = None,
        substrate_contributions: pd.DataFrame | None = None,
        caveats: tuple[ResultCaveat, ...] = (),
    ) -> None:
        self._init_kinase_workflow_result(
            dataset=dataset,
            references=references,
            scoring_result=scoring_result,
            prediction_result=prediction_result,
            eligibility_report=eligibility_report,
            site_attrition_summary=site_attrition_summary,
            attrition_provenance=attrition_provenance,
            activity_result=activity_result,
            provenance=provenance,
            substrate_contributions=substrate_contributions,
            caveats=caveats,
            assume_owned=False,
        )

    def _init_kinase_workflow_result(
        self,
        *,
        dataset: AnalysisReadyPhosphoDataset,
        references: ReferenceBundle,
        scoring_result: KinaseScoringResult,
        prediction_result: KinasePredictionResult,
        eligibility_report: KinaseEligibilityReport | None = None,
        site_attrition_summary: KinaseWorkflowSiteAttritionSummary | None = None,
        attrition_provenance: KinaseWorkflowAttritionProvenance | None = None,
        activity_result: KinaseActivityResult | None = None,
        provenance: RunProvenance | None = None,
        substrate_contributions: pd.DataFrame | None = None,
        caveats: tuple[ResultCaveat, ...] = (),
        assume_owned: bool,
    ) -> None:
        _require_instance(
            dataset,
            AnalysisReadyPhosphoDataset,
            field_name="kinase_workflow_result.dataset",
        )
        _require_instance(
            references,
            ReferenceBundle,
            field_name="kinase_workflow_result.references",
        )
        _require_instance(
            scoring_result,
            KinaseScoringResult,
            field_name="kinase_workflow_result.scoring_result",
        )
        _require_instance(
            prediction_result,
            KinasePredictionResult,
            field_name="kinase_workflow_result.prediction_result",
        )
        _require_optional_instance(
            eligibility_report,
            KinaseEligibilityReport,
            field_name="kinase_workflow_result.eligibility_report",
        )
        _require_optional_instance(
            site_attrition_summary,
            KinaseWorkflowSiteAttritionSummary,
            field_name="kinase_workflow_result.site_attrition_summary",
        )
        _require_optional_instance(
            activity_result,
            KinaseActivityResult,
            field_name="kinase_workflow_result.activity_result",
        )
        _require_optional_instance(
            provenance,
            RunProvenance,
            field_name="kinase_workflow_result.provenance",
        )
        object.__setattr__(self, "dataset", dataset)
        object.__setattr__(self, "references", references)
        object.__setattr__(self, "scoring_result", scoring_result)
        object.__setattr__(self, "prediction_result", prediction_result)
        object.__setattr__(self, "eligibility_report", eligibility_report)
        object.__setattr__(self, "site_attrition_summary", site_attrition_summary)
        object.__setattr__(
            self,
            "attrition_provenance",
            _own_attrition_provenance(attrition_provenance),
        )
        object.__setattr__(self, "activity_result", activity_result)
        object.__setattr__(self, "provenance", provenance)
        object.__setattr__(self, "caveats", _own_workflow_caveats(caveats))
        object.__setattr__(
            self,
            "_substrate_contributions",
            _own_optional_kinase_substrate_contributions(
                substrate_contributions,
                assume_owned=assume_owned,
            ),
        )

    @property
    def input_dataset_preprocessing_report(self) -> DatasetPreprocessingReport | None:
        """Return preprocessing provenance of the input analysis-ready dataset."""

        return self.dataset.preprocessing_report

    @property
    def substrate_contributions(self) -> pd.DataFrame | None:
        """Return optional substrate-level contribution rows."""

        return export_optional_dataframe(self._substrate_contributions)

    def substrate_contributions_dataframe(self) -> pd.DataFrame | None:
        """Return optional substrate-level contribution rows."""

        return export_optional_dataframe(self._substrate_contributions)

    @classmethod
    def _from_owned(
        cls,
        *,
        dataset: AnalysisReadyPhosphoDataset,
        references: ReferenceBundle,
        scoring_result: KinaseScoringResult,
        prediction_result: KinasePredictionResult,
        eligibility_report: KinaseEligibilityReport | None = None,
        site_attrition_summary: KinaseWorkflowSiteAttritionSummary | None = None,
        attrition_provenance: KinaseWorkflowAttritionProvenance | None = None,
        activity_result: KinaseActivityResult | None = None,
        provenance: RunProvenance | None = None,
        substrate_contributions: pd.DataFrame | None = None,
        caveats: tuple[ResultCaveat, ...] = (),
    ) -> KinaseWorkflowResult:
        result = object.__new__(cls)
        KinaseWorkflowResult._init_kinase_workflow_result(
            result,
            dataset=dataset,
            references=references,
            scoring_result=scoring_result,
            prediction_result=prediction_result,
            eligibility_report=eligibility_report,
            site_attrition_summary=site_attrition_summary,
            attrition_provenance=attrition_provenance,
            activity_result=activity_result,
            provenance=provenance,
            substrate_contributions=substrate_contributions,
            caveats=caveats,
            assume_owned=True,
        )
        return result


def _own_optional_kinase_substrate_contributions(
    table: pd.DataFrame | None,
    *,
    assume_owned: bool,
) -> pd.DataFrame | None:
    if table is None:
        return None
    return KinaseSubstrateContributionTable(
        frame=table,
        _assume_owned=assume_owned,
    ).frame


def _own_workflow_caveats(
    caveats: tuple[ResultCaveat, ...],
) -> tuple[ResultCaveat, ...]:
    return validate_result_caveats(
        caveats,
        field_name="kinase_workflow_result.caveats",
    )


def _own_attrition_provenance(
    value: KinaseWorkflowAttritionProvenance | None,
) -> KinaseWorkflowAttritionProvenance | None:
    if value is None or isinstance(value, KinaseWorkflowAttritionProvenance):
        return value
    raise PhosPyInputError(
        "kinase_workflow_result.attrition_provenance must be "
        "KinaseWorkflowAttritionProvenance or None"
    )


def _require_instance(
    value: object,
    expected_type: type[object],
    *,
    field_name: str,
) -> None:
    if not isinstance(value, expected_type):
        raise PhosPyInputError(f"{field_name} must be {expected_type.__name__}")


def _require_optional_instance(
    value: object | None,
    expected_type: type[object],
    *,
    field_name: str,
) -> None:
    if value is not None and not isinstance(value, expected_type):
        raise PhosPyInputError(f"{field_name} must be {expected_type.__name__} or None")


def _require_mapping(value: object, *, field_name: str) -> Mapping[str, object]:
    return freeze_json_mapping(value, field_name=field_name)


def _require_non_empty_text(value: object, *, field_name: str) -> str:
    if not isinstance(value, str) or value.strip() == "":
        raise PhosPyInputError(f"{field_name} must be a non-empty string")
    return value


def _require_policy_outcome(value: object) -> str:
    if isinstance(value, str) and value in _KINASE_ATTRITION_POLICY_OUTCOMES:
        return value
    allowed = ", ".join(sorted(_KINASE_ATTRITION_POLICY_OUTCOMES))
    raise PhosPyInputError(
        "kinase_workflow_result.attrition_provenance.policy_outcome must be "
        f"one of: {allowed}"
    )


__all__ = [
    "ActivityMethodDiagnostics",
    "KinaseActivityResult",
    "KinaseEligibilityReport",
    "KinasePredictionResult",
    "KinaseScoringResult",
    "KinaseWorkflowAttritionProvenance",
    "KinaseWorkflowCaveat",
    "KinaseWorkflowPreprocessingAttritionSummary",
    "KinaseWorkflowResult",
    "KinaseWorkflowScoringAttritionSummary",
    "KinaseWorkflowSiteAttritionSummary",
    "KseaZScoreActivityDiagnostics",
    "SsgseaSubstrateEnrichmentActivityDiagnostics",
    "WeightedSubstrateActivityDiagnostics",
]
