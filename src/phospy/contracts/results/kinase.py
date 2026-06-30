"""Public kinase workflow result contracts."""
# pyright: reportMissingTypeStubs=false

from __future__ import annotations

from dataclasses import dataclass, field

import pandas as pd

from phospy.frames.ownership import export_optional_dataframe
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
from phospy.tables.kinase import KinaseSubstrateContributionTable


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
    activity_result: KinaseActivityResult | None = None
    provenance: RunProvenance | None = None
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
        activity_result: KinaseActivityResult | None = None,
        provenance: RunProvenance | None = None,
        substrate_contributions: pd.DataFrame | None = None,
        *,
        _assume_owned: bool = False,
    ) -> None:
        object.__setattr__(self, "dataset", dataset)
        object.__setattr__(self, "references", references)
        object.__setattr__(self, "scoring_result", scoring_result)
        object.__setattr__(self, "prediction_result", prediction_result)
        object.__setattr__(self, "eligibility_report", eligibility_report)
        object.__setattr__(self, "site_attrition_summary", site_attrition_summary)
        object.__setattr__(self, "activity_result", activity_result)
        object.__setattr__(self, "provenance", provenance)
        object.__setattr__(
            self,
            "_substrate_contributions",
            _own_optional_kinase_substrate_contributions(
                substrate_contributions,
                assume_owned=_assume_owned,
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


__all__ = [
    "ActivityMethodDiagnostics",
    "KinaseActivityResult",
    "KinaseEligibilityReport",
    "KinasePredictionResult",
    "KinaseScoringResult",
    "KinaseWorkflowPreprocessingAttritionSummary",
    "KinaseWorkflowResult",
    "KinaseWorkflowScoringAttritionSummary",
    "KinaseWorkflowSiteAttritionSummary",
    "KseaZScoreActivityDiagnostics",
    "SsgseaSubstrateEnrichmentActivityDiagnostics",
    "WeightedSubstrateActivityDiagnostics",
]
