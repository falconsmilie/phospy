"""Public result models."""
# pyright: reportMissingTypeStubs=false, reportUnnecessaryIsInstance=false
# pandas has no bundled stubs here; runtime boundary guards are intentionally retained.

from __future__ import annotations

from dataclasses import dataclass, field

import pandas as pd

from phospy._frame_ownership import export_optional_dataframe, own_optional_dataframe
from phospy.activities.models import KinaseActivityResult
from phospy.datasets.models import AnalysisReadyPhosphoDataset
from phospy.errors.validation import WorkflowValidationError
from phospy.prediction.models import KinasePredictionResult, KinaseScoringResult
from phospy.provenance.models import RunProvenance
from phospy.references.models import ReferenceBundle
from phospy.signalomes.models import (
    KinaseNetwork,
    SignalomeAlignmentDiagnostics,
    SignalomeAssignments,
    SignalomeModules,
    SignalomeModuleSelectionDiagnostics,
    SignalomeScorePreconditioningDiagnostics,
    default_signalome_alignment_diagnostics,
    default_signalome_module_selection_diagnostics,
    default_signalome_score_preconditioning_diagnostics,
)
from phospy.tables.signalome import (
    SignalomeProteinSiteContext,
    SignalomeSiteContext,
)


@dataclass(frozen=True, slots=True)
class KinaseWorkflowPreprocessingAttritionSummary:
    """Preprocessing-owned site attrition counters composed into kinase results."""

    input_rows: int
    rows_removed_during_preprocessing: int
    rows_removed_invalid_or_missing_site_identifiers: int
    duplicate_sites_merged_or_resolved: int
    output_rows: int


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
class KinaseWorkflowResult:
    """Top-level public kinase workflow result."""

    dataset: AnalysisReadyPhosphoDataset
    references: ReferenceBundle
    scoring_result: KinaseScoringResult
    prediction_result: KinasePredictionResult
    site_attrition_summary: KinaseWorkflowSiteAttritionSummary | None = None
    activity_result: KinaseActivityResult | None = None
    provenance: RunProvenance | None = None


@dataclass(frozen=True, slots=True, init=False)
class SignalomeWorkflowResult:
    """Top-level public signalome workflow result.

    `expanded_signalome` is a flattened optional DataFrame contract populated by
    the supported signalome executor lane. It includes focal-kinase rows with
    linked-kinase metadata, regulated module IDs, and selected site-membership
    rows with stable `site_order`. `score_preconditioning_diagnostics` reports
    downstream-score row preconditioning counts and the active
    `SignalomeConfig.validation.score_preconditioning_policy`.
    `site_membership` and
    `protein_site_context` provide optional signalome provenance sidecars for
    site-level and protein-level phosphosite context.

    Provenance in this object describes owned internal state at creation time.
    Public export helpers return defensive snapshots; mutating exported
    DataFrames does not mutate this owning result.
    """

    dataset: AnalysisReadyPhosphoDataset
    kinase_result: KinaseWorkflowResult
    module_assignments: SignalomeAssignments
    signalome_modules: SignalomeModules
    kinase_network: KinaseNetwork
    module_selection_diagnostics: SignalomeModuleSelectionDiagnostics = field(
        default_factory=default_signalome_module_selection_diagnostics
    )
    score_preconditioning_diagnostics: SignalomeScorePreconditioningDiagnostics = field(
        default_factory=default_signalome_score_preconditioning_diagnostics
    )
    alignment_diagnostics: SignalomeAlignmentDiagnostics = field(
        default_factory=default_signalome_alignment_diagnostics
    )
    provenance: RunProvenance | None = None
    _expanded_signalome: pd.DataFrame | None = field(init=False, repr=False)
    _site_membership: pd.DataFrame | None = field(init=False, repr=False)
    _protein_site_context: pd.DataFrame | None = field(init=False, repr=False)
    _init_payload: (
        tuple[pd.DataFrame | None, pd.DataFrame | None, pd.DataFrame | None, bool]
        | None
    ) = field(init=False, repr=False, default=None)

    def __init__(
        self,
        dataset: AnalysisReadyPhosphoDataset,
        kinase_result: KinaseWorkflowResult,
        module_assignments: SignalomeAssignments,
        signalome_modules: SignalomeModules,
        kinase_network: KinaseNetwork,
        module_selection_diagnostics: SignalomeModuleSelectionDiagnostics | None = None,
        score_preconditioning_diagnostics: SignalomeScorePreconditioningDiagnostics
        | None = None,
        alignment_diagnostics: SignalomeAlignmentDiagnostics | None = None,
        expanded_signalome: pd.DataFrame | None = None,
        site_membership: pd.DataFrame | None = None,
        protein_site_context: pd.DataFrame | None = None,
        provenance: RunProvenance | None = None,
        _assume_owned: bool = False,
    ) -> None:
        object.__setattr__(self, "dataset", dataset)
        object.__setattr__(self, "kinase_result", kinase_result)
        object.__setattr__(self, "module_assignments", module_assignments)
        object.__setattr__(self, "signalome_modules", signalome_modules)
        object.__setattr__(self, "kinase_network", kinase_network)
        object.__setattr__(
            self,
            "module_selection_diagnostics",
            (
                default_signalome_module_selection_diagnostics()
                if module_selection_diagnostics is None
                else module_selection_diagnostics
            ),
        )
        object.__setattr__(
            self,
            "score_preconditioning_diagnostics",
            (
                default_signalome_score_preconditioning_diagnostics()
                if score_preconditioning_diagnostics is None
                else score_preconditioning_diagnostics
            ),
        )
        object.__setattr__(
            self,
            "alignment_diagnostics",
            (
                default_signalome_alignment_diagnostics()
                if alignment_diagnostics is None
                else alignment_diagnostics
            ),
        )
        object.__setattr__(self, "provenance", provenance)
        object.__setattr__(
            self,
            "_init_payload",
            (
                expanded_signalome,
                site_membership,
                protein_site_context,
                _assume_owned,
            ),
        )
        self.__post_init__()

    def __post_init__(self) -> None:
        payload = self._init_payload
        if payload is None:
            raise WorkflowValidationError(
                "signalome_result internal initialization payload missing"
            )
        expanded_signalome, site_membership, protein_site_context, _assume_owned = (
            payload
        )
        expanded_signalome = own_optional_dataframe(
            expanded_signalome,
            field_name="signalome_result.expanded_signalome",
            error_type=WorkflowValidationError,
            assume_owned=_assume_owned,
        )
        site_membership = own_optional_dataframe(
            site_membership,
            field_name="signalome_result.site_membership",
            error_type=WorkflowValidationError,
            assume_owned=_assume_owned,
        )
        protein_site_context = own_optional_dataframe(
            protein_site_context,
            field_name="signalome_result.protein_site_context",
            error_type=WorkflowValidationError,
            assume_owned=_assume_owned,
        )
        if site_membership is not None:
            site_membership = SignalomeSiteContext(
                frame=site_membership,
                _assume_owned=True,
            ).frame
        if protein_site_context is not None:
            protein_site_context = SignalomeProteinSiteContext(
                frame=protein_site_context,
                _assume_owned=True,
            ).frame
        if self.provenance is not None and not isinstance(
            self.provenance, RunProvenance
        ):
            raise WorkflowValidationError(
                "signalome_result.provenance must be RunProvenance or None"
            )
        object.__setattr__(self, "_expanded_signalome", expanded_signalome)
        object.__setattr__(self, "_site_membership", site_membership)
        object.__setattr__(self, "_protein_site_context", protein_site_context)
        object.__setattr__(self, "_init_payload", None)

    @property
    def expanded_signalome(self) -> pd.DataFrame | None:
        return export_optional_dataframe(self._expanded_signalome)

    @property
    def site_membership(self) -> pd.DataFrame | None:
        return export_optional_dataframe(self._site_membership)

    @property
    def protein_site_context(self) -> pd.DataFrame | None:
        return export_optional_dataframe(self._protein_site_context)

    @classmethod
    def _from_owned(
        cls,
        *,
        dataset: AnalysisReadyPhosphoDataset,
        kinase_result: KinaseWorkflowResult,
        module_assignments: SignalomeAssignments,
        signalome_modules: SignalomeModules,
        kinase_network: KinaseNetwork,
        module_selection_diagnostics: SignalomeModuleSelectionDiagnostics | None = None,
        score_preconditioning_diagnostics: SignalomeScorePreconditioningDiagnostics
        | None = None,
        alignment_diagnostics: SignalomeAlignmentDiagnostics | None = None,
        expanded_signalome: pd.DataFrame | None = None,
        site_membership: pd.DataFrame | None = None,
        protein_site_context: pd.DataFrame | None = None,
        provenance: RunProvenance | None = None,
    ) -> SignalomeWorkflowResult:
        return cls(
            dataset=dataset,
            kinase_result=kinase_result,
            module_assignments=module_assignments,
            signalome_modules=signalome_modules,
            kinase_network=kinase_network,
            module_selection_diagnostics=module_selection_diagnostics,
            score_preconditioning_diagnostics=score_preconditioning_diagnostics,
            alignment_diagnostics=alignment_diagnostics,
            expanded_signalome=expanded_signalome,
            site_membership=site_membership,
            protein_site_context=protein_site_context,
            provenance=provenance,
            _assume_owned=True,
        )

    def to_dataframe(self) -> pd.DataFrame | None:
        """Return an expanded-signalome snapshot when available."""

        return export_optional_dataframe(self._expanded_signalome)

    def site_membership_dataframe(self) -> pd.DataFrame | None:
        """Return a site-membership snapshot when available."""

        return export_optional_dataframe(self._site_membership)

    def protein_site_context_dataframe(self) -> pd.DataFrame | None:
        """Return a protein-site context snapshot when available."""

        return export_optional_dataframe(self._protein_site_context)


__all__ = [
    "KinaseActivityResult",
    "KinasePredictionResult",
    "KinaseScoringResult",
    "KinaseWorkflowResult",
    "SignalomeWorkflowResult",
]
