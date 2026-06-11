"""Public result models."""
# pyright: reportMissingTypeStubs=false, reportUnnecessaryIsInstance=false
# pandas has no bundled stubs here; runtime boundary guards are intentionally retained.

from __future__ import annotations

from dataclasses import dataclass, field

import pandas as pd

from phospy.errors.validation import WorkflowValidationError
from phospy.frames.ownership import export_optional_dataframe, own_optional_dataframe
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
from phospy.science.differential.models import DifferentialAnalysisResult
from phospy.science.prediction.models import KinasePredictionResult, KinaseScoringResult
from phospy.science.references.models import ReferenceBundle
from phospy.science.signalomes.constants import (
    DISPLAY_ID_COLUMN,
    EXPANDED_SIGNALOME_ROW_KIND_COLUMN,
    EXPANDED_SIGNALOME_ROW_KIND_SITE,
    SITE_KEY_COLUMN,
)
from phospy.science.signalomes.models import (
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
from phospy.science.sites.validation import require_site_key_series
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
    """Compact, user-facing kinase workflow eligibility counters."""

    total_dataset_sites: int
    sequence_complete_sites: int
    localisation_eligible_sites: int | None
    reference_overlap_sites: int
    excluded_no_reference_match: int
    excluded_low_localisation: int | None
    eligible_kinases: int
    excluded_kinases_below_min_substrates: int


@dataclass(frozen=True, slots=True)
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

    @property
    def input_dataset_preprocessing_report(self) -> DatasetPreprocessingReport | None:
        """Return preprocessing provenance of the input analysis-ready dataset."""

        return self.dataset.preprocessing_report


@dataclass(frozen=True, slots=True, init=False)
class SignalomeWorkflowResult:
    """Top-level public signalome workflow result.

    This is a workflow-owned result object. Direct construction is supported for
    bundle reconstruction and tests, but the signalome workflow is the
    recommended public construction path for scientific coherence across module
    assignment, module summary, kinase-network, and context sidecars.

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
    DataFrames does not mutate this owning result. Constructor validation is
    intentionally limited to ownership, provenance type, and public table
    identity contracts. It does not run clustering, mapping, scoring, reference
    resolution, or dataset repair.
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
        _validate_signalome_result_site_level_identity(
            dataset=self.dataset,
            module_assignments=self.module_assignments,
            expanded_signalome=expanded_signalome,
            site_membership=site_membership,
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

    @property
    def input_dataset_preprocessing_report(self) -> DatasetPreprocessingReport | None:
        """Return preprocessing provenance of the input analysis-ready dataset."""

        return self.dataset.preprocessing_report

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


def _validate_signalome_result_site_level_identity(
    *,
    dataset: AnalysisReadyPhosphoDataset,
    module_assignments: SignalomeAssignments,
    expanded_signalome: pd.DataFrame | None,
    site_membership: pd.DataFrame | None,
) -> None:
    dataset_identity = _signalome_dataset_identity_lookup(dataset)
    _validate_site_level_signalome_rows(
        table=module_assignments.table,
        field_name="signalome_result.module_assignments.table",
        dataset_identity=dataset_identity,
    )
    _validate_expanded_signalome_identity(
        expanded_signalome,
        dataset_identity=dataset_identity,
    )
    _validate_site_level_signalome_rows(
        table=site_membership,
        field_name="signalome_result.site_membership",
        dataset_identity=dataset_identity,
    )


def _validate_expanded_signalome_identity(
    expanded_signalome: pd.DataFrame | None,
    *,
    dataset_identity: dict[str, str],
) -> None:
    missing = [
        column
        for column in (SITE_KEY_COLUMN, DISPLAY_ID_COLUMN)
        if expanded_signalome is not None and column not in expanded_signalome.columns
    ]
    if missing:
        joined = ", ".join(missing)
        raise WorkflowValidationError(
            f"signalome_result.expanded_signalome is missing required columns: {joined}"
        )
    if expanded_signalome is None or expanded_signalome.empty:
        return
    site_rows = expanded_signalome
    if EXPANDED_SIGNALOME_ROW_KIND_COLUMN in expanded_signalome.columns:
        site_rows = expanded_signalome.loc[
            expanded_signalome.loc[:, EXPANDED_SIGNALOME_ROW_KIND_COLUMN].astype(str)
            == EXPANDED_SIGNALOME_ROW_KIND_SITE,
            :,
        ]
    _validate_site_level_signalome_rows(
        table=site_rows,
        field_name="signalome_result.expanded_signalome",
        dataset_identity=dataset_identity,
    )


def _validate_site_level_signalome_rows(
    *,
    table: pd.DataFrame | None,
    field_name: str,
    dataset_identity: dict[str, str],
) -> None:
    if table is None:
        return
    missing = [
        column for column in (SITE_KEY_COLUMN, DISPLAY_ID_COLUMN) if column not in table
    ]
    if missing:
        joined = ", ".join(missing)
        raise WorkflowValidationError(
            f"{field_name} is missing required columns: {joined}"
        )
    if table.empty:
        return
    for column_name in (SITE_KEY_COLUMN, DISPLAY_ID_COLUMN):
        invalid_count = int(
            sum(
                1
                for value in table.loc[:, column_name].tolist()
                if not isinstance(value, str) or value.strip() == ""
            )
        )
        if invalid_count:
            raise WorkflowValidationError(
                f"{field_name} site rows require non-empty "
                f"{column_name} values; invalid_count={invalid_count}"
            )
    site_keys = table.loc[:, SITE_KEY_COLUMN].astype(str)
    require_site_key_series(
        site_keys,
        field_name=f"{field_name}.{SITE_KEY_COLUMN}",
        error_type=WorkflowValidationError,
    )
    display_ids = table.loc[:, DISPLAY_ID_COLUMN].astype(str)
    missing_site_keys = [
        site_key
        for site_key in dict.fromkeys(site_keys.tolist())
        if site_key not in dataset_identity
    ]
    if missing_site_keys:
        preview = ", ".join(repr(value) for value in missing_site_keys[:5])
        suffix = "" if len(missing_site_keys) <= 5 else " ..."
        raise WorkflowValidationError(
            f"{field_name}.{SITE_KEY_COLUMN} values must align to "
            "signalome_result.dataset; missing_site_keys="
            f"{preview}{suffix}"
        )
    mismatches = [
        f"{site_key!r}: observed={display_id!r}, expected={dataset_identity[site_key]!r}"
        for site_key, display_id in zip(
            site_keys.tolist(),
            display_ids.tolist(),
            strict=True,
        )
        if dataset_identity[site_key] != display_id
    ]
    if mismatches:
        preview = "; ".join(mismatches[:5])
        suffix = "" if len(mismatches) <= 5 else " ..."
        raise WorkflowValidationError(
            f"{field_name}.{DISPLAY_ID_COLUMN} values must match "
            "signalome_result.dataset.site_metadata.display_id for each site_key; "
            f"mismatches={preview}{suffix}"
        )


def _signalome_dataset_identity_lookup(
    dataset: AnalysisReadyPhosphoDataset,
) -> dict[str, str]:
    site_metadata = dataset._borrow_site_metadata_frame()
    missing = [
        column
        for column in (SITE_KEY_COLUMN, DISPLAY_ID_COLUMN)
        if column not in site_metadata.columns
    ]
    if missing:
        joined = ", ".join(missing)
        raise WorkflowValidationError(
            "signalome_result.dataset.site_metadata is missing required columns: "
            f"{joined}"
        )
    site_keys = site_metadata.loc[:, SITE_KEY_COLUMN].astype(str)
    require_site_key_series(
        site_keys,
        field_name=f"signalome_result.dataset.site_metadata.{SITE_KEY_COLUMN}",
        error_type=WorkflowValidationError,
    )
    display_ids = site_metadata.loc[:, DISPLAY_ID_COLUMN]
    invalid_display_id_count = int(
        sum(
            1
            for value in display_ids.tolist()
            if not isinstance(value, str) or value.strip() == ""
        )
    )
    if invalid_display_id_count:
        raise WorkflowValidationError(
            "signalome_result.dataset.site_metadata.display_id values must be "
            "non-empty strings; "
            f"invalid_count={invalid_display_id_count}"
        )
    return {
        site_key: display_id
        for site_key, display_id in zip(
            site_keys.tolist(),
            display_ids.astype(str).tolist(),
            strict=True,
        )
    }


__all__ = [
    "ActivityMethodDiagnostics",
    "DifferentialAnalysisResult",
    "KinaseActivityResult",
    "KinasePredictionResult",
    "KinaseScoringResult",
    "KinaseWorkflowResult",
    "KseaZScoreActivityDiagnostics",
    "SsgseaSubstrateEnrichmentActivityDiagnostics",
    "SignalomeWorkflowResult",
    "WeightedSubstrateActivityDiagnostics",
]
