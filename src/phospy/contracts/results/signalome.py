"""Public signalome workflow result contracts."""
# pyright: reportMissingTypeStubs=false, reportUnnecessaryIsInstance=false

from __future__ import annotations

from dataclasses import dataclass, field

import pandas as pd

from phospy.contracts.result_caveats import ResultCaveat, validate_result_caveats
from phospy.contracts.results.kinase import KinaseWorkflowResult
from phospy.errors.validation import ContractValidationError, PhosPyValidationError
from phospy.frames.ownership import export_optional_dataframe, own_optional_dataframe
from phospy.provenance.models import RunProvenance
from phospy.science.datasets.models import (
    AnalysisReadyPhosphoDataset,
    DatasetPreprocessingReport,
)
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
from phospy.tables.signalome import (
    SignalomeProteinSiteContext,
    SignalomeSiteContext,
)


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
    intentionally limited to ownership, provenance type, and local public table
    shape contracts. It does not run clustering, mapping, scoring, reference
    resolution, dataset identity validation, dataset repair, file export,
    plotting, or report formatting.
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
    caveats: tuple[ResultCaveat, ...] = ()
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
        caveats: tuple[ResultCaveat, ...] = (),
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
            "caveats",
            validate_result_caveats(
                caveats,
                field_name="signalome_result.caveats",
                error_type=ContractValidationError,
            ),
        )
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
            raise ContractValidationError(
                "signalome_result internal initialization payload missing"
            )
        expanded_signalome, site_membership, protein_site_context, _assume_owned = (
            payload
        )
        expanded_signalome = own_optional_dataframe(
            expanded_signalome,
            field_name="signalome_result.expanded_signalome",
            error_type=ContractValidationError,
            assume_owned=_assume_owned,
        )
        site_membership = own_optional_dataframe(
            site_membership,
            field_name="signalome_result.site_membership",
            error_type=ContractValidationError,
            assume_owned=_assume_owned,
        )
        protein_site_context = own_optional_dataframe(
            protein_site_context,
            field_name="signalome_result.protein_site_context",
            error_type=ContractValidationError,
            assume_owned=_assume_owned,
        )
        _validate_expanded_signalome_shape(expanded_signalome)
        if site_membership is not None:
            try:
                site_membership = SignalomeSiteContext(
                    frame=site_membership,
                    _assume_owned=True,
                ).frame
            except PhosPyValidationError as exc:
                raise ContractValidationError(str(exc)) from exc
        if protein_site_context is not None:
            try:
                protein_site_context = SignalomeProteinSiteContext(
                    frame=protein_site_context,
                    _assume_owned=True,
                ).frame
            except PhosPyValidationError as exc:
                raise ContractValidationError(str(exc)) from exc
        if self.provenance is not None and not isinstance(
            self.provenance, RunProvenance
        ):
            raise ContractValidationError(
                "signalome_result.provenance must be RunProvenance or None"
            )
        object.__setattr__(self, "_expanded_signalome", expanded_signalome)
        object.__setattr__(self, "_site_membership", site_membership)
        object.__setattr__(self, "_protein_site_context", protein_site_context)
        object.__setattr__(self, "_init_payload", None)

    @property
    def expanded_signalome(self) -> pd.DataFrame | None:
        """Return an expanded-signalome snapshot when available."""

        return export_optional_dataframe(self._expanded_signalome)

    @property
    def site_membership(self) -> pd.DataFrame | None:
        """Return a site-membership snapshot when available."""

        return export_optional_dataframe(self._site_membership)

    @property
    def protein_site_context(self) -> pd.DataFrame | None:
        """Return a protein-site context snapshot when available."""

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
        caveats: tuple[ResultCaveat, ...] = (),
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
            caveats=caveats,
            _assume_owned=True,
        )

    def to_dataframe(self) -> pd.DataFrame | None:
        """Return an expanded-signalome snapshot, not an export."""

        return export_optional_dataframe(self._expanded_signalome)

    def site_membership_dataframe(self) -> pd.DataFrame | None:
        """Return a site-membership snapshot, not an export."""

        return export_optional_dataframe(self._site_membership)

    def protein_site_context_dataframe(self) -> pd.DataFrame | None:
        """Return a protein-site context snapshot, not an export."""

        return export_optional_dataframe(self._protein_site_context)


def _validate_expanded_signalome_shape(
    expanded_signalome: pd.DataFrame | None,
) -> None:
    missing = [
        column
        for column in (SITE_KEY_COLUMN, DISPLAY_ID_COLUMN)
        if expanded_signalome is not None and column not in expanded_signalome.columns
    ]
    if missing:
        joined = ", ".join(missing)
        raise ContractValidationError(
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
    _validate_non_empty_site_identity_columns(
        site_rows,
        field_name="signalome_result.expanded_signalome",
    )


def _validate_non_empty_site_identity_columns(
    table: pd.DataFrame,
    *,
    field_name: str,
) -> None:
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
            raise ContractValidationError(
                f"{field_name} site rows require non-empty "
                f"{column_name} values; invalid_count={invalid_count}"
            )


__all__ = [
    "SignalomeWorkflowResult",
]
