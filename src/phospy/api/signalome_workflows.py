from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Protocol, cast

import pandas as pd

from ..datasets.models import (
    SiteToProteinResolutionDiagnostics,
    SiteToProteinResolutionResult,
)
from ..prediction.results import KinasePredictionResult, PredMatResult
from ..prediction.scoring import KinaseScoringResult
from ..signalomes.analysis import execute_signalome_inputs
from ..signalomes.results import SignalomeResult
from ..validation.requests.signalome import SignalomeInputs, validate_signalome_request
from .contracts import SignalomeRunConfig

__all__ = ["SignalomeWorkflow"]


class _AnalysisReadyDatasetProtocol(Protocol):
    phospho_matrix: pd.DataFrame

    def resolve_site_to_protein_mapping(
        self,
        *,
        metadata_columns: Sequence[str] | None = None,
        fallback_policy: str = "strict",
        allow_gene_symbol_fallback: bool = False,
        allow_ambiguous_fallback: bool = False,
    ) -> pd.Series: ...


class _AnalysisReadyDatasetWithDiagnosticsProtocol(
    _AnalysisReadyDatasetProtocol, Protocol
):
    def resolve_site_to_protein_mapping_with_diagnostics(
        self,
        *,
        metadata_columns: Sequence[str] | None = None,
        fallback_policy: str = "strict",
        allow_gene_symbol_fallback: bool = False,
        allow_ambiguous_fallback: bool = False,
    ) -> SiteToProteinResolutionResult: ...


def _coerce_analysis_ready_dataset(
    dataset: object,
) -> _AnalysisReadyDatasetProtocol:
    if not hasattr(dataset, "phospho_matrix") or not hasattr(
        dataset,
        "resolve_site_to_protein_mapping",
    ):
        msg = (
            "dataset must be an AnalysisReadyPhosphoDataset for "
            "run_from_analysis_ready()."
        )
        raise TypeError(msg)
    return cast(_AnalysisReadyDatasetProtocol, dataset)


class SignalomeWorkflow:
    """Construct signalomes from validated scoring and prediction outputs."""

    def run(
        self,
        *,
        scoring_result: KinaseScoringResult,
        prediction_result: KinasePredictionResult | PredMatResult,
        expression_matrix: pd.DataFrame,
        kinases_of_interest: Sequence[str],
        site_to_protein: Mapping[str, str] | None = None,
        config: SignalomeRunConfig | None = None,
    ) -> SignalomeResult:
        resolved_config = SignalomeRunConfig.from_value(config)
        request = validate_signalome_request(
            scoring_result=scoring_result,
            prediction_result=prediction_result,
            expression_matrix=expression_matrix,
            kinases_of_interest=kinases_of_interest,
            site_to_protein=site_to_protein,
            kinase_network_threshold=resolved_config.kinase_network_threshold,
            kinase_network_policy=resolved_config.kinase_network_policy,
            assignment_policy=resolved_config.assignment_policy,
            signalome_cutoff=resolved_config.signalome_cutoff,
            module_count=resolved_config.module_count,
            min_kinase_module_share_percent=(
                resolved_config.min_kinase_module_share_percent
            ),
            module_selection_policy=resolved_config.module_selection_policy,
        )
        return execute_signalome_inputs(request)

    def run_from_analysis_ready(
        self,
        *,
        dataset: object,
        scoring_result: KinaseScoringResult,
        prediction_result: KinasePredictionResult | PredMatResult,
        kinases_of_interest: Sequence[str],
        site_to_protein: Mapping[str, str] | None = None,
        metadata_protein_columns: Sequence[str] | None = None,
        metadata_fallback_policy: str = "strict",
        allow_gene_symbol_fallback: bool = False,
        allow_ambiguous_metadata_mapping: bool = False,
        config: SignalomeRunConfig | None = None,
    ) -> SignalomeResult:
        """Run signalome analysis from an analysis-ready dataset boundary.

        Defaults to strict metadata resolution requiring a ``protein_id`` column.
        To opt in to metadata fallback columns, set
        ``metadata_fallback_policy="metadata"``.
        """
        analysis_ready_dataset = _coerce_analysis_ready_dataset(dataset)

        resolution_diagnostics = None
        if site_to_protein is not None:
            resolved_site_to_protein = dict(site_to_protein)
        elif hasattr(
            analysis_ready_dataset, "resolve_site_to_protein_mapping_with_diagnostics"
        ):
            dataset_with_diagnostics = cast(
                _AnalysisReadyDatasetWithDiagnosticsProtocol,
                analysis_ready_dataset,
            )
            resolution_result = dataset_with_diagnostics.resolve_site_to_protein_mapping_with_diagnostics(
                metadata_columns=metadata_protein_columns,
                fallback_policy=metadata_fallback_policy,
                allow_gene_symbol_fallback=allow_gene_symbol_fallback,
                allow_ambiguous_fallback=allow_ambiguous_metadata_mapping,
            )
            resolved_site_to_protein = dict(resolution_result.mapping.to_dict())
            resolution_diagnostics = resolution_result.diagnostics
        else:
            resolved_mapping = analysis_ready_dataset.resolve_site_to_protein_mapping(
                metadata_columns=metadata_protein_columns,
                fallback_policy=metadata_fallback_policy,
                allow_gene_symbol_fallback=allow_gene_symbol_fallback,
                allow_ambiguous_fallback=allow_ambiguous_metadata_mapping,
            )
            resolved_site_to_protein = dict(resolved_mapping.to_dict())
            diagnostics_attr = resolved_mapping.attrs.get(
                "site_to_protein_resolution_diagnostics"
            )
            if isinstance(diagnostics_attr, SiteToProteinResolutionDiagnostics):
                resolution_diagnostics = diagnostics_attr

        result = self.run(
            scoring_result=scoring_result,
            prediction_result=prediction_result,
            expression_matrix=analysis_ready_dataset.phospho_matrix,
            kinases_of_interest=kinases_of_interest,
            site_to_protein=resolved_site_to_protein,
            config=config,
        )
        if site_to_protein is None:
            result.attach_site_to_protein_resolution_diagnostics(resolution_diagnostics)
        return result

    def run_validated(
        self,
        request: SignalomeInputs,
    ) -> SignalomeResult:
        return execute_signalome_inputs(request)
