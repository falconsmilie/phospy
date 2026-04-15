from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Protocol, cast

import pandas as pd

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

        resolved_site_to_protein = (
            dict(site_to_protein)
            if site_to_protein is not None
            else analysis_ready_dataset.resolve_site_to_protein_mapping(
                metadata_columns=metadata_protein_columns,
                fallback_policy=metadata_fallback_policy,
                allow_gene_symbol_fallback=allow_gene_symbol_fallback,
                allow_ambiguous_fallback=allow_ambiguous_metadata_mapping,
            ).to_dict()
        )
        return self.run(
            scoring_result=scoring_result,
            prediction_result=prediction_result,
            expression_matrix=analysis_ready_dataset.phospho_matrix,
            kinases_of_interest=kinases_of_interest,
            site_to_protein=resolved_site_to_protein,
            config=config,
        )

    def run_validated(
        self,
        request: SignalomeInputs,
    ) -> SignalomeResult:
        return execute_signalome_inputs(request)
