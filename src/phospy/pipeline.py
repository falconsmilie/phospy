from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from .analysis import KinaseActivityAnalyzer, KinaseActivityResult
from .constants import ComparisonSpec
from .dataset import CoreProcessingResult, PhosphoDataset
from .io import load_phospho_table, load_pred_mat, load_total_table
from .validation.requests import CorePipelineRequest, KinaseActivityRequest
from .validation.tables import PredMatSchema


@dataclass(slots=True)
class CoreOutputs:
    core: CoreProcessingResult
    kinase_activity: KinaseActivityResult | None = None


class PhosRPipeline:
    """High-level orchestration around dataset processing and kinase analysis."""

    def __init__(
        self,
        dataset: PhosphoDataset,
        pred_mat: pd.DataFrame | None = None,
        localization_threshold: float = 0.75,
        min_observed: int = 4,
        max_unmatched_fraction: float = 0.0,
        kinase_activity_threshold: float = 0.6,
        kinase_activity_min_substrates: int = 3,
        kinase_activity_top_n_substrates: int = 20,
    ) -> None:
        activity_request = KinaseActivityRequest.validate_request(
            threshold=kinase_activity_threshold,
            min_substrates=kinase_activity_min_substrates,
            top_n_substrates=kinase_activity_top_n_substrates,
        )
        self.dataset = dataset
        self.pred_mat = (
            PredMatSchema.validate(pred_mat, context="pred_mat")
            if pred_mat is not None
            else None
        )
        self.localization_threshold = localization_threshold
        self.min_observed = min_observed
        self.max_unmatched_fraction = max_unmatched_fraction
        self.kinase_activity_threshold = activity_request.threshold
        self.kinase_activity_min_substrates = activity_request.min_substrates
        self.kinase_activity_top_n_substrates = activity_request.top_n_substrates

    @classmethod
    def from_request(cls, request: CorePipelineRequest) -> PhosRPipeline:
        dataset = PhosphoDataset(
            total_df=load_total_table(request.total_path),
            phospho_df=load_phospho_table(
                request.phospho_path,
                encoding=request.phospho_encoding,
            ),
            comparisons=request.comparisons,
        )
        pred_mat = (
            load_pred_mat(request.pred_mat_path)
            if request.pred_mat_path is not None
            else None
        )
        return cls(
            dataset=dataset,
            pred_mat=pred_mat,
            localization_threshold=request.localization_threshold,
            min_observed=request.min_observed,
            max_unmatched_fraction=request.max_unmatched_fraction,
            kinase_activity_threshold=request.kinase_activity_threshold,
            kinase_activity_min_substrates=request.kinase_activity_min_substrates,
            kinase_activity_top_n_substrates=request.kinase_activity_top_n_substrates,
        )

    @classmethod
    def from_files(
        cls,
        total_path: str | Path,
        phospho_path: str | Path,
        pred_mat_path: str | Path | None = None,
        comparisons: Sequence[ComparisonSpec] | None = None,
        phospho_encoding: str | None = None,
        localization_threshold: float = 0.75,
        min_observed: int = 4,
        max_unmatched_fraction: float = 0.0,
        kinase_activity_threshold: float = 0.6,
        kinase_activity_min_substrates: int = 3,
        kinase_activity_top_n_substrates: int = 20,
    ) -> PhosRPipeline:
        request = CorePipelineRequest.validate_request(
            total_path=Path(total_path),
            phospho_path=Path(phospho_path),
            pred_mat_path=Path(pred_mat_path) if pred_mat_path is not None else None,
            phospho_encoding=phospho_encoding,
            comparisons=tuple(comparisons) if comparisons is not None else None,
            localization_threshold=localization_threshold,
            min_observed=min_observed,
            max_unmatched_fraction=max_unmatched_fraction,
            kinase_activity_threshold=kinase_activity_threshold,
            kinase_activity_min_substrates=kinase_activity_min_substrates,
            kinase_activity_top_n_substrates=kinase_activity_top_n_substrates,
        )
        return cls.from_request(request)

    def run(self, outdir: str | Path | None = None) -> CoreOutputs:
        core = self.dataset.process_core(
            localization_threshold=self.localization_threshold,
            min_observed=self.min_observed,
            max_unmatched_fraction=self.max_unmatched_fraction,
        )
        if outdir is not None:
            self.dataset.write_core_outputs(core, outdir)

        kinase_activity = None
        if self.pred_mat is not None:
            analyzer = KinaseActivityAnalyzer.from_validated_pred_mat(self.pred_mat)
            activity_request = KinaseActivityRequest.validate_request(
                threshold=self.kinase_activity_threshold,
                min_substrates=self.kinase_activity_min_substrates,
                top_n_substrates=self.kinase_activity_top_n_substrates,
            )
            kinase_activity = analyzer.analyze(
                core.site_matrix.matrix,
                request=activity_request,
            )
            if outdir is not None:
                analyzer.write_outputs(kinase_activity, outdir)

        return CoreOutputs(core=core, kinase_activity=kinase_activity)
