from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from .analysis import KinaseActivityAnalyzer, KinaseActivityResult
from .constants import ComparisonSpec
from .dataset import CoreProcessingResult, PhosphoDataset
from .io import load_phospho_table, load_pred_mat, load_total_table
from .validation.requests import CorePipelineRequest


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
    ) -> None:
        self.dataset = dataset
        self.pred_mat = pred_mat.copy() if pred_mat is not None else None
        self.localization_threshold = localization_threshold
        self.min_observed = min_observed

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
    ) -> PhosRPipeline:
        request = CorePipelineRequest.validate_request(
            total_path=Path(total_path),
            phospho_path=Path(phospho_path),
            pred_mat_path=Path(pred_mat_path) if pred_mat_path is not None else None,
            phospho_encoding=phospho_encoding,
            comparisons=tuple(comparisons) if comparisons is not None else None,
            localization_threshold=localization_threshold,
            min_observed=min_observed,
        )
        return cls.from_request(request)

    def run(self, outdir: str | Path | None = None) -> CoreOutputs:
        core = self.dataset.process_core(
            localization_threshold=self.localization_threshold,
            min_observed=self.min_observed,
        )
        if outdir is not None:
            self.dataset.write_core_outputs(core, outdir)

        kinase_activity = None
        if self.pred_mat is not None:
            analyzer = KinaseActivityAnalyzer(self.pred_mat)
            kinase_activity = analyzer.analyze(core.site_matrix.matrix)
            if outdir is not None:
                analyzer.write_outputs(kinase_activity, outdir)

        return CoreOutputs(core=core, kinase_activity=kinase_activity)
