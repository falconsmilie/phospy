from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from .analysis import KinaseActivityResult, analyze_kinase_activity
from .constants import ComparisonSpec
from .core_processing import CorePreprocessingConfig, CoreProcessingResult
from .dataset import PhosphoDataset
from .dataset_loader import DatasetLoader
from .io import load_pred_mat
from .validation.requests import CorePipelineRequest
from .writers import CoreOutputWriter, KinaseActivityWriter


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
        total_sentinel: float = 10.0,
        phospho_sentinel: float = 12.0,
    ) -> None:
        self.dataset = dataset
        self.pred_mat = pred_mat.copy() if pred_mat is not None else None
        self.preprocessing_config = CorePreprocessingConfig(
            localization_threshold=localization_threshold,
            min_observed=min_observed,
            max_unmatched_fraction=max_unmatched_fraction,
            total_sentinel=total_sentinel,
            phospho_sentinel=phospho_sentinel,
        )
        self.localization_threshold = localization_threshold
        self.min_observed = min_observed
        self.max_unmatched_fraction = max_unmatched_fraction
        self.total_sentinel = total_sentinel
        self.phospho_sentinel = phospho_sentinel

    @classmethod
    def from_request(cls, request: CorePipelineRequest) -> PhosRPipeline:
        total_df, phospho_df = DatasetLoader().load(
            request.total_path,
            request.phospho_path,
            phospho_encoding=request.phospho_encoding,
        )
        dataset = PhosphoDataset._from_validated_frames(
            total_df=total_df,
            phospho_df=phospho_df,
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
            total_sentinel=request.total_sentinel,
            phospho_sentinel=request.phospho_sentinel,
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
        total_sentinel: float = 10.0,
        phospho_sentinel: float = 12.0,
    ) -> PhosRPipeline:
        request = CorePipelineRequest.validate_request(
            total_path=Path(total_path),
            phospho_path=Path(phospho_path),
            pred_mat_path=Path(pred_mat_path) if pred_mat_path is not None else None,
            phospho_encoding=phospho_encoding,
            comparisons=tuple(comparisons) if comparisons is not None else None,
            localization_threshold=localization_threshold,
            min_observed=min_observed,
            total_sentinel=total_sentinel,
            phospho_sentinel=phospho_sentinel,
            max_unmatched_fraction=max_unmatched_fraction,
        )
        return cls.from_request(request)

    def run(self, outdir: str | Path | None = None) -> CoreOutputs:
        core = self.dataset.process_core(config=self.preprocessing_config)
        if outdir is not None:
            CoreOutputWriter.write(core, outdir)

        kinase_activity = None
        if self.pred_mat is not None:
            kinase_activity = analyze_kinase_activity(
                self.pred_mat, core.site_matrix.matrix
            )
            if outdir is not None:
                KinaseActivityWriter.write(kinase_activity, outdir)

        return CoreOutputs(core=core, kinase_activity=kinase_activity)
