from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from .analysis import KinaseActivityAnalyzer, KinaseActivityResult
from .dataset import CoreProcessingResult, PhosphoDataset


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
    ) -> None:
        self.dataset = dataset
        self.pred_mat = pred_mat.copy() if pred_mat is not None else None

    @classmethod
    def from_files(
        cls,
        total_path: str | Path,
        phospho_path: str | Path,
        pred_mat_path: str | Path | None = None,
    ) -> "PhosRPipeline":
        dataset = PhosphoDataset.from_files(total_path=total_path, phospho_path=phospho_path)
        pred_mat = pd.read_csv(pred_mat_path, index_col=0) if pred_mat_path is not None else None
        return cls(dataset=dataset, pred_mat=pred_mat)

    def run(self, outdir: str | Path | None = None) -> CoreOutputs:
        core = self.dataset.process_core()
        if outdir is not None:
            self.dataset.write_core_outputs(core, outdir)

        kinase_activity = None
        if self.pred_mat is not None:
            analyzer = KinaseActivityAnalyzer(self.pred_mat)
            kinase_activity = analyzer.analyze(core.site_matrix.matrix)
            if outdir is not None:
                analyzer.write_outputs(kinase_activity, outdir)

        return CoreOutputs(core=core, kinase_activity=kinase_activity)


def run_core_pipeline(
    total_path: str | Path,
    phospho_path: str | Path,
    outdir: str | Path,
    pred_mat_path: str | Path | None = None,
) -> CoreOutputs:
    pipeline = PhosRPipeline.from_files(
        total_path=total_path,
        phospho_path=phospho_path,
        pred_mat_path=pred_mat_path,
    )
    return pipeline.run(outdir=outdir)
