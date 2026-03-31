from __future__ import annotations

import json
import shutil
from collections.abc import Sequence
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from tempfile import TemporaryDirectory
from uuid import uuid4

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


def _package_version() -> str:
    try:
        from importlib.metadata import version

        return version("phospy")
    except Exception:
        return "unknown"


class PhosRPipeline:
    """High-level orchestration around dataset processing and kinase analysis."""

    def __init__(
        self,
        dataset: PhosphoDataset,
        pred_mat: pd.DataFrame | None = None,
        preprocessing_config: CorePreprocessingConfig | None = None,
        localization_threshold: float = 0.75,
        min_observed: int = 4,
        max_unmatched_fraction: float = 0.0,
        total_sentinel: float = 10.0,
        phospho_sentinel: float = 12.0,
    ) -> None:
        self.dataset = dataset
        self.pred_mat = pred_mat.copy() if pred_mat is not None else None
        self.preprocessing_config = preprocessing_config or CorePreprocessingConfig(
            localization_threshold=localization_threshold,
            min_observed=min_observed,
            max_unmatched_fraction=max_unmatched_fraction,
            total_sentinel=total_sentinel,
            phospho_sentinel=phospho_sentinel,
        )

    @classmethod
    def from_request(cls, request: CorePipelineRequest) -> PhosRPipeline:
        total_df, phospho_df = DatasetLoader().load(
            request.total_path,
            request.phospho_path,
            phospho_encoding=request.phospho_encoding,
        )
        dataset = PhosphoDataset.from_validated_inputs(
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
            preprocessing_config=CorePreprocessingConfig(
                localization_threshold=request.localization_threshold,
                min_observed=request.min_observed,
                max_unmatched_fraction=request.max_unmatched_fraction,
                total_sentinel=request.total_sentinel,
                phospho_sentinel=request.phospho_sentinel,
            ),
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

        kinase_activity = None
        if self.pred_mat is not None:
            kinase_activity = analyze_kinase_activity(
                self.pred_mat,
                core.site_matrix.matrix,
            )

        if outdir is not None:
            self._write_outputs_atomically(
                outdir=outdir,
                core=core,
                kinase_activity=kinase_activity,
            )

        return CoreOutputs(core=core, kinase_activity=kinase_activity)

    def _write_outputs_atomically(
        self,
        *,
        outdir: str | Path,
        core: CoreProcessingResult,
        kinase_activity: KinaseActivityResult | None,
    ) -> None:
        target_dir = Path(outdir)
        parent_dir = target_dir.parent
        parent_dir.mkdir(parents=True, exist_ok=True)

        with TemporaryDirectory(
            dir=parent_dir,
            prefix=f".{target_dir.name}.tmp-",
        ) as staging_dir_str:
            staging_dir = Path(staging_dir_str)
            CoreOutputWriter.write(core, staging_dir)
            if kinase_activity is not None:
                KinaseActivityWriter.write(kinase_activity, staging_dir)
            self._write_run_manifest(
                outdir=staging_dir,
                core=core,
                kinase_activity=kinase_activity,
            )
            self._publish_output_directory(
                staging_dir=staging_dir,
                target_dir=target_dir,
            )

    def _write_run_manifest(
        self,
        *,
        outdir: Path,
        core: CoreProcessingResult,
        kinase_activity: KinaseActivityResult | None,
    ) -> None:
        manifest = {
            "status": "success",
            "generated_at_utc": datetime.now(timezone.utc).isoformat(),
            "package_version": _package_version(),
            "has_kinase_activity": kinase_activity is not None,
            "core_rows": {
                "total_unique": int(core.total_unique.shape[0]),
                "total_filtered": int(core.total_filtered.shape[0]),
                "phospho_filtered": int(core.phospho_filtered.shape[0]),
                "phospho_corrected": int(core.phospho_corrected.shape[0]),
                "site_matrix": int(core.site_matrix.matrix.shape[0]),
            },
            "preprocessing_config": asdict(self.preprocessing_config),
        }
        (outdir / "run_manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )

    @staticmethod
    def _publish_output_directory(*, staging_dir: Path, target_dir: Path) -> None:
        if not target_dir.exists():
            staging_dir.replace(target_dir)
            return

        backup_dir = target_dir.with_name(f".{target_dir.name}.backup-{uuid4().hex}")
        target_dir.replace(backup_dir)
        try:
            staging_dir.replace(target_dir)
        except Exception:
            backup_dir.replace(target_dir)
            raise
        else:
            shutil.rmtree(backup_dir)
