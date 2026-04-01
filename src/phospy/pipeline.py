from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd

from .analysis import KinaseActivityAnalyzer, KinaseActivityResult
from .constants import ComparisonSpec
from .core_processing import CorePreprocessingConfig, CoreProcessingResult
from .dataset import PhosphoDataset
from .dataset_loader import DatasetLoader
from .io import load_pred_mat
from .publishing import OutputPublisher, RunManifestWriter
from .validation.errors import RequestValidationError, TableSchemaError
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
        preprocessing_config: CorePreprocessingConfig | None = None,
        localization_threshold: float = 0.75,
        min_observed: int = 4,
        max_unmatched_fraction: float = 0.0,
        total_sentinel: float = 10.0,
        phospho_sentinel: float = 12.0,
        *,
        manifest_writer: RunManifestWriter | None = None,
        output_publisher: OutputPublisher | None = None,
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
        self.kinase_activity_analyzer = KinaseActivityAnalyzer()
        self.manifest_writer = manifest_writer or RunManifestWriter()
        self.output_publisher = output_publisher or OutputPublisher()

    @classmethod
    def from_request(cls, request: CorePipelineRequest) -> PhosRPipeline:
        validated_inputs = DatasetLoader().load(
            request.total_path,
            request.phospho_path,
            phospho_encoding=request.phospho_encoding,
        )
        dataset = PhosphoDataset.from_validated_inputs(
            validated_inputs,
            comparisons=request.comparisons,
        )
        pred_mat = None
        if request.pred_mat_path is not None:
            try:
                pred_mat = load_pred_mat(request.pred_mat_path)
            except TableSchemaError:
                raise
            except (
                OSError,
                UnicodeError,
                pd.errors.ParserError,
                pd.errors.EmptyDataError,
            ) as error:
                msg = (
                    f"Invalid core pipeline request: pred_mat_path: "
                    f"unable to read pred_mat ({request.pred_mat_path}): {error}"
                )
                raise RequestValidationError(msg) from error
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
        core = self.dataset.preprocessing.run(config=self.preprocessing_config)

        kinase_activity = None
        if self.pred_mat is not None:
            kinase_activity = self.kinase_activity_analyzer.analyze(
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
            CoreOutputWriter().write(core, staging_dir)
            if kinase_activity is not None:
                KinaseActivityWriter.write(kinase_activity, staging_dir)
            self.manifest_writer.write(
                outdir=staging_dir,
                core=core,
                kinase_activity=kinase_activity,
                preprocessing_config=self.preprocessing_config,
            )
            self.output_publisher.publish(
                staging_dir=staging_dir,
                target_dir=target_dir,
            )
