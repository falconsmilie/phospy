from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd

from .analysis import KinaseActivityAnalyzer, KinaseActivityResult
from .constants import ComparisonSpec
from .core_processing import CorePreprocessingConfig, CoreProcessingResult
from .dataset import PhosphoDataset
from .dataset_loader import DatasetLoader
from .dataset_schema import DatasetSchema
from .io import load_pred_mat
from .publishing import OutputPublisher, RunManifestWriter
from .validation.errors import RequestValidationError, TableSchemaError
from .validation.requests import CorePipelineRequest
from .writers import CoreOutputWriter, KinaseActivityWriter


@dataclass(slots=True)
class CoreOutputs:
    core: CoreProcessingResult
    kinase_activity: KinaseActivityResult | None = None


@dataclass(frozen=True, slots=True)
class PipelineConstructionInputs:
    dataset: PhosphoDataset
    pred_mat: pd.DataFrame | None
    preprocessing_config: CorePreprocessingConfig


class PipelineRequestLoader:
    def __init__(
        self,
        *,
        dataset_loader_factory: Callable[[DatasetSchema], DatasetLoader] | None = None,
        pred_mat_loader: Callable[[str | Path], pd.DataFrame] | None = None,
    ) -> None:
        self.dataset_loader_factory = dataset_loader_factory or DatasetLoader
        self.pred_mat_loader = (
            load_pred_mat if pred_mat_loader is None else pred_mat_loader
        )

    def load(self, request: CorePipelineRequest) -> PipelineConstructionInputs:
        validated_inputs = self.dataset_loader_factory(
            schema=request.dataset_schema
        ).load(
            request.total_path,
            request.phospho_path,
            phospho_encoding=request.phospho_encoding,
        )
        dataset = PhosphoDataset.from_validated_inputs(
            validated_inputs,
            comparisons=request.comparisons,
        )
        return PipelineConstructionInputs(
            dataset=dataset,
            pred_mat=self._load_pred_mat(request),
            preprocessing_config=CorePreprocessingConfig(
                localization_threshold=request.localization_threshold,
                min_observed=request.min_observed,
                max_unmatched_fraction=request.max_unmatched_fraction,
                total_sentinel=request.total_sentinel,
                phospho_sentinel=request.phospho_sentinel,
            ),
        )

    def _load_pred_mat(self, request: CorePipelineRequest) -> pd.DataFrame | None:
        if request.pred_mat_path is None:
            return None
        try:
            return self.pred_mat_loader(request.pred_mat_path)
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


class PipelineExecutionRunner:
    def __init__(self, *, kinase_activity_analyzer: KinaseActivityAnalyzer) -> None:
        self.kinase_activity_analyzer = kinase_activity_analyzer

    def run(
        self,
        *,
        dataset: PhosphoDataset,
        pred_mat: pd.DataFrame | None,
        preprocessing_config: CorePreprocessingConfig,
    ) -> CoreOutputs:
        core = dataset.preprocessing.run(config=preprocessing_config)

        kinase_activity = None
        if pred_mat is not None:
            kinase_activity = self.kinase_activity_analyzer.analyze(
                pred_mat,
                core.site_matrix.matrix,
            )

        return CoreOutputs(core=core, kinase_activity=kinase_activity)


class PipelineOutputCoordinator:
    def publish(
        self,
        *,
        outdir: str | Path,
        outputs: CoreOutputs,
        preprocessing_config: CorePreprocessingConfig,
        manifest_writer: RunManifestWriter,
        output_publisher: OutputPublisher,
    ) -> None:
        target_dir = Path(outdir)
        parent_dir = target_dir.parent
        parent_dir.mkdir(parents=True, exist_ok=True)

        with TemporaryDirectory(
            dir=parent_dir,
            prefix=f".{target_dir.name}.tmp-",
        ) as staging_dir_str:
            staging_dir = Path(staging_dir_str)
            CoreOutputWriter().write(outputs.core, staging_dir)
            if outputs.kinase_activity is not None:
                KinaseActivityWriter.write(outputs.kinase_activity, staging_dir)
            manifest_writer.write(
                outdir=staging_dir,
                core=outputs.core,
                kinase_activity=outputs.kinase_activity,
                preprocessing_config=preprocessing_config,
            )
            output_publisher.publish(
                staging_dir=staging_dir,
                target_dir=target_dir,
            )


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
        self.execution_runner = PipelineExecutionRunner(
            kinase_activity_analyzer=self.kinase_activity_analyzer
        )
        self.output_coordinator = PipelineOutputCoordinator()

    @classmethod
    def from_request(cls, request: CorePipelineRequest) -> PhosRPipeline:
        inputs = PipelineRequestLoader().load(request)
        return cls(
            dataset=inputs.dataset,
            pred_mat=inputs.pred_mat,
            preprocessing_config=inputs.preprocessing_config,
        )

    @classmethod
    def from_files(
        cls,
        total_path: str | Path,
        phospho_path: str | Path,
        pred_mat_path: str | Path | None = None,
        comparisons: Sequence[ComparisonSpec] | None = None,
        phospho_encoding: str | None = None,
        schema: DatasetSchema | None = None,
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
            schema=schema or DatasetSchema(),
            comparisons=tuple(comparisons) if comparisons is not None else None,
            localization_threshold=localization_threshold,
            min_observed=min_observed,
            total_sentinel=total_sentinel,
            phospho_sentinel=phospho_sentinel,
            max_unmatched_fraction=max_unmatched_fraction,
        )
        return cls.from_request(request)

    def run(self, outdir: str | Path | None = None) -> CoreOutputs:
        outputs = self.execution_runner.run(
            dataset=self.dataset,
            pred_mat=self.pred_mat,
            preprocessing_config=self.preprocessing_config,
        )

        if outdir is not None:
            self.output_coordinator.publish(
                outdir=outdir,
                outputs=outputs,
                preprocessing_config=self.preprocessing_config,
                manifest_writer=self.manifest_writer,
                output_publisher=self.output_publisher,
            )

        return outputs
