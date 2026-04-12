from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import TYPE_CHECKING

import pandas as pd

from .activities.analysis import KinaseActivityAnalyzer
from .activities.results import KinaseActivityResult
from .datasets.loaders import DatasetLoader
from .datasets.models import PhosphoDataset
from .datasets.schema import DatasetSchema
from .errors import RequestValidationError, TableSchemaError
from .internal.constants import ComparisonSpec
from .io import load_pred_mat
from .io.publishing import OutputPublisher, RunManifestWriter
from .io.writers import CoreOutputWriter, KinaseActivityWriter
from .preprocessing.core import CorePreprocessingConfig, CoreProcessingResult
from .validation.requests.pipeline import (
    CorePipelineRequest,
    PipelineInputs,
    build_pipeline_inputs,
    validate_pipeline_construction_request,
    validate_pipeline_runtime_compatibility,
)

__all__ = ["PhosRPipeline"]


if TYPE_CHECKING:
    from .prediction.results import PredMatResult


@dataclass(slots=True)
class CoreOutputs:
    """Outputs returned by ``PhosRPipeline.run()`` for one pipeline run."""

    core: CoreProcessingResult
    kinase_activity: KinaseActivityResult | None = None


class _PipelineRequestLoader:
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

    def load(self, request: CorePipelineRequest) -> PipelineInputs:
        validated_inputs = self.dataset_loader_factory(
            schema=request.dataset_schema
        ).load(
            request.total_path,
            request.phospho_path,
            phospho_encoding=request.phospho_encoding,
        )
        dataset = PhosphoDataset.from_loaded_inputs(
            validated_inputs,
            comparisons=request.comparisons,
        )
        return build_pipeline_inputs(
            dataset=dataset,
            pred_mat=self._load_pred_mat(request),
            localization_threshold=request.localization_threshold,
            min_observed=request.min_observed,
            max_unmatched_fraction=request.max_unmatched_fraction,
            total_sentinel=request.total_sentinel,
            phospho_sentinel=request.phospho_sentinel,
            kinase_activity_threshold=request.kinase_activity_threshold,
            kinase_activity_min_substrates=request.kinase_activity_min_substrates,
            kinase_activity_top_n_substrates=request.kinase_activity_top_n_substrates,
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


def _run_pipeline_request(
    *,
    request: PipelineInputs,
    kinase_activity_analyzer: KinaseActivityAnalyzer,
) -> CoreOutputs:
    core = request.dataset.preprocessing.run(config=request.preprocessing_config)

    kinase_activity = None
    kinase_activity_request = validate_pipeline_runtime_compatibility(
        request=request,
        site_matrix=core.site_matrix.matrix,
    )
    if kinase_activity_request is not None:
        kinase_activity = kinase_activity_analyzer.run_validated(
            kinase_activity_request
        )

    return CoreOutputs(core=core, kinase_activity=kinase_activity)


def _publish_pipeline_outputs(
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
            KinaseActivityWriter().write(outputs.kinase_activity, staging_dir)
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
        pred_mat: pd.DataFrame | PredMatResult | None = None,
        preprocessing_config: CorePreprocessingConfig | None = None,
        localization_threshold: float = 0.75,
        min_observed: int = 4,
        max_unmatched_fraction: float = 0.0,
        total_sentinel: float = 10.0,
        phospho_sentinel: float = 12.0,
        kinase_activity_threshold: float = 0.6,
        kinase_activity_min_substrates: int = 3,
        kinase_activity_top_n_substrates: int = 20,
        *,
        manifest_writer: RunManifestWriter | None = None,
        output_publisher: OutputPublisher | None = None,
    ) -> None:
        request = validate_pipeline_construction_request(
            dataset=dataset,
            pred_mat=pred_mat,
            preprocessing_config=preprocessing_config,
            localization_threshold=localization_threshold,
            min_observed=min_observed,
            max_unmatched_fraction=max_unmatched_fraction,
            total_sentinel=total_sentinel,
            phospho_sentinel=phospho_sentinel,
            kinase_activity_threshold=kinase_activity_threshold,
            kinase_activity_min_substrates=kinase_activity_min_substrates,
            kinase_activity_top_n_substrates=kinase_activity_top_n_substrates,
        )
        self._initialize_from_inputs(
            request=request,
            manifest_writer=manifest_writer,
            output_publisher=output_publisher,
        )

    def _initialize_from_inputs(
        self,
        *,
        request: PipelineInputs,
        manifest_writer: RunManifestWriter | None = None,
        output_publisher: OutputPublisher | None = None,
    ) -> None:
        object.__setattr__(self, "request", request)
        self.dataset = request.dataset
        self.pred_mat = request.pred_mat
        self.preprocessing_config = request.preprocessing_config
        self.kinase_activity_analyzer = KinaseActivityAnalyzer()
        self.manifest_writer = manifest_writer or RunManifestWriter()
        self.output_publisher = output_publisher or OutputPublisher()

    @classmethod
    def _from_request(
        cls,
        request: CorePipelineRequest,
        *,
        manifest_writer: RunManifestWriter | None = None,
        output_publisher: OutputPublisher | None = None,
    ) -> PhosRPipeline:
        inputs = _PipelineRequestLoader().load(request)
        return cls._from_inputs(
            inputs,
            manifest_writer=manifest_writer,
            output_publisher=output_publisher,
        )

    @classmethod
    def _from_inputs(
        cls,
        request: PipelineInputs,
        *,
        manifest_writer: RunManifestWriter | None = None,
        output_publisher: OutputPublisher | None = None,
    ) -> PhosRPipeline:
        instance = cls.__new__(cls)
        instance._initialize_from_inputs(
            request=request,
            manifest_writer=manifest_writer,
            output_publisher=output_publisher,
        )
        return instance

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
        kinase_activity_threshold: float = 0.6,
        kinase_activity_min_substrates: int = 3,
        kinase_activity_top_n_substrates: int = 20,
        *,
        manifest_writer: RunManifestWriter | None = None,
        output_publisher: OutputPublisher | None = None,
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
            kinase_activity_threshold=kinase_activity_threshold,
            kinase_activity_min_substrates=kinase_activity_min_substrates,
            kinase_activity_top_n_substrates=kinase_activity_top_n_substrates,
        )
        return cls._from_request(
            request,
            manifest_writer=manifest_writer,
            output_publisher=output_publisher,
        )

    def run(self, outdir: str | Path | None = None) -> CoreOutputs:
        outputs = _run_pipeline_request(
            request=self.request,
            kinase_activity_analyzer=self.kinase_activity_analyzer,
        )

        if outdir is not None:
            _publish_pipeline_outputs(
                outdir=outdir,
                outputs=outputs,
                preprocessing_config=self.preprocessing_config,
                manifest_writer=self.manifest_writer,
                output_publisher=self.output_publisher,
            )

        return outputs
