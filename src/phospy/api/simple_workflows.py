from __future__ import annotations

from pathlib import Path

import pandas as pd

from ..activities.analysis import KinaseActivityAnalyzer
from ..internal.defaults import DEFAULT_MOTIF_FLANK_SIZE
from ..internal.kinase_workflows import PredMatWorkflow
from ..internal.types import PREDICTION_SVM_MODE_DEFAULT, PredictionSvmMode
from ..preprocessing.core import CorePreprocessingConfig
from ..preprocessing.modes import AnalysisReadyDatasetBuilder
from ..references import (
    BundledReferenceProvider,
    ReferenceProvider,
)
from .contracts import (
    DatasetLoadOptions,
    KinaseActivityConfig,
    PredictionRunConfig,
)
from .workflow_results import SimpleKinaseWorkflowResult

__all__ = ["SimpleKinaseWorkflow"]


class SimpleKinaseExecutionService:
    """Validated execution path for the public simple kinase workflow."""

    def __init__(
        self,
        *,
        analysis_ready_builder: AnalysisReadyDatasetBuilder,
        reference_provider: ReferenceProvider,
        pred_mat_workflow: PredMatWorkflow,
        activity_analyzer: KinaseActivityAnalyzer,
    ) -> None:
        self._analysis_ready_builder = analysis_ready_builder
        self._reference_provider = reference_provider
        self._pred_mat_workflow = pred_mat_workflow
        self._activity_analyzer = activity_analyzer

    @classmethod
    def create_default(
        cls,
        *,
        flank_size: int = DEFAULT_MOTIF_FLANK_SIZE,
        kernel: str = "rbf",
        svm_mode: PredictionSvmMode = PREDICTION_SVM_MODE_DEFAULT,
        reference_provider: ReferenceProvider | None = None,
        activity_analyzer: KinaseActivityAnalyzer | None = None,
        analysis_ready_builder: AnalysisReadyDatasetBuilder | None = None,
    ) -> SimpleKinaseExecutionService:
        """Documented default factory path for runtime collaborator assembly."""
        return cls(
            analysis_ready_builder=(
                AnalysisReadyDatasetBuilder()
                if analysis_ready_builder is None
                else analysis_ready_builder
            ),
            reference_provider=(
                BundledReferenceProvider()
                if reference_provider is None
                else reference_provider
            ),
            pred_mat_workflow=PredMatWorkflow(
                flank_size=flank_size,
                kernel=kernel,
                svm_mode=svm_mode,
            ),
            activity_analyzer=(
                KinaseActivityAnalyzer()
                if activity_analyzer is None
                else activity_analyzer
            ),
        )

    def run(
        self,
        *,
        phospho: pd.DataFrame | str | Path,
        species: str,
        total: pd.DataFrame | str | Path | None,
        reference: str,
        dataset_options: DatasetLoadOptions | None,
        preprocessing_config: CorePreprocessingConfig | None,
        prediction_config: PredictionRunConfig | None,
        activity_config: KinaseActivityConfig | None,
    ) -> SimpleKinaseWorkflowResult:
        resolved_dataset_options = DatasetLoadOptions.from_value(dataset_options)
        resolved_preprocessing_config = (
            CorePreprocessingConfig()
            if preprocessing_config is None
            else preprocessing_config
        )
        resolved_prediction_config = PredictionRunConfig.from_value(prediction_config)
        resolved_activity_config = KinaseActivityConfig.from_value(activity_config)

        analysis_ready_dataset = self._analysis_ready_builder.build(
            phospho=phospho,
            total=total,
            phospho_encoding=resolved_dataset_options.phospho_encoding,
            schema=resolved_dataset_options.schema,
            comparisons=resolved_dataset_options.comparisons,
            preprocessing_config=resolved_preprocessing_config,
            source="simple kinase workflow",
            phospho_only_source="simple kinase workflow (phospho only)",
        )
        reference_bundle = self._reference_provider.resolve(
            species=species,
            reference=reference,
        )
        workflow_result = self._pred_mat_workflow.run(
            phospho_matrix=analysis_ready_dataset.phospho_matrix,
            site_sequences=analysis_ready_dataset.site_sequences,
            reference_bundle=reference_bundle,
            prediction_config=resolved_prediction_config,
        )
        kinase_activity_result = self._activity_analyzer.run(
            pred_mat=workflow_result.pred_mat_result,
            phospho_matrix=analysis_ready_dataset.phospho_matrix,
            threshold=resolved_activity_config.threshold,
            min_substrates=resolved_activity_config.min_substrates,
            top_n_substrates=resolved_activity_config.top_n_substrates,
        )
        return SimpleKinaseWorkflowResult(
            analysis_ready_dataset=analysis_ready_dataset,
            reference_bundle=reference_bundle,
            scoring_result=workflow_result.scoring_result,
            prediction_result=workflow_result.prediction_result,
            pred_mat_result=workflow_result.pred_mat_result,
            kinase_activity_result=kinase_activity_result,
        )


class SimpleKinaseWorkflow:
    """Run the common end-to-end kinase inference lane from user-shaped inputs."""

    def __init__(
        self,
        flank_size: int = DEFAULT_MOTIF_FLANK_SIZE,
        kernel: str = "rbf",
        svm_mode: PredictionSvmMode = PREDICTION_SVM_MODE_DEFAULT,
        *,
        execution_service: SimpleKinaseExecutionService | None = None,
    ) -> None:
        self._execution_service = (
            SimpleKinaseExecutionService.create_default(
                flank_size=flank_size,
                kernel=kernel,
                svm_mode=svm_mode,
            )
            if execution_service is None
            else execution_service
        )

    def run(
        self,
        *,
        phospho: pd.DataFrame | str | Path,
        species: str,
        total: pd.DataFrame | str | Path | None = None,
        reference: str = "auto",
        dataset_options: DatasetLoadOptions | None = None,
        preprocessing_config: CorePreprocessingConfig | None = None,
        prediction_config: PredictionRunConfig | None = None,
        activity_config: KinaseActivityConfig | None = None,
    ) -> SimpleKinaseWorkflowResult:
        return self._execution_service.run(
            phospho=phospho,
            species=species,
            total=total,
            reference=reference,
            prediction_config=prediction_config,
            dataset_options=dataset_options,
            preprocessing_config=preprocessing_config,
            activity_config=activity_config,
        )
