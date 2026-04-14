from __future__ import annotations

from pathlib import Path

import pandas as pd

from ..activities.analysis import KinaseActivityAnalyzer
from ..internal.defaults import DEFAULT_MOTIF_FLANK_SIZE
from ..internal.kinase_workflows import PredMatWorkflow
from ..internal.types import PREDICTION_SVM_MODE_DEFAULT, PredictionSvmMode
from ..orchestration import KinaseOrchestrationService
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


class SimpleKinaseWorkflow:
    """Run the common end-to-end kinase inference lane from user-shaped inputs."""

    def __init__(
        self,
        flank_size: int = DEFAULT_MOTIF_FLANK_SIZE,
        kernel: str = "rbf",
        svm_mode: PredictionSvmMode = PREDICTION_SVM_MODE_DEFAULT,
        *,
        reference_provider: ReferenceProvider | None = None,
        activity_analyzer: KinaseActivityAnalyzer | None = None,
        analysis_ready_builder: AnalysisReadyDatasetBuilder | None = None,
    ) -> None:
        self.pred_mat_workflow = PredMatWorkflow(
            flank_size=flank_size,
            kernel=kernel,
            svm_mode=svm_mode,
        )
        self.reference_provider = (
            BundledReferenceProvider()
            if reference_provider is None
            else reference_provider
        )
        self.activity_analyzer = (
            KinaseActivityAnalyzer() if activity_analyzer is None else activity_analyzer
        )
        self.analysis_ready_builder = (
            AnalysisReadyDatasetBuilder()
            if analysis_ready_builder is None
            else analysis_ready_builder
        )
        self._orchestration = KinaseOrchestrationService()

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
        return self._orchestration.run_simple_workflow(
            phospho=phospho,
            species=species,
            total=total,
            reference=reference,
            prediction_config=prediction_config,
            dataset_options=dataset_options,
            preprocessing_config=preprocessing_config,
            activity_config=activity_config,
            analysis_ready_builder=self.analysis_ready_builder,
            reference_provider=self.reference_provider,
            pred_mat_workflow=self.pred_mat_workflow,
            activity_analyzer=self.activity_analyzer,
        )
