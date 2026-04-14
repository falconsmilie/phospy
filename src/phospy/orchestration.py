from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd

from .activities.analysis import KinaseActivityAnalyzer
from .activities.results import KinaseActivityResult
from .preprocessing.core import CorePreprocessingConfig, CoreProcessingResult
from .validation.requests.pipeline import (
    PipelineInputs,
    validate_pipeline_runtime_compatibility,
)

if TYPE_CHECKING:
    from .api.contracts import (
        DatasetLoadOptions,
        KinaseActivityConfig,
        PredictionRunConfig,
    )
    from .api.workflow_results import SimpleKinaseWorkflowResult
    from .internal.kinase_workflows import PredMatWorkflow
    from .preprocessing.modes import AnalysisReadyDatasetBuilder
    from .references import ReferenceProvider


class KinaseOrchestrationService:
    """Shared orchestration lane for high-level kinase entrypoints."""

    def run_simple_workflow(
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
        analysis_ready_builder: AnalysisReadyDatasetBuilder,
        reference_provider: ReferenceProvider,
        pred_mat_workflow: PredMatWorkflow,
        activity_analyzer: KinaseActivityAnalyzer,
    ) -> SimpleKinaseWorkflowResult:
        from .api.contracts import (
            DatasetLoadOptions,
            KinaseActivityConfig,
            PredictionRunConfig,
        )
        from .api.workflow_results import SimpleKinaseWorkflowResult

        resolved_dataset_options = DatasetLoadOptions.from_value(dataset_options)
        resolved_preprocessing_config = (
            CorePreprocessingConfig()
            if preprocessing_config is None
            else preprocessing_config
        )
        resolved_prediction_config = PredictionRunConfig.from_value(prediction_config)
        resolved_activity_config = KinaseActivityConfig.from_value(activity_config)

        analysis_ready_dataset = analysis_ready_builder.build(
            phospho=phospho,
            total=total,
            phospho_encoding=resolved_dataset_options.phospho_encoding,
            schema=resolved_dataset_options.schema,
            comparisons=resolved_dataset_options.comparisons,
            preprocessing_config=resolved_preprocessing_config,
            source="simple kinase workflow",
            phospho_only_source="simple kinase workflow (phospho only)",
        )
        reference_bundle = reference_provider.resolve(
            species=species,
            reference=reference,
        )
        workflow_result = pred_mat_workflow.run(
            phospho_matrix=analysis_ready_dataset.phospho_matrix,
            site_sequences=analysis_ready_dataset.site_sequences,
            reference_bundle=reference_bundle,
            prediction_config=resolved_prediction_config,
        )
        kinase_activity_result = activity_analyzer.run(
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

    def run_pipeline_runtime(
        self,
        *,
        request: PipelineInputs,
        kinase_activity_analyzer: KinaseActivityAnalyzer,
    ) -> tuple[CoreProcessingResult, KinaseActivityResult | None]:
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
        return core, kinase_activity
