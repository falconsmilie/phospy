from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Protocol

import pandas as pd

from ..activities.analysis import KinaseActivityAnalyzer
from ..internal.defaults import DEFAULT_MOTIF_FLANK_SIZE
from ..internal.types import PREDICTION_SVM_MODE_DEFAULT, PredictionSvmMode
from ..prediction.engines import (
    KinaseWorkflowExecutionResult,
    KinaseWorkflowExecutor,
)
from ..preprocessing.core import CorePreprocessingConfig
from ..preprocessing.modes import AnalysisReadyDatasetBuilder
from ..references import ReferenceBundle, ReferenceProvider
from ..validation.requests.workflow import WorkflowInputs
from .contracts import (
    DatasetLoadOptions,
    KinaseActivityConfig,
    PredictionRunConfig,
)
from .simple_workflow_composition import (
    SimpleKinaseExecutionGraph,
    create_default_simple_kinase_execution_graph,
)
from .workflow_results import SimpleKinaseWorkflowResult

__all__ = ["SimpleKinaseWorkflow"]


@dataclass(frozen=True, slots=True)
class _ResolvedSimpleKinaseConfigs:
    dataset_options: DatasetLoadOptions
    preprocessing_config: CorePreprocessingConfig
    prediction_config: PredictionRunConfig
    activity_config: KinaseActivityConfig


def _resolve_simple_kinase_configs(
    *,
    dataset_options: DatasetLoadOptions | None,
    preprocessing_config: CorePreprocessingConfig | None,
    prediction_config: PredictionRunConfig | None,
    activity_config: KinaseActivityConfig | None,
) -> _ResolvedSimpleKinaseConfigs:
    return _ResolvedSimpleKinaseConfigs(
        dataset_options=DatasetLoadOptions.from_value(dataset_options),
        preprocessing_config=(
            CorePreprocessingConfig()
            if preprocessing_config is None
            else preprocessing_config
        ),
        prediction_config=PredictionRunConfig.from_value(prediction_config),
        activity_config=KinaseActivityConfig.from_value(activity_config),
    )


def _validate_prediction_request(
    *,
    workflow_executor: KinaseWorkflowExecutor,
    phospho_matrix: pd.DataFrame,
    site_sequences: Mapping[str, str] | pd.Series,
    reference_bundle: ReferenceBundle,
    prediction_config: PredictionRunConfig,
) -> WorkflowInputs:
    return workflow_executor.validate_request(
        phospho_matrix=phospho_matrix,
        site_sequences=site_sequences,
        reference_bundle=reference_bundle,
        min_substrates=prediction_config.min_substrates,
        min_motif_size=prediction_config.min_motif_size,
        allow_profile_only_fallback=prediction_config.allow_profile_only_fallback,
        ensemble_size=prediction_config.ensemble_size,
        top=prediction_config.top,
        score_threshold=prediction_config.score_threshold,
        inclusion=prediction_config.inclusion,
        n_iterations=prediction_config.n_iterations,
        random_state=prediction_config.random_state,
        svm_mode=prediction_config.svm_mode,
        profile_policy=prediction_config.profile_policy,
    )


class _SimpleKinaseExecutionService(Protocol):
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
    ) -> SimpleKinaseWorkflowResult: ...


class SimpleKinaseWorkflow:
    """Run the common end-to-end kinase inference lane from user-shaped inputs."""

    def __init__(
        self,
        flank_size: int = DEFAULT_MOTIF_FLANK_SIZE,
        kernel: str = "rbf",
        svm_mode: PredictionSvmMode = PREDICTION_SVM_MODE_DEFAULT,
        *,
        execution_service: _SimpleKinaseExecutionService | None = None,
        execution_graph: SimpleKinaseExecutionGraph | None = None,
        analysis_ready_builder: AnalysisReadyDatasetBuilder | None = None,
        reference_provider: ReferenceProvider | None = None,
        activity_analyzer: KinaseActivityAnalyzer | None = None,
        workflow_executor: KinaseWorkflowExecutor | None = None,
    ) -> None:
        has_overrides = (
            analysis_ready_builder is not None
            or reference_provider is not None
            or activity_analyzer is not None
            or workflow_executor is not None
        )
        if execution_service is not None:
            if execution_graph is not None or has_overrides:
                msg = (
                    "execution_service cannot be combined with execution_graph "
                    "or collaborator overrides."
                )
                raise ValueError(msg)
            self._execution_service = execution_service
            self._analysis_ready_builder = None
            self._reference_provider = None
            self._activity_analyzer = None
            self._workflow_executor = None
            return

        if execution_graph is not None and has_overrides:
            msg = "execution_graph cannot be combined with collaborator overrides."
            raise ValueError(msg)
        if execution_graph is None:
            execution_graph = create_default_simple_kinase_execution_graph(
                flank_size=flank_size,
                kernel=kernel,
                svm_mode=svm_mode,
                analysis_ready_builder=analysis_ready_builder,
                reference_provider=reference_provider,
                activity_analyzer=activity_analyzer,
                workflow_executor=workflow_executor,
            )

        self._execution_service = None
        self._analysis_ready_builder = execution_graph.analysis_ready_builder
        self._reference_provider = execution_graph.reference_provider
        self._activity_analyzer = execution_graph.activity_analyzer
        self._workflow_executor = execution_graph.workflow_executor

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
        if self._execution_service is not None:
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

        resolved_configs = _resolve_simple_kinase_configs(
            dataset_options=dataset_options,
            preprocessing_config=preprocessing_config,
            prediction_config=prediction_config,
            activity_config=activity_config,
        )
        assert self._analysis_ready_builder is not None
        assert self._reference_provider is not None
        assert self._activity_analyzer is not None
        assert self._workflow_executor is not None

        analysis_ready_dataset = self._analysis_ready_builder.build(
            phospho=phospho,
            total=total,
            phospho_encoding=resolved_configs.dataset_options.phospho_encoding,
            schema=resolved_configs.dataset_options.schema,
            comparisons=resolved_configs.dataset_options.comparisons,
            preprocessing_config=resolved_configs.preprocessing_config,
            source="simple kinase workflow",
            phospho_only_source="simple kinase workflow (phospho only)",
        )
        reference_bundle = self._reference_provider.resolve(
            species=species,
            reference=reference,
        )
        request = _validate_prediction_request(
            workflow_executor=self._workflow_executor,
            phospho_matrix=analysis_ready_dataset.phospho_matrix,
            site_sequences=analysis_ready_dataset.site_sequences,
            reference_bundle=reference_bundle,
            prediction_config=resolved_configs.prediction_config,
        )
        workflow_result: KinaseWorkflowExecutionResult = (
            self._workflow_executor.execute_validated_request(request)
        )
        pred_mat_result = workflow_result.prediction_result.pred_mat_result
        kinase_activity_result = self._activity_analyzer.run(
            pred_mat=pred_mat_result,
            phospho_matrix=analysis_ready_dataset.phospho_matrix,
            threshold=resolved_configs.activity_config.threshold,
            min_substrates=resolved_configs.activity_config.min_substrates,
            top_n_substrates=resolved_configs.activity_config.top_n_substrates,
        )
        return SimpleKinaseWorkflowResult(
            analysis_ready_dataset=analysis_ready_dataset,
            reference_bundle=reference_bundle,
            scoring_result=workflow_result.scoring_result,
            prediction_result=workflow_result.prediction_result,
            pred_mat_result=pred_mat_result,
            kinase_activity_result=kinase_activity_result,
        )
