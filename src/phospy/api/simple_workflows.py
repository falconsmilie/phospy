from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import pandas as pd

from ..activities.analysis import KinaseActivityAnalyzer
from ..datasets.schema import DatasetSchema
from ..internal.constants import (
    DEFAULT_PHOSPHO_SENTINEL,
    DEFAULT_TOTAL_SENTINEL,
    ComparisonSpec,
)
from ..internal.types import PredictionSvmMode
from ..preprocessing.core import CorePreprocessingConfig
from ..preprocessing.modes import AnalysisReadyDatasetBuilder
from ..profiles import KinaseProfilePolicy
from ..references import (
    BundledReferenceProvider,
    ReferenceProvider,
)
from .kinase_workflows import PredMatWorkflow
from .workflow_results import SimpleKinaseWorkflowResult

__all__ = ["SimpleKinaseWorkflow"]


class SimpleKinaseWorkflow:
    """Run the common end-to-end kinase inference lane from user-shaped inputs."""

    def __init__(
        self,
        flank_size: int = 7,
        kernel: str = "rbf",
        svm_mode: PredictionSvmMode = "default",
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

    def run(
        self,
        *,
        phospho: pd.DataFrame | str | Path,
        species: str,
        total: pd.DataFrame | str | Path | None = None,
        reference: str = "auto",
        phospho_encoding: str | None = None,
        comparisons: Sequence[ComparisonSpec] | None = None,
        schema: DatasetSchema | None = None,
        preprocessing_config: CorePreprocessingConfig | None = None,
        localization_threshold: float = 0.75,
        min_observed: int = 4,
        max_unmatched_fraction: float = 0.0,
        total_sentinel: float | int = DEFAULT_TOTAL_SENTINEL,
        phospho_sentinel: float | int = DEFAULT_PHOSPHO_SENTINEL,
        min_substrates: int = 1,
        min_motif_size: int = 1,
        allow_profile_only_fallback: bool = False,
        ensemble_size: int = 10,
        top: int = 50,
        score_threshold: float = 0.8,
        inclusion: int = 20,
        n_iterations: int = 5,
        random_state: int | None = None,
        svm_mode: PredictionSvmMode | None = None,
        profile_policy: KinaseProfilePolicy | None = None,
        kinase_activity_threshold: float = 0.6,
        kinase_activity_min_substrates: int = 3,
        kinase_activity_top_n_substrates: int = 20,
    ) -> SimpleKinaseWorkflowResult:
        resolved_schema = schema or DatasetSchema()
        analysis_ready_dataset = self.analysis_ready_builder.build(
            phospho=phospho,
            total=total,
            phospho_encoding=phospho_encoding,
            schema=resolved_schema,
            comparisons=comparisons,
            preprocessing_config=preprocessing_config,
            localization_threshold=localization_threshold,
            min_observed=min_observed,
            max_unmatched_fraction=max_unmatched_fraction,
            total_sentinel=total_sentinel,
            phospho_sentinel=phospho_sentinel,
            source="simple kinase workflow",
            phospho_only_source="simple kinase workflow (phospho only)",
        )
        reference_bundle = self.reference_provider.resolve(
            species=species,
            reference=reference,
        )
        workflow_result = self.pred_mat_workflow.run(
            phospho_matrix=analysis_ready_dataset.phospho_matrix,
            site_sequences=analysis_ready_dataset.site_sequences,
            reference_bundle=reference_bundle,
            min_substrates=min_substrates,
            min_motif_size=min_motif_size,
            allow_profile_only_fallback=allow_profile_only_fallback,
            ensemble_size=ensemble_size,
            top=top,
            score_threshold=score_threshold,
            inclusion=inclusion,
            n_iterations=n_iterations,
            random_state=random_state,
            svm_mode=svm_mode,
            profile_policy=profile_policy,
        )
        kinase_activity_result = self.activity_analyzer.run(
            pred_mat=workflow_result.pred_mat_result,
            phospho_matrix=analysis_ready_dataset.phospho_matrix,
            threshold=kinase_activity_threshold,
            min_substrates=kinase_activity_min_substrates,
            top_n_substrates=kinase_activity_top_n_substrates,
        )
        return SimpleKinaseWorkflowResult(
            analysis_ready_dataset=analysis_ready_dataset,
            reference_bundle=reference_bundle,
            workflow_result=workflow_result,
            kinase_activity_result=kinase_activity_result,
        )
