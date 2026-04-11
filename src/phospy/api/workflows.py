from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from ..activities import KinaseActivityAnalyzer
from ..activities.results import KinaseActivityResult
from ..constants import (
    DEFAULT_PHOSPHO_SENTINEL,
    DEFAULT_TOTAL_SENTINEL,
    ComparisonSpec,
)
from ..datasets import AnalysisReadyPhosphoDataset, DatasetSchema
from ..motifs import MotifScoringResult
from ..prediction import (
    KinasePredictionResult,
    KinaseScoringResult,
    KinaseWorkflowExecutor,
    PredMatResult,
)
from ..preprocessing.core import CorePreprocessingConfig
from ..preprocessing.modes import AnalysisReadyDatasetBuilder
from ..profiles import KinaseProfileResult
from ..references import (
    BundledReferenceProvider,
    ReferenceBundle,
    ReferenceProvider,
)
from ..signalomes import SignalomeResult, execute_validated_signalome_request
from ..types import PredictionSvmMode
from ..validation.requests import (
    ValidatedSignalomeRequest,
    ValidatedWorkflowRequest,
    validate_signalome_request,
)

__all__ = [
    "KinaseWorkflow",
    "KinaseWorkflowResult",
    "PredMatWorkflow",
    "PredMatWorkflowResult",
    "SignalomeWorkflow",
    "SimpleKinaseWorkflow",
    "SimpleKinaseWorkflowResult",
]


@dataclass(slots=True)
class KinaseWorkflowResult:
    """Workflow outputs for a single native scoring and prediction run."""

    profile_result: KinaseProfileResult
    motif_result: MotifScoringResult | None
    scoring_result: KinaseScoringResult
    prediction_result: KinasePredictionResult


@dataclass(slots=True)
class PredMatWorkflowResult:
    """Stable result bundle for one public predMat generation run.

    The recommended predMat contract is exposed through ``pred_mat_result``.
    """

    scoring_result: KinaseScoringResult
    prediction_result: KinasePredictionResult
    pred_mat_result: PredMatResult

    def close(self) -> None:
        """Release owned trace resources, if any are attached downstream."""

        self.prediction_result.close()

    def __enter__(self) -> PredMatWorkflowResult:
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc: BaseException | None,
        tb: object | None,
    ) -> None:
        self.close()


@dataclass(slots=True)
class SimpleKinaseWorkflowResult:
    """Owned result bundle for the high-level common kinase workflow."""

    analysis_ready_dataset: AnalysisReadyPhosphoDataset
    reference_bundle: ReferenceBundle
    workflow_result: PredMatWorkflowResult
    kinase_activity_result: KinaseActivityResult

    @property
    def pred_mat_result(self) -> PredMatResult:
        return self.workflow_result.pred_mat_result

    @property
    def scoring_result(self) -> KinaseScoringResult:
        return self.workflow_result.scoring_result

    @property
    def prediction_result(self) -> KinasePredictionResult:
        return self.workflow_result.prediction_result

    def close(self) -> None:
        self.workflow_result.close()

    def __enter__(self) -> SimpleKinaseWorkflowResult:
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc: BaseException | None,
        tb: object | None,
    ) -> None:
        self.close()


class KinaseWorkflow:
    """Run the native kinase scoring and prediction workflow end to end."""

    def __init__(
        self,
        flank_size: int = 7,
        kernel: str = "rbf",
        svm_mode: PredictionSvmMode = "default",
    ) -> None:
        self._executor = KinaseWorkflowExecutor(
            flank_size=flank_size,
            kernel=kernel,
            svm_mode=svm_mode,
        )

    def _validate_request(
        self,
        *,
        phospho_matrix: pd.DataFrame,
        substrate_map: Mapping[str, Sequence[str]] | None = None,
        site_sequences: Mapping[str, str] | pd.Series | None = None,
        motif_sequences: Mapping[str, Sequence[str]] | None = None,
        reference_bundle: ReferenceBundle | None = None,
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
    ) -> ValidatedWorkflowRequest:
        return self._executor.validate_request(
            phospho_matrix=phospho_matrix,
            substrate_map=substrate_map,
            site_sequences=site_sequences,
            motif_sequences=motif_sequences,
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
        )

    def run(
        self,
        phospho_matrix: pd.DataFrame,
        substrate_map: Mapping[str, Sequence[str]] | None = None,
        site_sequences: Mapping[str, str] | pd.Series | None = None,
        motif_sequences: Mapping[str, Sequence[str]] | None = None,
        reference_bundle: ReferenceBundle | None = None,
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
    ) -> KinaseWorkflowResult:
        request = self._validate_request(
            phospho_matrix=phospho_matrix,
            substrate_map=substrate_map,
            site_sequences=site_sequences,
            motif_sequences=motif_sequences,
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
        )
        return self.run_validated(request)

    def run_validated(
        self,
        request: ValidatedWorkflowRequest,
    ) -> KinaseWorkflowResult:
        result = self._executor.execute_validated_request(request)
        return KinaseWorkflowResult(
            profile_result=result.profile_result,
            motif_result=result.motif_result,
            scoring_result=result.scoring_result,
            prediction_result=result.prediction_result,
        )


class PredMatWorkflow:
    """Generate a predMat from phosphosite and sequence inputs."""

    def __init__(
        self,
        flank_size: int = 7,
        kernel: str = "rbf",
        svm_mode: PredictionSvmMode = "default",
    ) -> None:
        self._executor = KinaseWorkflowExecutor(
            flank_size=flank_size,
            kernel=kernel,
            svm_mode=svm_mode,
        )

    def _validate_request(
        self,
        *,
        phospho_matrix: pd.DataFrame,
        substrate_map: Mapping[str, Sequence[str]] | None = None,
        site_sequences: Mapping[str, str] | pd.Series | None = None,
        motif_sequences: Mapping[str, Sequence[str]] | None = None,
        reference_bundle: ReferenceBundle | None = None,
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
    ) -> ValidatedWorkflowRequest:
        return self._executor.validate_request(
            phospho_matrix=phospho_matrix,
            substrate_map=substrate_map,
            site_sequences=site_sequences,
            motif_sequences=motif_sequences,
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
        )

    def run(
        self,
        phospho_matrix: pd.DataFrame,
        substrate_map: Mapping[str, Sequence[str]] | None = None,
        site_sequences: Mapping[str, str] | pd.Series | None = None,
        motif_sequences: Mapping[str, Sequence[str]] | None = None,
        reference_bundle: ReferenceBundle | None = None,
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
    ) -> PredMatWorkflowResult:
        request = self._validate_request(
            phospho_matrix=phospho_matrix,
            substrate_map=substrate_map,
            site_sequences=site_sequences,
            motif_sequences=motif_sequences,
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
        )
        return self.run_validated(request)

    def run_validated(
        self,
        request: ValidatedWorkflowRequest,
    ) -> PredMatWorkflowResult:
        result = self._executor.execute_validated_request(request)
        return PredMatWorkflowResult(
            scoring_result=result.scoring_result,
            prediction_result=result.prediction_result,
            pred_mat_result=result.prediction_result.pred_mat_result,
        )


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


class SignalomeWorkflow:
    """Construct signalomes from validated scoring and prediction outputs."""

    def _validate_request(
        self,
        *,
        scoring_result: KinaseScoringResult,
        prediction_result: KinasePredictionResult | PredMatResult,
        expression_matrix: pd.DataFrame,
        kinases_of_interest: Sequence[str],
        site_to_protein: Mapping[str, str] | None = None,
        kinase_network_threshold: float = 0.9,
        signalome_cutoff: float = 0.5,
        module_count: int | None = None,
        min_kinase_module_share_percent: float = 1.0,
    ) -> ValidatedSignalomeRequest:
        return validate_signalome_request(
            scoring_result=scoring_result,
            prediction_result=prediction_result,
            expression_matrix=expression_matrix,
            kinases_of_interest=kinases_of_interest,
            site_to_protein=site_to_protein,
            kinase_network_threshold=kinase_network_threshold,
            signalome_cutoff=signalome_cutoff,
            module_count=module_count,
            min_kinase_module_share_percent=min_kinase_module_share_percent,
        )

    def run(
        self,
        *,
        scoring_result: KinaseScoringResult,
        prediction_result: KinasePredictionResult | PredMatResult,
        expression_matrix: pd.DataFrame,
        kinases_of_interest: Sequence[str],
        site_to_protein: Mapping[str, str] | None = None,
        kinase_network_threshold: float = 0.9,
        signalome_cutoff: float = 0.5,
        module_count: int | None = None,
        min_kinase_module_share_percent: float = 1.0,
    ) -> SignalomeResult:
        request = self._validate_request(
            scoring_result=scoring_result,
            prediction_result=prediction_result,
            expression_matrix=expression_matrix,
            kinases_of_interest=kinases_of_interest,
            site_to_protein=site_to_protein,
            kinase_network_threshold=kinase_network_threshold,
            signalome_cutoff=signalome_cutoff,
            module_count=module_count,
            min_kinase_module_share_percent=min_kinase_module_share_percent,
        )
        return self.run_validated(request)

    def run_validated(
        self,
        request: ValidatedSignalomeRequest,
    ) -> SignalomeResult:
        return execute_validated_signalome_request(request)
