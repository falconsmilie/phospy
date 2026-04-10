from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from .analysis import KinaseActivityAnalyzer, KinaseActivityResult
from .constants import (
    DEFAULT_PHOSPHO_SENTINEL,
    DEFAULT_TOTAL_SENTINEL,
    TOTAL_GENE_COLUMN,
    ComparisonSpec,
)
from .core_processing import (
    CorePreprocessingConfig,
    CoreProcessingResult,
    resolve_core_preprocessing_config,
)
from .dataset import AnalysisReadyPhosphoDataset, PhosphoDataset
from .dataset_loader import DatasetLoader
from .dataset_schema import DatasetSchema
from .dataset_site_matrix import DatasetSiteMatrix
from .motifs import (
    BundledReferenceProvider,
    MotifScoringResult,
    ReferenceBundle,
    ReferenceProvider,
)
from .prediction import KinasePredictionResult, KinasePredictor, PredMatResult
from .preprocessing_services import PhosphoPreprocessor, ProteinCorrectionService
from .profiles import KinaseProfileResult, build_kinase_substrate_profiles
from .scoring import KinaseScorer, KinaseScoringResult
from .signalome_construction import execute_validated_signalome_request
from .signalomes import SignalomeResult
from .types import PredictionSvmMode
from .validation.requests import (
    ValidatedSignalomeRequest,
    ValidatedWorkflowRequest,
    validate_signalome_request,
    validate_workflow_request,
)

__all__ = [
    "KinaseWorkflow",
    "KinaseWorkflowResult",
    "PredMatWorkflow",
    "PredMatWorkflowResult",
    "SignalomeWorkflow",
    "SimpleKinaseWorkflow",
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


class _WorkflowBase:
    def __init__(
        self,
        flank_size: int = 7,
        kernel: str = "rbf",
        svm_mode: PredictionSvmMode = "default",
    ) -> None:
        self.flank_size = flank_size
        self.kernel = kernel
        self.svm_mode = svm_mode

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
        return validate_workflow_request(
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
            flank_size=self.flank_size,
            default_svm_mode=self.svm_mode,
        )

    def _execute_validated_request(
        self,
        request: ValidatedWorkflowRequest,
    ) -> KinaseWorkflowResult:
        raw_request = request.request
        phospho_matrix = request.phospho_matrix

        profile_result = build_kinase_substrate_profiles(
            substrate_map=raw_request.substrate_map,
            phospho_matrix=phospho_matrix,
            min_substrates=raw_request.min_substrates,
        )

        scorer = KinaseScorer(profile_result.profile_matrix)
        motif_result: MotifScoringResult | None = None
        scoring_matrix = phospho_matrix

        if request.motif_scorer is not None:
            scoring_matrix = phospho_matrix.loc[list(request.scoring_site_index)]
            motif_result = request.motif_scorer.score_sequences(
                seqs=raw_request.site_sequences,
                site_index=request.scoring_site_index,
                min_motif_size=raw_request.min_motif_size,
            )
            scoring_result = scorer.score(
                phospho_matrix=scoring_matrix,
                motif_scores=motif_result.motif_scores,
                motif_sizes=motif_result.motif_sizes,
                profile_sizes=profile_result.substrate_counts.astype(float),
                allow_profile_only_fallback=raw_request.allow_profile_only_fallback,
            )
        else:
            scoring_result = scorer.score(phospho_matrix=scoring_matrix)

        predictor = KinasePredictor(
            kernel=self.kernel,
            svm_mode=request.predictor_svm_mode,
        )
        prediction_result = predictor.predict_from_scoring_result(
            scoring_result=scoring_result,
            ensemble_size=raw_request.ensemble_size,
            top=raw_request.top,
            score_threshold=raw_request.score_threshold,
            inclusion=raw_request.inclusion,
            n_iterations=raw_request.n_iterations,
            random_state=raw_request.random_state,
            allow_profile_only_fallback=raw_request.allow_profile_only_fallback,
            svm_mode=raw_request.svm_mode,
        )

        return KinaseWorkflowResult(
            profile_result=profile_result,
            motif_result=motif_result,
            scoring_result=scoring_result,
            prediction_result=prediction_result,
        )


class KinaseWorkflow(_WorkflowBase):
    """Run the native kinase scoring and prediction workflow end to end."""

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
        return self._execute_validated_request(request)


class PredMatWorkflow(_WorkflowBase):
    """Generate a predMat from phosphosite and sequence inputs."""

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
        result = self._execute_validated_request(request)
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
        analysis_ready_dataset = self._build_analysis_ready_dataset(
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

    def _build_analysis_ready_dataset(
        self,
        *,
        phospho: pd.DataFrame | str | Path,
        total: pd.DataFrame | str | Path | None,
        phospho_encoding: str | None,
        schema: DatasetSchema,
        comparisons: Sequence[ComparisonSpec] | None,
        preprocessing_config: CorePreprocessingConfig | None,
        localization_threshold: float,
        min_observed: int,
        max_unmatched_fraction: float,
        total_sentinel: float | int,
        phospho_sentinel: float | int,
    ) -> AnalysisReadyPhosphoDataset:
        if total is None:
            return self._build_phospho_only_analysis_ready_dataset(
                phospho=phospho,
                phospho_encoding=phospho_encoding,
                schema=schema,
                comparisons=comparisons,
                preprocessing_config=preprocessing_config,
                localization_threshold=localization_threshold,
                min_observed=min_observed,
                phospho_sentinel=phospho_sentinel,
            )

        dataset_loader = DatasetLoader(schema=schema)
        phospho_df = self._resolve_phospho_input(
            phospho,
            dataset_loader=dataset_loader,
            phospho_encoding=phospho_encoding,
        )
        total_df = self._resolve_total_input(total, dataset_loader=dataset_loader)
        dataset = PhosphoDataset.from_loaded_inputs(
            dataset_loader.validate_inputs(total_df=total_df, phospho_df=phospho_df),
            comparisons=comparisons,
        )
        return dataset.run_analysis_ready(
            localization_threshold=localization_threshold,
            min_observed=min_observed,
            max_unmatched_fraction=max_unmatched_fraction,
            total_sentinel=total_sentinel,
            phospho_sentinel=phospho_sentinel,
            config=preprocessing_config,
            source="simple kinase workflow",
        )

    def _build_phospho_only_analysis_ready_dataset(
        self,
        *,
        phospho: pd.DataFrame | str | Path,
        phospho_encoding: str | None,
        schema: DatasetSchema,
        comparisons: Sequence[ComparisonSpec] | None,
        preprocessing_config: CorePreprocessingConfig | None,
        localization_threshold: float,
        min_observed: int,
        phospho_sentinel: float | int,
    ) -> AnalysisReadyPhosphoDataset:
        dataset_loader = DatasetLoader(schema=schema)
        phospho_df = self._resolve_phospho_input(
            phospho,
            dataset_loader=dataset_loader,
            phospho_encoding=phospho_encoding,
        )
        resolved_config = resolve_core_preprocessing_config(
            config=preprocessing_config,
            localization_threshold=localization_threshold,
            min_observed=min_observed,
            phospho_sentinel=phospho_sentinel,
            context="SimpleKinaseWorkflow.run()",
            config_param_name="preprocessing_config",
        )
        phospho_filtered = PhosphoPreprocessor(schema=schema).prepare(
            phospho_df,
            localization_threshold=resolved_config.localization_threshold,
            sentinel=resolved_config.phospho_sentinel,
            min_observed=resolved_config.min_observed,
        )
        phospho_corrected = phospho_filtered.rename(
            columns=dict(zip(schema.phospho_cols, schema.corrected_cols, strict=True))
        )
        phospho_corrected = ProteinCorrectionService(
            schema=schema,
            comparisons=comparisons,
        ).add_pairwise_comparisons(phospho_corrected)
        site_matrix = DatasetSiteMatrix(schema=schema).build(phospho_corrected)
        core_result = CoreProcessingResult(
            total_unique=pd.DataFrame(columns=[TOTAL_GENE_COLUMN, *schema.total_cols]),
            total_filtered=pd.DataFrame(
                columns=[TOTAL_GENE_COLUMN, *schema.total_cols]
            ),
            phospho_filtered=phospho_filtered,
            phospho_corrected=phospho_corrected,
            site_matrix=site_matrix,
        )
        return AnalysisReadyPhosphoDataset.from_core_processing_result(
            core_result,
            schema=schema,
            comparisons=comparisons,
            source="simple kinase workflow (phospho only)",
        )

    @staticmethod
    def _resolve_total_input(
        total: pd.DataFrame | str | Path,
        *,
        dataset_loader: DatasetLoader,
    ) -> pd.DataFrame:
        if isinstance(total, pd.DataFrame):
            return dataset_loader.validate_total(total)
        return dataset_loader.load_total(total)

    @staticmethod
    def _resolve_phospho_input(
        phospho: pd.DataFrame | str | Path,
        *,
        dataset_loader: DatasetLoader,
        phospho_encoding: str | None,
    ) -> pd.DataFrame:
        if isinstance(phospho, pd.DataFrame):
            return dataset_loader.validate_phospho(phospho)
        return dataset_loader.load_phospho(phospho, encoding=phospho_encoding)


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
