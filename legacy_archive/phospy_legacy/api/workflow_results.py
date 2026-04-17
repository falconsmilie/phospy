from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

from ..activities.results import KinaseActivityResult
from ..datasets.models import AnalysisReadyPhosphoDataset
from ..internal.pandas_copy import detached_frame_copy, detached_series_copy
from ..prediction.results import KinasePredictionResult, PredMatResult
from ..prediction.scoring import KinaseScoringResult
from ..references import ReferenceBundle

if TYPE_CHECKING:
    from ..io.readers import SimpleKinaseWorkflowOutputBundle
    from ..preprocessing import CorePreprocessingConfig
    from .contracts import (
        DatasetLoadOptions,
        KinaseActivityConfig,
        PredictionRunConfig,
        SimpleKinaseWorkflowBundleMetadata,
        SimpleKinaseWorkflowConfigSnapshot,
    )

__all__ = ["SimpleKinaseWorkflowResult"]


@dataclass(slots=True, init=False)
class SimpleKinaseWorkflowResult:
    """Owned result bundle for the high-level common kinase workflow.

    The primary public contract is nested:
    - ``scoring_result`` for ``profile_scores``, ``combined_scores``, and ``weights``
    - ``prediction_result`` for ``pred_matrix`` and ``substrate_list``
    - ``kinase_activity_result`` for downstream activity tables

    ``pred_mat_result`` remains as a convenience seam for the canonical predMat table.

    Ownership convention:
    - ``to_*`` methods return detached safe copies
    - ``to_owned_*`` methods return cheap shared owned state
    - ``to_mutable_*_unsafe`` methods return explicit mutable shared state
    """

    _analysis_ready_dataset: AnalysisReadyPhosphoDataset
    _reference_bundle: ReferenceBundle
    _scoring_result: KinaseScoringResult
    _prediction_result: KinasePredictionResult
    _kinase_activity_result: KinaseActivityResult

    def __init__(
        self,
        *,
        analysis_ready_dataset: AnalysisReadyPhosphoDataset,
        reference_bundle: ReferenceBundle,
        scoring_result: KinaseScoringResult,
        prediction_result: KinasePredictionResult,
        kinase_activity_result: KinaseActivityResult,
    ) -> None:
        self._analysis_ready_dataset = analysis_ready_dataset
        self._reference_bundle = reference_bundle
        self._scoring_result = scoring_result
        self._prediction_result = prediction_result
        self._kinase_activity_result = kinase_activity_result

    @property
    def analysis_ready_dataset(self) -> AnalysisReadyPhosphoDataset:
        """Return a detached analysis-ready dataset boundary."""

        if not isinstance(self._analysis_ready_dataset, AnalysisReadyPhosphoDataset):
            return self._analysis_ready_dataset
        return self.to_analysis_ready_dataset()

    @property
    def reference_bundle(self) -> ReferenceBundle:
        """Return the reference bundle associated with this workflow run."""

        return self._reference_bundle

    @property
    def scoring_result(self) -> KinaseScoringResult:
        """Return detached scoring tables wrapped as a scoring result object."""

        if not isinstance(self._scoring_result, KinaseScoringResult):
            return self._scoring_result
        return self.to_scoring_result()

    @property
    def prediction_result(self) -> KinasePredictionResult:
        """Return a detached prediction result object."""

        if not isinstance(self._prediction_result, KinasePredictionResult):
            return self._prediction_result
        return self.to_prediction_result()

    @property
    def kinase_activity_result(self) -> KinaseActivityResult:
        """Return detached kinase activity tables wrapped as a result object."""

        if not isinstance(self._kinase_activity_result, KinaseActivityResult):
            return self._kinase_activity_result
        return self.to_kinase_activity_result()

    @property
    def pred_mat_result(self) -> PredMatResult:
        """Canonical predMat output for this run."""

        return self._prediction_result.pred_mat_result

    def to_analysis_ready_dataset(self) -> AnalysisReadyPhosphoDataset:
        """Return a detached analysis-ready dataset boundary."""

        dataset = self._analysis_ready_dataset
        if not isinstance(dataset, AnalysisReadyPhosphoDataset):
            return dataset
        return AnalysisReadyPhosphoDataset.from_external(
            phospho_matrix=dataset.to_owned_phospho_matrix(),
            site_metadata=dataset.to_owned_site_metadata(),
            site_sequences=dataset.to_owned_site_sequences(),
            phospho_corrected=dataset.to_owned_phospho_corrected(),
            provenance=dataset.provenance,
        )

    def to_owned_analysis_ready_dataset(self) -> AnalysisReadyPhosphoDataset:
        """Return cheap shared owned analysis-ready dataset state."""

        return self._analysis_ready_dataset

    def to_mutable_analysis_ready_dataset_unsafe(self) -> AnalysisReadyPhosphoDataset:
        """Return explicit mutable shared analysis-ready dataset state."""

        return self._analysis_ready_dataset

    def to_scoring_result(self) -> KinaseScoringResult:
        """Return detached scoring tables wrapped as a scoring result object."""

        if not isinstance(self._scoring_result, KinaseScoringResult):
            return self._scoring_result
        combined_scores = self._scoring_result.combined_scores
        weights = self._scoring_result.weights
        return KinaseScoringResult(
            profile_scores=detached_frame_copy(self._scoring_result.profile_scores),
            combined_scores=(
                None
                if combined_scores is None
                else detached_frame_copy(combined_scores)
            ),
            weights=None if weights is None else detached_frame_copy(weights),
        )

    def to_owned_scoring_result(self) -> KinaseScoringResult:
        """Return cheap shared owned scoring-result state."""

        return self._scoring_result

    def to_mutable_scoring_result_unsafe(self) -> KinaseScoringResult:
        """Return explicit mutable shared scoring-result state."""

        return self._scoring_result

    def to_prediction_result(self) -> KinasePredictionResult:
        """Return a detached prediction result object."""

        if not isinstance(self._prediction_result, KinasePredictionResult):
            return self._prediction_result
        return KinasePredictionResult(
            pred_matrix=self._prediction_result.to_pred_matrix(),
            substrate_list=self._prediction_result.to_substrate_list(),
            debug_traces=self._prediction_result.debug_traces,
            trace_level=self._prediction_result.trace_level,
            trace_sink=None,
            owns_trace_sink=False,
        )

    def to_owned_prediction_result(self) -> KinasePredictionResult:
        """Return cheap shared owned prediction-result state."""

        return self._prediction_result

    def to_mutable_prediction_result_unsafe(self) -> KinasePredictionResult:
        """Return explicit mutable shared prediction-result state."""

        return self._prediction_result

    def to_kinase_activity_result(self) -> KinaseActivityResult:
        """Return detached kinase activity tables wrapped as a result object."""

        activity = self._kinase_activity_result
        if not isinstance(activity, KinaseActivityResult):
            return activity
        return KinaseActivityResult(
            weighted_activity=detached_frame_copy(activity.weighted_activity),
            ksea_scores=detached_frame_copy(activity.ksea_scores),
            ksea_counts=detached_series_copy(activity.ksea_counts),
            target_counts=detached_series_copy(activity.target_counts),
            target_table=detached_frame_copy(activity.target_table),
            overlap_summary=activity.overlap_summary,
        )

    def to_owned_kinase_activity_result(self) -> KinaseActivityResult:
        """Return cheap shared owned kinase-activity-result state."""

        return self._kinase_activity_result

    def to_mutable_kinase_activity_result_unsafe(self) -> KinaseActivityResult:
        """Return explicit mutable shared kinase-activity-result state."""

        return self._kinase_activity_result

    def save_output_bundle(
        self,
        outdir: str | Path,
        *,
        config_snapshot: SimpleKinaseWorkflowConfigSnapshot
        | Mapping[str, object]
        | None = None,
        dataset_options: DatasetLoadOptions | Mapping[str, object] | None = None,
        preprocessing_config: CorePreprocessingConfig | None = None,
        prediction_config: PredictionRunConfig | Mapping[str, object] | None = None,
        activity_config: KinaseActivityConfig | Mapping[str, object] | None = None,
    ) -> Path:
        """Persist this workflow result to an explicit output-bundle directory."""

        from ..io.publishing import save_simple_kinase_workflow_output_bundle

        return save_simple_kinase_workflow_output_bundle(
            result=self,
            outdir=outdir,
            config_snapshot=config_snapshot,
            dataset_options=dataset_options,
            preprocessing_config=preprocessing_config,
            prediction_config=prediction_config,
            activity_config=activity_config,
        )

    @staticmethod
    def load_output_bundle_metadata(
        bundle_dir: str | Path,
    ) -> SimpleKinaseWorkflowBundleMetadata:
        """Load metadata from a previously saved workflow output bundle."""

        from ..io.publishing import load_simple_kinase_workflow_output_bundle_metadata

        return load_simple_kinase_workflow_output_bundle_metadata(bundle_dir)

    @staticmethod
    def load_output_bundle(
        bundle_dir: str | Path,
        *,
        table_ids: Sequence[str] | None = None,
    ) -> SimpleKinaseWorkflowOutputBundle:
        """Load metadata and selected tables from a saved output bundle."""

        from ..io.publishing import load_simple_kinase_workflow_output_bundle

        return load_simple_kinase_workflow_output_bundle(
            bundle_dir,
            table_ids=tuple(table_ids) if table_ids is not None else None,
        )

    def close(self) -> None:
        self._prediction_result.close()

    def __enter__(self) -> SimpleKinaseWorkflowResult:
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc: BaseException | None,
        tb: object | None,
    ) -> None:
        self.close()
