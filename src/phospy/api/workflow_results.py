from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd

from ..activities.results import KinaseActivityResult
from ..datasets.models import AnalysisReadyPhosphoDataset
from ..internal.pandas_copy import detached_frame_copy
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


@dataclass(slots=True)
class SimpleKinaseWorkflowResult:
    """Owned result bundle for the high-level common kinase workflow.

    Ownership convention:
    - ``to_*`` methods return detached safe copies
    - ``to_owned_*`` methods return cheap shared owned state
    - ``to_mutable_*_unsafe`` methods return explicit mutable shared state
    """

    analysis_ready_dataset: AnalysisReadyPhosphoDataset
    reference_bundle: ReferenceBundle
    scoring_result: KinaseScoringResult
    prediction_result: KinasePredictionResult
    kinase_activity_result: KinaseActivityResult

    @property
    def pred_mat_result(self) -> PredMatResult:
        """Canonical predMat output for this run."""

        return self.prediction_result.pred_mat_result

    @property
    def profile_scores(self) -> pd.DataFrame:
        """Return a detached profile-based scoring table from the scoring stage."""

        return self.to_profile_scores()

    @property
    def combined_scores(self) -> pd.DataFrame | None:
        """Return detached combined motif/profile scores when available."""

        return self.to_combined_scores()

    @property
    def weights(self) -> pd.DataFrame | None:
        """Return detached score-combination weights when available."""

        return self.to_weights()

    @property
    def substrate_list(self) -> dict[str, list[str]]:
        """Return detached predicted substrate memberships keyed by kinase."""

        return self.to_substrate_list()

    def to_profile_scores(self) -> pd.DataFrame:
        """Return detached profile-based scoring outputs."""

        return detached_frame_copy(self.scoring_result.profile_scores)

    def to_owned_profile_scores(self) -> pd.DataFrame:
        """Return cheap shared owned profile-based scoring outputs (no copy)."""

        return self.scoring_result.profile_scores

    def to_mutable_profile_scores_unsafe(self) -> pd.DataFrame:
        """Return explicit mutable shared profile-based scoring outputs."""

        return self.scoring_result.profile_scores

    def to_combined_scores(self) -> pd.DataFrame | None:
        """Return detached combined motif/profile outputs when available."""

        combined = self.scoring_result.combined_scores
        if combined is None:
            return None
        return detached_frame_copy(combined)

    def to_owned_combined_scores(self) -> pd.DataFrame | None:
        """Return cheap shared owned combined motif/profile outputs."""

        return self.scoring_result.combined_scores

    def to_mutable_combined_scores_unsafe(self) -> pd.DataFrame | None:
        """Return explicit mutable shared combined motif/profile outputs."""

        return self.scoring_result.combined_scores

    def to_weights(self) -> pd.DataFrame | None:
        """Return detached score-combination weights when available."""

        weights = self.scoring_result.weights
        if weights is None:
            return None
        return detached_frame_copy(weights)

    def to_owned_weights(self) -> pd.DataFrame | None:
        """Return cheap shared owned score-combination weights."""

        return self.scoring_result.weights

    def to_mutable_weights_unsafe(self) -> pd.DataFrame | None:
        """Return explicit mutable shared score-combination weights."""

        return self.scoring_result.weights

    def to_substrate_list(self) -> dict[str, list[str]]:
        """Return detached predicted substrate memberships keyed by kinase."""

        if hasattr(self.prediction_result, "to_substrate_list"):
            return self.prediction_result.to_substrate_list()
        return {
            kinase: list(substrates)
            for kinase, substrates in self.prediction_result.substrate_list.items()
        }

    def to_owned_substrate_list(self) -> dict[str, list[str]]:
        """Return cheap shared owned substrate memberships (no copy)."""

        if hasattr(self.prediction_result, "to_owned_substrate_list"):
            return self.prediction_result.to_owned_substrate_list()
        return self.prediction_result.substrate_list

    def to_mutable_substrate_list_unsafe(self) -> dict[str, list[str]]:
        """Return explicit mutable shared substrate memberships."""

        if hasattr(self.prediction_result, "to_mutable_substrate_list_unsafe"):
            return self.prediction_result.to_mutable_substrate_list_unsafe()
        return self.prediction_result.substrate_list

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
        self.prediction_result.close()

    def __enter__(self) -> SimpleKinaseWorkflowResult:
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc: BaseException | None,
        tb: object | None,
    ) -> None:
        self.close()
