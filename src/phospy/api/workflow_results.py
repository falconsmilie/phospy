from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from ..activities.results import KinaseActivityResult
from ..datasets.models import AnalysisReadyPhosphoDataset
from ..prediction.results import KinasePredictionResult, PredMatResult
from ..prediction.scoring import KinaseScoringResult
from ..references import ReferenceBundle

__all__ = ["SimpleKinaseWorkflowResult"]


@dataclass(slots=True)
class SimpleKinaseWorkflowResult:
    """Owned result bundle for the high-level common kinase workflow."""

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
        """Profile-based scoring table from the scoring stage."""

        return self.scoring_result.profile_scores

    @property
    def combined_scores(self) -> pd.DataFrame | None:
        """Combined motif/profile scores when motif scoring is available."""

        return self.scoring_result.combined_scores

    @property
    def weights(self) -> pd.DataFrame | None:
        """Score-combination weights when motif scoring is available."""

        return self.scoring_result.weights

    @property
    def substrate_list(self) -> dict[str, list[str]]:
        """Predicted substrate memberships keyed by kinase."""

        return self.prediction_result.substrate_list

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
