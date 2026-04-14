from __future__ import annotations

from dataclasses import dataclass

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
    pred_mat_result: PredMatResult
    kinase_activity_result: KinaseActivityResult

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
