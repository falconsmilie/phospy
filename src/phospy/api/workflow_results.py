from __future__ import annotations

from dataclasses import dataclass

from ..activities.results import KinaseActivityResult
from ..datasets import AnalysisReadyPhosphoDataset
from ..motifs import MotifScoringResult
from ..prediction import KinasePredictionResult, KinaseScoringResult, PredMatResult
from ..profiles import KinaseProfileResult
from ..references import ReferenceBundle

__all__ = [
    "KinaseWorkflowResult",
    "PredMatWorkflowResult",
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
