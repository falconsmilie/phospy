from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Generic, TypeVar

import pandas as pd

from ..api.contracts import PredictionRunConfig
from ..internal.defaults import DEFAULT_MOTIF_FLANK_SIZE
from ..internal.types import PredictionSvmMode
from ..motifs import MotifScoringResult
from ..prediction.engines import (
    KinaseWorkflowExecutionResult,
    KinaseWorkflowExecutor,
)
from ..prediction.results import KinasePredictionResult, PredMatResult
from ..prediction.scoring import KinaseScoringResult
from ..profiles import KinaseProfileResult
from ..references import ReferenceBundle
from ..validation.requests.workflow import WorkflowInputs

__all__ = ["KinaseWorkflow", "PredMatWorkflow"]


@dataclass(slots=True)
class KinaseWorkflowResult:
    profile_result: KinaseProfileResult
    motif_result: MotifScoringResult | None
    scoring_result: KinaseScoringResult
    prediction_result: KinasePredictionResult


@dataclass(slots=True)
class PredMatWorkflowResult:
    scoring_result: KinaseScoringResult
    prediction_result: KinasePredictionResult
    pred_mat_result: PredMatResult

    def close(self) -> None:
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


_WorkflowResultT = TypeVar(
    "_WorkflowResultT",
    KinaseWorkflowResult,
    PredMatWorkflowResult,
)


class _BaseKinaseWorkflow(Generic[_WorkflowResultT]):
    """Shared validated execution path for internal kinase workflow adapters."""

    def __init__(
        self,
        flank_size: int = DEFAULT_MOTIF_FLANK_SIZE,
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
        prediction_config: PredictionRunConfig | None = None,
    ) -> WorkflowInputs:
        resolved_config = PredictionRunConfig.from_value(prediction_config)
        return self._executor.validate_request(
            phospho_matrix=phospho_matrix,
            substrate_map=substrate_map,
            site_sequences=site_sequences,
            motif_sequences=motif_sequences,
            reference_bundle=reference_bundle,
            min_substrates=resolved_config.min_substrates,
            min_motif_size=resolved_config.min_motif_size,
            allow_profile_only_fallback=resolved_config.allow_profile_only_fallback,
            ensemble_size=resolved_config.ensemble_size,
            top=resolved_config.top,
            score_threshold=resolved_config.score_threshold,
            inclusion=resolved_config.inclusion,
            n_iterations=resolved_config.n_iterations,
            random_state=resolved_config.random_state,
            svm_mode=resolved_config.svm_mode,
            profile_policy=resolved_config.profile_policy,
        )

    def run(
        self,
        phospho_matrix: pd.DataFrame,
        substrate_map: Mapping[str, Sequence[str]] | None = None,
        site_sequences: Mapping[str, str] | pd.Series | None = None,
        motif_sequences: Mapping[str, Sequence[str]] | None = None,
        reference_bundle: ReferenceBundle | None = None,
        prediction_config: PredictionRunConfig | None = None,
    ) -> _WorkflowResultT:
        request = self._validate_request(
            phospho_matrix=phospho_matrix,
            substrate_map=substrate_map,
            site_sequences=site_sequences,
            motif_sequences=motif_sequences,
            reference_bundle=reference_bundle,
            prediction_config=prediction_config,
        )
        return self.run_validated(request)

    def run_validated(
        self,
        request: WorkflowInputs,
    ) -> _WorkflowResultT:
        result = self._executor.execute_validated_request(request)
        return self._package_result(result)

    def _package_result(
        self,
        result: KinaseWorkflowExecutionResult,
    ) -> _WorkflowResultT:
        raise NotImplementedError


class KinaseWorkflow(_BaseKinaseWorkflow[KinaseWorkflowResult]):
    """Run the native kinase scoring and prediction workflow end to end."""

    def _package_result(
        self,
        result: KinaseWorkflowExecutionResult,
    ) -> KinaseWorkflowResult:
        return KinaseWorkflowResult(
            profile_result=result.profile_result,
            motif_result=result.motif_result,
            scoring_result=result.scoring_result,
            prediction_result=result.prediction_result,
        )


class PredMatWorkflow(_BaseKinaseWorkflow[PredMatWorkflowResult]):
    """Generate a predMat from phosphosite and sequence inputs."""

    def _package_result(
        self,
        result: KinaseWorkflowExecutionResult,
    ) -> PredMatWorkflowResult:
        return PredMatWorkflowResult(
            scoring_result=result.scoring_result,
            prediction_result=result.prediction_result,
            pred_mat_result=result.prediction_result.pred_mat_result,
        )
