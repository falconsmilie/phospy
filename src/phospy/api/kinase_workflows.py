from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Generic, TypeVar

import pandas as pd

from ..internal.types import PredictionSvmMode
from ..prediction.engines import (
    KinaseWorkflowExecutionResult,
    KinaseWorkflowExecutor,
)
from ..profiles import KinaseProfilePolicy
from ..references import ReferenceBundle
from ..validation.requests.workflow import WorkflowInputs
from .workflow_results import KinaseWorkflowResult, PredMatWorkflowResult

__all__ = ["KinaseWorkflow", "PredMatWorkflow"]

_WorkflowResultT = TypeVar(
    "_WorkflowResultT",
    KinaseWorkflowResult,
    PredMatWorkflowResult,
)


class _BaseKinaseWorkflow(Generic[_WorkflowResultT]):
    """Shared validated execution path for native kinase workflow adapters."""

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
        profile_policy: KinaseProfilePolicy | None = None,
    ) -> WorkflowInputs:
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
            profile_policy=profile_policy,
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
        profile_policy: KinaseProfilePolicy | None = None,
    ) -> _WorkflowResultT:
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
            profile_policy=profile_policy,
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
