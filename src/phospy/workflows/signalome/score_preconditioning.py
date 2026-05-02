"""Score preconditioning for signalome interpretation."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from phospy.api.configs import (
    SIGNALOME_SCORE_PRECONDITIONING_POLICY_ALLOW_AND_REPORT,
    SIGNALOME_SCORE_PRECONDITIONING_POLICY_ERROR_ON_DROP,
    SignalomeScorePreconditioningPolicy,
)
from phospy.signalomes.models import SignalomeScorePreconditioningDiagnostics
from phospy.workflows.signalome.boundary_errors import (
    raise_signalome_boundary_error,
    raise_wrapped_signalome_boundary_error,
)
from phospy.workflows.signalome.constants import (
    SIGNALOME_INTERPRETER_SCORE_PRECONDITIONING_SEAM,
)


@dataclass(frozen=True, slots=True)
class SignalomeScorePreconditioningResult:
    downstream_score_matrix: pd.DataFrame
    diagnostics: SignalomeScorePreconditioningDiagnostics


class SignalomeScorePreconditioner:
    """Prepare aligned downstream scores for score-driven signalome stages."""

    def run(
        self,
        *,
        score_matrix: pd.DataFrame,
        policy: SignalomeScorePreconditioningPolicy,
    ) -> SignalomeScorePreconditioningResult:
        if score_matrix.empty:
            return SignalomeScorePreconditioningResult(
                downstream_score_matrix=score_matrix,
                diagnostics=self._diagnostics(
                    input_row_count=0,
                    dropped_all_missing_row_count=0,
                    retained_row_count=0,
                    policy=policy,
                ),
            )
        try:
            score_values = score_matrix.to_numpy(dtype=float, copy=False)
            infinite_mask = np.isinf(score_values)
        except (TypeError, ValueError) as exc:
            raise_wrapped_signalome_boundary_error(
                stage_name="signalome.score_preconditioning_conversion",
                seam="signalome.interpreter.score_preconditioning_conversion",
                field_name=(
                    "signalome workflow request kinase_result.scoring_result."
                    "downstream_score_matrix"
                ),
                operation="converting downstream score matrix to float for finite checks",
                next_action=(
                    "ensure kinase scoring outputs contain numeric finite values "
                    "before running SignalomeWorkflow"
                ),
                original_error=exc,
                aligned_score_sites=int(score_matrix.shape[0]),
                aligned_score_kinases=int(score_matrix.shape[1]),
            )
        if infinite_mask.any():
            raise_signalome_boundary_error(
                seam=SIGNALOME_INTERPRETER_SCORE_PRECONDITIONING_SEAM,
                next_action=(
                    "rerun kinase workflow and ensure scoring outputs contain "
                    "finite values only"
                ),
                aligned_score_sites=int(score_matrix.shape[0]),
                aligned_score_kinases=int(score_matrix.shape[1]),
                infinite_score_entries=int(infinite_mask.sum()),
            )
        supported_rows = score_matrix.notna().any(axis=1)
        supported_row_mask = supported_rows.to_numpy(dtype=bool, copy=False)
        input_row_count = int(score_matrix.shape[0])
        retained_row_count = int(supported_row_mask.sum())
        dropped_all_missing_row_count = int(input_row_count - retained_row_count)
        diagnostics = self._diagnostics(
            input_row_count=input_row_count,
            dropped_all_missing_row_count=dropped_all_missing_row_count,
            retained_row_count=retained_row_count,
            policy=policy,
        )
        if policy == SIGNALOME_SCORE_PRECONDITIONING_POLICY_ERROR_ON_DROP:
            if dropped_all_missing_row_count > 0:
                raise_signalome_boundary_error(
                    seam=SIGNALOME_INTERPRETER_SCORE_PRECONDITIONING_SEAM,
                    next_action=(
                        "set "
                        "config.validation.score_preconditioning_policy="
                        "'allow_and_report' "
                        "to proceed with explicit row dropping, or ensure upstream "
                        "downstream scores contain non-missing support for every "
                        "interpreted site"
                    ),
                    aligned_score_sites=input_row_count,
                    aligned_score_kinases=int(score_matrix.shape[1]),
                    dropped_all_missing_row_count=dropped_all_missing_row_count,
                    retained_row_count=retained_row_count,
                    score_preconditioning_policy=policy,
                )
            return SignalomeScorePreconditioningResult(
                downstream_score_matrix=score_matrix,
                diagnostics=diagnostics,
            )
        if policy != SIGNALOME_SCORE_PRECONDITIONING_POLICY_ALLOW_AND_REPORT:
            raise_signalome_boundary_error(
                seam=SIGNALOME_INTERPRETER_SCORE_PRECONDITIONING_SEAM,
                next_action=(
                    "use a supported score preconditioning policy from "
                    "SignalomeConfig.validation.score_preconditioning_policy"
                ),
                score_preconditioning_policy=policy,
            )
        if supported_row_mask.all():
            return SignalomeScorePreconditioningResult(
                downstream_score_matrix=score_matrix,
                diagnostics=diagnostics,
            )
        return SignalomeScorePreconditioningResult(
            downstream_score_matrix=score_matrix.iloc[supported_row_mask, :],
            diagnostics=diagnostics,
        )

    @staticmethod
    def _diagnostics(
        *,
        input_row_count: int,
        dropped_all_missing_row_count: int,
        retained_row_count: int,
        policy: SignalomeScorePreconditioningPolicy,
    ) -> SignalomeScorePreconditioningDiagnostics:
        return SignalomeScorePreconditioningDiagnostics(
            input_row_count=int(input_row_count),
            dropped_all_missing_row_count=int(dropped_all_missing_row_count),
            retained_row_count=int(retained_row_count),
            policy=policy,
        )


__all__ = [
    "SignalomeScorePreconditioner",
    "SignalomeScorePreconditioningResult",
]
