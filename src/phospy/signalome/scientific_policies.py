"""Scientific policy records for signalome score preconditioning."""

from __future__ import annotations

from dataclasses import dataclass

from phospy.provenance.scientific_policy_models import (
    ScientificPolicyId,
    ScientificPolicyRecord,
)


@dataclass(frozen=True, slots=True)
class ScorePreconditioningPolicy:
    """Executable metadata policy for score preconditioning behavior."""

    policy: str
    input_row_count: int
    dropped_all_missing_row_count: int
    retained_row_count: int
    row_retention_rule: str = "drop_rows_with_all_scores_missing"
    retained_partial_missing_rows: bool = True

    @property
    def record(self) -> ScientificPolicyRecord:
        return build_score_preconditioning_policy(
            policy=self.policy,
            input_row_count=self.input_row_count,
            dropped_all_missing_row_count=self.dropped_all_missing_row_count,
            retained_row_count=self.retained_row_count,
            row_retention_rule=self.row_retention_rule,
            retained_partial_missing_rows=self.retained_partial_missing_rows,
        )


def build_score_preconditioning_policy(
    *,
    policy: str,
    input_row_count: int,
    dropped_all_missing_row_count: int,
    retained_row_count: int,
    row_retention_rule: str = "drop_rows_with_all_scores_missing",
    retained_partial_missing_rows: bool = True,
) -> ScientificPolicyRecord:
    base_policy = resolve_score_preconditioning_policy(policy=policy)
    return ScientificPolicyRecord(
        id=ScientificPolicyId.SIGNALOME_SCORE_PRECONDITIONING,
        name=base_policy.name,
        version="1",
        description=(
            "Preconditions aligned downstream score rows before signalome "
            "construction by handling unsupported all-missing rows explicitly."
        ),
        parameters={
            "policy": str(policy),
            "row_retention_rule": str(row_retention_rule),
            "retained_partial_missing_rows": bool(retained_partial_missing_rows),
            "input_row_count": int(input_row_count),
            "dropped_all_missing_row_count": int(dropped_all_missing_row_count),
            "retained_row_count": int(retained_row_count),
        },
        assumptions=(
            "All-missing score rows are scientifically unsupported for score-driven "
            "signalome construction.",
            "Preconditioning policy determines whether row dropping is allowed or "
            "treated as a boundary error.",
            "Row retention changes site coverage and therefore can change final "
            "score-derived signalome assignments and module summaries.",
        ),
        output_scale="Retained downstream score matrix rows for signalome execution.",
        quantitative_meaning="retained_signalome_score_rows",
    )


SCORE_PRECONDITIONING_ERROR_ON_DROP_POLICY = ScientificPolicyRecord(
    id=ScientificPolicyId.SIGNALOME_SCORE_PRECONDITIONING,
    name="score_preconditioning_error_on_drop_v1",
    version="1",
    description=(
        "Treat dropped all-missing downstream score rows as a hard workflow boundary "
        "error."
    ),
    parameters={
        "policy": "error_on_drop",
        "row_retention_rule": "drop_rows_with_all_scores_missing",
        "retained_partial_missing_rows": True,
    },
    assumptions=(
        "All interpreted sites must retain at least one finite downstream score.",
        "Any unsupported all-missing row invalidates signalome construction.",
    ),
    output_scale="Validation-only policy for preconditioning gate behavior.",
    quantitative_meaning="preconditioning_validation_rule",
)


SCORE_PRECONDITIONING_ALLOW_AND_REPORT_POLICY = ScientificPolicyRecord(
    id=ScientificPolicyId.SIGNALOME_SCORE_PRECONDITIONING,
    name="score_preconditioning_allow_and_report_v1",
    version="1",
    description=(
        "Allow dropping all-missing downstream score rows and report the drop "
        "diagnostics in provenance."
    ),
    parameters={
        "policy": "allow_and_report",
        "row_retention_rule": "drop_rows_with_all_scores_missing",
        "retained_partial_missing_rows": True,
    },
    assumptions=(
        "All-missing rows provide no usable downstream score evidence.",
        "Retained-site coverage can change as unsupported rows are removed.",
    ),
    output_scale="Validation-and-row-retention policy for score preconditioning.",
    quantitative_meaning="preconditioning_row_retention_rule",
)


def resolve_score_preconditioning_policy(*, policy: str) -> ScientificPolicyRecord:
    if policy == "error_on_drop":
        return SCORE_PRECONDITIONING_ERROR_ON_DROP_POLICY
    if policy == "allow_and_report":
        return SCORE_PRECONDITIONING_ALLOW_AND_REPORT_POLICY
    return ScientificPolicyRecord(
        id=ScientificPolicyId.SIGNALOME_SCORE_PRECONDITIONING,
        name=f"score_preconditioning_{policy}_v1",
        version="1",
        description="Configured score preconditioning policy.",
        parameters={"policy": str(policy)},
        assumptions=(),
        output_scale="Validation-and-row-retention policy for score preconditioning.",
        quantitative_meaning="preconditioning_row_retention_rule",
    )


__all__ = [
    "SCORE_PRECONDITIONING_ALLOW_AND_REPORT_POLICY",
    "SCORE_PRECONDITIONING_ERROR_ON_DROP_POLICY",
    "ScorePreconditioningPolicy",
    "build_score_preconditioning_policy",
    "resolve_score_preconditioning_policy",
]
