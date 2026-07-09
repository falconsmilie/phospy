"""Structured caveats for signalome workflow results."""

from __future__ import annotations

from phospy.contracts.result_caveats import ResultCaveat
from phospy.science.scoring.policy_models import DownstreamScoreSource
from phospy.validation.identity_contracts import (
    validate_reference_context_compatibility,
)
from phospy.workflows.intensity_scale_evidence import (
    build_declared_input_intensity_scale_caveat,
)
from phospy.workflows.kinase.caveats import (
    KINASE_NON_DEFAULT_REFERENCE_SOURCE_CAVEAT_CODE,
    KINASE_REFERENCE_SCORE_FALLBACK_CAVEAT_CODE,
    KINASE_SCORING_LIMITATION_CAVEAT_CODE,
)
from phospy.workflows.result_caveat_helpers import (
    build_localisation_policy_details,
    build_reference_context_compatibility_caveat,
    deduplicate_caveats,
    is_permissive_localisation_requirement,
)
from phospy.workflows.signalome.contracts import ResolvedSignalomeWorkflowRequest

SIGNALOME_UPSTREAM_KINASE_ATTRITION_CAVEAT_CODE = (
    "signalome_upstream_kinase_score_attrition"
)
SIGNALOME_PREDICTION_REFERENCE_LIMITATION_CAVEAT_CODE = (
    "signalome_prediction_reference_limitations"
)
SIGNALOME_PERMISSIVE_LOCALISATION_POLICY_CAVEAT_CODE = (
    "signalome_permissive_localisation_policy"
)
SIGNALOME_PROTEIN_ID_GROUPING_ASSUMPTION_CAVEAT_CODE = (
    "signalome_protein_id_grouping_assumption"
)


def build_signalome_result_caveats(
    *,
    request: ResolvedSignalomeWorkflowRequest,
) -> tuple[ResultCaveat, ...]:
    """Build compact machine-readable caveats for signalome workflow results."""

    caveats: list[ResultCaveat] = []
    declared_input_scale = build_declared_input_intensity_scale_caveat(
        dataset=request.dataset,
        workflow_scope="signalome",
    )
    if declared_input_scale is not None:
        caveats.append(declared_input_scale)
    upstream_attrition = _upstream_kinase_attrition_caveat(request)
    if upstream_attrition is not None:
        caveats.append(upstream_attrition)
    prediction_reference = _prediction_reference_limitation_caveat(request)
    if prediction_reference is not None:
        caveats.append(prediction_reference)
    caveats.extend(_reference_context_unknown_caveats(request))
    localisation = _permissive_localisation_caveat(request)
    if localisation is not None:
        caveats.append(localisation)
    caveats.append(_protein_id_grouping_assumption_caveat(request))
    return deduplicate_caveats(caveats)


def _upstream_kinase_attrition_caveat(
    request: ResolvedSignalomeWorkflowRequest,
) -> ResultCaveat | None:
    attrition = request.kinase_result.attrition_provenance
    if attrition is None:
        return None
    metrics = dict(attrition.metrics)
    scored_fraction = _optional_float(metrics.get("scored_fraction"))
    if scored_fraction is None and attrition.policy_outcome == "passed":
        return None
    if scored_fraction is not None and scored_fraction >= 1.0:
        if attrition.policy_outcome == "passed":
            return None
    details: dict[str, object] = {
        "upstream_policy_outcome": attrition.policy_outcome,
        "downstream_score_source": request.downstream_score_source.value,
        "signalome_score_site_count": int(request.downstream_score_matrix.shape[0]),
        "signalome_score_kinase_count": int(request.downstream_score_matrix.shape[1]),
        "score_preconditioning_policy": (
            request.score_preconditioning_diagnostics.policy
        ),
        "score_preconditioning_retained_row_count": int(
            request.score_preconditioning_diagnostics.retained_row_count
        ),
        "score_preconditioning_dropped_all_missing_row_count": int(
            request.score_preconditioning_diagnostics.dropped_all_missing_row_count
        ),
    }
    for key in (
        "total_dataset_sites",
        "reference_overlap_sites",
        "sequence_supported_sites",
        "scored_sites",
        "reference_overlap_fraction",
        "sequence_supported_fraction",
        "scored_fraction",
    ):
        if key in metrics:
            details[key] = metrics[key]
    return ResultCaveat(
        code=SIGNALOME_UPSTREAM_KINASE_ATTRITION_CAVEAT_CODE,
        severity="warning",
        message=(
            "Signalome workflow inherited upstream kinase score attrition; "
            "module and network outputs only reflect sites retained in the "
            "upstream score matrix."
        ),
        details=details,
    )


def _prediction_reference_limitation_caveat(
    request: ResolvedSignalomeWorkflowRequest,
) -> ResultCaveat | None:
    upstream_codes = {caveat.code for caveat in request.kinase_result.caveats}
    source = request.downstream_score_source
    reference_provenance = request.kinase_result.references.provenance
    source_type = (
        None if reference_provenance is None else reference_provenance.source_type
    )
    should_emit = (
        source is not DownstreamScoreSource.RANK_WEIGHTED_FUSION_SCORES
        or source_type != "bundled"
        or bool(
            upstream_codes
            & {
                KINASE_NON_DEFAULT_REFERENCE_SOURCE_CAVEAT_CODE,
                KINASE_REFERENCE_SCORE_FALLBACK_CAVEAT_CODE,
                KINASE_SCORING_LIMITATION_CAVEAT_CODE,
            }
        )
    )
    if not should_emit:
        return None
    details = {
        "downstream_score_source": source.value,
        "upstream_scoring_mode": request.kinase_result.scoring_result.scoring_mode,
        "upstream_score_scale": request.kinase_result.scoring_result.score_scale,
        "upstream_reference_source_type": source_type or "unknown",
        "upstream_reference_bundle_id": (
            None if reference_provenance is None else reference_provenance.bundle_id
        ),
        "upstream_caveat_codes": sorted(upstream_codes),
    }
    return ResultCaveat(
        code=SIGNALOME_PREDICTION_REFERENCE_LIMITATION_CAVEAT_CODE,
        severity="info",
        message=(
            "Signalome workflow consumes upstream kinase prediction and score "
            "matrices; reference and scoring limitations from the kinase result "
            "carry into signalome summaries."
        ),
        details=details,
    )


def _reference_context_unknown_caveats(
    request: ResolvedSignalomeWorkflowRequest,
) -> tuple[ResultCaveat, ...]:
    caveats: list[ResultCaveat] = []
    dataset_context = request.dataset.reference_context
    result_provenance = request.kinase_result.provenance
    for operation, right_context in (
        (
            "signalome workflow result dataset/upstream kinase result",
            None if result_provenance is None else result_provenance.reference_context,
        ),
        (
            "signalome workflow result dataset/upstream kinase reference",
            None
            if request.kinase_result.references.provenance is None
            else request.kinase_result.references.provenance.reference_context,
        ),
    ):
        warning = validate_reference_context_compatibility(
            dataset_context,
            right_context,
            operation=operation,
            allow_unknown=True,
        )
        if warning is None:
            continue
        caveats.append(
            build_reference_context_compatibility_caveat(
                warning,
                workflow_scope="signalome",
            )
        )
    return tuple(caveats)


def _permissive_localisation_caveat(
    request: ResolvedSignalomeWorkflowRequest,
) -> ResultCaveat | None:
    requirement = request.execution_config.localisation_requirement
    if not is_permissive_localisation_requirement(requirement):
        return None
    details = build_localisation_policy_details(
        site_metadata=request.dataset.site_metadata,
        requirement=requirement,
        workflow_scope="signalome",
    )
    details.update(
        {
            "policy_is_permissive": True,
            "retained_signalome_site_count": int(request.site_to_protein.shape[0]),
        }
    )
    return ResultCaveat(
        code=SIGNALOME_PERMISSIVE_LOCALISATION_POLICY_CAVEAT_CODE,
        severity="warning",
        message=(
            "Signalome workflow localisation policy does not enforce a minimum "
            "localisation probability; unknown or low-confidence phosphosites can "
            "contribute to modules and kinase score-profile associations."
        ),
        details=details,
    )


def _protein_id_grouping_assumption_caveat(
    request: ResolvedSignalomeWorkflowRequest,
) -> ResultCaveat:
    protein_ids = request.site_to_protein.astype(str)
    counts = protein_ids.value_counts(sort=False)
    multi_site_counts = counts.loc[counts > 1]
    details = {
        "grouping_column": "protein_id",
        "grouping_source": "dataset.site_metadata.protein_id",
        "retained_signalome_site_count": int(protein_ids.shape[0]),
        "unique_protein_group_count": int(counts.shape[0]),
        "multi_site_protein_group_count": int(multi_site_counts.shape[0]),
        "max_sites_per_protein_group": (0 if counts.empty else int(counts.max())),
        "site_to_protein_index_matches_score_matrix": bool(
            request.site_to_protein.index.equals(request.downstream_score_matrix.index)
        ),
        "protein_id_grouping_is_explicit": True,
    }
    return ResultCaveat(
        code=SIGNALOME_PROTEIN_ID_GROUPING_ASSUMPTION_CAVEAT_CODE,
        severity="info",
        message=(
            "Signalome workflow groups phosphosites by dataset.site_metadata."
            "protein_id; module context assumes those protein_id values are the "
            "intended protein grouping labels."
        ),
        details=details,
    )


def _optional_float(value: object) -> float | None:
    if isinstance(value, bool):
        return None
    if isinstance(value, (int, float)):
        return float(value)
    return None


__all__ = [
    "SIGNALOME_PERMISSIVE_LOCALISATION_POLICY_CAVEAT_CODE",
    "SIGNALOME_PREDICTION_REFERENCE_LIMITATION_CAVEAT_CODE",
    "SIGNALOME_PROTEIN_ID_GROUPING_ASSUMPTION_CAVEAT_CODE",
    "SIGNALOME_UPSTREAM_KINASE_ATTRITION_CAVEAT_CODE",
    "build_signalome_result_caveats",
]
