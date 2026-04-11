from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import TYPE_CHECKING

import pandas as pd

from ...signalome_site_ids import resolve_signalome_site_to_protein
from ..errors import NoCandidateKinasesError, RequestValidationError
from ..schema.tables import PredMatSchema

if TYPE_CHECKING:
    from ...prediction.results import KinasePredictionResult, PredMatResult
    from ...prediction.scoring import KinaseScoringResult


def validate_signalome_site_grouping(
    *,
    site_ids: Sequence[str],
    site_to_protein: Mapping[str, str] | None,
) -> pd.Series:
    """Validate or derive the site-to-protein grouping used for signalomes."""

    return resolve_signalome_site_to_protein(
        site_ids=site_ids,
        site_to_protein=site_to_protein,
        missing_mapping_context=(
            "site_to_protein must define a protein ID for every aligned phosphosite "
            "row. Missing mappings for"
        ),
        invalid_mapping_context=(
            "site_to_protein must map aligned phosphosite rows to non-empty protein "
            "IDs. Invalid mappings for"
        ),
        invalid_site_id_context=(
            "Signalome construction requires either an explicit site_to_protein "
            "mapping or phosphosite identifiers in the supported 'PROTEIN;SITE;...' "
            "format. Invalid aligned site IDs"
        ),
    )


def resolve_scoring_matrix(scoring_result: KinaseScoringResult) -> pd.DataFrame:
    from ...prediction.scoring import KinaseScoringResult

    if not isinstance(scoring_result, KinaseScoringResult):
        msg = "scoring_result must be a KinaseScoringResult"
        raise RequestValidationError(msg)

    if scoring_result.combined_scores is not None:
        return scoring_result.combined_scores
    return scoring_result.profile_scores


def validate_prediction_result_pred_mat(
    prediction_result: KinasePredictionResult | PredMatResult,
) -> pd.DataFrame:
    pred_mat = resolve_pred_mat(prediction_result)
    if pred_mat.shape[1] == 0:
        msg = (
            "prediction_result does not contain any kinase columns because no "
            "candidate kinases qualified for prediction. Regenerate predMat with "
            "less restrictive top, score_threshold, or inclusion settings."
        )
        raise NoCandidateKinasesError(msg)
    return PredMatSchema.validate(pred_mat, context="prediction_result")


def resolve_pred_mat(
    prediction_result: KinasePredictionResult | PredMatResult,
) -> pd.DataFrame:
    from ...prediction.results import KinasePredictionResult, PredMatResult

    if isinstance(prediction_result, KinasePredictionResult):
        return prediction_result.pred_mat_result.to_frame(copy=False)
    if isinstance(prediction_result, PredMatResult):
        return prediction_result.to_frame(copy=False)
    msg = "prediction_result must be a KinasePredictionResult or PredMatResult"
    raise RequestValidationError(msg)


__all__ = [
    "resolve_pred_mat",
    "resolve_scoring_matrix",
    "validate_prediction_result_pred_mat",
    "validate_signalome_site_grouping",
]
