from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING

import pandas as pd
from pydantic import Field, ValidationError, field_validator

from ...errors import RequestValidationError
from ...internal.types import (
    SignalomeAssignmentPolicy,
    SignalomeKinaseNetworkPolicy,
)
from ...signalomes.clustering import SignalomeModuleSelectionPolicy
from ..compatibility import validate_signalome_alignment
from ..domain import (
    resolve_scoring_matrix,
    validate_prediction_result_pred_mat,
    validate_signalome_site_grouping,
)
from ..values.collections import (
    normalize_site_to_protein_mapping,
    normalize_string_sequence,
)
from ..values.enums import (
    validate_kinase_network_policy,
    validate_signalome_assignment_policy,
)
from .shared import PhospyRequestModel

if TYPE_CHECKING:
    from ...prediction.results import KinasePredictionResult, PredMatResult
    from ...prediction.scoring import KinaseScoringResult


class SignalomeRequest(PhospyRequestModel):
    """Raw boundary request for public signalome construction."""

    kinases_of_interest: tuple[str, ...]
    site_to_protein: dict[str, str] | None = None
    kinase_network_threshold: float = Field(default=0.9, ge=0.0, le=1.0)
    kinase_network_policy: SignalomeKinaseNetworkPolicy = "positive_only"
    assignment_policy: SignalomeAssignmentPolicy = "cutoff_binary"
    signalome_cutoff: float = Field(default=0.5, ge=0.0, le=1.0)
    module_count: int | None = Field(default=None, ge=1)
    min_kinase_module_share_percent: float = Field(default=1.0, ge=0.0)
    module_selection_policy: SignalomeModuleSelectionPolicy = Field(
        default_factory=SignalomeModuleSelectionPolicy
    )

    @field_validator("kinases_of_interest", mode="before")
    @classmethod
    def normalize_kinases_of_interest(
        cls,
        value: Sequence[str],
    ) -> tuple[str, ...]:
        return normalize_string_sequence(
            value,
            field_name="kinases_of_interest",
            empty_message="kinases_of_interest must contain at least one kinase name",
            invalid_message=(
                "kinases_of_interest must be provided as a sequence of kinase names"
            ),
            deduplicate=True,
        )

    @field_validator("site_to_protein", mode="before")
    @classmethod
    def normalize_site_to_protein(
        cls,
        value: object,
    ) -> dict[str, str] | None:
        return normalize_site_to_protein_mapping(value)

    @field_validator("module_selection_policy", mode="before")
    @classmethod
    def normalize_module_selection_policy(
        cls,
        value: object,
    ) -> SignalomeModuleSelectionPolicy:
        return SignalomeModuleSelectionPolicy.from_value(value)

    @field_validator("kinase_network_policy", mode="before")
    @classmethod
    def normalize_kinase_network_policy(
        cls,
        value: object,
    ) -> SignalomeKinaseNetworkPolicy:
        if not isinstance(value, str):
            msg = (
                "kinase_network_policy must be one of: 'positive_only', "
                "'absolute_threshold', 'signed'"
            )
            raise TypeError(msg)
        return validate_kinase_network_policy(value)  # type: ignore[arg-type]

    @field_validator("assignment_policy", mode="before")
    @classmethod
    def normalize_assignment_policy(
        cls,
        value: object,
    ) -> SignalomeAssignmentPolicy:
        if not isinstance(value, str):
            msg = "assignment_policy must be one of: 'cutoff_binary', 'weighted_top'"
            raise TypeError(msg)
        return validate_signalome_assignment_policy(value)  # type: ignore[arg-type]

    @classmethod
    def validate_request(cls, **data: object) -> SignalomeRequest:
        try:
            return cls.model_validate(data)
        except ValidationError as error:
            raise RequestValidationError.from_pydantic(
                context="Invalid signalome request",
                error=error,
            ) from error


@dataclass(slots=True)
class SignalomeInputs:
    """Trusted aligned signalome inputs owned by the signalome boundary."""

    scoring_matrix: pd.DataFrame
    pred_mat: pd.DataFrame
    expression_matrix: pd.DataFrame
    site_to_protein: pd.Series
    kinases_of_interest: tuple[str, ...]
    kinase_network_threshold: float
    kinase_network_policy: SignalomeKinaseNetworkPolicy
    assignment_policy: SignalomeAssignmentPolicy
    signalome_cutoff: float
    module_count: int | None
    min_kinase_module_share_percent: float
    module_selection_policy: SignalomeModuleSelectionPolicy

    @classmethod
    def from_trusted_inputs(
        cls,
        *,
        scoring_matrix: pd.DataFrame,
        pred_mat: pd.DataFrame,
        expression_matrix: pd.DataFrame,
        kinases_of_interest: Sequence[str],
        site_to_protein: Mapping[str, str] | pd.Series | None = None,
        kinase_network_threshold: float,
        kinase_network_policy: SignalomeKinaseNetworkPolicy,
        assignment_policy: SignalomeAssignmentPolicy,
        signalome_cutoff: float,
        module_count: int | None,
        min_kinase_module_share_percent: float,
        module_selection_policy: SignalomeModuleSelectionPolicy,
        scoring_context: str,
        pred_mat_context: str,
        expression_context: str,
    ) -> SignalomeInputs:
        (
            validated_scoring_matrix,
            validated_pred_mat,
            validated_expression_matrix,
            common_sites,
        ) = validate_signalome_alignment(
            scoring_matrix=scoring_matrix,
            pred_mat=pred_mat,
            expression_matrix=expression_matrix,
            kinases_of_interest=kinases_of_interest,
            module_count=module_count,
            scoring_context=scoring_context,
            pred_mat_context=pred_mat_context,
            expression_context=expression_context,
        )

        validated_site_to_protein = validate_signalome_site_grouping(
            site_ids=common_sites,
            site_to_protein=site_to_protein,
        )

        return cls(
            scoring_matrix=validated_scoring_matrix,
            pred_mat=validated_pred_mat,
            expression_matrix=validated_expression_matrix,
            site_to_protein=validated_site_to_protein,
            kinases_of_interest=tuple(kinases_of_interest),
            kinase_network_threshold=kinase_network_threshold,
            kinase_network_policy=kinase_network_policy,
            assignment_policy=assignment_policy,
            signalome_cutoff=signalome_cutoff,
            module_count=module_count,
            min_kinase_module_share_percent=min_kinase_module_share_percent,
            module_selection_policy=module_selection_policy,
        )


def validate_signalome_request(
    *,
    scoring_result: KinaseScoringResult,
    prediction_result: KinasePredictionResult | PredMatResult,
    expression_matrix: pd.DataFrame,
    kinases_of_interest: Sequence[str],
    site_to_protein: Mapping[str, str] | None = None,
    kinase_network_threshold: float = 0.9,
    kinase_network_policy: SignalomeKinaseNetworkPolicy = "positive_only",
    assignment_policy: SignalomeAssignmentPolicy = "cutoff_binary",
    signalome_cutoff: float = 0.5,
    module_count: int | None = None,
    min_kinase_module_share_percent: float = 1.0,
    module_selection_policy: SignalomeModuleSelectionPolicy | None = None,
) -> SignalomeInputs:
    """Validate raw signalome inputs and return trusted aligned inputs."""

    request = SignalomeRequest.validate_request(
        kinases_of_interest=kinases_of_interest,
        site_to_protein=site_to_protein,
        kinase_network_threshold=kinase_network_threshold,
        kinase_network_policy=kinase_network_policy,
        assignment_policy=assignment_policy,
        signalome_cutoff=signalome_cutoff,
        module_count=module_count,
        min_kinase_module_share_percent=min_kinase_module_share_percent,
        module_selection_policy=module_selection_policy,
    )

    return SignalomeInputs.from_trusted_inputs(
        scoring_matrix=resolve_scoring_matrix(scoring_result),
        pred_mat=validate_prediction_result_pred_mat(prediction_result),
        expression_matrix=expression_matrix,
        kinases_of_interest=request.kinases_of_interest,
        site_to_protein=request.site_to_protein,
        kinase_network_threshold=request.kinase_network_threshold,
        kinase_network_policy=request.kinase_network_policy,
        assignment_policy=request.assignment_policy,
        signalome_cutoff=request.signalome_cutoff,
        module_count=request.module_count,
        min_kinase_module_share_percent=request.min_kinase_module_share_percent,
        module_selection_policy=request.module_selection_policy,
        scoring_context="scoring_result",
        pred_mat_context="prediction_result",
        expression_context="expression_matrix",
    )
