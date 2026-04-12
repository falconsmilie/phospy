from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import pandas as pd
from pydantic import Field, ValidationError, field_validator, model_validator

from ...errors import RequestValidationError
from ...internal.types import PredictionSvmMode
from ...motifs import KinaseMotifScorer
from ...references import ReferenceBundle
from ..compatibility import validate_workflow_matrix_inputs
from ..domain import resolve_reference_bundle_inputs
from ..values.collections import (
    normalize_sequence_mapping,
    normalize_site_sequence_series,
)
from .shared import PhospyRequestModel


class KinaseWorkflowRequest(PhospyRequestModel):
    """Raw boundary request for native kinase workflow execution."""

    phospho_matrix: pd.DataFrame
    substrate_map: dict[str, tuple[str, ...]]
    site_sequences: pd.Series | None = None
    motif_sequences: dict[str, tuple[str, ...]] | None = None
    min_substrates: int = Field(default=1, ge=1)
    min_motif_size: int = Field(default=1, ge=1)
    allow_profile_only_fallback: bool = False
    ensemble_size: int = Field(default=10, ge=1)
    top: int = Field(default=50, ge=1)
    score_threshold: float = Field(default=0.8, ge=0.0, le=1.0)
    inclusion: int = Field(default=20, ge=1)
    n_iterations: int = Field(default=5, ge=1)
    random_state: int | None = None
    svm_mode: PredictionSvmMode | None = None

    @field_validator("substrate_map", mode="before")
    @classmethod
    def validate_substrate_map(
        cls,
        value: Mapping[str, Sequence[str]],
    ) -> dict[str, tuple[str, ...]]:
        return normalize_sequence_mapping(
            value,
            field_name="substrate_map",
            empty_message="substrate_map must not be empty",
        )

    @field_validator("site_sequences", mode="before")
    @classmethod
    def validate_site_sequences(
        cls,
        value: Mapping[str, str] | pd.Series | None,
    ) -> pd.Series | None:
        return normalize_site_sequence_series(value)

    @field_validator("motif_sequences", mode="before")
    @classmethod
    def validate_motif_sequences(
        cls,
        value: Mapping[str, Sequence[str]] | None,
    ) -> dict[str, tuple[str, ...]] | None:
        if value is None:
            return None
        return normalize_sequence_mapping(
            value,
            field_name="motif_sequences",
            empty_message=(
                "motif_sequences must not be empty; pass None and set "
                "allow_profile_only_fallback=True for profile-only prediction"
            ),
        )

    @model_validator(mode="after")
    def validate_cross_field_requirements(self) -> KinaseWorkflowRequest:
        if self.motif_sequences is None and not self.allow_profile_only_fallback:
            msg = (
                "motif_sequences are required for end-to-end prediction unless "
                "allow_profile_only_fallback=True"
            )
            raise ValueError(msg)

        if self.motif_sequences is not None and self.site_sequences is None:
            msg = "site_sequences are required when motif_sequences are provided"
            raise ValueError(msg)

        return self

    @classmethod
    def validate_request(cls, **data: object) -> KinaseWorkflowRequest:
        try:
            return cls.model_validate(data)
        except ValidationError as error:
            raise RequestValidationError.from_pydantic(
                context="Invalid kinase workflow request",
                error=error,
            ) from error


@dataclass(slots=True)
class WorkflowInputs:
    """Trusted workflow inputs owned by the workflow boundary."""

    request: KinaseWorkflowRequest
    phospho_matrix: pd.DataFrame
    scoring_site_index: tuple[str, ...]
    motif_scorer: KinaseMotifScorer | None
    predictor_svm_mode: PredictionSvmMode


def build_workflow_inputs(
    request: KinaseWorkflowRequest,
    *,
    flank_size: int,
    default_svm_mode: PredictionSvmMode,
    context: str = "Kinase workflow inputs",
) -> WorkflowInputs:
    """Build trusted workflow inputs from a validated raw request."""

    validated_matrix, scoring_site_index = validate_workflow_matrix_inputs(
        request.phospho_matrix,
        request.substrate_map,
        request.site_sequences,
        require_site_sequences_for_prediction=request.motif_sequences is not None,
        context=context,
    )
    owned_request = _copy_workflow_request_owned_state(
        request,
        phospho_matrix=validated_matrix,
    )
    motif_scorer = (
        None
        if request.motif_sequences is None
        else KinaseMotifScorer.from_substrate_sequences(
            motif_sequences=request.motif_sequences,
            flank_size=flank_size,
        )
    )
    return WorkflowInputs(
        request=owned_request,
        phospho_matrix=validated_matrix,
        scoring_site_index=scoring_site_index,
        motif_scorer=motif_scorer,
        predictor_svm_mode=(
            default_svm_mode if request.svm_mode is None else request.svm_mode
        ),
    )


def validate_workflow_request(
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
    flank_size: int = 7,
    default_svm_mode: PredictionSvmMode = "default",
    context: str = "Kinase workflow inputs",
) -> WorkflowInputs:
    """Validate raw workflow inputs and return trusted workflow inputs."""

    resolved_substrate_map, resolved_motif_sequences = resolve_reference_bundle_inputs(
        substrate_map=substrate_map,
        motif_sequences=motif_sequences,
        reference_bundle=reference_bundle,
    )
    request = KinaseWorkflowRequest.validate_request(
        phospho_matrix=phospho_matrix,
        substrate_map=resolved_substrate_map,
        site_sequences=site_sequences,
        motif_sequences=resolved_motif_sequences,
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
    )
    return build_workflow_inputs(
        request,
        flank_size=flank_size,
        default_svm_mode=default_svm_mode,
        context=context,
    )


def validate_workflow_inputs(
    phospho_matrix: pd.DataFrame,
    substrate_map: Mapping[str, Sequence[str]] | None = None,
    site_sequences: Mapping[str, str] | pd.Series | None = None,
    motif_sequences: Mapping[str, Sequence[str]] | None = None,
    *,
    reference_bundle: ReferenceBundle | None = None,
    flank_size: int = 7,
    context: str = "Kinase workflow inputs",
) -> pd.DataFrame:
    """Validate workflow matrix compatibility without building runtime inputs."""

    resolved_substrate_map, resolved_motif_sequences = resolve_reference_bundle_inputs(
        substrate_map=substrate_map,
        motif_sequences=motif_sequences,
        reference_bundle=reference_bundle,
    )
    validated_matrix, _ = validate_workflow_matrix_inputs(
        phospho_matrix,
        resolved_substrate_map,
        site_sequences,
        require_site_sequences_for_prediction=resolved_motif_sequences is not None,
        context=context,
    )
    if resolved_motif_sequences is not None:
        KinaseMotifScorer.from_substrate_sequences(
            motif_sequences=resolved_motif_sequences,
            flank_size=flank_size,
        )
    return validated_matrix


def _copy_workflow_request_owned_state(
    request: KinaseWorkflowRequest,
    *,
    phospho_matrix: pd.DataFrame,
) -> KinaseWorkflowRequest:
    site_sequences = request.site_sequences
    if site_sequences is not None:
        site_sequences = site_sequences.copy(deep=True)

    return request.model_copy(
        update={
            "phospho_matrix": phospho_matrix,
            "site_sequences": site_sequences,
        }
    )
