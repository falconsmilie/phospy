from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import pandas as pd
from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    ValidationError,
    field_validator,
    model_validator,
)

from ...errors import RequestValidationError
from ...internal.defaults import (
    DEFAULT_MOTIF_FLANK_SIZE,
    DEFAULT_PREDICTION_ENSEMBLE_SIZE,
    DEFAULT_PREDICTION_INCLUSION,
    DEFAULT_PREDICTION_N_ITERATIONS,
    DEFAULT_PREDICTION_SCORE_THRESHOLD,
    DEFAULT_PREDICTION_TOP,
)
from ...internal.types import PREDICTION_SVM_MODE_DEFAULT, PredictionSvmMode
from ...prediction.motif_scoring import KinaseMotifScorer
from ...prediction.profiles import KinaseProfilePolicy
from ...references import ReferenceBundle
from ..compatibility import validate_workflow_matrix_inputs
from ..domain import resolve_reference_bundle_inputs
from ..values.collections import (
    normalize_sequence_mapping,
    normalize_site_sequence_series,
)


class KinaseWorkflowRequest(BaseModel):
    """Raw boundary request for native kinase workflow execution."""

    model_config = ConfigDict(arbitrary_types_allowed=True, frozen=True)

    phospho_matrix: pd.DataFrame
    substrate_map: dict[str, tuple[str, ...]]
    site_sequences: pd.Series | None = None
    motif_sequences: dict[str, tuple[str, ...]] | None = None
    min_substrates: int = Field(default=1, ge=1)
    min_motif_size: int = Field(default=1, ge=1)
    allow_profile_only_fallback: bool = False
    ensemble_size: int = Field(default=DEFAULT_PREDICTION_ENSEMBLE_SIZE, ge=1)
    top: int = Field(default=DEFAULT_PREDICTION_TOP, ge=1)
    score_threshold: float = Field(
        default=DEFAULT_PREDICTION_SCORE_THRESHOLD, ge=0.0, le=1.0
    )
    inclusion: int = Field(default=DEFAULT_PREDICTION_INCLUSION, ge=1)
    n_iterations: int = Field(default=DEFAULT_PREDICTION_N_ITERATIONS, ge=1)
    random_state: int | None = None
    svm_mode: PredictionSvmMode | None = None
    profile_policy: KinaseProfilePolicy = Field(default_factory=KinaseProfilePolicy)

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

    @field_validator("profile_policy", mode="before")
    @classmethod
    def validate_profile_policy(cls, value: object) -> KinaseProfilePolicy:
        return KinaseProfilePolicy.from_value(value)

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

    phospho_matrix: pd.DataFrame
    substrate_map: dict[str, tuple[str, ...]]
    site_sequences: pd.Series | None
    scoring_site_index: tuple[str, ...]
    motif_scorer: KinaseMotifScorer | None
    min_substrates: int
    min_motif_size: int
    allow_profile_only_fallback: bool
    ensemble_size: int
    top: int
    score_threshold: float
    inclusion: int
    n_iterations: int
    random_state: int | None
    svm_mode: PredictionSvmMode | None
    predictor_svm_mode: PredictionSvmMode
    profile_policy: KinaseProfilePolicy


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
    owned_site_sequences = request.site_sequences
    if owned_site_sequences is not None:
        owned_site_sequences = owned_site_sequences.copy(deep=True)
    motif_scorer = (
        None
        if request.motif_sequences is None
        else KinaseMotifScorer.from_substrate_sequences(
            motif_sequences=request.motif_sequences,
            flank_size=flank_size,
        )
    )
    return WorkflowInputs(
        phospho_matrix=validated_matrix,
        substrate_map=dict(request.substrate_map),
        site_sequences=owned_site_sequences,
        scoring_site_index=scoring_site_index,
        motif_scorer=motif_scorer,
        min_substrates=request.min_substrates,
        min_motif_size=request.min_motif_size,
        allow_profile_only_fallback=request.allow_profile_only_fallback,
        ensemble_size=request.ensemble_size,
        top=request.top,
        score_threshold=request.score_threshold,
        inclusion=request.inclusion,
        n_iterations=request.n_iterations,
        random_state=request.random_state,
        svm_mode=request.svm_mode,
        predictor_svm_mode=(
            default_svm_mode if request.svm_mode is None else request.svm_mode
        ),
        profile_policy=request.profile_policy,
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
    ensemble_size: int = DEFAULT_PREDICTION_ENSEMBLE_SIZE,
    top: int = DEFAULT_PREDICTION_TOP,
    score_threshold: float = DEFAULT_PREDICTION_SCORE_THRESHOLD,
    inclusion: int = DEFAULT_PREDICTION_INCLUSION,
    n_iterations: int = DEFAULT_PREDICTION_N_ITERATIONS,
    random_state: int | None = None,
    svm_mode: PredictionSvmMode | None = None,
    profile_policy: KinaseProfilePolicy | None = None,
    flank_size: int = DEFAULT_MOTIF_FLANK_SIZE,
    default_svm_mode: PredictionSvmMode = PREDICTION_SVM_MODE_DEFAULT,
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
        profile_policy=profile_policy,
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
    flank_size: int = DEFAULT_MOTIF_FLANK_SIZE,
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
