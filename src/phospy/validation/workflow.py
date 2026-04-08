from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import pandas as pd
from pydantic import Field, ValidationError, field_validator, model_validator

from ..motifs import KinaseMotifScorer
from ..types import PredictionSvmMode
from ._models import PhospyRequestModel
from .errors import InputCompatibilityError, RequestValidationError
from .tables import SiteMatrixSchema


class KinaseWorkflowRequest(PhospyRequestModel):
    """Raw boundary request for native kinase workflow execution.

    The validated request model owns all non-scalar workflow inputs except the
    phosphosite matrix, which is copied later at the trusted workflow boundary.
    Mapping-backed sequence collections are normalized into detached concrete
    containers so caller mutation cannot change validated request state.
    """

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
        return _freeze_sequence_mapping(
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
        if value is None:
            return None
        if isinstance(value, pd.Series):
            normalized = value.copy(deep=True).astype(object)
            normalized.index = normalized.index.map(str)
            normalized[:] = normalized.map(str)
            return normalized
        if isinstance(value, Mapping):
            normalized = {
                str(site_id): str(sequence) for site_id, sequence in value.items()
            }
            return pd.Series(normalized, dtype=object)
        msg = (
            "site_sequences must be provided as a mapping keyed by phosphosite ID "
            "or as a pandas Series with an explicit phosphosite index; plain "
            "sequences are not supported"
        )
        raise ValueError(msg)

    @field_validator("motif_sequences", mode="before")
    @classmethod
    def validate_motif_sequences(
        cls,
        value: Mapping[str, Sequence[str]] | None,
    ) -> dict[str, tuple[str, ...]] | None:
        if value is None:
            return None
        return _freeze_sequence_mapping(
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
class ValidatedWorkflowRequest:
    """Trusted workflow inputs owned by the workflow boundary."""

    request: KinaseWorkflowRequest
    phospho_matrix: pd.DataFrame
    scoring_site_index: tuple[str, ...]
    motif_scorer: KinaseMotifScorer | None
    predictor_svm_mode: PredictionSvmMode


# Backward-compatible alias for older internal names.
ValidatedKinaseWorkflowInputs = ValidatedWorkflowRequest


def build_validated_workflow_request(
    request: KinaseWorkflowRequest,
    *,
    flank_size: int,
    default_svm_mode: PredictionSvmMode,
    context: str = "Kinase workflow inputs",
) -> ValidatedWorkflowRequest:
    """Build a trusted workflow request from validated raw options."""
    validated_matrix, scoring_site_index = _validate_workflow_matrix_inputs(
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
    return ValidatedWorkflowRequest(
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
    substrate_map: Mapping[str, Sequence[str]],
    site_sequences: Mapping[str, str] | pd.Series | None = None,
    motif_sequences: Mapping[str, Sequence[str]] | None = None,
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
) -> ValidatedWorkflowRequest:
    request = KinaseWorkflowRequest.validate_request(
        phospho_matrix=phospho_matrix,
        substrate_map=substrate_map,
        site_sequences=site_sequences,
        motif_sequences=motif_sequences,
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
    return build_validated_workflow_request(
        request,
        flank_size=flank_size,
        default_svm_mode=default_svm_mode,
        context=context,
    )


def build_workflow_request_inputs(
    request: KinaseWorkflowRequest,
    *,
    flank_size: int,
    default_svm_mode: PredictionSvmMode = "default",
    context: str = "Kinase workflow inputs",
) -> ValidatedWorkflowRequest:
    """Compatibility wrapper around :func:`build_validated_workflow_request`."""

    return build_validated_workflow_request(
        request,
        flank_size=flank_size,
        default_svm_mode=default_svm_mode,
        context=context,
    )


def validate_workflow_inputs(
    phospho_matrix: pd.DataFrame,
    substrate_map: Mapping[str, Sequence[str]],
    site_sequences: Mapping[str, str] | pd.Series | None,
    motif_sequences: Mapping[str, Sequence[str]] | None,
    *,
    flank_size: int = 7,
    context: str = "Kinase workflow inputs",
) -> pd.DataFrame:
    validated_matrix, _ = _validate_workflow_matrix_inputs(
        phospho_matrix,
        substrate_map,
        site_sequences,
        require_site_sequences_for_prediction=motif_sequences is not None,
        context=context,
    )
    if motif_sequences is not None:
        KinaseMotifScorer.from_substrate_sequences(
            motif_sequences=motif_sequences,
            flank_size=flank_size,
        )
    return validated_matrix


def _copy_workflow_request_owned_state(
    request: KinaseWorkflowRequest,
    *,
    phospho_matrix: pd.DataFrame,
) -> KinaseWorkflowRequest:
    """Return a trusted request model that fully owns workflow boundary state.

    ``KinaseWorkflowRequest`` already normalizes mapping-backed collections into
    detached concrete containers. This helper completes the ownership transfer by
    replacing the caller-managed phosphosite matrix and cloning the normalized
    site-sequence series so the trusted request bundle has no shared mutable
    boundary state left.
    """

    site_sequences = request.site_sequences
    if site_sequences is not None:
        site_sequences = site_sequences.copy(deep=True)

    return request.model_copy(
        update={
            "phospho_matrix": phospho_matrix,
            "site_sequences": site_sequences,
        }
    )


def _validate_workflow_matrix_inputs(
    phospho_matrix: pd.DataFrame,
    substrate_map: Mapping[str, Sequence[str]],
    site_sequences: Mapping[str, str] | pd.Series | None,
    *,
    require_site_sequences_for_prediction: bool,
    context: str,
) -> tuple[pd.DataFrame, tuple[str, ...]]:
    validated_matrix = SiteMatrixSchema.validate(
        phospho_matrix,
        context="phospho_matrix",
    )

    overlapping_sites = {
        site
        for sites in substrate_map.values()
        for site in sites
        if site in validated_matrix.index
    }
    if not overlapping_sites:
        msg = f"{context} contain no overlap between substrate_map and phospho_matrix"
        raise InputCompatibilityError(msg)

    scoring_site_index = tuple(str(site) for site in validated_matrix.index)

    if require_site_sequences_for_prediction:
        if site_sequences is None:
            msg = "site_sequences are required when motif_sequences are provided"
            raise InputCompatibilityError(msg)

        sequence_index = _extract_sequence_index(site_sequences)
        scoring_site_index = tuple(
            site for site in validated_matrix.index if str(site) in sequence_index
        )
        if not scoring_site_index:
            msg = (
                f"{context} contain no phosphosites with sequence coverage "
                "required for scoring and prediction"
            )
            raise InputCompatibilityError(msg)

    return validated_matrix, scoring_site_index


def _extract_sequence_index(
    site_sequences: Mapping[str, str] | pd.Series,
) -> set[str]:
    if isinstance(site_sequences, pd.Series):
        return {str(value) for value in site_sequences.index}
    if isinstance(site_sequences, Mapping):
        return {str(value) for value in site_sequences}
    msg = (
        "site_sequences must be provided as a mapping keyed by phosphosite ID "
        "or as a pandas Series with an explicit phosphosite index; plain "
        "sequences are not supported"
    )
    raise InputCompatibilityError(msg)


def _freeze_sequence_mapping(
    value: Mapping[str, Sequence[str]],
    *,
    field_name: str,
    empty_message: str,
) -> dict[str, tuple[str, ...]]:
    if not isinstance(value, Mapping):
        msg = f"{field_name} must be provided as a mapping"
        raise ValueError(msg)
    if not value:
        raise ValueError(empty_message)

    normalized: dict[str, tuple[str, ...]] = {}
    for key, raw_sequences in value.items():
        if isinstance(raw_sequences, (str, bytes)):
            msg = (
                f"{field_name}[{key!r}] must be a sequence of values, "
                "not a plain string"
            )
            raise ValueError(msg)
        normalized[str(key)] = tuple(str(item) for item in raw_sequences)
    return normalized


__all__ = [
    "KinaseWorkflowRequest",
    "ValidatedKinaseWorkflowInputs",
    "ValidatedWorkflowRequest",
    "build_validated_workflow_request",
    "build_workflow_request_inputs",
    "validate_workflow_inputs",
    "validate_workflow_request",
]
