from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any

import pandas as pd
from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    TypeAdapter,
    ValidationError,
    field_validator,
    model_validator,
)

from ..constants import ComparisonSpec
from ..core_processing import (
    CorePreprocessingConfig,
    resolve_core_preprocessing_config,
)
from ..dataset_schema import DatasetSchema
from ..motifs import KinaseMotifScorer
from ..signalome_site_ids import resolve_signalome_site_to_protein
from ..types import (
    PredictionSvmMode,
    PredictionTraceFormat,
    PredictionTraceLevel,
)
from .collections import (
    normalize_sequence_mapping,
    normalize_site_sequence_series,
    normalize_site_to_protein_mapping,
    normalize_string_sequence,
)
from .compatibility import (
    validate_pred_mat_overlap,
    validate_signalome_alignment,
    validate_workflow_matrix_inputs,
)
from .errors import (
    InputCompatibilityError,
    NoCandidateKinasesError,
    RequestValidationError,
)
from .identifiers import validate_existing_file_path
from .schemas import (
    PhosphoInputSchema,
    PredictionScoreMatrixSchema,
    PredMatSchema,
    SiteMatrixSchema,
    TotalInputSchema,
    normalize_pred_mat_input,
)

if TYPE_CHECKING:
    from ..dataset import PhosphoDataset
    from ..prediction.models import KinasePredictionResult, PredMatResult
    from ..scoring import KinaseScoringResult


_PREDICTION_SVM_MODE_ADAPTER = TypeAdapter(PredictionSvmMode)
_PREDICTION_TRACE_FORMAT_ADAPTER = TypeAdapter(PredictionTraceFormat)
_PREDICTION_TRACE_LEVEL_ADAPTER = TypeAdapter(PredictionTraceLevel)


class PhospyRequestModel(BaseModel):
    """Base request model with shared validation behaviour."""

    model_config = ConfigDict(arbitrary_types_allowed=True, frozen=True)


def validate_adapter_value(
    *,
    value: object,
    adapter: TypeAdapter[object],
    field_name: str,
    context: str,
) -> object:
    """Validate an adapter-backed field and raise a package-level request error."""

    try:
        return adapter.validate_python(value)
    except ValidationError as error:
        details = error.errors(include_url=False)
        message = str(details[0].get("msg", "Invalid value")) if details else str(error)
        raise RequestValidationError(f"{context}: {field_name}: {message}") from error


@dataclass(frozen=True, slots=True)
class ValidatedDatasetPaths:
    """Validated file-backed dataset input paths."""

    total_path: Path
    phospho_path: Path


@dataclass(slots=True)
class ValidatedDatasetInputs:
    """Trusted validated bundle for the public :class:`phospy.PhosphoDataset` boundary."""

    schema: DatasetSchema
    total_df: pd.DataFrame
    phospho_df: pd.DataFrame
    comparisons: tuple[ComparisonSpec, ...] | None = None


def validate_dataset_file_paths(
    total_path: str | Path,
    phospho_path: str | Path,
) -> ValidatedDatasetPaths:
    """Validate dataset file paths before table loading."""

    return ValidatedDatasetPaths(
        total_path=validate_existing_file_path(
            total_path,
            context="total input table path",
        ),
        phospho_path=validate_existing_file_path(
            phospho_path,
            context="phospho input table path",
        ),
    )


def validate_dataset_frames(
    *,
    total_df: pd.DataFrame,
    phospho_df: pd.DataFrame,
    schema: DatasetSchema,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Validate in-memory dataset tables against the configured schema."""

    validated_total = TotalInputSchema.validate(
        total_df,
        total_cols=schema.total_cols,
    )
    validated_phospho = PhosphoInputSchema.validate(
        phospho_df,
        phospho_cols=schema.phospho_cols,
    )
    return validated_total, validated_phospho


def validate_dataset_request(
    *,
    total_df: pd.DataFrame,
    phospho_df: pd.DataFrame,
    schema: DatasetSchema | None = None,
    comparisons: Sequence[ComparisonSpec] | None = None,
    context: str = "PhosphoDataset",
) -> ValidatedDatasetInputs:
    """Validate raw dataset inputs for the public dataset boundary."""

    resolved_schema = schema or DatasetSchema()
    validated_total, validated_phospho = validate_dataset_frames(
        total_df=total_df,
        phospho_df=phospho_df,
        schema=resolved_schema,
    )
    validated_comparisons = _validate_dataset_comparisons(
        schema=resolved_schema,
        comparisons=comparisons,
        context=context,
    )
    return ValidatedDatasetInputs(
        schema=resolved_schema,
        total_df=validated_total,
        phospho_df=validated_phospho,
        comparisons=validated_comparisons,
    )


def build_validated_dataset_inputs(
    *,
    schema: DatasetSchema,
    total_df: pd.DataFrame,
    phospho_df: pd.DataFrame,
    comparisons: Sequence[ComparisonSpec] | None = None,
    context: str = "PhosphoDataset",
) -> ValidatedDatasetInputs:
    """Build a validated dataset request from already validated frames."""

    return ValidatedDatasetInputs(
        schema=schema,
        total_df=total_df,
        phospho_df=phospho_df,
        comparisons=_validate_dataset_comparisons(
            schema=schema,
            comparisons=comparisons,
            context=context,
        ),
    )


class KinaseActivityRequest(PhospyRequestModel):
    """Raw boundary options for downstream kinase activity analysis."""

    threshold: float = Field(default=0.6, ge=0.0, le=1.0)
    min_substrates: int = Field(default=3, ge=1)
    top_n_substrates: int = Field(default=20, ge=1)

    @classmethod
    def validate_request(cls, **data: object) -> KinaseActivityRequest:
        try:
            return cls.model_validate(data)
        except ValidationError as error:
            raise RequestValidationError.from_pydantic(
                context="Invalid kinase activity request",
                error=error,
            ) from error


@dataclass(slots=True)
class ValidatedAnalysisRequest:
    """Trusted validated bundle for the public :class:`phospy.KinaseActivityAnalyzer` API."""

    request: KinaseActivityRequest
    pred_mat: pd.DataFrame
    phospho_matrix: pd.DataFrame

    @classmethod
    def from_trusted_inputs(
        cls,
        *,
        request: KinaseActivityRequest,
        pred_mat: pd.DataFrame,
        phospho_matrix: pd.DataFrame,
        pred_context: str = "pred_mat",
        matrix_context: str = "phospho_matrix",
        min_overlap: int = 1,
        min_fraction: float = 0.1,
    ) -> ValidatedAnalysisRequest:
        """Build a validated analysis request from already-owned validated matrices."""

        validate_pred_mat_overlap(
            pred_mat,
            phospho_matrix,
            pred_context=pred_context,
            matrix_context=matrix_context,
            min_overlap=min_overlap,
            min_fraction=min_fraction,
        )
        return cls(
            request=request,
            pred_mat=pred_mat,
            phospho_matrix=phospho_matrix,
        )


def validate_analysis_request(
    *,
    pred_mat: pd.DataFrame | PredMatResult,
    phospho_matrix: pd.DataFrame,
    threshold: float = 0.6,
    min_substrates: int = 3,
    top_n_substrates: int = 20,
    pred_context: str = "pred_mat",
    matrix_context: str = "phospho_matrix",
    min_overlap: int = 1,
    min_fraction: float = 0.1,
) -> ValidatedAnalysisRequest:
    """Validate raw analysis inputs and return a trusted analysis request."""

    request = KinaseActivityRequest.validate_request(
        threshold=threshold,
        min_substrates=min_substrates,
        top_n_substrates=top_n_substrates,
    )
    normalized_pred_mat = normalize_pred_mat_input(pred_mat)
    if normalized_pred_mat is None:
        msg = f"{pred_context} must be provided"
        raise RequestValidationError(msg)
    if normalized_pred_mat.shape[1] == 0:
        msg = (
            f"{pred_context} does not contain any kinase columns because no "
            "candidate kinases qualified for prediction. Regenerate predMat with "
            "less restrictive top, score_threshold, or inclusion settings."
        )
        raise NoCandidateKinasesError(msg)
    validated_pred_mat = PredMatSchema.validate(
        normalized_pred_mat,
        context=pred_context,
    )
    validated_matrix = SiteMatrixSchema.validate(phospho_matrix, context=matrix_context)
    return ValidatedAnalysisRequest.from_trusted_inputs(
        request=request,
        pred_mat=validated_pred_mat,
        phospho_matrix=validated_matrix,
        pred_context=pred_context,
        matrix_context=matrix_context,
        min_overlap=min_overlap,
        min_fraction=min_fraction,
    )


class PredictionRequest(PhospyRequestModel):
    """Validated boundary request for prediction execution."""

    combined_scores: pd.DataFrame
    ensemble_size: int = Field(ge=1)
    top: int = Field(ge=1)
    score_threshold: float = Field(ge=0.0, le=1.0)
    inclusion: int = Field(ge=1)
    n_iterations: int = Field(ge=1)
    random_state: int | None = None
    debug_kinases: tuple[str, ...] | None = None
    debug_top_n: int = Field(default=10, ge=1)
    svm_mode: PredictionSvmMode
    sampling_trace: Any | None = None
    trace_level: PredictionTraceLevel = "none"
    trace_sink_format: PredictionTraceFormat = "csv"
    trace_sink: Any | None = None

    @field_validator("combined_scores")
    @classmethod
    def validate_combined_scores(cls, value: pd.DataFrame) -> pd.DataFrame:
        return PredictionScoreMatrixSchema.validate(
            value,
            context="combined_scores",
        )

    @field_validator("debug_kinases", mode="before")
    @classmethod
    def normalize_debug_kinases(
        cls,
        value: Sequence[str] | None,
    ) -> tuple[str, ...] | None:
        if value is None:
            return None
        return tuple(value)

    @model_validator(mode="after")
    def validate_trace_sink_requirements(self) -> PredictionRequest:
        if self.trace_level != "full" and self.trace_sink is not None:
            msg = "trace_sink may only be provided when trace_level='full'"
            raise ValueError(msg)
        return self

    @classmethod
    def validate_request(
        cls,
        *,
        default_svm_mode: PredictionSvmMode,
        capture_debug_trace: bool = False,
        trace_sink_format: PredictionTraceFormat = "csv",
        **data: object,
    ) -> PredictionRequest:
        resolved_trace_level = validate_adapter_value(
            value=(
                "summary"
                if data.get("trace_level") is None and capture_debug_trace
                else data.get("trace_level") or "none"
            ),
            adapter=_PREDICTION_TRACE_LEVEL_ADAPTER,
            field_name="trace_level",
            context="Invalid prediction request",
        )
        resolved_trace_format = validate_adapter_value(
            value=trace_sink_format,
            adapter=_PREDICTION_TRACE_FORMAT_ADAPTER,
            field_name="trace_sink_format",
            context="Invalid prediction request",
        )
        resolved_svm_mode = validate_adapter_value(
            value=(
                default_svm_mode if data.get("svm_mode") is None else data["svm_mode"]
            ),
            adapter=_PREDICTION_SVM_MODE_ADAPTER,
            field_name="svm_mode",
            context="Invalid prediction request",
        )

        request_data = dict(data)
        request_data["svm_mode"] = resolved_svm_mode
        request_data["trace_level"] = resolved_trace_level
        request_data["trace_sink_format"] = resolved_trace_format

        try:
            return cls.model_validate(request_data)
        except ValidationError as error:
            raise RequestValidationError.from_pydantic(
                context="Invalid prediction request",
                error=error,
            ) from error


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
class ValidatedWorkflowRequest:
    """Trusted workflow inputs owned by the workflow boundary."""

    request: KinaseWorkflowRequest
    phospho_matrix: pd.DataFrame
    scoring_site_index: tuple[str, ...]
    motif_scorer: KinaseMotifScorer | None
    predictor_svm_mode: PredictionSvmMode


ValidatedKinaseWorkflowInputs = ValidatedWorkflowRequest


def build_validated_workflow_request(
    request: KinaseWorkflowRequest,
    *,
    flank_size: int,
    default_svm_mode: PredictionSvmMode,
    context: str = "Kinase workflow inputs",
) -> ValidatedWorkflowRequest:
    """Build a trusted workflow request from validated raw options."""

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
    validated_matrix, _ = validate_workflow_matrix_inputs(
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


class SignalomeRequest(PhospyRequestModel):
    """Raw boundary request for public signalome construction."""

    kinases_of_interest: tuple[str, ...]
    site_to_protein: dict[str, str] | None = None
    kinase_network_threshold: float = Field(default=0.9, ge=0.0, le=1.0)
    signalome_cutoff: float = Field(default=0.5, ge=0.0, le=1.0)
    module_count: int | None = Field(default=None, ge=1)
    min_kinase_module_share_percent: float = Field(default=1.0, ge=0.0)

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
class ValidatedSignalomeRequest:
    """Trusted aligned inputs for signalome construction."""

    request: SignalomeRequest
    scoring_matrix: pd.DataFrame
    pred_mat: pd.DataFrame
    expression_matrix: pd.DataFrame
    site_to_protein: pd.Series


def validate_signalome_request(
    *,
    scoring_result: KinaseScoringResult,
    prediction_result: KinasePredictionResult | PredMatResult,
    expression_matrix: pd.DataFrame,
    kinases_of_interest: Sequence[str],
    site_to_protein: Mapping[str, str] | None = None,
    kinase_network_threshold: float = 0.9,
    signalome_cutoff: float = 0.5,
    module_count: int | None = None,
    min_kinase_module_share_percent: float = 1.0,
) -> ValidatedSignalomeRequest:
    """Validate raw signalome inputs and return a trusted aligned bundle."""

    request = SignalomeRequest.validate_request(
        kinases_of_interest=kinases_of_interest,
        site_to_protein=site_to_protein,
        kinase_network_threshold=kinase_network_threshold,
        signalome_cutoff=signalome_cutoff,
        module_count=module_count,
        min_kinase_module_share_percent=min_kinase_module_share_percent,
    )

    return _build_validated_signalome_request(
        request=request,
        scoring_matrix=_resolve_scoring_matrix(scoring_result),
        pred_mat=_validate_prediction_result_pred_mat(prediction_result),
        expression_matrix=expression_matrix,
        scoring_context="scoring_result",
        pred_mat_context="prediction_result",
        expression_context="expression_matrix",
    )


class CorePipelineRequest(PhospyRequestModel):
    """Validated file-backed boundary request for pipeline construction."""

    total_path: Path
    phospho_path: Path
    pred_mat_path: Path | None = None
    phospho_encoding: str | None = None
    dataset_schema: DatasetSchema = Field(default_factory=DatasetSchema, alias="schema")
    comparisons: tuple[ComparisonSpec, ...] | None = None
    localization_threshold: float = Field(default=0.75, ge=0.0, le=1.0)
    min_observed: int = Field(default=4, ge=1)
    total_sentinel: float = 10.0
    phospho_sentinel: float = 12.0
    max_unmatched_fraction: float = Field(default=0.0, ge=0.0, le=1.0)
    kinase_activity_threshold: float = Field(default=0.6, ge=0.0, le=1.0)
    kinase_activity_min_substrates: int = Field(default=3, ge=1)
    kinase_activity_top_n_substrates: int = Field(default=20, ge=1)

    @field_validator("total_path", "phospho_path", "pred_mat_path")
    @classmethod
    def validate_existing_path(cls, value: Path | None) -> Path | None:
        if value is None:
            return None
        try:
            return validate_existing_file_path(value, context="core pipeline file path")
        except RequestValidationError as error:
            raise ValueError(str(error)) from error

    @model_validator(mode="after")
    def validate_comparisons(self) -> CorePipelineRequest:
        try:
            validated = self.dataset_schema.validate_comparisons(
                self.comparisons,
                context="Core pipeline request",
            )
        except (InputCompatibilityError, TypeError, ValueError) as error:
            raise ValueError(str(error)) from error
        object.__setattr__(self, "comparisons", validated)
        return self

    @classmethod
    def validate_request(cls, **data: object) -> CorePipelineRequest:
        try:
            return cls.model_validate(data)
        except ValidationError as error:
            raise RequestValidationError.from_pydantic(
                context="Invalid core pipeline request",
                error=error,
            ) from error


@dataclass(slots=True)
class ValidatedPipelineRequest:
    """Trusted pipeline inputs owned by the pipeline boundary."""

    dataset: PhosphoDataset
    pred_mat: pd.DataFrame | None
    preprocessing_config: CorePreprocessingConfig
    kinase_activity_request: KinaseActivityRequest | None


def build_pipeline_request(
    *,
    dataset: PhosphoDataset,
    validated_pred_mat: pd.DataFrame | None = None,
    preprocessing_config: CorePreprocessingConfig | None = None,
    localization_threshold: float = 0.75,
    min_observed: int = 4,
    max_unmatched_fraction: float = 0.0,
    total_sentinel: float = 10.0,
    phospho_sentinel: float = 12.0,
    kinase_activity_threshold: float = 0.6,
    kinase_activity_min_substrates: int = 3,
    kinase_activity_top_n_substrates: int = 20,
) -> ValidatedPipelineRequest:
    """Build a trusted pipeline request from already-owned inputs."""

    from ..dataset import PhosphoDataset

    if not isinstance(dataset, PhosphoDataset):
        msg = (
            "Invalid pipeline construction request: dataset must be a "
            "PhosphoDataset instance"
        )
        raise RequestValidationError(msg)

    try:
        resolved_config = resolve_core_preprocessing_config(
            config=preprocessing_config,
            localization_threshold=localization_threshold,
            min_observed=min_observed,
            max_unmatched_fraction=max_unmatched_fraction,
            total_sentinel=total_sentinel,
            phospho_sentinel=phospho_sentinel,
            context="Invalid pipeline construction request",
            config_param_name="preprocessing_config",
        )
    except (TypeError, ValueError) as error:
        raise RequestValidationError(str(error)) from error

    if validated_pred_mat is not None and not isinstance(
        validated_pred_mat, pd.DataFrame
    ):
        msg = (
            "Invalid pipeline construction request: pred_mat must be a "
            "pandas DataFrame when provided"
        )
        raise RequestValidationError(msg)

    kinase_activity_request = None
    if validated_pred_mat is not None:
        kinase_activity_request = KinaseActivityRequest.validate_request(
            threshold=kinase_activity_threshold,
            min_substrates=kinase_activity_min_substrates,
            top_n_substrates=kinase_activity_top_n_substrates,
        )

    return ValidatedPipelineRequest(
        dataset=dataset,
        pred_mat=validated_pred_mat,
        preprocessing_config=resolved_config,
        kinase_activity_request=kinase_activity_request,
    )


def validate_pipeline_construction_request(
    *,
    dataset: PhosphoDataset,
    pred_mat: pd.DataFrame | PredMatResult | None = None,
    preprocessing_config: CorePreprocessingConfig | None = None,
    localization_threshold: float = 0.75,
    min_observed: int = 4,
    max_unmatched_fraction: float = 0.0,
    total_sentinel: float = 10.0,
    phospho_sentinel: float = 12.0,
    kinase_activity_threshold: float = 0.6,
    kinase_activity_min_substrates: int = 3,
    kinase_activity_top_n_substrates: int = 20,
) -> ValidatedPipelineRequest:
    """Validate raw in-memory inputs for pipeline construction only."""

    normalized_pred_mat = normalize_pred_mat_input(pred_mat)
    validated_pred_mat = None
    if normalized_pred_mat is not None:
        validated_pred_mat = PredMatSchema.validate(
            normalized_pred_mat,
            context="pipeline pred_mat",
        )

    return build_pipeline_request(
        dataset=dataset,
        validated_pred_mat=validated_pred_mat,
        preprocessing_config=preprocessing_config,
        localization_threshold=localization_threshold,
        min_observed=min_observed,
        max_unmatched_fraction=max_unmatched_fraction,
        total_sentinel=total_sentinel,
        phospho_sentinel=phospho_sentinel,
        kinase_activity_threshold=kinase_activity_threshold,
        kinase_activity_min_substrates=kinase_activity_min_substrates,
        kinase_activity_top_n_substrates=kinase_activity_top_n_substrates,
    )


def validate_pipeline_runtime_compatibility(
    *,
    request: ValidatedPipelineRequest,
    site_matrix: pd.DataFrame,
) -> ValidatedAnalysisRequest | None:
    """Validate post-preprocessing overlap before kinase analysis runs."""

    if request.pred_mat is None or request.kinase_activity_request is None:
        return None

    try:
        return ValidatedAnalysisRequest.from_trusted_inputs(
            request=request.kinase_activity_request,
            pred_mat=request.pred_mat,
            phospho_matrix=site_matrix,
            pred_context="pipeline pred_mat",
            matrix_context="preprocessed site matrix",
        )
    except InputCompatibilityError as error:
        raise InputCompatibilityError(
            f"Pipeline runtime compatibility failed after preprocessing: {error}"
        ) from error


def _validate_dataset_comparisons(
    *,
    schema: DatasetSchema,
    comparisons: Sequence[ComparisonSpec] | None,
    context: str,
) -> tuple[ComparisonSpec, ...] | None:
    try:
        return schema.validate_comparisons(comparisons, context=context)
    except (InputCompatibilityError, TypeError, ValueError):
        raise


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


def _build_validated_signalome_request(
    *,
    request: SignalomeRequest,
    scoring_matrix: pd.DataFrame,
    pred_mat: pd.DataFrame,
    expression_matrix: pd.DataFrame,
    scoring_context: str,
    pred_mat_context: str,
    expression_context: str,
) -> ValidatedSignalomeRequest:
    (
        validated_scoring_matrix,
        validated_pred_mat,
        validated_expression_matrix,
        common_sites,
    ) = validate_signalome_alignment(
        scoring_matrix=scoring_matrix,
        pred_mat=pred_mat,
        expression_matrix=expression_matrix,
        kinases_of_interest=request.kinases_of_interest,
        module_count=request.module_count,
        scoring_context=scoring_context,
        pred_mat_context=pred_mat_context,
        expression_context=expression_context,
    )

    validated_site_to_protein = _validate_signalome_site_grouping(
        site_ids=common_sites,
        site_to_protein=request.site_to_protein,
    )

    return ValidatedSignalomeRequest(
        request=request,
        scoring_matrix=validated_scoring_matrix,
        pred_mat=validated_pred_mat,
        expression_matrix=validated_expression_matrix,
        site_to_protein=validated_site_to_protein,
    )


def _validate_signalome_site_grouping(
    *,
    site_ids: Sequence[str],
    site_to_protein: Mapping[str, str] | None,
) -> pd.Series:
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


def _resolve_scoring_matrix(scoring_result: KinaseScoringResult) -> pd.DataFrame:
    from ..scoring import KinaseScoringResult

    if not isinstance(scoring_result, KinaseScoringResult):
        msg = "scoring_result must be a KinaseScoringResult"
        raise RequestValidationError(msg)

    if scoring_result.combined_scores is not None:
        return scoring_result.combined_scores
    return scoring_result.profile_scores


def _validate_prediction_result_pred_mat(
    prediction_result: KinasePredictionResult | PredMatResult,
) -> pd.DataFrame:
    pred_mat = _resolve_pred_mat(prediction_result)
    if pred_mat.shape[1] == 0:
        msg = (
            "prediction_result does not contain any kinase columns because no "
            "candidate kinases qualified for prediction. Regenerate predMat with "
            "less restrictive top, score_threshold, or inclusion settings."
        )
        raise NoCandidateKinasesError(msg)
    return PredMatSchema.validate(pred_mat, context="prediction_result")


def _resolve_pred_mat(
    prediction_result: KinasePredictionResult | PredMatResult,
) -> pd.DataFrame:
    from ..prediction.models import KinasePredictionResult, PredMatResult

    if isinstance(prediction_result, KinasePredictionResult):
        return prediction_result.pred_mat_result.to_frame(copy=False)
    if isinstance(prediction_result, PredMatResult):
        return prediction_result.to_frame(copy=False)
    msg = "prediction_result must be a KinasePredictionResult or PredMatResult"
    raise RequestValidationError(msg)


__all__ = [
    "CorePipelineRequest",
    "KinaseActivityRequest",
    "KinaseWorkflowRequest",
    "PhospyRequestModel",
    "PredictionRequest",
    "SignalomeRequest",
    "ValidatedAnalysisRequest",
    "ValidatedDatasetInputs",
    "ValidatedDatasetPaths",
    "ValidatedKinaseWorkflowInputs",
    "ValidatedPipelineRequest",
    "ValidatedSignalomeRequest",
    "ValidatedWorkflowRequest",
    "build_pipeline_request",
    "build_validated_dataset_inputs",
    "build_validated_workflow_request",
    "build_workflow_request_inputs",
    "validate_adapter_value",
    "validate_analysis_request",
    "validate_dataset_file_paths",
    "validate_dataset_frames",
    "validate_dataset_request",
    "validate_pipeline_construction_request",
    "validate_pipeline_runtime_compatibility",
    "validate_signalome_request",
    "validate_workflow_inputs",
    "validate_workflow_request",
    "validate_existing_file_path",
]
