"""Internal intensity-scale-state establishment for dataset builder execution."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

import numpy as np
import pandas as pd

from phospy.errors.transformations import (
    PhosPyTransformationError,
    TransformationStateEstablishmentError,
    TransformerExecutionError,
)
from phospy.errors.validation import TransformationValidationError
from phospy.science.transformations._authority import (
    _dataset_resolver_establishment_authority,
)
from phospy.science.transformations.contracts import Transformer
from phospy.science.transformations.models import (
    DeclaredIntensityScaleDiagnosticPolicy,
    IntensityScaleEstablishmentMode,
    IntensityScaleEvidenceLevel,
    IntensityScaleKind,
    IntensityScaleState,
    establish_intensity_scale_state,
)
from phospy.validation.transformations.state import IntensityScaleStateValidator


@dataclass(frozen=True, slots=True)
class ResolvedIntensityScale:
    """Quantitative matrices paired with an established intensity scale state."""

    phospho: pd.DataFrame
    total: pd.DataFrame | None
    intensity_scale_state: IntensityScaleState


class DatasetIntensityScaleResolver:
    """Resolve intensity scale state for a dataset build request.

    Supported establishment path:
    A configured transformer that establishes state from matrices.

    Public builder policy wires this resolver with identity pass-through by
    default. Pass-through transformers cannot establish scale from unknown
    inputs; they can only preserve explicitly declared input scale state.
    """

    def __init__(self, *, transformer: Transformer | None = None) -> None:
        self._transformer = transformer
        self._state_validator = IntensityScaleStateValidator()

    def run(
        self,
        *,
        phospho: pd.DataFrame,
        total: pd.DataFrame | None,
        expected_scale_kind: IntensityScaleKind | None = None,
        declared_input_scale_state: IntensityScaleState | None = None,
        declared_input_establishment_mode: IntensityScaleEstablishmentMode | None = (
            None
        ),
        declared_scale_diagnostic_policy: DeclaredIntensityScaleDiagnosticPolicy
        | str = DeclaredIntensityScaleDiagnosticPolicy.WARN,
        input_declaration_source: str | None = None,
        scale_establishment_parameters: Mapping[str, object] | None = None,
        scale_establishment_evidence_level: IntensityScaleEvidenceLevel | None = None,
        establishment_transformer_name: str | None = None,
        establishment_trace_id: str | None = None,
    ) -> ResolvedIntensityScale:
        diagnostic_policy = _normalize_declared_scale_diagnostic_policy(
            declared_scale_diagnostic_policy
        )
        if self._transformer is None:
            raise TransformationStateEstablishmentError(
                "unable to establish intensity scale state with confidence: "
                "no supported intensity-scale establisher is configured. "
                "Configure AnalysisReadyDatasetBuilder("
                "executor=DatasetBuildExecutor(transformer=...))."
            )
        preserves_input_scale_state = self._capability_enabled(
            self._transformer,
            "preserves_input_scale_state",
        )
        changes_numeric_values = self._capability_enabled(
            self._transformer,
            "changes_numeric_values",
            default=True,
        )
        identity_like_transformer = (
            preserves_input_scale_state and not changes_numeric_values
        )

        if declared_input_scale_state is not None and not preserves_input_scale_state:
            raise TransformationStateEstablishmentError(
                "unsupported identity state establishment: declared input intensity "
                "scale state can only be preserved by a transformer with "
                "preserves_input_scale_state=True. Use a preserving transformer for "
                "explicit input declarations, or remove the declaration and establish "
                "state via a scale-changing transformer."
            )

        try:
            transformed = self._transformer.run(phospho=phospho, total=total)
        except PhosPyTransformationError:
            raise
        except (TypeError, ValueError) as exc:
            raise TransformerExecutionError(
                "configured transformer failed while establishing intensity scale "
                "state from dataset matrices"
            ) from exc

        if (total is None) is not (transformed.total is None):
            raise TransformationStateEstablishmentError(
                "configured transformer changed total-matrix presence while "
                "establishing intensity scale state; this is unsupported. "
                "Use a transformer that preserves phospho/total matrix presence."
            )

        state = transformed.state
        if declared_input_scale_state is not None and preserves_input_scale_state:
            self._validate_state(
                declared_input_scale_state,
                has_total_matrix=transformed.total is not None,
                source="declared input intensity scale state",
            )
            state = declared_input_scale_state
        elif (
            preserves_input_scale_state
            and expected_scale_kind is IntensityScaleKind.LOG2
            and state.phospho.kind is not IntensityScaleKind.LOG2
        ):
            raise TransformationStateEstablishmentError(
                "missing intensity state evidence for expected 'log2': "
                "configured transformer preserves declared input state and cannot "
                "establish log2 from unknown/raw input without an explicit trusted "
                "declaration."
            )
        if declared_input_scale_state is None and identity_like_transformer:
            raise TransformationStateEstablishmentError(
                "unable to establish intensity scale state with confidence: "
                "pass-through/identity transformer cannot establish scientific "
                "input scale from values alone. Provide an explicit "
                "input_intensity_scale declaration or apply a scale-changing "
                "transformation."
            )
        if (
            expected_scale_kind is not None
            and state.phospho.kind is not expected_scale_kind
        ):
            source_label = (
                "declared input intensity state"
                if declared_input_scale_state is not None
                else "configured transformer output"
            )
            raise TransformationStateEstablishmentError(
                "mismatched expected intensity state: "
                f"expected '{expected_scale_kind.value}' but {source_label} is "
                f"'{state.phospho.kind.value}'"
            )

        self._validate_state(
            state,
            has_total_matrix=transformed.total is not None,
            source="configured transformer",
        )
        establishment_mode = self._resolve_establishment_mode(
            declared_input_scale_state=declared_input_scale_state,
            declared_input_establishment_mode=declared_input_establishment_mode,
            identity_like_transformer=identity_like_transformer,
        )
        diagnostic_warnings = ()
        if establishment_mode is IntensityScaleEstablishmentMode.DECLARED:
            diagnostic_warnings = _run_declared_scale_sanity_diagnostics(
                declared_scale_kind=state.phospho.kind,
                phospho=transformed.phospho,
                total=transformed.total,
            )
            if (
                diagnostic_warnings
                and diagnostic_policy is DeclaredIntensityScaleDiagnosticPolicy.ERROR
            ):
                raise _declared_scale_diagnostic_error(
                    declared_scale_kind=state.phospho.kind,
                    diagnostic_warnings=diagnostic_warnings,
                )
        transformer_source = (
            f"{self._transformer.__class__.__module__}."
            f"{self._transformer.__class__.__qualname__}"
        )
        transformer_name = self._resolve_transformer_name(
            establishment_mode=establishment_mode,
            default_transformer_name=transformer_source,
            explicit_transformer_name=establishment_transformer_name,
        )
        evidence_level = self._resolve_evidence_level(
            establishment_mode=establishment_mode,
            explicit_evidence_level=scale_establishment_evidence_level,
        )
        established_state = establish_intensity_scale_state(
            state,
            established_via=transformer_source,
            establishment_mode=establishment_mode,
            evidence_level=evidence_level,
            transformer_name=transformer_name,
            input_declaration_source=input_declaration_source,
            parameters=(
                {}
                if scale_establishment_parameters is None
                else scale_establishment_parameters
            ),
            trace_id=establishment_trace_id,
            diagnostic_warnings=diagnostic_warnings,
            _authority=_dataset_resolver_establishment_authority(),
        )
        return ResolvedIntensityScale(
            phospho=transformed.phospho,
            total=transformed.total,
            intensity_scale_state=established_state,
        )

    @staticmethod
    def _resolve_establishment_mode(
        *,
        declared_input_scale_state: IntensityScaleState | None,
        declared_input_establishment_mode: IntensityScaleEstablishmentMode | None,
        identity_like_transformer: bool,
    ) -> IntensityScaleEstablishmentMode:
        if declared_input_scale_state is not None:
            if declared_input_establishment_mode is None:
                return IntensityScaleEstablishmentMode.DECLARED
            if isinstance(
                declared_input_establishment_mode,
                IntensityScaleEstablishmentMode,
            ):
                return declared_input_establishment_mode
            return IntensityScaleEstablishmentMode(
                str(declared_input_establishment_mode)
            )
        if identity_like_transformer:
            return IntensityScaleEstablishmentMode.IDENTITY
        return IntensityScaleEstablishmentMode.DERIVED

    @staticmethod
    def _capability_enabled(
        transformer: Transformer,
        capability_name: str,
        *,
        default: bool = False,
    ) -> bool:
        value = getattr(transformer, capability_name, default)
        if value is None:
            return default
        return bool(value)

    @staticmethod
    def _resolve_transformer_name(
        *,
        establishment_mode: IntensityScaleEstablishmentMode,
        default_transformer_name: str,
        explicit_transformer_name: str | None,
    ) -> str | None:
        if explicit_transformer_name is not None:
            normalized = str(explicit_transformer_name).strip()
            if normalized:
                return normalized
        if establishment_mode is IntensityScaleEstablishmentMode.DECLARED:
            return None
        return default_transformer_name

    @staticmethod
    def _resolve_evidence_level(
        *,
        establishment_mode: IntensityScaleEstablishmentMode,
        explicit_evidence_level: IntensityScaleEvidenceLevel | None,
    ) -> IntensityScaleEvidenceLevel:
        if explicit_evidence_level is not None:
            if isinstance(explicit_evidence_level, IntensityScaleEvidenceLevel):
                return explicit_evidence_level
            return IntensityScaleEvidenceLevel(str(explicit_evidence_level))
        if establishment_mode is IntensityScaleEstablishmentMode.DECLARED:
            return IntensityScaleEvidenceLevel.DECLARED_BY_USER
        return IntensityScaleEvidenceLevel.UNKNOWN

    def _validate_state(
        self,
        state: IntensityScaleState,
        *,
        has_total_matrix: bool,
        source: str,
    ) -> None:
        try:
            self._state_validator.run(
                intensity_scale_state=state,
                has_total_matrix=has_total_matrix,
            )
        except TransformationValidationError as exc:
            raise TransformationStateEstablishmentError(
                f"{source} produced an invalid intensity scale state: {exc}"
            ) from exc


def _normalize_declared_scale_diagnostic_policy(
    policy: DeclaredIntensityScaleDiagnosticPolicy | str,
) -> DeclaredIntensityScaleDiagnosticPolicy:
    if isinstance(policy, DeclaredIntensityScaleDiagnosticPolicy):
        return policy
    try:
        return DeclaredIntensityScaleDiagnosticPolicy(str(policy).strip())
    except ValueError as exc:
        supported = ", ".join(
            item.value for item in DeclaredIntensityScaleDiagnosticPolicy
        )
        raise TransformationStateEstablishmentError(
            "unsupported declared intensity-scale diagnostic policy "
            f"{policy!r}; supported: {supported}"
        ) from exc


def _declared_scale_diagnostic_error(
    *,
    declared_scale_kind: IntensityScaleKind,
    diagnostic_warnings: tuple[str, ...],
) -> TransformationStateEstablishmentError:
    warning_count = len(diagnostic_warnings)
    warning_label = "warning" if warning_count == 1 else "warnings"
    return TransformationStateEstablishmentError(
        f"declared {declared_scale_kind.value} intensity scale produced "
        f"{warning_count} diagnostic {warning_label}; first warning: "
        f"{diagnostic_warnings[0]}. Suggested fix: correct input_intensity_scale, "
        "apply a PhosPy transform, or explicitly opt into suspicious declaration "
        "override with "
        "allow_suspicious_declared_input_intensity_scale=True on "
        "DatasetBuildRequest."
    )


def _run_declared_scale_sanity_diagnostics(
    *,
    declared_scale_kind: IntensityScaleKind,
    phospho: pd.DataFrame,
    total: pd.DataFrame | None,
) -> tuple[str, ...]:
    warnings: list[str] = []
    for matrix_label, matrix in _iter_numeric_matrices(phospho=phospho, total=total):
        finite_values = _finite_numeric_values(matrix)
        if finite_values.size == 0:
            continue
        median_value = float(np.median(finite_values))
        min_value = float(np.min(finite_values))
        max_value = float(np.max(finite_values))
        dynamic_range = max_value - min_value
        integer_like_fraction = _fraction_integer_like(finite_values)
        near_zero_signed_fraction = _fraction_near_zero_signed_values(finite_values)
        likely_log2_fraction = _fraction_within_log2_like_range(finite_values)
        if declared_scale_kind is IntensityScaleKind.LOG2:
            large_value_fraction = float(np.mean(finite_values >= 25.0))
            extreme_value_fraction = float(np.mean(finite_values >= 35.0))
            if (
                min_value >= 0.0
                and median_value >= 20.0
                and large_value_fraction >= 0.8
            ):
                warnings.append(
                    f"{matrix_label}: declared log2 scale is suspicious; "
                    f"median={median_value:.3f}, min={min_value:.3f}, "
                    f"fraction_ge_25={large_value_fraction:.3f}"
                )
            if max_value >= 35.0 or extreme_value_fraction >= 0.2:
                warnings.append(
                    f"{matrix_label}: declared log2 scale has very large positive values; "
                    f"max={max_value:.3f}, fraction_ge_35={extreme_value_fraction:.3f}"
                )
            if min_value <= -40.0 or dynamic_range >= 80.0:
                warnings.append(
                    f"{matrix_label}: declared log2 scale has highly suspicious range; "
                    f"min={min_value:.3f}, max={max_value:.3f}, range={dynamic_range:.3f}"
                )
            if median_value >= 1000.0 and integer_like_fraction >= 0.95:
                warnings.append(
                    f"{matrix_label}: declared log2 scale strongly resembles raw linear intensities; "
                    f"median={median_value:.3f}, integer_like_fraction={integer_like_fraction:.3f}"
                )
            continue
        if declared_scale_kind is IntensityScaleKind.LINEAR:
            negative_fraction = float(np.mean(finite_values < 0.0))
            has_positive = bool(np.any(finite_values > 0.0))
            has_negative = bool(np.any(finite_values < 0.0))
            if negative_fraction > 0.0:
                warnings.append(
                    f"{matrix_label}: declared linear scale contains negative values; "
                    f"fraction_negative={negative_fraction:.3f}"
                )
            if near_zero_signed_fraction >= 0.7 and has_positive and has_negative:
                warnings.append(
                    f"{matrix_label}: declared linear scale shows many small signed values "
                    f"consistent with log-ratio style data; "
                    f"fraction_abs_le_5_with_sign={near_zero_signed_fraction:.3f}"
                )
            if (
                likely_log2_fraction >= 0.95
                and integer_like_fraction <= 0.8
                and max_value <= 20.0
                and median_value <= 8.0
                and (has_negative or near_zero_signed_fraction >= 0.7)
            ):
                warnings.append(
                    f"{matrix_label}: declared linear scale strongly resembles "
                    f"already log-transformed data; "
                    f"fraction_in_log2_range={likely_log2_fraction:.3f}, "
                    f"integer_like_fraction={integer_like_fraction:.3f}"
                )
    return tuple(dict.fromkeys(warnings))


def _iter_numeric_matrices(
    *,
    phospho: pd.DataFrame,
    total: pd.DataFrame | None,
) -> tuple[tuple[str, pd.DataFrame], ...]:
    matrices: list[tuple[str, pd.DataFrame]] = [("phospho", phospho)]
    if total is not None:
        matrices.append(("total", total))
    return tuple(matrices)


def _finite_numeric_values(frame: pd.DataFrame) -> np.ndarray:
    numeric = frame.select_dtypes(include="number")
    if numeric.empty:
        return np.array([], dtype=float)
    values = numeric.to_numpy(dtype=float, copy=False).reshape(-1)
    finite_mask = np.isfinite(values)
    if not bool(finite_mask.any()):
        return np.array([], dtype=float)
    return values[finite_mask]


def _fraction_integer_like(values: np.ndarray) -> float:
    return float(np.mean(np.isclose(values, np.round(values), atol=1e-9)))


def _fraction_near_zero_signed_values(values: np.ndarray) -> float:
    near_zero_mask = np.abs(values) <= 5.0
    if not bool(near_zero_mask.any()):
        return 0.0
    signed_mask = values != 0.0
    return float(np.mean(near_zero_mask & signed_mask))


def _fraction_within_log2_like_range(values: np.ndarray) -> float:
    return float(np.mean((values >= -10.0) & (values <= 20.0)))
