"""Intensity-scale-state establishment for Dataset construction."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

import pandas as pd

from phospy.errors.transformations import (
    PhosPyTransformationError,
    TransformationStateEstablishmentError,
    TransformerExecutionError,
)
from phospy.errors.validation import TransformationValidationError
from phospy.quantitative._authority import (
    _dataset_resolver_establishment_authority,
)
from phospy.quantitative.contracts import Transformer
from phospy.quantitative.declared_scale_diagnostics import (
    _run_declared_scale_sanity_diagnostics,
)
from phospy.quantitative.models import (
    DeclaredIntensityScaleDiagnosticPolicy,
    IntensityScaleEstablishmentMode,
    IntensityScaleEvidenceLevel,
    IntensityScaleKind,
    IntensityScaleState,
    establish_intensity_scale_state,
)
from phospy.quantitative.state_coherence import (
    require_intensity_scale_state_coherence,
)
from phospy.quantitative.transformers import IdentityTransformer, Log2Transformer
from phospy.quantitative.transformers._identifiers import (
    resolve_transformer_identifier,
)

_BUILT_IN_TRANSFORMER_IDENTIFIERS = {
    IdentityTransformer: IdentityTransformer.stable_identifier,
    Log2Transformer: Log2Transformer.stable_identifier,
}


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
            if (
                diagnostic_warnings
                and state.phospho.kind is IntensityScaleKind.LOG2
                and diagnostic_policy is DeclaredIntensityScaleDiagnosticPolicy.WARN
            ):
                diagnostic_warnings = (
                    *diagnostic_warnings,
                    _dataset_build_override_scope_warning(),
                )
        transformer_source = resolve_transformer_identifier(
            self._transformer,
            built_in_identifiers=_BUILT_IN_TRANSFORMER_IDENTIFIERS,
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
            require_intensity_scale_state_coherence(
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


def _dataset_build_override_scope_warning() -> str:
    return (
        "The suspicious declared input intensity scale was accepted for dataset "
        "construction only. Downstream workflows may independently validate this "
        "declaration and may require their own explicit override. Differential "
        "analysis requires "
        "DifferentialAnalysisConfig.allow_suspicious_declared_input_scale=True."
    )
