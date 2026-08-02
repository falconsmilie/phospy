"""Processing-state data models."""

from __future__ import annotations

from dataclasses import dataclass, field

from phospy.errors.input import DatasetProcessingStateError, PhosPyInputError
from phospy.science.datasets.preprocessing.imputation_scale_policy import (
    IMPUTATION_OPERATION_ORDERS,
)
from phospy.science.datasets.preprocessing.policy_models import (
    MissingDataPolicy,
    NormalisationPolicy,
    SiteSequenceConflictPolicy,
    SiteSequenceResolutionMode,
    TotalProteinCorrectionPolicy,
)
from phospy.science.transformations.models import (
    IntensityScaleState,
    QuantitativeMeaning,
)

from .missing_data import MissingDataDiagnostics
from .total_protein import TotalProteinCorrectionDiagnostics


@dataclass(frozen=True, slots=True)
class MissingDataState:
    """Missing-data policy state at the analysis-ready dataset boundary."""

    policy: MissingDataPolicy
    min_observed_values: int | None
    complete_matrix: bool
    imputed: bool
    diagnostics: MissingDataDiagnostics | None = None
    has_missing_values: bool | None = None
    missing_value_count: int | None = None
    imputation_input_scale: str | None = None
    imputation_operation_order: str | None = None

    def __post_init__(self) -> None:
        try:
            policy = MissingDataPolicy.parse(
                self.policy,
                field_name="dataset processing state missing_data.policy",
            )
        except PhosPyInputError as exc:
            raise DatasetProcessingStateError(str(exc)) from exc
        object.__setattr__(self, "policy", policy)
        complete_matrix = _require_bool(
            self.complete_matrix,
            field_name="dataset processing state missing_data.complete_matrix",
        )
        imputed = _require_bool(
            self.imputed,
            field_name="dataset processing state missing_data.imputed",
        )
        min_observed_values = require_optional_non_negative_int(
            self.min_observed_values,
            field_name="dataset processing state missing_data.min_observed_values",
        )
        has_missing_values = require_optional_bool(
            self.has_missing_values,
            field_name="dataset processing state missing_data.has_missing_values",
        )
        missing_value_count = require_optional_non_negative_int(
            self.missing_value_count,
            field_name="dataset processing state missing_data.missing_value_count",
        )
        imputation_input_scale = _require_optional_imputation_input_scale_state(
            self.imputation_input_scale,
        )
        imputation_operation_order = _require_optional_imputation_operation_order_state(
            self.imputation_operation_order,
        )
        diagnostics = self.diagnostics
        if diagnostics is not None and not isinstance(
            diagnostics, MissingDataDiagnostics
        ):
            diagnostics = MissingDataDiagnostics.from_payload(
                diagnostics,
                field_name="dataset processing state missing_data.diagnostics",
            )

        if diagnostics is not None:
            diagnostic_missing_count = _diagnostic_int(
                diagnostics,
                key="output_missing_cell_count",
                field_name=(
                    "dataset processing state missing_data.diagnostics."
                    "output_missing_cell_count"
                ),
            )
            if missing_value_count is None:
                missing_value_count = diagnostic_missing_count
            elif missing_value_count != diagnostic_missing_count:
                raise DatasetProcessingStateError(
                    "dataset processing state missing_data.missing_value_count "
                    "must match missing_data.diagnostics.output_missing_cell_count"
                )
            diagnostic_has_missing_values = diagnostic_missing_count > 0
            if has_missing_values is None:
                has_missing_values = diagnostic_has_missing_values
            elif has_missing_values != diagnostic_has_missing_values:
                raise DatasetProcessingStateError(
                    "dataset processing state missing_data.has_missing_values "
                    "must match missing_data.diagnostics.output_missing_cell_count"
                )
            imputation_input_scale = _resolve_diagnostic_optional_str(
                diagnostics=diagnostics,
                key="imputation_input_scale",
                current=imputation_input_scale,
                field_name=(
                    "dataset processing state missing_data.imputation_input_scale"
                ),
            )
            imputation_operation_order = _resolve_diagnostic_optional_str(
                diagnostics=diagnostics,
                key="imputation_operation_order",
                current=imputation_operation_order,
                field_name=(
                    "dataset processing state missing_data.imputation_operation_order"
                ),
            )

        if complete_matrix and missing_value_count is None:
            missing_value_count = 0
        if has_missing_values is None and missing_value_count is not None:
            has_missing_values = missing_value_count > 0
        if has_missing_values is None and complete_matrix:
            has_missing_values = False

        if complete_matrix and (missing_value_count or 0) > 0:
            raise DatasetProcessingStateError(
                "dataset processing state missing_data.complete_matrix cannot be "
                "True when missing_value_count is greater than zero"
            )
        if complete_matrix and has_missing_values:
            raise DatasetProcessingStateError(
                "dataset processing state missing_data.complete_matrix cannot be "
                "True when has_missing_values is True"
            )
        if has_missing_values is False and (missing_value_count or 0) > 0:
            raise DatasetProcessingStateError(
                "dataset processing state missing_data.has_missing_values cannot be "
                "False when missing_value_count is greater than zero"
            )

        imputed_cell_count = (
            None
            if diagnostics is None
            else _diagnostic_int(
                diagnostics,
                key="imputed_cell_count",
                field_name=(
                    "dataset processing state missing_data.diagnostics."
                    "imputed_cell_count"
                ),
            )
        )
        if imputed:
            if self.policy is MissingDataPolicy.FORBID:
                raise DatasetProcessingStateError(
                    "dataset processing state missing_data.imputed cannot be True "
                    "when missing_data.policy is 'forbid'"
                )
            if diagnostics is None:
                raise DatasetProcessingStateError(
                    "dataset processing state missing_data.imputed requires "
                    "missing_data.diagnostics with imputation provenance"
                )
            imputation_method_id = diagnostics.get("imputation_method_id")
            if not isinstance(imputation_method_id, str) or not imputation_method_id:
                raise DatasetProcessingStateError(
                    "dataset processing state missing_data.imputed requires "
                    "missing_data.diagnostics.imputation_method_id"
                )
            if imputation_method_id == "forbid":
                raise DatasetProcessingStateError(
                    "dataset processing state missing_data.imputed requires an "
                    "imputation method, not 'forbid'"
                )
            if imputed_cell_count is None or imputed_cell_count <= 0:
                raise DatasetProcessingStateError(
                    "dataset processing state missing_data.imputed requires "
                    "missing_data.diagnostics.imputed_cell_count greater than zero"
                )
        elif imputed_cell_count is not None and imputed_cell_count > 0:
            raise DatasetProcessingStateError(
                "dataset processing state missing_data.imputed cannot be False "
                "when diagnostics record imputed cells"
            )

        object.__setattr__(self, "complete_matrix", complete_matrix)
        object.__setattr__(self, "imputed", imputed)
        object.__setattr__(self, "min_observed_values", min_observed_values)
        object.__setattr__(self, "diagnostics", diagnostics)
        object.__setattr__(self, "has_missing_values", has_missing_values)
        object.__setattr__(self, "missing_value_count", missing_value_count)
        object.__setattr__(self, "imputation_input_scale", imputation_input_scale)
        object.__setattr__(
            self,
            "imputation_operation_order",
            imputation_operation_order,
        )


@dataclass(frozen=True, slots=True)
class NormalisationState:
    """Normalisation policy state at the analysis-ready dataset boundary."""

    policy: str

    def __post_init__(self) -> None:
        try:
            policy = NormalisationPolicy.parse(
                self.policy,
                field_name="dataset processing state normalisation.policy",
            )
        except PhosPyInputError as exc:
            raise DatasetProcessingStateError(str(exc)) from exc
        object.__setattr__(self, "policy", policy)


@dataclass(frozen=True, slots=True)
class SiteSequenceResolutionRowDiagnostic:
    """Durable per-row site-sequence resolution diagnostic record."""

    row_index: int
    row_id: str
    site_id: str | None
    status: str
    existing_site_sequence: str | None
    fasta_site_sequence: str | None
    resolved_site_sequence: str | None
    action: str
    reason: str | None
    conflict_policy: str | None
    resolver_version: str | None
    fasta_source_path: str | None
    fasta_sha256: str | None


@dataclass(frozen=True, slots=True)
class SiteSequenceResolutionState:
    """Dataset site-sequence FASTA-resolution state at preprocessing boundary."""

    configured: bool
    mode: str | None
    flank_size: int | None
    fasta_source_path: str | None
    fasta_source_label: str | None
    fasta_sha256: str | None
    resolver_version: str | None
    resolved_site_count: int
    unresolved_site_count: int
    unresolved_counts_by_reason: dict[str, int]
    filled_missing_count: int
    replaced_existing_count: int
    preserved_existing_count: int
    existing_sequence_conflict_count: int
    conflict_policy: str | None = None
    row_diagnostics: tuple[SiteSequenceResolutionRowDiagnostic, ...] = ()

    def __post_init__(self) -> None:
        configured = _require_bool(
            self.configured,
            field_name="dataset processing state site_sequence_resolution.configured",
        )
        mode = None
        if self.mode is not None:
            try:
                mode = SiteSequenceResolutionMode.parse(
                    self.mode,
                    field_name="dataset processing state site_sequence_resolution.mode",
                ).value
            except PhosPyInputError as exc:
                raise DatasetProcessingStateError(str(exc)) from exc
        flank_size = require_optional_non_negative_int(
            self.flank_size,
            field_name="dataset processing state site_sequence_resolution.flank_size",
        )
        fasta_source_path = require_optional_str(
            self.fasta_source_path,
            field_name=(
                "dataset processing state site_sequence_resolution.fasta_source_path"
            ),
        )
        fasta_source_label = require_optional_str(
            self.fasta_source_label,
            field_name=(
                "dataset processing state site_sequence_resolution.fasta_source_label"
            ),
        )
        fasta_sha256 = require_optional_str(
            self.fasta_sha256,
            field_name="dataset processing state site_sequence_resolution.fasta_sha256",
        )
        resolver_version = require_optional_str(
            self.resolver_version,
            field_name=(
                "dataset processing state site_sequence_resolution.resolver_version"
            ),
        )
        resolved_site_count = _require_non_negative_int(
            self.resolved_site_count,
            field_name=(
                "dataset processing state site_sequence_resolution.resolved_site_count"
            ),
        )
        unresolved_site_count = _require_non_negative_int(
            self.unresolved_site_count,
            field_name=(
                "dataset processing state site_sequence_resolution."
                "unresolved_site_count"
            ),
        )
        filled_missing_count = _require_non_negative_int(
            self.filled_missing_count,
            field_name=(
                "dataset processing state site_sequence_resolution.filled_missing_count"
            ),
        )
        replaced_existing_count = _require_non_negative_int(
            self.replaced_existing_count,
            field_name=(
                "dataset processing state site_sequence_resolution."
                "replaced_existing_count"
            ),
        )
        preserved_existing_count = _require_non_negative_int(
            self.preserved_existing_count,
            field_name=(
                "dataset processing state site_sequence_resolution."
                "preserved_existing_count"
            ),
        )
        existing_sequence_conflict_count = _require_non_negative_int(
            self.existing_sequence_conflict_count,
            field_name=(
                "dataset processing state site_sequence_resolution."
                "existing_sequence_conflict_count"
            ),
        )
        unresolved_counts_by_reason = _require_string_non_negative_int_mapping(
            self.unresolved_counts_by_reason,
            field_name=(
                "dataset processing state site_sequence_resolution."
                "unresolved_counts_by_reason"
            ),
        )
        row_diagnostics = _require_site_sequence_row_diagnostics(self.row_diagnostics)
        row_conflict_count = _count_site_sequence_row_conflicts(row_diagnostics)
        if row_conflict_count > 0 and existing_sequence_conflict_count == 0:
            raise DatasetProcessingStateError(
                "dataset processing state site_sequence_resolution records "
                "conflict row diagnostics but existing_sequence_conflict_count is zero"
            )
        if (
            existing_sequence_conflict_count > 0
            and row_conflict_count > existing_sequence_conflict_count
        ):
            raise DatasetProcessingStateError(
                "dataset processing state site_sequence_resolution."
                "existing_sequence_conflict_count must cover conflict row diagnostics"
            )
        if existing_sequence_conflict_count > 0 and not row_diagnostics:
            raise DatasetProcessingStateError(
                "dataset processing state site_sequence_resolution says conflicts "
                "exist but has no row diagnostics"
            )

        conflict_policy = None
        if self.conflict_policy is not None:
            try:
                conflict_policy = SiteSequenceConflictPolicy.parse(
                    self.conflict_policy,
                    field_name=(
                        "dataset processing state site_sequence_resolution."
                        "conflict_policy"
                    ),
                ).value
            except PhosPyInputError as exc:
                raise DatasetProcessingStateError(str(exc)) from exc
        if configured and conflict_policy is None:
            raise DatasetProcessingStateError(
                "dataset processing state site_sequence_resolution.conflict_policy "
                "is required when site-sequence resolution is configured"
            )
        if existing_sequence_conflict_count > 0 and conflict_policy is None:
            raise DatasetProcessingStateError(
                "dataset processing state site_sequence_resolution.conflict_policy "
                "is required when conflicts are recorded"
            )
        if configured and not (fasta_source_path or fasta_source_label or fasta_sha256):
            raise DatasetProcessingStateError(
                "dataset processing state site_sequence_resolution requires a "
                "non-empty FASTA/source identifier when configured"
            )
        if (
            configured
            and (resolved_site_count > 0 or filled_missing_count > 0)
            and not (fasta_sha256 or resolver_version)
        ):
            raise DatasetProcessingStateError(
                "dataset processing state site_sequence_resolution records "
                "reference-resolved sequences but no FASTA hash or resolver version"
            )
        if (
            conflict_policy == SiteSequenceConflictPolicy.ERROR.value
            and existing_sequence_conflict_count > 0
        ):
            raise DatasetProcessingStateError(
                "dataset processing state site_sequence_resolution.conflict_policy "
                "cannot be 'error' when unresolved conflicts are recorded as tolerated"
            )
        if conflict_policy == SiteSequenceConflictPolicy.ERROR.value:
            tolerated_conflicts = [
                item
                for item in row_diagnostics
                if _site_sequence_row_records_conflict(item)
                and item.action in {"preserve_existing", "replace_existing"}
            ]
            if tolerated_conflicts:
                raise DatasetProcessingStateError(
                    "dataset processing state site_sequence_resolution.conflict_policy "
                    "cannot be 'error' when conflict diagnostics record tolerated "
                    "actions"
                )

        object.__setattr__(self, "configured", configured)
        object.__setattr__(self, "mode", mode)
        object.__setattr__(self, "flank_size", flank_size)
        object.__setattr__(self, "fasta_source_path", fasta_source_path)
        object.__setattr__(self, "fasta_source_label", fasta_source_label)
        object.__setattr__(self, "fasta_sha256", fasta_sha256)
        object.__setattr__(self, "resolver_version", resolver_version)
        object.__setattr__(self, "resolved_site_count", resolved_site_count)
        object.__setattr__(self, "unresolved_site_count", unresolved_site_count)
        object.__setattr__(
            self, "unresolved_counts_by_reason", unresolved_counts_by_reason
        )
        object.__setattr__(self, "filled_missing_count", filled_missing_count)
        object.__setattr__(self, "replaced_existing_count", replaced_existing_count)
        object.__setattr__(self, "preserved_existing_count", preserved_existing_count)
        object.__setattr__(
            self,
            "existing_sequence_conflict_count",
            existing_sequence_conflict_count,
        )
        object.__setattr__(self, "conflict_policy", conflict_policy)
        object.__setattr__(self, "row_diagnostics", row_diagnostics)


@dataclass(frozen=True, slots=True)
class TotalProteinCorrectionState:
    """Total-protein correction state at the analysis-ready dataset boundary."""

    policy: TotalProteinCorrectionPolicy
    applied: bool
    formula: str | None = None
    requires_log_scale: bool | None = False
    input_scale: str | None = None
    output_scale: str | None = None
    quantitative_meaning: QuantitativeMeaning | None = None
    diagnostics: TotalProteinCorrectionDiagnostics | None = None

    def __post_init__(self) -> None:
        try:
            policy = TotalProteinCorrectionPolicy.parse(
                self.policy,
                field_name="dataset processing state total_protein_correction.policy",
            )
        except PhosPyInputError as exc:
            raise DatasetProcessingStateError(str(exc)) from exc
        object.__setattr__(self, "policy", policy)
        applied = _require_bool(
            self.applied,
            field_name="dataset processing state total_protein_correction.applied",
        )
        formula = require_optional_str(
            self.formula,
            field_name="dataset processing state total_protein_correction.formula",
        )
        requires_log_scale = require_optional_bool(
            self.requires_log_scale,
            field_name=(
                "dataset processing state total_protein_correction.requires_log_scale"
            ),
        )
        input_scale = require_optional_str(
            self.input_scale,
            field_name="dataset processing state total_protein_correction.input_scale",
        )
        output_scale = require_optional_str(
            self.output_scale,
            field_name="dataset processing state total_protein_correction.output_scale",
        )
        quantitative_meaning = self.quantitative_meaning
        if quantitative_meaning is not None and not isinstance(
            quantitative_meaning, QuantitativeMeaning
        ):
            try:
                quantitative_meaning = QuantitativeMeaning(str(quantitative_meaning))
            except ValueError as exc:
                supported = ", ".join(member.value for member in QuantitativeMeaning)
                raise PhosPyInputError(
                    "dataset processing state total_protein_correction."
                    "quantitative_meaning must be one of: "
                    f"{supported}"
                ) from exc
            object.__setattr__(self, "quantitative_meaning", quantitative_meaning)
        diagnostics = self.diagnostics
        if diagnostics is not None and not isinstance(
            diagnostics, TotalProteinCorrectionDiagnostics
        ):
            diagnostics = TotalProteinCorrectionDiagnostics.from_payload(
                diagnostics,
                field_name=(
                    "dataset processing state total_protein_correction.diagnostics"
                ),
            )
        if applied:
            if self.policy is TotalProteinCorrectionPolicy.NONE:
                raise DatasetProcessingStateError(
                    "dataset processing state total_protein_correction.applied "
                    "cannot be True when policy is 'none'"
                )
            if formula is None:
                raise DatasetProcessingStateError(
                    "dataset processing state total_protein_correction.applied "
                    "requires a correction formula/method"
                )
            if diagnostics is None:
                raise DatasetProcessingStateError(
                    "dataset processing state total_protein_correction.applied "
                    "requires correction diagnostics/provenance"
                )
            _require_total_correction_diagnostics_policy(
                diagnostics,
                expected_policy=self.policy,
                applied=applied,
            )
        else:
            if self.policy is not TotalProteinCorrectionPolicy.NONE:
                raise DatasetProcessingStateError(
                    "dataset processing state total_protein_correction.applied "
                    "cannot be False when policy records a correction method"
                )
            if formula is not None:
                raise DatasetProcessingStateError(
                    "dataset processing state total_protein_correction.formula "
                    "must be None when correction was not applied"
                )
            if (
                requires_log_scale is True
                or input_scale is not None
                or output_scale is not None
            ):
                raise DatasetProcessingStateError(
                    "dataset processing state total_protein_correction scale "
                    "provenance must be absent when correction was not applied"
                )
            if diagnostics is not None:
                _require_noop_total_correction_diagnostics(diagnostics)
        object.__setattr__(self, "applied", applied)
        object.__setattr__(self, "formula", formula)
        object.__setattr__(self, "requires_log_scale", requires_log_scale)
        object.__setattr__(self, "input_scale", input_scale)
        object.__setattr__(self, "output_scale", output_scale)
        object.__setattr__(self, "diagnostics", diagnostics)


def _require_bool(value: object, *, field_name: str) -> bool:
    if not isinstance(value, bool):
        raise DatasetProcessingStateError(f"{field_name} must be a bool")
    return value


def require_optional_bool(value: object, *, field_name: str) -> bool | None:
    if value is None:
        return None
    return _require_bool(value, field_name=field_name)


def _require_non_negative_int(value: object, *, field_name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise DatasetProcessingStateError(f"{field_name} must be an int")
    if value < 0:
        raise DatasetProcessingStateError(f"{field_name} must be >= 0")
    return value


def require_optional_non_negative_int(
    value: object,
    *,
    field_name: str,
) -> int | None:
    if value is None:
        return None
    return _require_non_negative_int(value, field_name=field_name)


def require_optional_str(value: object, *, field_name: str) -> str | None:
    if value is None:
        return None
    if not isinstance(value, str):
        raise DatasetProcessingStateError(f"{field_name} must be a string")
    normalized = value.strip()
    if not normalized:
        raise DatasetProcessingStateError(f"{field_name} must be a non-empty string")
    return normalized


def _diagnostic_int(
    diagnostics: MissingDataDiagnostics | TotalProteinCorrectionDiagnostics,
    *,
    key: str,
    field_name: str,
) -> int:
    value = diagnostics.get(key)
    if isinstance(value, bool) or not isinstance(value, int):
        raise DatasetProcessingStateError(f"{field_name} must be an int")
    return value


def _resolve_diagnostic_optional_str(
    *,
    diagnostics: MissingDataDiagnostics,
    key: str,
    current: str | None,
    field_name: str,
) -> str | None:
    diagnostic_value = diagnostics.get(key)
    if diagnostic_value is None:
        return current
    diagnostic_text = require_optional_str(
        diagnostic_value,
        field_name=f"dataset processing state missing_data.diagnostics.{key}",
    )
    if current is None:
        return diagnostic_text
    if current == diagnostic_text:
        return current
    raise DatasetProcessingStateError(
        f"{field_name} must match missing_data.diagnostics.{key}"
    )


def _require_optional_imputation_input_scale_state(value: object) -> str | None:
    parsed = require_optional_str(
        value,
        field_name="dataset processing state missing_data.imputation_input_scale",
    )
    if parsed is None or parsed in {"linear", "log2"}:
        return parsed
    raise DatasetProcessingStateError(
        "dataset processing state missing_data.imputation_input_scale must be "
        "one of: linear, log2"
    )


def _require_optional_imputation_operation_order_state(value: object) -> str | None:
    parsed = require_optional_str(
        value,
        field_name="dataset processing state missing_data.imputation_operation_order",
    )
    if parsed is None or parsed in IMPUTATION_OPERATION_ORDERS:
        return parsed
    supported = ", ".join(sorted(IMPUTATION_OPERATION_ORDERS))
    raise DatasetProcessingStateError(
        "dataset processing state missing_data.imputation_operation_order must be "
        f"one of: {supported}"
    )


def _require_string_non_negative_int_mapping(
    value: object,
    *,
    field_name: str,
) -> dict[str, int]:
    if not isinstance(value, dict):
        raise DatasetProcessingStateError(f"{field_name} must be a dict")
    parsed: dict[str, int] = {}
    for raw_key, raw_value in value.items():
        key = require_optional_str(raw_key, field_name=f"{field_name}.<key>")
        if key is None:  # pragma: no cover - require_optional_str handles this
            raise DatasetProcessingStateError(f"{field_name}.<key> is required")
        parsed[key] = _require_non_negative_int(
            raw_value,
            field_name=f"{field_name}.{key}",
        )
    return parsed


def _require_site_sequence_row_diagnostics(
    value: object,
) -> tuple[SiteSequenceResolutionRowDiagnostic, ...]:
    if not isinstance(value, tuple):
        raise DatasetProcessingStateError(
            "dataset processing state site_sequence_resolution.row_diagnostics "
            "must be a tuple"
        )
    for item in value:
        if not isinstance(item, SiteSequenceResolutionRowDiagnostic):
            raise DatasetProcessingStateError(
                "dataset processing state site_sequence_resolution.row_diagnostics "
                "must contain SiteSequenceResolutionRowDiagnostic instances"
            )
    return value


def _site_sequence_row_records_conflict(
    item: SiteSequenceResolutionRowDiagnostic,
) -> bool:
    return "conflict" in item.status.lower() or (
        item.reason is not None and "conflict" in item.reason.lower()
    )


def _count_site_sequence_row_conflicts(
    row_diagnostics: tuple[SiteSequenceResolutionRowDiagnostic, ...],
) -> int:
    return sum(
        1 for item in row_diagnostics if _site_sequence_row_records_conflict(item)
    )


def _require_total_correction_diagnostics_policy(
    diagnostics: TotalProteinCorrectionDiagnostics,
    *,
    expected_policy: TotalProteinCorrectionPolicy,
    applied: bool,
) -> None:
    if (
        applied
        and diagnostics.get("policy") is None
        and diagnostics.get("resolved_policy") is None
    ):
        raise DatasetProcessingStateError(
            "dataset processing state total_protein_correction.diagnostics "
            "must record policy or resolved_policy when correction was applied"
        )
    for key in ("policy", "resolved_policy"):
        value = diagnostics.get(key)
        if value is None:
            continue
        parsed = TotalProteinCorrectionPolicy.parse(
            value,
            field_name=(
                f"dataset processing state total_protein_correction.diagnostics.{key}"
            ),
        )
        if parsed is TotalProteinCorrectionPolicy.NONE and applied:
            raise DatasetProcessingStateError(
                "dataset processing state total_protein_correction.diagnostics "
                "cannot record policy 'none' when correction was applied"
            )
        if parsed is not expected_policy:
            raise DatasetProcessingStateError(
                "dataset processing state total_protein_correction.diagnostics "
                f"{key} must match total_protein_correction.policy"
            )


def _require_noop_total_correction_diagnostics(
    diagnostics: TotalProteinCorrectionDiagnostics,
) -> None:
    for key in ("policy", "requested_policy", "resolved_policy"):
        value = diagnostics.get(key)
        if value is None:
            continue
        parsed = TotalProteinCorrectionPolicy.parse(
            value,
            field_name=(
                f"dataset processing state total_protein_correction.diagnostics.{key}"
            ),
        )
        if parsed is not TotalProteinCorrectionPolicy.NONE:
            raise DatasetProcessingStateError(
                "dataset processing state total_protein_correction.diagnostics "
                "must not record correction provenance when correction was not applied"
            )
    if diagnostics.get("formula") is not None:
        raise DatasetProcessingStateError(
            "dataset processing state total_protein_correction.diagnostics.formula "
            "must be absent when correction was not applied"
        )
    corrected_row_count = diagnostics.get("corrected_row_count")
    if isinstance(corrected_row_count, int) and corrected_row_count > 0:
        raise DatasetProcessingStateError(
            "dataset processing state total_protein_correction.diagnostics "
            "must not record corrected rows when correction was not applied"
        )
    correction_provenance_keys = (
        "requires_log_scale",
        "input_scale",
        "output_scale",
        "matched_rows",
        "identity_mode",
        "identity_matching_policy",
        "phosphosite_key",
        "total_protein_key",
        "mapping_phosphosite_key",
        "mapping_total_protein_key",
        "mapping_table_fingerprint",
        "duplicate_policy",
        "unmatched_policy",
        "phosphosite_row_count",
        "total_protein_row_count",
        "corrected_row_count",
        "uncorrected_row_count",
        "unused_total_protein_row_count",
        "total_rows_used_by_multiple_phosphosites",
        "corrected_phosphosite_row_ids",
        "corrected_phosphosite_to_total_protein_row_id",
        "unmatched_phosphosite_row_ids",
        "uncorrected_phosphosite_row_reasons",
        "unused_total_protein_row_ids",
        "gene_symbol_matching_used",
        "gene_symbol_identity_warning",
        "total_table_hash",
        "input_phospho_hash",
        "output_phospho_hash",
    )
    for key in correction_provenance_keys:
        if diagnostics.get(key) is not None:
            raise DatasetProcessingStateError(
                "dataset processing state total_protein_correction.diagnostics "
                "must not record correction provenance when correction was not applied"
            )


@dataclass(frozen=True, slots=True)
class SiteMatrixState:
    """Site-matrix construction state at the analysis-ready dataset boundary."""

    policy: str
    constructed: bool
    missing_data_policy: str
    minimum_observed_values: int | None
    duplicate_site_policy: str


@dataclass(frozen=True, slots=True)
class ComparisonState:
    """Comparison-building state at the analysis-ready dataset boundary."""

    policy: str
    sample_group_column: str
    pairs: tuple[tuple[str, str], ...] | None


@dataclass(frozen=True, slots=True)
class RuvReadinessState:
    """Report-only RUV-readiness metadata state."""

    enabled: bool
    ready: bool
    reasons: tuple[str, ...]
    control_feature_column: str
    replicate_group_column: str
    batch_column: str | None
    control_feature_count: int
    replicate_group_count: int
    batch_count: int | None
    requires_complete_matrix: bool
    matrix_complete: bool
    imputation_method_id: str | None
    missingness_mask_preserved: bool


def default_ruv_readiness_state() -> RuvReadinessState:
    """Return the default disabled readiness state."""

    return RuvReadinessState(
        enabled=False,
        ready=False,
        reasons=("not configured",),
        control_feature_column="is_control_feature",
        replicate_group_column="replicate_group",
        batch_column="batch",
        control_feature_count=0,
        replicate_group_count=0,
        batch_count=0,
        requires_complete_matrix=True,
        matrix_complete=False,
        imputation_method_id=None,
        missingness_mask_preserved=False,
    )


@dataclass(frozen=True, slots=True)
class DatasetProcessingState:
    """Compact summary of preprocessing state at the analysis-ready boundary."""

    intensity_scale: IntensityScaleState
    site_sequence_resolution: SiteSequenceResolutionState
    missing_data: MissingDataState
    normalisation: NormalisationState
    total_protein_correction: TotalProteinCorrectionState
    site_matrix: SiteMatrixState
    comparisons: ComparisonState
    ruv_readiness: RuvReadinessState = field(
        default_factory=default_ruv_readiness_state
    )

    def __post_init__(self) -> None:
        _require_processing_state_instance(
            self.intensity_scale,
            expected_type=IntensityScaleState,
            field_name="dataset processing state intensity_scale",
        )
        _require_processing_state_instance(
            self.site_sequence_resolution,
            expected_type=SiteSequenceResolutionState,
            field_name="dataset processing state site_sequence_resolution",
        )
        _require_processing_state_instance(
            self.missing_data,
            expected_type=MissingDataState,
            field_name="dataset processing state missing_data",
        )
        _require_processing_state_instance(
            self.normalisation,
            expected_type=NormalisationState,
            field_name="dataset processing state normalisation",
        )
        _require_processing_state_instance(
            self.total_protein_correction,
            expected_type=TotalProteinCorrectionState,
            field_name="dataset processing state total_protein_correction",
        )
        _require_processing_state_instance(
            self.site_matrix,
            expected_type=SiteMatrixState,
            field_name="dataset processing state site_matrix",
        )
        _require_processing_state_instance(
            self.comparisons,
            expected_type=ComparisonState,
            field_name="dataset processing state comparisons",
        )
        _require_processing_state_instance(
            self.ruv_readiness,
            expected_type=RuvReadinessState,
            field_name="dataset processing state ruv_readiness",
        )


def _require_processing_state_instance(
    value: object,
    *,
    expected_type: type[object],
    field_name: str,
) -> None:
    if not isinstance(value, expected_type):
        raise DatasetProcessingStateError(
            f"{field_name} must be a {expected_type.__name__} instance"
        )
