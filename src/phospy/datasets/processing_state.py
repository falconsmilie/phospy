"""Dataset preprocessing-state summary models."""

from __future__ import annotations

from collections.abc import Iterator, Mapping
from dataclasses import dataclass, field
from typing import TypeAlias

from phospy.errors.input import PhosPyInputError
from phospy.transformations.models import IntensityScaleState

JsonPrimitive: TypeAlias = None | str | bool | int | float
JsonValue: TypeAlias = JsonPrimitive | list["JsonValue"] | dict[str, "JsonValue"]

MISSING_DATA_DIAGNOSTICS_SCHEMA_VERSION_V1 = 1
TOTAL_PROTEIN_CORRECTION_DIAGNOSTICS_SCHEMA_VERSION_V1 = 1
_V1_KNOWN_MISSING_DATA_DIAGNOSTICS_FIELDS = frozenset(
    (
        "diagnostics_schema_version",
        "missing_data_policy",
        "imputation_method_id",
        "imputation_method_family",
        "input_missing_cell_count",
        "output_missing_cell_count",
        "imputed_cell_count",
        "affected_row_count",
        "affected_column_count",
        "affected_row_ids",
        "affected_column_ids",
        "imputed_row_ids",
        "imputed_column_ids",
        "dropped_row_ids",
        "random_seed",
        "method_parameters",
        "matrix_scale_requirement",
        "stage_order",
        "missingness_mask_hash",
        "left_censored_assumption",
        "rows_not_imputable",
    )
)
_V1_KNOWN_DIAGNOSTICS_FIELDS = frozenset(
    (
        "diagnostics_schema_version",
        "policy",
        "requested_policy",
        "resolved_policy",
        "formula",
        "requires_log_scale",
        "input_scale",
        "output_scale",
        "quantitative_meaning",
        "matched_rows",
        "identity_mode",
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
)


class TotalProteinCorrectionDiagnostics(Mapping[str, JsonValue]):
    """Typed diagnostics contract for total-protein correction state."""

    diagnostics_schema_version: int

    @classmethod
    def from_payload(
        cls,
        payload: object,
        *,
        field_name: str,
    ) -> TotalProteinCorrectionDiagnostics:
        mapping = _require_mapping(payload, field_name=field_name)
        return TotalProteinCorrectionDiagnosticsV1.from_mapping(
            mapping,
            field_name=field_name,
        )

    def to_payload(self) -> dict[str, JsonValue]:
        """Return normalized diagnostics payload suitable for bundle JSON."""

        raise NotImplementedError


class MissingDataDiagnostics(Mapping[str, JsonValue]):
    """Typed diagnostics contract for missing-data preprocessing state."""

    diagnostics_schema_version: int

    @classmethod
    def from_payload(
        cls,
        payload: object,
        *,
        field_name: str,
    ) -> MissingDataDiagnostics:
        mapping = _require_mapping(payload, field_name=field_name)
        return MissingDataDiagnosticsV1.from_mapping(mapping, field_name=field_name)

    def to_payload(self) -> dict[str, JsonValue]:
        """Return normalized diagnostics payload suitable for bundle JSON."""

        raise NotImplementedError


@dataclass(frozen=True, slots=True, eq=False)
class MissingDataDiagnosticsV1(MissingDataDiagnostics):
    """Versioned missing-data diagnostics payload (schema v1)."""

    missing_data_policy: str
    input_missing_cell_count: int
    output_missing_cell_count: int
    imputed_cell_count: int
    affected_row_count: int
    affected_column_count: int
    affected_row_ids: tuple[str, ...]
    affected_column_ids: tuple[str, ...]
    imputed_row_ids: tuple[str, ...]
    imputed_column_ids: tuple[str, ...]
    dropped_row_ids: tuple[str, ...]
    method_parameters: dict[str, JsonValue]
    stage_order: tuple[str, ...]
    missingness_mask_hash: str
    rows_not_imputable: tuple[str, ...]
    imputation_method_id: str | None = None
    imputation_method_family: str | None = None
    random_seed: int | None = None
    matrix_scale_requirement: str | None = None
    left_censored_assumption: bool | None = None
    diagnostics_schema_version: int = MISSING_DATA_DIAGNOSTICS_SCHEMA_VERSION_V1
    _payload: dict[str, JsonValue] = field(init=False, repr=False, compare=False)

    @classmethod
    def from_mapping(
        cls,
        payload: Mapping[str, object],
        *,
        field_name: str,
    ) -> MissingDataDiagnosticsV1:
        if "diagnostics_schema_version" not in payload:
            raise PhosPyInputError(
                f"{field_name}.diagnostics_schema_version is required"
            )
        return cls._from_versioned_payload(payload, field_name=field_name)

    @classmethod
    def _from_versioned_payload(
        cls,
        payload: Mapping[str, object],
        *,
        field_name: str,
    ) -> MissingDataDiagnosticsV1:
        _require_string_keys(payload, field_name=field_name)
        version = _require_int(
            payload.get("diagnostics_schema_version"),
            field_name=f"{field_name}.diagnostics_schema_version",
        )
        if version != MISSING_DATA_DIAGNOSTICS_SCHEMA_VERSION_V1:
            raise PhosPyInputError(
                f"{field_name}.diagnostics_schema_version={version!r} is unsupported; "
                f"expected {MISSING_DATA_DIAGNOSTICS_SCHEMA_VERSION_V1}"
            )
        unknown_fields = sorted(
            key
            for key in payload
            if key not in _V1_KNOWN_MISSING_DATA_DIAGNOSTICS_FIELDS
        )
        if unknown_fields:
            raise PhosPyInputError(
                f"{field_name} contains unsupported field(s): "
                + ", ".join(unknown_fields)
            )
        return cls(
            diagnostics_schema_version=version,
            missing_data_policy=_require_required_str(
                payload.get("missing_data_policy"),
                field_name=f"{field_name}.missing_data_policy",
            ),
            imputation_method_id=_require_optional_str(
                payload.get("imputation_method_id"),
                field_name=f"{field_name}.imputation_method_id",
            ),
            imputation_method_family=_require_optional_str(
                payload.get("imputation_method_family"),
                field_name=f"{field_name}.imputation_method_family",
            ),
            input_missing_cell_count=_require_required_non_negative_int(
                payload.get("input_missing_cell_count"),
                field_name=f"{field_name}.input_missing_cell_count",
            ),
            output_missing_cell_count=_require_required_non_negative_int(
                payload.get("output_missing_cell_count"),
                field_name=f"{field_name}.output_missing_cell_count",
            ),
            imputed_cell_count=_require_required_non_negative_int(
                payload.get("imputed_cell_count"),
                field_name=f"{field_name}.imputed_cell_count",
            ),
            affected_row_count=_require_required_non_negative_int(
                payload.get("affected_row_count"),
                field_name=f"{field_name}.affected_row_count",
            ),
            affected_column_count=_require_required_non_negative_int(
                payload.get("affected_column_count"),
                field_name=f"{field_name}.affected_column_count",
            ),
            affected_row_ids=_require_required_string_tuple(
                payload.get("affected_row_ids"),
                field_name=f"{field_name}.affected_row_ids",
            ),
            affected_column_ids=_require_required_string_tuple(
                payload.get("affected_column_ids"),
                field_name=f"{field_name}.affected_column_ids",
            ),
            imputed_row_ids=_require_required_string_tuple(
                payload.get("imputed_row_ids"),
                field_name=f"{field_name}.imputed_row_ids",
            ),
            imputed_column_ids=_require_required_string_tuple(
                payload.get("imputed_column_ids"),
                field_name=f"{field_name}.imputed_column_ids",
            ),
            dropped_row_ids=_require_required_string_tuple(
                payload.get("dropped_row_ids"),
                field_name=f"{field_name}.dropped_row_ids",
            ),
            random_seed=_require_optional_int(
                payload.get("random_seed"),
                field_name=f"{field_name}.random_seed",
            ),
            method_parameters=_require_json_mapping(
                payload.get("method_parameters"),
                field_name=f"{field_name}.method_parameters",
            ),
            matrix_scale_requirement=_require_optional_str(
                payload.get("matrix_scale_requirement"),
                field_name=f"{field_name}.matrix_scale_requirement",
            ),
            stage_order=_require_required_string_tuple(
                payload.get("stage_order"),
                field_name=f"{field_name}.stage_order",
            ),
            missingness_mask_hash=_require_required_str(
                payload.get("missingness_mask_hash"),
                field_name=f"{field_name}.missingness_mask_hash",
            ),
            left_censored_assumption=_require_optional_bool(
                payload.get("left_censored_assumption"),
                field_name=f"{field_name}.left_censored_assumption",
            ),
            rows_not_imputable=_require_required_string_tuple(
                payload.get("rows_not_imputable"),
                field_name=f"{field_name}.rows_not_imputable",
            ),
        )

    def __post_init__(self) -> None:
        if (
            self.diagnostics_schema_version
            != MISSING_DATA_DIAGNOSTICS_SCHEMA_VERSION_V1
        ):
            raise PhosPyInputError(
                "dataset processing state missing_data diagnostics schema version "
                f"must be {MISSING_DATA_DIAGNOSTICS_SCHEMA_VERSION_V1}"
            )
        missing_data_policy = _require_required_str(
            self.missing_data_policy,
            field_name=(
                "dataset processing state missing_data.diagnostics.missing_data_policy"
            ),
        )
        imputation_method_id = _require_optional_str(
            self.imputation_method_id,
            field_name=(
                "dataset processing state missing_data.diagnostics.imputation_method_id"
            ),
        )
        imputation_method_family = _require_optional_str(
            self.imputation_method_family,
            field_name=(
                "dataset processing state missing_data.diagnostics."
                "imputation_method_family"
            ),
        )
        input_missing_cell_count = _require_required_non_negative_int(
            self.input_missing_cell_count,
            field_name=(
                "dataset processing state missing_data.diagnostics."
                "input_missing_cell_count"
            ),
        )
        output_missing_cell_count = _require_required_non_negative_int(
            self.output_missing_cell_count,
            field_name=(
                "dataset processing state missing_data.diagnostics."
                "output_missing_cell_count"
            ),
        )
        imputed_cell_count = _require_required_non_negative_int(
            self.imputed_cell_count,
            field_name=(
                "dataset processing state missing_data.diagnostics.imputed_cell_count"
            ),
        )
        affected_row_count = _require_required_non_negative_int(
            self.affected_row_count,
            field_name=(
                "dataset processing state missing_data.diagnostics.affected_row_count"
            ),
        )
        affected_column_count = _require_required_non_negative_int(
            self.affected_column_count,
            field_name=(
                "dataset processing state missing_data.diagnostics."
                "affected_column_count"
            ),
        )
        affected_row_ids = _require_required_string_tuple(
            self.affected_row_ids,
            field_name=(
                "dataset processing state missing_data.diagnostics.affected_row_ids"
            ),
        )
        affected_column_ids = _require_required_string_tuple(
            self.affected_column_ids,
            field_name=(
                "dataset processing state missing_data.diagnostics.affected_column_ids"
            ),
        )
        imputed_row_ids = _require_required_string_tuple(
            self.imputed_row_ids,
            field_name=(
                "dataset processing state missing_data.diagnostics.imputed_row_ids"
            ),
        )
        imputed_column_ids = _require_required_string_tuple(
            self.imputed_column_ids,
            field_name=(
                "dataset processing state missing_data.diagnostics.imputed_column_ids"
            ),
        )
        dropped_row_ids = _require_required_string_tuple(
            self.dropped_row_ids,
            field_name=(
                "dataset processing state missing_data.diagnostics.dropped_row_ids"
            ),
        )
        random_seed = _require_optional_int(
            self.random_seed,
            field_name="dataset processing state missing_data.diagnostics.random_seed",
        )
        method_parameters = _require_json_mapping(
            self.method_parameters,
            field_name=(
                "dataset processing state missing_data.diagnostics.method_parameters"
            ),
        )
        matrix_scale_requirement = _require_optional_str(
            self.matrix_scale_requirement,
            field_name=(
                "dataset processing state missing_data.diagnostics."
                "matrix_scale_requirement"
            ),
        )
        stage_order = _require_required_string_tuple(
            self.stage_order,
            field_name="dataset processing state missing_data.diagnostics.stage_order",
        )
        missingness_mask_hash = _require_required_str(
            self.missingness_mask_hash,
            field_name=(
                "dataset processing state missing_data.diagnostics.missingness_mask_hash"
            ),
        )
        left_censored_assumption = _require_optional_bool(
            self.left_censored_assumption,
            field_name=(
                "dataset processing state missing_data.diagnostics."
                "left_censored_assumption"
            ),
        )
        rows_not_imputable = _require_required_string_tuple(
            self.rows_not_imputable,
            field_name=(
                "dataset processing state missing_data.diagnostics.rows_not_imputable"
            ),
        )
        payload: dict[str, JsonValue] = {
            "diagnostics_schema_version": self.diagnostics_schema_version,
            "missing_data_policy": missing_data_policy,
            "input_missing_cell_count": input_missing_cell_count,
            "output_missing_cell_count": output_missing_cell_count,
            "imputed_cell_count": imputed_cell_count,
            "affected_row_count": affected_row_count,
            "affected_column_count": affected_column_count,
            "affected_row_ids": list(affected_row_ids),
            "affected_column_ids": list(affected_column_ids),
            "imputed_row_ids": list(imputed_row_ids),
            "imputed_column_ids": list(imputed_column_ids),
            "dropped_row_ids": list(dropped_row_ids),
            "method_parameters": dict(method_parameters),
            "stage_order": list(stage_order),
            "missingness_mask_hash": missingness_mask_hash,
            "rows_not_imputable": list(rows_not_imputable),
        }
        _set_optional_payload_value(
            payload, "imputation_method_id", imputation_method_id
        )
        _set_optional_payload_value(
            payload, "imputation_method_family", imputation_method_family
        )
        _set_optional_payload_value(payload, "random_seed", random_seed)
        _set_optional_payload_value(
            payload, "matrix_scale_requirement", matrix_scale_requirement
        )
        _set_optional_payload_value(
            payload, "left_censored_assumption", left_censored_assumption
        )

        object.__setattr__(self, "missing_data_policy", missing_data_policy)
        object.__setattr__(self, "imputation_method_id", imputation_method_id)
        object.__setattr__(self, "imputation_method_family", imputation_method_family)
        object.__setattr__(self, "input_missing_cell_count", input_missing_cell_count)
        object.__setattr__(self, "output_missing_cell_count", output_missing_cell_count)
        object.__setattr__(self, "imputed_cell_count", imputed_cell_count)
        object.__setattr__(self, "affected_row_count", affected_row_count)
        object.__setattr__(self, "affected_column_count", affected_column_count)
        object.__setattr__(self, "affected_row_ids", affected_row_ids)
        object.__setattr__(self, "affected_column_ids", affected_column_ids)
        object.__setattr__(self, "imputed_row_ids", imputed_row_ids)
        object.__setattr__(self, "imputed_column_ids", imputed_column_ids)
        object.__setattr__(self, "dropped_row_ids", dropped_row_ids)
        object.__setattr__(self, "random_seed", random_seed)
        object.__setattr__(self, "method_parameters", dict(method_parameters))
        object.__setattr__(self, "matrix_scale_requirement", matrix_scale_requirement)
        object.__setattr__(self, "stage_order", stage_order)
        object.__setattr__(self, "missingness_mask_hash", missingness_mask_hash)
        object.__setattr__(self, "left_censored_assumption", left_censored_assumption)
        object.__setattr__(self, "rows_not_imputable", rows_not_imputable)
        object.__setattr__(self, "_payload", payload)

    def __getitem__(self, key: str) -> JsonValue:
        return self._payload[key]

    def __iter__(self) -> Iterator[str]:
        return iter(self._payload)

    def __len__(self) -> int:
        return len(self._payload)

    def to_payload(self) -> dict[str, JsonValue]:
        return dict(self._payload)


@dataclass(frozen=True, slots=True, eq=False)
class TotalProteinCorrectionDiagnosticsV1(TotalProteinCorrectionDiagnostics):
    """Versioned total-protein correction diagnostics payload (schema v1)."""

    policy: str | None = None
    requested_policy: str | None = None
    resolved_policy: str | None = None
    formula: str | None = None
    requires_log_scale: bool | None = None
    input_scale: str | None = None
    output_scale: str | None = None
    quantitative_meaning: str | None = None
    matched_rows: int | None = None
    identity_mode: str | None = None
    phosphosite_key: str | None = None
    total_protein_key: str | None = None
    mapping_phosphosite_key: str | None = None
    mapping_total_protein_key: str | None = None
    mapping_table_fingerprint: str | None = None
    duplicate_policy: str | None = None
    unmatched_policy: str | None = None
    phosphosite_row_count: int | None = None
    total_protein_row_count: int | None = None
    corrected_row_count: int | None = None
    uncorrected_row_count: int | None = None
    unused_total_protein_row_count: int | None = None
    total_rows_used_by_multiple_phosphosites: int | None = None
    corrected_phosphosite_row_ids: tuple[str, ...] | None = None
    corrected_phosphosite_to_total_protein_row_id: dict[str, str] | None = None
    unmatched_phosphosite_row_ids: tuple[str, ...] | None = None
    uncorrected_phosphosite_row_reasons: dict[str, str] | None = None
    unused_total_protein_row_ids: tuple[str, ...] | None = None
    gene_symbol_matching_used: bool | None = None
    gene_symbol_identity_warning: str | None = None
    total_table_hash: str | None = None
    input_phospho_hash: str | None = None
    output_phospho_hash: str | None = None
    diagnostics_schema_version: int = (
        TOTAL_PROTEIN_CORRECTION_DIAGNOSTICS_SCHEMA_VERSION_V1
    )
    _payload: dict[str, JsonValue] = field(init=False, repr=False, compare=False)

    @classmethod
    def from_mapping(
        cls,
        payload: Mapping[str, object],
        *,
        field_name: str,
    ) -> TotalProteinCorrectionDiagnosticsV1:
        if "diagnostics_schema_version" not in payload:
            raise PhosPyInputError(
                f"{field_name}.diagnostics_schema_version is required"
            )
        return cls._from_versioned_payload(payload, field_name=field_name)

    @classmethod
    def _from_versioned_payload(
        cls,
        payload: Mapping[str, object],
        *,
        field_name: str,
    ) -> TotalProteinCorrectionDiagnosticsV1:
        _require_string_keys(payload, field_name=field_name)
        version = _require_int(
            payload.get("diagnostics_schema_version"),
            field_name=f"{field_name}.diagnostics_schema_version",
        )
        if version != TOTAL_PROTEIN_CORRECTION_DIAGNOSTICS_SCHEMA_VERSION_V1:
            raise PhosPyInputError(
                f"{field_name}.diagnostics_schema_version={version!r} is unsupported; "
                f"expected {TOTAL_PROTEIN_CORRECTION_DIAGNOSTICS_SCHEMA_VERSION_V1}"
            )
        unknown_fields = sorted(
            key for key in payload if key not in _V1_KNOWN_DIAGNOSTICS_FIELDS
        )
        if unknown_fields:
            raise PhosPyInputError(
                f"{field_name} contains unsupported field(s): "
                + ", ".join(unknown_fields)
            )
        return cls(
            policy=_require_optional_str(
                payload.get("policy"),
                field_name=f"{field_name}.policy",
            ),
            requested_policy=_require_optional_str(
                payload.get("requested_policy"),
                field_name=f"{field_name}.requested_policy",
            ),
            resolved_policy=_require_optional_str(
                payload.get("resolved_policy"),
                field_name=f"{field_name}.resolved_policy",
            ),
            formula=_require_optional_str(
                payload.get("formula"),
                field_name=f"{field_name}.formula",
            ),
            requires_log_scale=_require_optional_bool(
                payload.get("requires_log_scale"),
                field_name=f"{field_name}.requires_log_scale",
            ),
            input_scale=_require_optional_str(
                payload.get("input_scale"),
                field_name=f"{field_name}.input_scale",
            ),
            output_scale=_require_optional_str(
                payload.get("output_scale"),
                field_name=f"{field_name}.output_scale",
            ),
            quantitative_meaning=_require_required_str(
                payload.get("quantitative_meaning"),
                field_name=f"{field_name}.quantitative_meaning",
            ),
            matched_rows=_require_optional_non_negative_int(
                payload.get("matched_rows"),
                field_name=f"{field_name}.matched_rows",
            ),
            identity_mode=_require_optional_str(
                payload.get("identity_mode"),
                field_name=f"{field_name}.identity_mode",
            ),
            phosphosite_key=_require_optional_str(
                payload.get("phosphosite_key"),
                field_name=f"{field_name}.phosphosite_key",
            ),
            total_protein_key=_require_optional_str(
                payload.get("total_protein_key"),
                field_name=f"{field_name}.total_protein_key",
            ),
            mapping_phosphosite_key=_require_optional_str(
                payload.get("mapping_phosphosite_key"),
                field_name=f"{field_name}.mapping_phosphosite_key",
            ),
            mapping_total_protein_key=_require_optional_str(
                payload.get("mapping_total_protein_key"),
                field_name=f"{field_name}.mapping_total_protein_key",
            ),
            mapping_table_fingerprint=_require_optional_str(
                payload.get("mapping_table_fingerprint"),
                field_name=f"{field_name}.mapping_table_fingerprint",
            ),
            duplicate_policy=_require_optional_str(
                payload.get("duplicate_policy"),
                field_name=f"{field_name}.duplicate_policy",
            ),
            unmatched_policy=_require_optional_str(
                payload.get("unmatched_policy"),
                field_name=f"{field_name}.unmatched_policy",
            ),
            phosphosite_row_count=_require_optional_non_negative_int(
                payload.get("phosphosite_row_count"),
                field_name=f"{field_name}.phosphosite_row_count",
            ),
            total_protein_row_count=_require_optional_non_negative_int(
                payload.get("total_protein_row_count"),
                field_name=f"{field_name}.total_protein_row_count",
            ),
            corrected_row_count=_require_optional_non_negative_int(
                payload.get("corrected_row_count"),
                field_name=f"{field_name}.corrected_row_count",
            ),
            uncorrected_row_count=_require_optional_non_negative_int(
                payload.get("uncorrected_row_count"),
                field_name=f"{field_name}.uncorrected_row_count",
            ),
            unused_total_protein_row_count=_require_optional_non_negative_int(
                payload.get("unused_total_protein_row_count"),
                field_name=f"{field_name}.unused_total_protein_row_count",
            ),
            total_rows_used_by_multiple_phosphosites=_require_optional_non_negative_int(
                payload.get("total_rows_used_by_multiple_phosphosites"),
                field_name=f"{field_name}.total_rows_used_by_multiple_phosphosites",
            ),
            corrected_phosphosite_row_ids=_require_optional_string_tuple(
                payload.get("corrected_phosphosite_row_ids"),
                field_name=f"{field_name}.corrected_phosphosite_row_ids",
            ),
            corrected_phosphosite_to_total_protein_row_id=_require_optional_string_to_string_mapping(
                payload.get("corrected_phosphosite_to_total_protein_row_id"),
                field_name=(
                    f"{field_name}.corrected_phosphosite_to_total_protein_row_id"
                ),
            ),
            unmatched_phosphosite_row_ids=_require_optional_string_tuple(
                payload.get("unmatched_phosphosite_row_ids"),
                field_name=f"{field_name}.unmatched_phosphosite_row_ids",
            ),
            uncorrected_phosphosite_row_reasons=_require_optional_string_to_string_mapping(
                payload.get("uncorrected_phosphosite_row_reasons"),
                field_name=f"{field_name}.uncorrected_phosphosite_row_reasons",
            ),
            unused_total_protein_row_ids=_require_optional_string_tuple(
                payload.get("unused_total_protein_row_ids"),
                field_name=f"{field_name}.unused_total_protein_row_ids",
            ),
            gene_symbol_matching_used=_require_optional_bool(
                payload.get("gene_symbol_matching_used"),
                field_name=f"{field_name}.gene_symbol_matching_used",
            ),
            gene_symbol_identity_warning=_require_optional_str(
                payload.get("gene_symbol_identity_warning"),
                field_name=f"{field_name}.gene_symbol_identity_warning",
            ),
            total_table_hash=_require_optional_str(
                payload.get("total_table_hash"),
                field_name=f"{field_name}.total_table_hash",
            ),
            input_phospho_hash=_require_optional_str(
                payload.get("input_phospho_hash"),
                field_name=f"{field_name}.input_phospho_hash",
            ),
            output_phospho_hash=_require_optional_str(
                payload.get("output_phospho_hash"),
                field_name=f"{field_name}.output_phospho_hash",
            ),
            diagnostics_schema_version=version,
        )

    def __post_init__(self) -> None:
        if (
            self.diagnostics_schema_version
            != TOTAL_PROTEIN_CORRECTION_DIAGNOSTICS_SCHEMA_VERSION_V1
        ):
            raise PhosPyInputError(
                "dataset processing state total_protein_correction diagnostics "
                f"schema version must be {TOTAL_PROTEIN_CORRECTION_DIAGNOSTICS_SCHEMA_VERSION_V1}"
            )
        policy = _require_optional_str(
            self.policy,
            field_name="dataset processing state total_protein_correction.diagnostics.policy",
        )
        requested_policy = _require_optional_str(
            self.requested_policy,
            field_name=(
                "dataset processing state total_protein_correction.diagnostics."
                "requested_policy"
            ),
        )
        resolved_policy = _require_optional_str(
            self.resolved_policy,
            field_name=(
                "dataset processing state total_protein_correction.diagnostics."
                "resolved_policy"
            ),
        )
        formula = _require_optional_str(
            self.formula,
            field_name="dataset processing state total_protein_correction.diagnostics.formula",
        )
        requires_log_scale = _require_optional_bool(
            self.requires_log_scale,
            field_name=(
                "dataset processing state total_protein_correction.diagnostics."
                "requires_log_scale"
            ),
        )
        input_scale = _require_optional_str(
            self.input_scale,
            field_name=(
                "dataset processing state total_protein_correction.diagnostics."
                "input_scale"
            ),
        )
        output_scale = _require_optional_str(
            self.output_scale,
            field_name=(
                "dataset processing state total_protein_correction.diagnostics."
                "output_scale"
            ),
        )
        quantitative_meaning = _require_required_str(
            self.quantitative_meaning,
            field_name=(
                "dataset processing state total_protein_correction.diagnostics."
                "quantitative_meaning"
            ),
        )
        matched_rows = _require_optional_non_negative_int(
            self.matched_rows,
            field_name=(
                "dataset processing state total_protein_correction.diagnostics."
                "matched_rows"
            ),
        )
        identity_mode = _require_optional_str(
            self.identity_mode,
            field_name=(
                "dataset processing state total_protein_correction.diagnostics."
                "identity_mode"
            ),
        )
        phosphosite_key = _require_optional_str(
            self.phosphosite_key,
            field_name=(
                "dataset processing state total_protein_correction.diagnostics."
                "phosphosite_key"
            ),
        )
        total_protein_key = _require_optional_str(
            self.total_protein_key,
            field_name=(
                "dataset processing state total_protein_correction.diagnostics."
                "total_protein_key"
            ),
        )
        mapping_phosphosite_key = _require_optional_str(
            self.mapping_phosphosite_key,
            field_name=(
                "dataset processing state total_protein_correction.diagnostics."
                "mapping_phosphosite_key"
            ),
        )
        mapping_total_protein_key = _require_optional_str(
            self.mapping_total_protein_key,
            field_name=(
                "dataset processing state total_protein_correction.diagnostics."
                "mapping_total_protein_key"
            ),
        )
        mapping_table_fingerprint = _require_optional_str(
            self.mapping_table_fingerprint,
            field_name=(
                "dataset processing state total_protein_correction.diagnostics."
                "mapping_table_fingerprint"
            ),
        )
        duplicate_policy = _require_optional_str(
            self.duplicate_policy,
            field_name=(
                "dataset processing state total_protein_correction.diagnostics."
                "duplicate_policy"
            ),
        )
        unmatched_policy = _require_optional_str(
            self.unmatched_policy,
            field_name=(
                "dataset processing state total_protein_correction.diagnostics."
                "unmatched_policy"
            ),
        )
        phosphosite_row_count = _require_optional_non_negative_int(
            self.phosphosite_row_count,
            field_name=(
                "dataset processing state total_protein_correction.diagnostics."
                "phosphosite_row_count"
            ),
        )
        total_protein_row_count = _require_optional_non_negative_int(
            self.total_protein_row_count,
            field_name=(
                "dataset processing state total_protein_correction.diagnostics."
                "total_protein_row_count"
            ),
        )
        corrected_row_count = _require_optional_non_negative_int(
            self.corrected_row_count,
            field_name=(
                "dataset processing state total_protein_correction.diagnostics."
                "corrected_row_count"
            ),
        )
        uncorrected_row_count = _require_optional_non_negative_int(
            self.uncorrected_row_count,
            field_name=(
                "dataset processing state total_protein_correction.diagnostics."
                "uncorrected_row_count"
            ),
        )
        unused_total_protein_row_count = _require_optional_non_negative_int(
            self.unused_total_protein_row_count,
            field_name=(
                "dataset processing state total_protein_correction.diagnostics."
                "unused_total_protein_row_count"
            ),
        )
        total_rows_used_by_multiple_phosphosites = _require_optional_non_negative_int(
            self.total_rows_used_by_multiple_phosphosites,
            field_name=(
                "dataset processing state total_protein_correction.diagnostics."
                "total_rows_used_by_multiple_phosphosites"
            ),
        )
        corrected_phosphosite_row_ids = _require_optional_string_tuple(
            self.corrected_phosphosite_row_ids,
            field_name=(
                "dataset processing state total_protein_correction.diagnostics."
                "corrected_phosphosite_row_ids"
            ),
        )
        corrected_phosphosite_to_total_protein_row_id = (
            _require_optional_string_to_string_mapping(
                self.corrected_phosphosite_to_total_protein_row_id,
                field_name=(
                    "dataset processing state total_protein_correction.diagnostics."
                    "corrected_phosphosite_to_total_protein_row_id"
                ),
            )
        )
        unmatched_phosphosite_row_ids = _require_optional_string_tuple(
            self.unmatched_phosphosite_row_ids,
            field_name=(
                "dataset processing state total_protein_correction.diagnostics."
                "unmatched_phosphosite_row_ids"
            ),
        )
        uncorrected_phosphosite_row_reasons = (
            _require_optional_string_to_string_mapping(
                self.uncorrected_phosphosite_row_reasons,
                field_name=(
                    "dataset processing state total_protein_correction.diagnostics."
                    "uncorrected_phosphosite_row_reasons"
                ),
            )
        )
        unused_total_protein_row_ids = _require_optional_string_tuple(
            self.unused_total_protein_row_ids,
            field_name=(
                "dataset processing state total_protein_correction.diagnostics."
                "unused_total_protein_row_ids"
            ),
        )
        gene_symbol_matching_used = _require_optional_bool(
            self.gene_symbol_matching_used,
            field_name=(
                "dataset processing state total_protein_correction.diagnostics."
                "gene_symbol_matching_used"
            ),
        )
        gene_symbol_identity_warning = _require_optional_str(
            self.gene_symbol_identity_warning,
            field_name=(
                "dataset processing state total_protein_correction.diagnostics."
                "gene_symbol_identity_warning"
            ),
        )
        total_table_hash = _require_optional_str(
            self.total_table_hash,
            field_name=(
                "dataset processing state total_protein_correction.diagnostics."
                "total_table_hash"
            ),
        )
        input_phospho_hash = _require_optional_str(
            self.input_phospho_hash,
            field_name=(
                "dataset processing state total_protein_correction.diagnostics."
                "input_phospho_hash"
            ),
        )
        output_phospho_hash = _require_optional_str(
            self.output_phospho_hash,
            field_name=(
                "dataset processing state total_protein_correction.diagnostics."
                "output_phospho_hash"
            ),
        )

        payload: dict[str, JsonValue] = {
            "diagnostics_schema_version": self.diagnostics_schema_version
        }
        _set_optional_payload_value(payload, "policy", policy)
        _set_optional_payload_value(payload, "requested_policy", requested_policy)
        _set_optional_payload_value(payload, "resolved_policy", resolved_policy)
        _set_optional_payload_value(payload, "formula", formula)
        _set_optional_payload_value(payload, "requires_log_scale", requires_log_scale)
        _set_optional_payload_value(payload, "input_scale", input_scale)
        _set_optional_payload_value(payload, "output_scale", output_scale)
        payload["quantitative_meaning"] = quantitative_meaning
        _set_optional_payload_value(payload, "matched_rows", matched_rows)
        _set_optional_payload_value(payload, "identity_mode", identity_mode)
        _set_optional_payload_value(payload, "phosphosite_key", phosphosite_key)
        _set_optional_payload_value(payload, "total_protein_key", total_protein_key)
        _set_optional_payload_value(
            payload, "mapping_phosphosite_key", mapping_phosphosite_key
        )
        _set_optional_payload_value(
            payload, "mapping_total_protein_key", mapping_total_protein_key
        )
        _set_optional_payload_value(
            payload, "mapping_table_fingerprint", mapping_table_fingerprint
        )
        _set_optional_payload_value(payload, "duplicate_policy", duplicate_policy)
        _set_optional_payload_value(payload, "unmatched_policy", unmatched_policy)
        _set_optional_payload_value(
            payload, "phosphosite_row_count", phosphosite_row_count
        )
        _set_optional_payload_value(
            payload, "total_protein_row_count", total_protein_row_count
        )
        _set_optional_payload_value(payload, "corrected_row_count", corrected_row_count)
        _set_optional_payload_value(
            payload, "uncorrected_row_count", uncorrected_row_count
        )
        _set_optional_payload_value(
            payload, "unused_total_protein_row_count", unused_total_protein_row_count
        )
        _set_optional_payload_value(
            payload,
            "total_rows_used_by_multiple_phosphosites",
            total_rows_used_by_multiple_phosphosites,
        )
        if corrected_phosphosite_row_ids is not None:
            payload["corrected_phosphosite_row_ids"] = list(
                corrected_phosphosite_row_ids
            )
        if corrected_phosphosite_to_total_protein_row_id is not None:
            payload["corrected_phosphosite_to_total_protein_row_id"] = dict(
                corrected_phosphosite_to_total_protein_row_id
            )
        if unmatched_phosphosite_row_ids is not None:
            payload["unmatched_phosphosite_row_ids"] = list(
                unmatched_phosphosite_row_ids
            )
        if uncorrected_phosphosite_row_reasons is not None:
            payload["uncorrected_phosphosite_row_reasons"] = dict(
                uncorrected_phosphosite_row_reasons
            )
        if unused_total_protein_row_ids is not None:
            payload["unused_total_protein_row_ids"] = list(unused_total_protein_row_ids)
        _set_optional_payload_value(
            payload, "gene_symbol_matching_used", gene_symbol_matching_used
        )
        _set_optional_payload_value(
            payload, "gene_symbol_identity_warning", gene_symbol_identity_warning
        )
        _set_optional_payload_value(payload, "total_table_hash", total_table_hash)
        _set_optional_payload_value(payload, "input_phospho_hash", input_phospho_hash)
        _set_optional_payload_value(payload, "output_phospho_hash", output_phospho_hash)

        object.__setattr__(self, "policy", policy)
        object.__setattr__(self, "requested_policy", requested_policy)
        object.__setattr__(self, "resolved_policy", resolved_policy)
        object.__setattr__(self, "formula", formula)
        object.__setattr__(self, "requires_log_scale", requires_log_scale)
        object.__setattr__(self, "input_scale", input_scale)
        object.__setattr__(self, "output_scale", output_scale)
        object.__setattr__(self, "quantitative_meaning", quantitative_meaning)
        object.__setattr__(self, "matched_rows", matched_rows)
        object.__setattr__(self, "identity_mode", identity_mode)
        object.__setattr__(self, "phosphosite_key", phosphosite_key)
        object.__setattr__(self, "total_protein_key", total_protein_key)
        object.__setattr__(self, "mapping_phosphosite_key", mapping_phosphosite_key)
        object.__setattr__(self, "mapping_total_protein_key", mapping_total_protein_key)
        object.__setattr__(self, "mapping_table_fingerprint", mapping_table_fingerprint)
        object.__setattr__(self, "duplicate_policy", duplicate_policy)
        object.__setattr__(self, "unmatched_policy", unmatched_policy)
        object.__setattr__(self, "phosphosite_row_count", phosphosite_row_count)
        object.__setattr__(self, "total_protein_row_count", total_protein_row_count)
        object.__setattr__(self, "corrected_row_count", corrected_row_count)
        object.__setattr__(self, "uncorrected_row_count", uncorrected_row_count)
        object.__setattr__(
            self, "unused_total_protein_row_count", unused_total_protein_row_count
        )
        object.__setattr__(
            self,
            "total_rows_used_by_multiple_phosphosites",
            total_rows_used_by_multiple_phosphosites,
        )
        object.__setattr__(
            self,
            "corrected_phosphosite_row_ids",
            corrected_phosphosite_row_ids,
        )
        object.__setattr__(
            self,
            "corrected_phosphosite_to_total_protein_row_id",
            corrected_phosphosite_to_total_protein_row_id,
        )
        object.__setattr__(
            self, "unmatched_phosphosite_row_ids", unmatched_phosphosite_row_ids
        )
        object.__setattr__(
            self,
            "uncorrected_phosphosite_row_reasons",
            uncorrected_phosphosite_row_reasons,
        )
        object.__setattr__(
            self, "unused_total_protein_row_ids", unused_total_protein_row_ids
        )
        object.__setattr__(self, "gene_symbol_matching_used", gene_symbol_matching_used)
        object.__setattr__(
            self, "gene_symbol_identity_warning", gene_symbol_identity_warning
        )
        object.__setattr__(self, "total_table_hash", total_table_hash)
        object.__setattr__(self, "input_phospho_hash", input_phospho_hash)
        object.__setattr__(self, "output_phospho_hash", output_phospho_hash)
        object.__setattr__(self, "_payload", payload)

    def __getitem__(self, key: str) -> JsonValue:
        return self._payload[key]

    def __iter__(self) -> Iterator[str]:
        return iter(self._payload)

    def __len__(self) -> int:
        return len(self._payload)

    def to_payload(self) -> dict[str, JsonValue]:
        return dict(self._payload)


@dataclass(frozen=True, slots=True)
class MissingDataState:
    """Missing-data policy state at the analysis-ready dataset boundary."""

    policy: str
    min_observed_values: int | None
    complete_matrix: bool
    imputed: bool
    diagnostics: MissingDataDiagnostics | None = None

    def __post_init__(self) -> None:
        if self.diagnostics is None:
            return
        if isinstance(self.diagnostics, MissingDataDiagnostics):
            return
        normalized = MissingDataDiagnostics.from_payload(
            self.diagnostics,
            field_name="dataset processing state missing_data.diagnostics",
        )
        object.__setattr__(self, "diagnostics", normalized)


@dataclass(frozen=True, slots=True)
class NormalisationState:
    """Normalisation policy state at the analysis-ready dataset boundary."""

    policy: str


@dataclass(frozen=True, slots=True)
class SiteSequenceResolutionState:
    """Dataset site-sequence FASTA-resolution state at preprocessing boundary."""

    configured: bool
    mode: str | None
    flank_size: int | None
    fasta_sha256: str | None
    resolved_site_count: int
    unresolved_site_count: int
    unresolved_counts_by_reason: dict[str, int]


@dataclass(frozen=True, slots=True)
class TotalProteinCorrectionState:
    """Total-protein correction state at the analysis-ready dataset boundary."""

    policy: str
    applied: bool
    formula: str | None = None
    requires_log_scale: bool | None = False
    input_scale: str | None = None
    output_scale: str | None = None
    quantitative_meaning: str | None = None
    diagnostics: TotalProteinCorrectionDiagnostics | None = None

    def __post_init__(self) -> None:
        if self.diagnostics is None:
            return
        if isinstance(self.diagnostics, TotalProteinCorrectionDiagnostics):
            return
        normalized = TotalProteinCorrectionDiagnostics.from_payload(
            self.diagnostics,
            field_name="dataset processing state total_protein_correction.diagnostics",
        )
        object.__setattr__(self, "diagnostics", normalized)


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
class DatasetProcessingState:
    """Compact summary of preprocessing state at the analysis-ready boundary."""

    intensity_scale: IntensityScaleState
    site_sequence_resolution: SiteSequenceResolutionState
    missing_data: MissingDataState
    normalisation: NormalisationState
    total_protein_correction: TotalProteinCorrectionState
    site_matrix: SiteMatrixState
    comparisons: ComparisonState


def _require_mapping(value: object, *, field_name: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise PhosPyInputError(f"{field_name} must be an object")
    return value


def _require_string_keys(value: Mapping[str, object], *, field_name: str) -> None:
    for key in value:
        if not isinstance(key, str):
            raise PhosPyInputError(
                f"{field_name} must contain only string keys; got key "
                f"{key!r} ({type(key).__name__})"
            )


def _require_optional_str(value: object, *, field_name: str) -> str | None:
    if value is None:
        return None
    if not isinstance(value, str):
        raise PhosPyInputError(f"{field_name} must be a string")
    normalized = value.strip()
    if not normalized:
        raise PhosPyInputError(f"{field_name} must be a non-empty string")
    return normalized


def _require_required_str(value: object, *, field_name: str) -> str:
    if value is None:
        raise PhosPyInputError(f"{field_name} is required")
    parsed = _require_optional_str(value, field_name=field_name)
    if parsed is None:  # pragma: no cover - defensive guard
        raise PhosPyInputError(f"{field_name} is required")
    return parsed


def _require_optional_bool(value: object, *, field_name: str) -> bool | None:
    if value is None:
        return None
    if not isinstance(value, bool):
        raise PhosPyInputError(f"{field_name} must be a bool")
    return value


def _require_int(value: object, *, field_name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise PhosPyInputError(f"{field_name} must be an int")
    return value


def _require_optional_int(value: object, *, field_name: str) -> int | None:
    if value is None:
        return None
    return _require_int(value, field_name=field_name)


def _require_optional_non_negative_int(value: object, *, field_name: str) -> int | None:
    if value is None:
        return None
    parsed = _require_int(value, field_name=field_name)
    if parsed < 0:
        raise PhosPyInputError(f"{field_name} must be >= 0")
    return parsed


def _require_required_non_negative_int(value: object, *, field_name: str) -> int:
    if value is None:
        raise PhosPyInputError(f"{field_name} is required")
    parsed = _require_optional_non_negative_int(value, field_name=field_name)
    if parsed is None:  # pragma: no cover - defensive guard
        raise PhosPyInputError(f"{field_name} is required")
    return parsed


def _require_optional_string_tuple(
    value: object, *, field_name: str
) -> tuple[str, ...] | None:
    if value is None:
        return None
    if not isinstance(value, (list, tuple)):
        raise PhosPyInputError(f"{field_name} must be an array of strings")
    parsed: list[str] = []
    for position, item in enumerate(value):
        parsed.append(
            _require_required_str(
                item,
                field_name=f"{field_name}[{position}]",
            )
        )
    return tuple(parsed)


def _require_required_string_tuple(
    value: object, *, field_name: str
) -> tuple[str, ...]:
    if value is None:
        raise PhosPyInputError(f"{field_name} is required")
    parsed = _require_optional_string_tuple(value, field_name=field_name)
    if parsed is None:  # pragma: no cover - defensive guard
        raise PhosPyInputError(f"{field_name} is required")
    return parsed


def _require_optional_string_to_string_mapping(
    value: object,
    *,
    field_name: str,
) -> dict[str, str] | None:
    if value is None:
        return None
    if not isinstance(value, Mapping):
        raise PhosPyInputError(f"{field_name} must be an object of string mappings")
    parsed: dict[str, str] = {}
    for raw_key, raw_value in value.items():
        key = _require_required_str(
            raw_key,
            field_name=f"{field_name}.<key>",
        )
        parsed[key] = _require_required_str(
            raw_value,
            field_name=f"{field_name}.{key}",
        )
    return parsed


def _require_json_value(value: object, *, field_name: str) -> JsonValue:
    if value is None or isinstance(value, (str, bool, int, float)):
        return value
    if isinstance(value, list):
        return [
            _require_json_value(item, field_name=f"{field_name}[]") for item in value
        ]
    if isinstance(value, tuple):
        return [
            _require_json_value(item, field_name=f"{field_name}[]") for item in value
        ]
    if isinstance(value, Mapping):
        normalized: dict[str, JsonValue] = {}
        for raw_key, raw_value in value.items():
            key = _require_required_str(raw_key, field_name=f"{field_name}.<key>")
            normalized[key] = _require_json_value(
                raw_value, field_name=f"{field_name}.{key}"
            )
        return normalized
    raise PhosPyInputError(
        f"{field_name} must be JSON-compatible (null, bool, int, float, string, array, or object)"
    )


def _require_json_mapping(
    value: object,
    *,
    field_name: str,
) -> dict[str, JsonValue]:
    mapping = _require_mapping(value, field_name=field_name)
    normalized: dict[str, JsonValue] = {}
    for raw_key, raw_value in mapping.items():
        key = _require_required_str(raw_key, field_name=f"{field_name}.<key>")
        normalized[key] = _require_json_value(
            raw_value, field_name=f"{field_name}.{key}"
        )
    return normalized


def _set_optional_payload_value(
    payload: dict[str, JsonValue],
    key: str,
    value: JsonValue,
) -> None:
    if value is not None:
        payload[key] = value
