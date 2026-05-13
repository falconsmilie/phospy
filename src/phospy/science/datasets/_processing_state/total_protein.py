"""Total-protein correction diagnostics models and parsing contracts."""

from __future__ import annotations

from collections.abc import Iterator, Mapping
from dataclasses import dataclass, field

from phospy.errors.input import PhosPyInputError
from phospy.science.datasets.preprocessing.policy_models import (
    TotalProteinCorrectionIdentityMatchingPolicy,
    TotalProteinCorrectionPolicy,
)
from phospy.science.transformations.models import QuantitativeMeaning

from .json_contracts import (
    TOTAL_PROTEIN_CORRECTION_DIAGNOSTICS_SCHEMA_VERSION_V1,
    V1_KNOWN_TOTAL_PROTEIN_DIAGNOSTICS_FIELDS,
    JsonValue,
    require_int,
    require_mapping,
    require_optional_bool,
    require_optional_non_negative_int,
    require_optional_str,
    require_optional_string_to_string_mapping,
    require_optional_string_tuple,
    require_required_str,
    require_string_keys,
    set_optional_payload_value,
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
        mapping = require_mapping(payload, field_name=field_name)
        return TotalProteinCorrectionDiagnosticsV1.from_mapping(
            mapping,
            field_name=field_name,
        )

    def to_payload(self) -> dict[str, JsonValue]:
        """Return normalised diagnostics payload suitable for bundle JSON."""

        raise NotImplementedError


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
    identity_matching_policy: str | None = None
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
        require_string_keys(payload, field_name=field_name)
        version = require_int(
            payload.get("diagnostics_schema_version"),
            field_name=f"{field_name}.diagnostics_schema_version",
        )
        if version != TOTAL_PROTEIN_CORRECTION_DIAGNOSTICS_SCHEMA_VERSION_V1:
            raise PhosPyInputError(
                f"{field_name}.diagnostics_schema_version={version!r} is unsupported; "
                f"expected {TOTAL_PROTEIN_CORRECTION_DIAGNOSTICS_SCHEMA_VERSION_V1}"
            )
        unknown_fields = sorted(
            key
            for key in payload
            if key not in V1_KNOWN_TOTAL_PROTEIN_DIAGNOSTICS_FIELDS
        )
        if unknown_fields:
            raise PhosPyInputError(
                f"{field_name} contains unsupported field(s): "
                + ", ".join(unknown_fields)
            )
        return cls(
            policy=require_optional_str(
                payload.get("policy"),
                field_name=f"{field_name}.policy",
            ),
            requested_policy=require_optional_str(
                payload.get("requested_policy"),
                field_name=f"{field_name}.requested_policy",
            ),
            resolved_policy=require_optional_str(
                payload.get("resolved_policy"),
                field_name=f"{field_name}.resolved_policy",
            ),
            formula=require_optional_str(
                payload.get("formula"),
                field_name=f"{field_name}.formula",
            ),
            requires_log_scale=require_optional_bool(
                payload.get("requires_log_scale"),
                field_name=f"{field_name}.requires_log_scale",
            ),
            input_scale=require_optional_str(
                payload.get("input_scale"),
                field_name=f"{field_name}.input_scale",
            ),
            output_scale=require_optional_str(
                payload.get("output_scale"),
                field_name=f"{field_name}.output_scale",
            ),
            quantitative_meaning=require_required_str(
                payload.get("quantitative_meaning"),
                field_name=f"{field_name}.quantitative_meaning",
            ),
            matched_rows=require_optional_non_negative_int(
                payload.get("matched_rows"),
                field_name=f"{field_name}.matched_rows",
            ),
            identity_mode=require_optional_str(
                payload.get("identity_mode"),
                field_name=f"{field_name}.identity_mode",
            ),
            identity_matching_policy=require_optional_str(
                payload.get("identity_matching_policy"),
                field_name=f"{field_name}.identity_matching_policy",
            ),
            phosphosite_key=require_optional_str(
                payload.get("phosphosite_key"),
                field_name=f"{field_name}.phosphosite_key",
            ),
            total_protein_key=require_optional_str(
                payload.get("total_protein_key"),
                field_name=f"{field_name}.total_protein_key",
            ),
            mapping_phosphosite_key=require_optional_str(
                payload.get("mapping_phosphosite_key"),
                field_name=f"{field_name}.mapping_phosphosite_key",
            ),
            mapping_total_protein_key=require_optional_str(
                payload.get("mapping_total_protein_key"),
                field_name=f"{field_name}.mapping_total_protein_key",
            ),
            mapping_table_fingerprint=require_optional_str(
                payload.get("mapping_table_fingerprint"),
                field_name=f"{field_name}.mapping_table_fingerprint",
            ),
            duplicate_policy=require_optional_str(
                payload.get("duplicate_policy"),
                field_name=f"{field_name}.duplicate_policy",
            ),
            unmatched_policy=require_optional_str(
                payload.get("unmatched_policy"),
                field_name=f"{field_name}.unmatched_policy",
            ),
            phosphosite_row_count=require_optional_non_negative_int(
                payload.get("phosphosite_row_count"),
                field_name=f"{field_name}.phosphosite_row_count",
            ),
            total_protein_row_count=require_optional_non_negative_int(
                payload.get("total_protein_row_count"),
                field_name=f"{field_name}.total_protein_row_count",
            ),
            corrected_row_count=require_optional_non_negative_int(
                payload.get("corrected_row_count"),
                field_name=f"{field_name}.corrected_row_count",
            ),
            uncorrected_row_count=require_optional_non_negative_int(
                payload.get("uncorrected_row_count"),
                field_name=f"{field_name}.uncorrected_row_count",
            ),
            unused_total_protein_row_count=require_optional_non_negative_int(
                payload.get("unused_total_protein_row_count"),
                field_name=f"{field_name}.unused_total_protein_row_count",
            ),
            total_rows_used_by_multiple_phosphosites=require_optional_non_negative_int(
                payload.get("total_rows_used_by_multiple_phosphosites"),
                field_name=f"{field_name}.total_rows_used_by_multiple_phosphosites",
            ),
            corrected_phosphosite_row_ids=require_optional_string_tuple(
                payload.get("corrected_phosphosite_row_ids"),
                field_name=f"{field_name}.corrected_phosphosite_row_ids",
            ),
            corrected_phosphosite_to_total_protein_row_id=require_optional_string_to_string_mapping(
                payload.get("corrected_phosphosite_to_total_protein_row_id"),
                field_name=f"{field_name}.corrected_phosphosite_to_total_protein_row_id",
            ),
            unmatched_phosphosite_row_ids=require_optional_string_tuple(
                payload.get("unmatched_phosphosite_row_ids"),
                field_name=f"{field_name}.unmatched_phosphosite_row_ids",
            ),
            uncorrected_phosphosite_row_reasons=require_optional_string_to_string_mapping(
                payload.get("uncorrected_phosphosite_row_reasons"),
                field_name=f"{field_name}.uncorrected_phosphosite_row_reasons",
            ),
            unused_total_protein_row_ids=require_optional_string_tuple(
                payload.get("unused_total_protein_row_ids"),
                field_name=f"{field_name}.unused_total_protein_row_ids",
            ),
            gene_symbol_matching_used=require_optional_bool(
                payload.get("gene_symbol_matching_used"),
                field_name=f"{field_name}.gene_symbol_matching_used",
            ),
            gene_symbol_identity_warning=require_optional_str(
                payload.get("gene_symbol_identity_warning"),
                field_name=f"{field_name}.gene_symbol_identity_warning",
            ),
            total_table_hash=require_optional_str(
                payload.get("total_table_hash"),
                field_name=f"{field_name}.total_table_hash",
            ),
            input_phospho_hash=require_optional_str(
                payload.get("input_phospho_hash"),
                field_name=f"{field_name}.input_phospho_hash",
            ),
            output_phospho_hash=require_optional_str(
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
        policy = require_optional_str(
            self.policy,
            field_name="dataset processing state total_protein_correction.diagnostics.policy",
        )
        if policy is not None:
            policy = TotalProteinCorrectionPolicy.parse(
                policy,
                field_name="dataset processing state total_protein_correction.diagnostics.policy",
            ).value
        requested_policy = require_optional_str(
            self.requested_policy,
            field_name="dataset processing state total_protein_correction.diagnostics.requested_policy",
        )
        if requested_policy is not None:
            requested_policy = TotalProteinCorrectionPolicy.parse(
                requested_policy,
                field_name=(
                    "dataset processing state total_protein_correction."
                    "diagnostics.requested_policy"
                ),
            ).value
        resolved_policy = require_optional_str(
            self.resolved_policy,
            field_name="dataset processing state total_protein_correction.diagnostics.resolved_policy",
        )
        if resolved_policy is not None:
            resolved_policy = TotalProteinCorrectionPolicy.parse(
                resolved_policy,
                field_name=(
                    "dataset processing state total_protein_correction."
                    "diagnostics.resolved_policy"
                ),
            ).value
        formula = require_optional_str(
            self.formula,
            field_name="dataset processing state total_protein_correction.diagnostics.formula",
        )
        requires_log_scale = require_optional_bool(
            self.requires_log_scale,
            field_name="dataset processing state total_protein_correction.diagnostics.requires_log_scale",
        )
        input_scale = require_optional_str(
            self.input_scale,
            field_name="dataset processing state total_protein_correction.diagnostics.input_scale",
        )
        output_scale = require_optional_str(
            self.output_scale,
            field_name="dataset processing state total_protein_correction.diagnostics.output_scale",
        )
        quantitative_meaning_raw = require_required_str(
            self.quantitative_meaning,
            field_name=(
                "dataset processing state total_protein_correction.diagnostics."
                "quantitative_meaning"
            ),
        )
        try:
            quantitative_meaning = QuantitativeMeaning(quantitative_meaning_raw).value
        except ValueError as exc:
            supported = ", ".join(member.value for member in QuantitativeMeaning)
            raise PhosPyInputError(
                "dataset processing state total_protein_correction.diagnostics."
                "quantitative_meaning must be one of: "
                f"{supported}"
            ) from exc
        matched_rows = require_optional_non_negative_int(
            self.matched_rows,
            field_name="dataset processing state total_protein_correction.diagnostics.matched_rows",
        )
        identity_mode = require_optional_str(
            self.identity_mode,
            field_name="dataset processing state total_protein_correction.diagnostics.identity_mode",
        )
        identity_matching_policy = require_optional_str(
            self.identity_matching_policy,
            field_name=(
                "dataset processing state total_protein_correction.diagnostics."
                "identity_matching_policy"
            ),
        )
        if identity_matching_policy is not None:
            identity_matching_policy = (
                TotalProteinCorrectionIdentityMatchingPolicy.parse(
                    identity_matching_policy,
                    field_name=(
                        "dataset processing state total_protein_correction."
                        "diagnostics.identity_matching_policy"
                    ),
                ).value
            )
        phosphosite_key = require_optional_str(
            self.phosphosite_key,
            field_name="dataset processing state total_protein_correction.diagnostics.phosphosite_key",
        )
        total_protein_key = require_optional_str(
            self.total_protein_key,
            field_name="dataset processing state total_protein_correction.diagnostics.total_protein_key",
        )
        mapping_phosphosite_key = require_optional_str(
            self.mapping_phosphosite_key,
            field_name="dataset processing state total_protein_correction.diagnostics.mapping_phosphosite_key",
        )
        mapping_total_protein_key = require_optional_str(
            self.mapping_total_protein_key,
            field_name="dataset processing state total_protein_correction.diagnostics.mapping_total_protein_key",
        )
        mapping_table_fingerprint = require_optional_str(
            self.mapping_table_fingerprint,
            field_name="dataset processing state total_protein_correction.diagnostics.mapping_table_fingerprint",
        )
        duplicate_policy = require_optional_str(
            self.duplicate_policy,
            field_name="dataset processing state total_protein_correction.diagnostics.duplicate_policy",
        )
        unmatched_policy = require_optional_str(
            self.unmatched_policy,
            field_name="dataset processing state total_protein_correction.diagnostics.unmatched_policy",
        )
        phosphosite_row_count = require_optional_non_negative_int(
            self.phosphosite_row_count,
            field_name="dataset processing state total_protein_correction.diagnostics.phosphosite_row_count",
        )
        total_protein_row_count = require_optional_non_negative_int(
            self.total_protein_row_count,
            field_name="dataset processing state total_protein_correction.diagnostics.total_protein_row_count",
        )
        corrected_row_count = require_optional_non_negative_int(
            self.corrected_row_count,
            field_name="dataset processing state total_protein_correction.diagnostics.corrected_row_count",
        )
        uncorrected_row_count = require_optional_non_negative_int(
            self.uncorrected_row_count,
            field_name="dataset processing state total_protein_correction.diagnostics.uncorrected_row_count",
        )
        unused_total_protein_row_count = require_optional_non_negative_int(
            self.unused_total_protein_row_count,
            field_name=(
                "dataset processing state total_protein_correction.diagnostics."
                "unused_total_protein_row_count"
            ),
        )
        total_rows_used_by_multiple_phosphosites = require_optional_non_negative_int(
            self.total_rows_used_by_multiple_phosphosites,
            field_name=(
                "dataset processing state total_protein_correction.diagnostics."
                "total_rows_used_by_multiple_phosphosites"
            ),
        )
        corrected_phosphosite_row_ids = require_optional_string_tuple(
            self.corrected_phosphosite_row_ids,
            field_name=(
                "dataset processing state total_protein_correction.diagnostics."
                "corrected_phosphosite_row_ids"
            ),
        )
        corrected_phosphosite_to_total_protein_row_id = (
            require_optional_string_to_string_mapping(
                self.corrected_phosphosite_to_total_protein_row_id,
                field_name=(
                    "dataset processing state total_protein_correction.diagnostics."
                    "corrected_phosphosite_to_total_protein_row_id"
                ),
            )
        )
        unmatched_phosphosite_row_ids = require_optional_string_tuple(
            self.unmatched_phosphosite_row_ids,
            field_name=(
                "dataset processing state total_protein_correction.diagnostics."
                "unmatched_phosphosite_row_ids"
            ),
        )
        uncorrected_phosphosite_row_reasons = require_optional_string_to_string_mapping(
            self.uncorrected_phosphosite_row_reasons,
            field_name=(
                "dataset processing state total_protein_correction.diagnostics."
                "uncorrected_phosphosite_row_reasons"
            ),
        )
        unused_total_protein_row_ids = require_optional_string_tuple(
            self.unused_total_protein_row_ids,
            field_name=(
                "dataset processing state total_protein_correction.diagnostics."
                "unused_total_protein_row_ids"
            ),
        )
        gene_symbol_matching_used = require_optional_bool(
            self.gene_symbol_matching_used,
            field_name=(
                "dataset processing state total_protein_correction.diagnostics."
                "gene_symbol_matching_used"
            ),
        )
        gene_symbol_identity_warning = require_optional_str(
            self.gene_symbol_identity_warning,
            field_name=(
                "dataset processing state total_protein_correction.diagnostics."
                "gene_symbol_identity_warning"
            ),
        )
        total_table_hash = require_optional_str(
            self.total_table_hash,
            field_name=(
                "dataset processing state total_protein_correction.diagnostics."
                "total_table_hash"
            ),
        )
        input_phospho_hash = require_optional_str(
            self.input_phospho_hash,
            field_name=(
                "dataset processing state total_protein_correction.diagnostics."
                "input_phospho_hash"
            ),
        )
        output_phospho_hash = require_optional_str(
            self.output_phospho_hash,
            field_name=(
                "dataset processing state total_protein_correction.diagnostics."
                "output_phospho_hash"
            ),
        )

        payload: dict[str, JsonValue] = {
            "diagnostics_schema_version": self.diagnostics_schema_version
        }
        set_optional_payload_value(payload, "policy", policy)
        set_optional_payload_value(payload, "requested_policy", requested_policy)
        set_optional_payload_value(payload, "resolved_policy", resolved_policy)
        set_optional_payload_value(payload, "formula", formula)
        set_optional_payload_value(payload, "requires_log_scale", requires_log_scale)
        set_optional_payload_value(payload, "input_scale", input_scale)
        set_optional_payload_value(payload, "output_scale", output_scale)
        payload["quantitative_meaning"] = quantitative_meaning
        set_optional_payload_value(payload, "matched_rows", matched_rows)
        set_optional_payload_value(payload, "identity_mode", identity_mode)
        set_optional_payload_value(
            payload,
            "identity_matching_policy",
            identity_matching_policy,
        )
        set_optional_payload_value(payload, "phosphosite_key", phosphosite_key)
        set_optional_payload_value(payload, "total_protein_key", total_protein_key)
        set_optional_payload_value(
            payload, "mapping_phosphosite_key", mapping_phosphosite_key
        )
        set_optional_payload_value(
            payload, "mapping_total_protein_key", mapping_total_protein_key
        )
        set_optional_payload_value(
            payload, "mapping_table_fingerprint", mapping_table_fingerprint
        )
        set_optional_payload_value(payload, "duplicate_policy", duplicate_policy)
        set_optional_payload_value(payload, "unmatched_policy", unmatched_policy)
        set_optional_payload_value(
            payload, "phosphosite_row_count", phosphosite_row_count
        )
        set_optional_payload_value(
            payload, "total_protein_row_count", total_protein_row_count
        )
        set_optional_payload_value(payload, "corrected_row_count", corrected_row_count)
        set_optional_payload_value(
            payload, "uncorrected_row_count", uncorrected_row_count
        )
        set_optional_payload_value(
            payload, "unused_total_protein_row_count", unused_total_protein_row_count
        )
        set_optional_payload_value(
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
        set_optional_payload_value(
            payload, "gene_symbol_matching_used", gene_symbol_matching_used
        )
        set_optional_payload_value(
            payload, "gene_symbol_identity_warning", gene_symbol_identity_warning
        )
        set_optional_payload_value(payload, "total_table_hash", total_table_hash)
        set_optional_payload_value(payload, "input_phospho_hash", input_phospho_hash)
        set_optional_payload_value(payload, "output_phospho_hash", output_phospho_hash)

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
        object.__setattr__(self, "identity_matching_policy", identity_matching_policy)
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
