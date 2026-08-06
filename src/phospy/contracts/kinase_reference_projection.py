"""Typed kinase reference-projection provenance contract."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

from phospy.errors.workflows import WorkflowBoundaryError

KINASE_REFERENCE_PROJECTION_INTERPRETER_VERSION = (
    "phospy.workflows.kinase.reference_projector.v1"
)
KINASE_REFERENCE_PROJECTION_SCHEMA_VERSION = 1
KINASE_REFERENCE_PROJECTION_UNMATCHED_IDENTIFIER_EXAMPLE_LIMIT = 5
KINASE_REFERENCE_PROJECTION_SOURCE_IDENTIFIER_NAMESPACE = (
    "references.kinase_substrate_map.substrate_site"
)
KINASE_REFERENCE_PROJECTION_OUTPUT_IDENTIFIER_NAMESPACE = "dataset.site_key"
KINASE_REFERENCE_PROJECTION_SOURCE_IDENTITY_SEMANTICS = (
    "source reference substrate identifiers before dataset projection; values may "
    "be dataset site_key values, dataset display_id values, or unmatched reference "
    "substrate identifiers"
)
KINASE_REFERENCE_PROJECTION_OUTPUT_IDENTITY_SEMANTICS = (
    "dataset site_key rows produced by projecting source reference substrate "
    "identifiers through dataset site_key/display_id identity"
)
KINASE_REFERENCE_PROJECTION_SOURCE_IDENTIFIER_KIND_DATASET_SITE_KEY = "dataset_site_key"
KINASE_REFERENCE_PROJECTION_SOURCE_IDENTIFIER_KIND_DATASET_DISPLAY_ID = (
    "dataset_display_id"
)
KINASE_REFERENCE_PROJECTION_SOURCE_IDENTIFIER_KIND_UNMATCHED_REFERENCE_SUBSTRATE = (
    "unmatched_reference_substrate_identifier"
)
KINASE_REFERENCE_PROJECTION_SUPPORTED_SOURCE_IDENTIFIER_KINDS = frozenset(
    {
        KINASE_REFERENCE_PROJECTION_SOURCE_IDENTIFIER_KIND_DATASET_SITE_KEY,
        KINASE_REFERENCE_PROJECTION_SOURCE_IDENTIFIER_KIND_DATASET_DISPLAY_ID,
        (
            KINASE_REFERENCE_PROJECTION_SOURCE_IDENTIFIER_KIND_UNMATCHED_REFERENCE_SUBSTRATE
        ),
    }
)
KINASE_REFERENCE_PROJECTION_ONE_TO_MANY_DIAGNOSTICS = (
    "display_reference_matching.one_to_many_display_reference_matches"
)

_MATCHED_SOURCE_IDENTIFIER_KINDS = frozenset(
    {
        KINASE_REFERENCE_PROJECTION_SOURCE_IDENTIFIER_KIND_DATASET_SITE_KEY,
        KINASE_REFERENCE_PROJECTION_SOURCE_IDENTIFIER_KIND_DATASET_DISPLAY_ID,
    }
)
_PROJECTION_SUMMARY_PAYLOAD_FIELDS = frozenset(
    {
        "schema_version",
        "source_reference_row_count",
        "unique_source_substrate_identifier_count",
        "matched_source_substrate_identifier_count",
        "unmatched_source_substrate_identifier_count",
        "matched_source_substrate_identifiers",
        "unmatched_source_substrate_identifiers",
        "unmatched_source_substrate_identifier_examples",
        "unmatched_identifier_example_limit",
        "projected_dataset_site_key_count",
        "source_identifier_namespace",
        "output_identifier_namespace",
        "source_identifier_kinds",
        "source_identity_semantics",
        "output_identity_semantics",
        "one_to_many_display_reference_match_count",
        "one_to_many_display_reference_site_key_rows",
        "one_to_many_projection_diagnostics",
        "interpreter_version",
    }
)


@dataclass(frozen=True, slots=True)
class KinaseReferenceProjectionSummary:
    """Canonical reference-substrate projection summary before unmatched rows vanish."""

    source_reference_row_count: int
    matched_source_substrate_identifiers: tuple[str, ...]
    unmatched_source_substrate_identifiers: tuple[str, ...]
    projected_dataset_site_key_count: int
    source_identifier_kinds: tuple[str, ...]
    one_to_many_display_reference_match_count: int
    one_to_many_display_reference_site_key_rows: int
    source_identifier_namespace: str = (
        KINASE_REFERENCE_PROJECTION_SOURCE_IDENTIFIER_NAMESPACE
    )
    output_identifier_namespace: str = (
        KINASE_REFERENCE_PROJECTION_OUTPUT_IDENTIFIER_NAMESPACE
    )
    source_identity_semantics: str = (
        KINASE_REFERENCE_PROJECTION_SOURCE_IDENTITY_SEMANTICS
    )
    output_identity_semantics: str = (
        KINASE_REFERENCE_PROJECTION_OUTPUT_IDENTITY_SEMANTICS
    )
    unmatched_identifier_example_limit: int = (
        KINASE_REFERENCE_PROJECTION_UNMATCHED_IDENTIFIER_EXAMPLE_LIMIT
    )
    schema_version: int = KINASE_REFERENCE_PROJECTION_SCHEMA_VERSION
    interpreter_version: str = KINASE_REFERENCE_PROJECTION_INTERPRETER_VERSION

    def __post_init__(self) -> None:
        source_reference_row_count = _require_non_negative_int(
            self.source_reference_row_count,
            field_name="kinase_reference_projection.source_reference_row_count",
        )
        projected_dataset_site_key_count = _require_non_negative_int(
            self.projected_dataset_site_key_count,
            field_name="kinase_reference_projection.projected_dataset_site_key_count",
        )
        one_to_many_display_reference_match_count = _require_non_negative_int(
            self.one_to_many_display_reference_match_count,
            field_name=(
                "kinase_reference_projection.one_to_many_display_reference_match_count"
            ),
        )
        one_to_many_display_reference_site_key_rows = _require_non_negative_int(
            self.one_to_many_display_reference_site_key_rows,
            field_name=(
                "kinase_reference_projection."
                "one_to_many_display_reference_site_key_rows"
            ),
        )
        example_limit = _require_non_negative_int(
            self.unmatched_identifier_example_limit,
            field_name=(
                "kinase_reference_projection.unmatched_identifier_example_limit"
            ),
        )
        schema_version = _require_supported_int(
            self.schema_version,
            field_name="kinase_reference_projection.schema_version",
            expected=KINASE_REFERENCE_PROJECTION_SCHEMA_VERSION,
        )
        matched = _unique_sorted_text_tuple(
            self.matched_source_substrate_identifiers,
            field_name=(
                "kinase_reference_projection.matched_source_substrate_identifiers"
            ),
        )
        unmatched = _unique_sorted_text_tuple(
            self.unmatched_source_substrate_identifiers,
            field_name=(
                "kinase_reference_projection.unmatched_source_substrate_identifiers"
            ),
        )
        overlap = set(matched).intersection(unmatched)
        if overlap:
            examples = sorted(overlap)[
                :KINASE_REFERENCE_PROJECTION_UNMATCHED_IDENTIFIER_EXAMPLE_LIMIT
            ]
            raise WorkflowBoundaryError(
                "kinase reference projection summary matched and unmatched "
                f"identifier sets must be disjoint; overlap_examples={examples}"
            )
        source_identifier_kinds = _validate_source_identifier_kinds(
            self.source_identifier_kinds,
            field_name="kinase_reference_projection.source_identifier_kinds",
            matched_source_identifier_count=len(matched),
            unmatched_source_identifier_count=len(unmatched),
        )
        _validate_projection_summary_count_invariants(
            source_reference_row_count=source_reference_row_count,
            matched_source_identifier_count=len(matched),
            unmatched_source_identifier_count=len(unmatched),
            projected_dataset_site_key_count=projected_dataset_site_key_count,
            one_to_many_display_reference_match_count=(
                one_to_many_display_reference_match_count
            ),
            one_to_many_display_reference_site_key_rows=(
                one_to_many_display_reference_site_key_rows
            ),
        )
        object.__setattr__(
            self,
            "source_reference_row_count",
            source_reference_row_count,
        )
        object.__setattr__(self, "matched_source_substrate_identifiers", matched)
        object.__setattr__(self, "unmatched_source_substrate_identifiers", unmatched)
        object.__setattr__(
            self,
            "projected_dataset_site_key_count",
            projected_dataset_site_key_count,
        )
        object.__setattr__(self, "source_identifier_kinds", source_identifier_kinds)
        object.__setattr__(
            self,
            "one_to_many_display_reference_match_count",
            one_to_many_display_reference_match_count,
        )
        object.__setattr__(
            self,
            "one_to_many_display_reference_site_key_rows",
            one_to_many_display_reference_site_key_rows,
        )
        object.__setattr__(
            self,
            "source_identifier_namespace",
            _require_supported_text(
                self.source_identifier_namespace,
                field_name="kinase_reference_projection.source_identifier_namespace",
                expected=KINASE_REFERENCE_PROJECTION_SOURCE_IDENTIFIER_NAMESPACE,
            ),
        )
        object.__setattr__(
            self,
            "output_identifier_namespace",
            _require_supported_text(
                self.output_identifier_namespace,
                field_name="kinase_reference_projection.output_identifier_namespace",
                expected=KINASE_REFERENCE_PROJECTION_OUTPUT_IDENTIFIER_NAMESPACE,
            ),
        )
        object.__setattr__(
            self,
            "source_identity_semantics",
            _require_supported_text(
                self.source_identity_semantics,
                field_name="kinase_reference_projection.source_identity_semantics",
                expected=KINASE_REFERENCE_PROJECTION_SOURCE_IDENTITY_SEMANTICS,
            ),
        )
        object.__setattr__(
            self,
            "output_identity_semantics",
            _require_supported_text(
                self.output_identity_semantics,
                field_name="kinase_reference_projection.output_identity_semantics",
                expected=KINASE_REFERENCE_PROJECTION_OUTPUT_IDENTITY_SEMANTICS,
            ),
        )
        object.__setattr__(
            self,
            "unmatched_identifier_example_limit",
            _require_supported_int(
                example_limit,
                field_name=(
                    "kinase_reference_projection.unmatched_identifier_example_limit"
                ),
                expected=KINASE_REFERENCE_PROJECTION_UNMATCHED_IDENTIFIER_EXAMPLE_LIMIT,
            ),
        )
        object.__setattr__(self, "schema_version", schema_version)
        object.__setattr__(
            self,
            "interpreter_version",
            _require_supported_text(
                self.interpreter_version,
                field_name="kinase_reference_projection.interpreter_version",
                expected=KINASE_REFERENCE_PROJECTION_INTERPRETER_VERSION,
            ),
        )

    @property
    def unique_source_substrate_identifier_count(self) -> int:
        return int(
            len(self.matched_source_substrate_identifiers)
            + len(self.unmatched_source_substrate_identifiers)
        )

    @property
    def matched_source_substrate_identifier_count(self) -> int:
        return len(self.matched_source_substrate_identifiers)

    @property
    def unmatched_source_substrate_identifier_count(self) -> int:
        return len(self.unmatched_source_substrate_identifiers)

    @property
    def unmatched_source_substrate_identifier_examples(self) -> tuple[str, ...]:
        return self.unmatched_source_substrate_identifiers[
            : self.unmatched_identifier_example_limit
        ]

    def to_payload(self) -> dict[str, object]:
        return {
            "schema_version": int(self.schema_version),
            "source_reference_row_count": int(self.source_reference_row_count),
            "unique_source_substrate_identifier_count": int(
                self.unique_source_substrate_identifier_count
            ),
            "matched_source_substrate_identifier_count": int(
                self.matched_source_substrate_identifier_count
            ),
            "unmatched_source_substrate_identifier_count": int(
                self.unmatched_source_substrate_identifier_count
            ),
            "matched_source_substrate_identifiers": list(
                self.matched_source_substrate_identifiers
            ),
            "unmatched_source_substrate_identifiers": list(
                self.unmatched_source_substrate_identifiers
            ),
            "unmatched_source_substrate_identifier_examples": list(
                self.unmatched_source_substrate_identifier_examples
            ),
            "unmatched_identifier_example_limit": int(
                self.unmatched_identifier_example_limit
            ),
            "projected_dataset_site_key_count": int(
                self.projected_dataset_site_key_count
            ),
            "source_identifier_namespace": self.source_identifier_namespace,
            "output_identifier_namespace": self.output_identifier_namespace,
            "source_identifier_kinds": list(self.source_identifier_kinds),
            "source_identity_semantics": self.source_identity_semantics,
            "output_identity_semantics": self.output_identity_semantics,
            "one_to_many_display_reference_match_count": int(
                self.one_to_many_display_reference_match_count
            ),
            "one_to_many_display_reference_site_key_rows": int(
                self.one_to_many_display_reference_site_key_rows
            ),
            "one_to_many_projection_diagnostics": (
                KINASE_REFERENCE_PROJECTION_ONE_TO_MANY_DIAGNOSTICS
            ),
            "interpreter_version": self.interpreter_version,
        }

    @classmethod
    def from_payload(
        cls,
        payload: Mapping[str, object],
    ) -> KinaseReferenceProjectionSummary:
        payload = _require_payload_mapping(payload)
        schema_version = _payload_int(payload, "schema_version")
        if schema_version != KINASE_REFERENCE_PROJECTION_SCHEMA_VERSION:
            raise WorkflowBoundaryError(
                "kinase_reference_projection.schema_version is unsupported; "
                f"expected={KINASE_REFERENCE_PROJECTION_SCHEMA_VERSION!r}, "
                f"got={schema_version!r}"
            )
        _reject_unsupported_payload_fields(payload)
        _require_payload_fields(payload)
        _require_payload_exact_text(
            payload,
            "source_identifier_namespace",
            expected=KINASE_REFERENCE_PROJECTION_SOURCE_IDENTIFIER_NAMESPACE,
        )
        _require_payload_exact_text(
            payload,
            "output_identifier_namespace",
            expected=KINASE_REFERENCE_PROJECTION_OUTPUT_IDENTIFIER_NAMESPACE,
        )
        _require_payload_exact_text(
            payload,
            "source_identity_semantics",
            expected=KINASE_REFERENCE_PROJECTION_SOURCE_IDENTITY_SEMANTICS,
        )
        _require_payload_exact_text(
            payload,
            "output_identity_semantics",
            expected=KINASE_REFERENCE_PROJECTION_OUTPUT_IDENTITY_SEMANTICS,
        )
        _require_payload_exact_int(
            payload,
            "unmatched_identifier_example_limit",
            expected=KINASE_REFERENCE_PROJECTION_UNMATCHED_IDENTIFIER_EXAMPLE_LIMIT,
        )
        _require_payload_exact_text(
            payload,
            "one_to_many_projection_diagnostics",
            expected=KINASE_REFERENCE_PROJECTION_ONE_TO_MANY_DIAGNOSTICS,
        )
        _require_payload_exact_text(
            payload,
            "interpreter_version",
            expected=KINASE_REFERENCE_PROJECTION_INTERPRETER_VERSION,
        )
        summary = cls(
            source_reference_row_count=_payload_int(
                payload,
                "source_reference_row_count",
            ),
            matched_source_substrate_identifiers=_payload_text_tuple(
                payload,
                "matched_source_substrate_identifiers",
            ),
            unmatched_source_substrate_identifiers=_payload_text_tuple(
                payload,
                "unmatched_source_substrate_identifiers",
            ),
            projected_dataset_site_key_count=_payload_int(
                payload,
                "projected_dataset_site_key_count",
            ),
            source_identifier_kinds=_payload_text_tuple(
                payload,
                "source_identifier_kinds",
            ),
            one_to_many_display_reference_match_count=_payload_int(
                payload,
                "one_to_many_display_reference_match_count",
            ),
            one_to_many_display_reference_site_key_rows=_payload_int(
                payload,
                "one_to_many_display_reference_site_key_rows",
            ),
            source_identifier_namespace=_payload_text(
                payload,
                "source_identifier_namespace",
            ),
            output_identifier_namespace=_payload_text(
                payload,
                "output_identifier_namespace",
            ),
            source_identity_semantics=_payload_text(
                payload,
                "source_identity_semantics",
            ),
            output_identity_semantics=_payload_text(
                payload,
                "output_identity_semantics",
            ),
            unmatched_identifier_example_limit=_payload_int(
                payload,
                "unmatched_identifier_example_limit",
            ),
            schema_version=_payload_int(payload, "schema_version"),
            interpreter_version=_payload_text(payload, "interpreter_version"),
        )
        _require_payload_matches_summary(payload, summary)
        return summary


def _unique_sorted_text_tuple(
    values: tuple[str, ...],
    *,
    field_name: str,
) -> tuple[str, ...]:
    return tuple(
        sorted(
            {_require_non_empty_text(value, field_name=field_name) for value in values}
        )
    )


def _validate_source_identifier_kinds(
    values: tuple[str, ...],
    *,
    field_name: str,
    matched_source_identifier_count: int,
    unmatched_source_identifier_count: int,
) -> tuple[str, ...]:
    kinds = _unique_sorted_text_tuple(values, field_name=field_name)
    unsupported = [
        value
        for value in kinds
        if value not in KINASE_REFERENCE_PROJECTION_SUPPORTED_SOURCE_IDENTIFIER_KINDS
    ]
    if unsupported:
        supported = ", ".join(
            sorted(KINASE_REFERENCE_PROJECTION_SUPPORTED_SOURCE_IDENTIFIER_KINDS)
        )
        raise WorkflowBoundaryError(
            f"{field_name} contains unsupported identifier kind(s): "
            f"{unsupported}; supported={supported}"
        )
    kind_set = set(kinds)
    has_matched_kind = bool(kind_set.intersection(_MATCHED_SOURCE_IDENTIFIER_KINDS))
    has_unmatched_kind = (
        KINASE_REFERENCE_PROJECTION_SOURCE_IDENTIFIER_KIND_UNMATCHED_REFERENCE_SUBSTRATE
        in kind_set
    )
    if unmatched_source_identifier_count > 0 and not has_unmatched_kind:
        raise WorkflowBoundaryError(
            f"{field_name} must include "
            f"{KINASE_REFERENCE_PROJECTION_SOURCE_IDENTIFIER_KIND_UNMATCHED_REFERENCE_SUBSTRATE!r} "
            "when unmatched source identifiers are recorded"
        )
    if unmatched_source_identifier_count == 0 and has_unmatched_kind:
        raise WorkflowBoundaryError(
            f"{field_name} must not include "
            f"{KINASE_REFERENCE_PROJECTION_SOURCE_IDENTIFIER_KIND_UNMATCHED_REFERENCE_SUBSTRATE!r} "
            "when no unmatched source identifiers are recorded"
        )
    if matched_source_identifier_count > 0 and not has_matched_kind:
        supported_matched = ", ".join(sorted(_MATCHED_SOURCE_IDENTIFIER_KINDS))
        raise WorkflowBoundaryError(
            f"{field_name} must include at least one matched identifier kind when "
            f"matched source identifiers are recorded; supported_matched={supported_matched}"
        )
    if matched_source_identifier_count == 0 and has_matched_kind:
        raise WorkflowBoundaryError(
            f"{field_name} must not include matched identifier kinds when no "
            "matched source identifiers are recorded"
        )
    if matched_source_identifier_count == 0 and unmatched_source_identifier_count == 0:
        if kinds:
            raise WorkflowBoundaryError(
                f"{field_name} must be empty when no source identifiers are recorded"
            )
    return kinds


def _validate_projection_summary_count_invariants(
    *,
    source_reference_row_count: int,
    matched_source_identifier_count: int,
    unmatched_source_identifier_count: int,
    projected_dataset_site_key_count: int,
    one_to_many_display_reference_match_count: int,
    one_to_many_display_reference_site_key_rows: int,
) -> None:
    unique_source_identifier_count = (
        matched_source_identifier_count + unmatched_source_identifier_count
    )
    if source_reference_row_count < unique_source_identifier_count:
        raise WorkflowBoundaryError(
            "kinase_reference_projection.source_reference_row_count must be "
            "greater than or equal to the unique source substrate-identifier count"
        )
    if unique_source_identifier_count > 0 and source_reference_row_count == 0:
        raise WorkflowBoundaryError(
            "kinase_reference_projection.source_reference_row_count must be "
            "non-zero when unique source substrate identifiers are recorded"
        )
    if matched_source_identifier_count > 0 and projected_dataset_site_key_count == 0:
        raise WorkflowBoundaryError(
            "kinase_reference_projection.projected_dataset_site_key_count must be "
            "non-zero when matched source identifiers are recorded"
        )
    if matched_source_identifier_count == 0 and projected_dataset_site_key_count != 0:
        raise WorkflowBoundaryError(
            "kinase_reference_projection.projected_dataset_site_key_count must be "
            "zero when no matched source identifiers are recorded"
        )
    if (
        one_to_many_display_reference_match_count == 0
        and one_to_many_display_reference_site_key_rows != 0
    ) or (
        one_to_many_display_reference_match_count != 0
        and one_to_many_display_reference_site_key_rows == 0
    ):
        raise WorkflowBoundaryError(
            "kinase_reference_projection one-to-many match count and projected "
            "site_key row count must either both be zero or both be non-zero"
        )
    if (
        one_to_many_display_reference_match_count > 0
        and one_to_many_display_reference_site_key_rows
        < one_to_many_display_reference_match_count * 2
    ):
        raise WorkflowBoundaryError(
            "kinase_reference_projection each one-to-many display reference match "
            "must represent at least two projected dataset site_key rows"
        )
    if one_to_many_display_reference_match_count > matched_source_identifier_count:
        raise WorkflowBoundaryError(
            "kinase_reference_projection.one_to_many_display_reference_match_count "
            "cannot exceed the matched source substrate-identifier count"
        )
    if one_to_many_display_reference_site_key_rows > projected_dataset_site_key_count:
        raise WorkflowBoundaryError(
            "kinase_reference_projection.one_to_many_display_reference_site_key_rows "
            "cannot exceed projected_dataset_site_key_count"
        )


def _require_non_empty_text(value: object, *, field_name: str) -> str:
    text = str(value).strip()
    if not text:
        raise WorkflowBoundaryError(f"{field_name} must be a non-empty string")
    return text


def _require_non_negative_int(value: object, *, field_name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise WorkflowBoundaryError(f"{field_name} must be a non-negative integer")
    normalized = int(value)
    if normalized < 0:
        raise WorkflowBoundaryError(f"{field_name} must be a non-negative integer")
    return normalized


def _require_supported_text(value: object, *, field_name: str, expected: str) -> str:
    text = _require_non_empty_text(value, field_name=field_name)
    if text != expected:
        raise WorkflowBoundaryError(
            f"{field_name} has unsupported value; expected={expected!r}, got={text!r}"
        )
    return expected


def _require_supported_int(value: object, *, field_name: str, expected: int) -> int:
    observed = _require_non_negative_int(value, field_name=field_name)
    if observed != expected:
        raise WorkflowBoundaryError(
            f"{field_name} has unsupported value; expected={expected!r}, "
            f"got={observed!r}"
        )
    return expected


def _require_payload_mapping(value: object) -> Mapping[str, object]:
    if isinstance(value, Mapping):
        return value
    raise WorkflowBoundaryError("kinase reference projection payload must be a mapping")


def _require_payload_fields(payload: Mapping[str, object]) -> None:
    missing = sorted(
        key for key in _PROJECTION_SUMMARY_PAYLOAD_FIELDS if key not in payload
    )
    if missing:
        raise WorkflowBoundaryError(
            "kinase reference projection payload is missing required key(s): "
            + ", ".join(missing)
        )


def _reject_unsupported_payload_fields(payload: Mapping[str, object]) -> None:
    unsupported = sorted(
        str(key)
        for key in payload.keys()
        if str(key) not in _PROJECTION_SUMMARY_PAYLOAD_FIELDS
    )
    if unsupported:
        raise WorkflowBoundaryError(
            "kinase reference projection payload contains unsupported key(s): "
            + ", ".join(unsupported)
        )


def _payload_text(payload: Mapping[str, object], key: str) -> str:
    if key not in payload:
        raise WorkflowBoundaryError(
            f"kinase reference projection payload is missing required key: {key}"
        )
    return _require_non_empty_text(
        payload[key],
        field_name=f"kinase_reference_projection.{key}",
    )


def _payload_int(payload: Mapping[str, object], key: str) -> int:
    if key not in payload:
        raise WorkflowBoundaryError(
            f"kinase reference projection payload is missing required key: {key}"
        )
    return _require_non_negative_int(
        payload[key],
        field_name=f"kinase_reference_projection.{key}",
    )


def _require_payload_exact_text(
    payload: Mapping[str, object],
    key: str,
    *,
    expected: str,
) -> None:
    observed = _payload_text(payload, key)
    if observed != expected:
        raise WorkflowBoundaryError(
            f"kinase_reference_projection.{key} has unsupported value; "
            f"expected={expected!r}, got={observed!r}"
        )


def _require_payload_exact_int(
    payload: Mapping[str, object],
    key: str,
    *,
    expected: int,
) -> None:
    observed = _payload_int(payload, key)
    if observed != expected:
        raise WorkflowBoundaryError(
            f"kinase_reference_projection.{key} has unsupported value; "
            f"expected={expected!r}, got={observed!r}"
        )


def _payload_text_tuple(payload: Mapping[str, object], key: str) -> tuple[str, ...]:
    if key not in payload:
        raise WorkflowBoundaryError(
            f"kinase reference projection payload is missing required key: {key}"
        )
    value = payload[key]
    if isinstance(value, str) or not isinstance(value, (list, tuple)):
        raise WorkflowBoundaryError(
            f"kinase_reference_projection.{key} must be a sequence of strings"
        )
    return tuple(
        _require_non_empty_text(
            item,
            field_name=f"kinase_reference_projection.{key}[]",
        )
        for item in value
    )


def _require_payload_matches_summary(
    payload: Mapping[str, object],
    summary: KinaseReferenceProjectionSummary,
) -> None:
    canonical_payload = summary.to_payload()
    for key in sorted(_PROJECTION_SUMMARY_PAYLOAD_FIELDS):
        observed = payload.get(key)
        expected = canonical_payload[key]
        if observed != expected:
            raise WorkflowBoundaryError(
                "kinase reference projection payload field disagrees with the "
                f"canonical projection summary: {key}; expected={expected!r}, "
                f"got={observed!r}"
            )


__all__ = [
    "KINASE_REFERENCE_PROJECTION_INTERPRETER_VERSION",
    "KINASE_REFERENCE_PROJECTION_ONE_TO_MANY_DIAGNOSTICS",
    "KINASE_REFERENCE_PROJECTION_OUTPUT_IDENTIFIER_NAMESPACE",
    "KINASE_REFERENCE_PROJECTION_OUTPUT_IDENTITY_SEMANTICS",
    "KINASE_REFERENCE_PROJECTION_SCHEMA_VERSION",
    "KINASE_REFERENCE_PROJECTION_SOURCE_IDENTIFIER_KIND_DATASET_DISPLAY_ID",
    "KINASE_REFERENCE_PROJECTION_SOURCE_IDENTIFIER_KIND_DATASET_SITE_KEY",
    "KINASE_REFERENCE_PROJECTION_SOURCE_IDENTIFIER_KIND_UNMATCHED_REFERENCE_SUBSTRATE",
    "KINASE_REFERENCE_PROJECTION_SOURCE_IDENTIFIER_NAMESPACE",
    "KINASE_REFERENCE_PROJECTION_SOURCE_IDENTITY_SEMANTICS",
    "KINASE_REFERENCE_PROJECTION_SUPPORTED_SOURCE_IDENTIFIER_KINDS",
    "KINASE_REFERENCE_PROJECTION_UNMATCHED_IDENTIFIER_EXAMPLE_LIMIT",
    "KinaseReferenceProjectionSummary",
]
