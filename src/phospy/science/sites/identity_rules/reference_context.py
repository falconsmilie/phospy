"""Reference-context compatibility rules for phosphosite identity."""

from __future__ import annotations

from phospy.errors.validation import PhosPyValidationError
from phospy.provenance.models import ReferenceContextProtocol
from phospy.science.sites.identity_rules.contracts import (
    REFERENCE_CONTEXT_IDENTITY_FIELDS,
    ReferenceContextCompatibilityWarning,
)


def validate_reference_context_compatibility(
    left: ReferenceContextProtocol | None,
    right: ReferenceContextProtocol | None,
    *,
    operation: str,
    allow_unknown: bool = False,
    error_type: type[Exception] = PhosPyValidationError,
) -> ReferenceContextCompatibilityWarning | None:
    """Require two biological reference contexts to describe the same identity."""

    resolved_operation = _required_operation_text(operation)
    if left is not None and right is not None:
        if left == right:
            return None
        raise error_type(
            _reference_context_mismatch_message(
                left=left,
                right=right,
                operation=resolved_operation,
            )
        )

    missing_contexts = _missing_reference_context_sides(left=left, right=right)
    if not missing_contexts:
        return None
    if allow_unknown:
        return ReferenceContextCompatibilityWarning(
            operation=resolved_operation,
            missing_contexts=missing_contexts,
            left_reference_context_id=_reference_context_id(left),
            right_reference_context_id=_reference_context_id(right),
        )
    raise error_type(
        "reference-context compatibility failed for "
        f"operation={resolved_operation!r}: unknown reference context for "
        f"{', '.join(missing_contexts)}; configure an explicit workflow policy "
        "only when unknown reference context is scientifically acceptable"
    )


def _required_operation_text(value: object) -> str:
    if not isinstance(value, str) or value.strip() == "":
        raise PhosPyValidationError(
            "reference-context compatibility operation must be a non-empty string"
        )
    return value.strip()


def _missing_reference_context_sides(
    *,
    left: ReferenceContextProtocol | None,
    right: ReferenceContextProtocol | None,
) -> tuple[str, ...]:
    missing: list[str] = []
    if left is None:
        missing.append("left")
    if right is None:
        missing.append("right")
    return tuple(missing)


def _reference_context_id(context: ReferenceContextProtocol | None) -> str | None:
    return None if context is None else context.reference_context_id


def _reference_context_mismatch_message(
    *,
    left: ReferenceContextProtocol,
    right: ReferenceContextProtocol,
    operation: str,
) -> str:
    mismatched_fields = [
        field_name
        for field_name in REFERENCE_CONTEXT_IDENTITY_FIELDS
        if getattr(left, field_name) != getattr(right, field_name)
    ]
    field_text = ", ".join(mismatched_fields) or "unknown"
    return (
        "reference-context compatibility failed for "
        f"operation={operation!r}: incompatible reference contexts; "
        f"mismatched_fields={field_text}; "
        f"left={_reference_context_summary(left)}; "
        f"right={_reference_context_summary(right)}"
    )


def _reference_context_summary(context: ReferenceContextProtocol) -> str:
    parts = [
        f"reference_context_id={context.reference_context_id!r}",
        f"organism={context.organism!r}",
        f"protein_namespace={context.protein_namespace!r}",
        f"source_name={context.source_name!r}",
        f"source_version={context.source_version!r}",
        f"proteome_version={context.proteome_version!r}",
        f"reference_table_sha256={context.reference_table_sha256!r}",
    ]
    return "{" + ", ".join(parts) + "}"


__all__ = ["validate_reference_context_compatibility"]
