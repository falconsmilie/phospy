"""Reference-context compatibility rules for reference consumers."""

from __future__ import annotations

from dataclasses import dataclass

from phospy.errors.validation import PhosPyValidationError
from phospy.provenance.result_caveats import ResultCaveat
from phospy.references.provenance.context import ReferenceContext

REFERENCE_CONTEXT_UNKNOWN_CAVEAT_CODE = "reference_context_unknown"
REFERENCE_CONTEXT_IDENTITY_FIELDS = (
    "organism",
    "protein_namespace",
    "source_name",
    "source_version",
    "proteome_version",
    "reference_table_sha256",
)


@dataclass(frozen=True, slots=True)
class ReferenceContextCompatibilityWarning:
    """Typed warning returned when unknown reference context is explicitly allowed."""

    operation: str
    missing_contexts: tuple[str, ...]
    left_reference_context_id: str | None
    right_reference_context_id: str | None
    code: str = REFERENCE_CONTEXT_UNKNOWN_CAVEAT_CODE
    severity: str = "warning"

    @property
    def message(self) -> str:
        missing = ", ".join(self.missing_contexts)
        noun = "contexts are" if len(self.missing_contexts) > 1 else "context is"
        return (
            "Reference-context compatibility could not be proven because "
            f"{missing} {noun} unknown for operation={self.operation!r}."
        )

    def to_payload(self) -> dict[str, object]:
        return {
            "code": self.code,
            "severity": self.severity,
            "message": self.message,
            "operation": self.operation,
            "missing_contexts": list(self.missing_contexts),
            "left_reference_context_id": self.left_reference_context_id,
            "right_reference_context_id": self.right_reference_context_id,
        }


def validate_reference_context_compatibility(
    left: ReferenceContext | None,
    right: ReferenceContext | None,
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


def build_reference_context_compatibility_caveat(
    warning: ReferenceContextCompatibilityWarning,
    *,
    policy: object,
    workflow_scope: str,
) -> ResultCaveat:
    """Convert a compatibility warning into a result caveat."""

    details = warning.to_payload()
    details["policy"] = str(policy)
    details["workflow_scope"] = workflow_scope
    return ResultCaveat(
        code=warning.code,
        severity="warning",
        message=warning.message,
        details=details,
    )


def _required_operation_text(value: object) -> str:
    if not isinstance(value, str) or value.strip() == "":
        raise PhosPyValidationError(
            "reference-context compatibility operation must be a non-empty string"
        )
    return value.strip()


def _missing_reference_context_sides(
    *,
    left: ReferenceContext | None,
    right: ReferenceContext | None,
) -> tuple[str, ...]:
    missing: list[str] = []
    if left is None:
        missing.append("left")
    if right is None:
        missing.append("right")
    return tuple(missing)


def _reference_context_id(context: ReferenceContext | None) -> str | None:
    return None if context is None else context.reference_context_id


def _reference_context_mismatch_message(
    *,
    left: ReferenceContext,
    right: ReferenceContext,
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


def _reference_context_summary(context: ReferenceContext) -> str:
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


__all__ = [
    "REFERENCE_CONTEXT_IDENTITY_FIELDS",
    "REFERENCE_CONTEXT_UNKNOWN_CAVEAT_CODE",
    "ReferenceContextCompatibilityWarning",
    "build_reference_context_compatibility_caveat",
    "validate_reference_context_compatibility",
]
