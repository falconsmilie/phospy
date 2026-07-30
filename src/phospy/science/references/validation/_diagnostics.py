"""Reference validation diagnostic formatting helpers."""

from __future__ import annotations

from dataclasses import dataclass

from phospy.science.references.errors import ReferenceManifestError
from phospy.science.references.manifest import (
    RedistributionEvidenceType,
    RedistributionStatus,
    ReferenceManifest,
)

_MISSING = object()


@dataclass(frozen=True, slots=True)
class _ReferenceDiagnosticContext:
    reference_id: object = _MISSING
    display_name: object = _MISSING
    organism: object = _MISSING
    namespace: object = _MISSING
    redistribution_status: object = _MISSING

    @classmethod
    def from_manifest(cls, manifest: ReferenceManifest) -> _ReferenceDiagnosticContext:
        return cls(
            reference_id=manifest.reference_id,
            display_name=manifest.display_name,
            organism=manifest.organism,
            namespace=manifest.protein_namespace,
            redistribution_status=manifest.redistribution_status,
        )

    @classmethod
    def from_raw_payload(
        cls,
        payload: dict[str, object],
    ) -> _ReferenceDiagnosticContext:
        return cls(
            reference_id=_raw_context_value(payload, "reference_id"),
            display_name=_raw_context_value(payload, "display_name"),
            organism=_raw_context_value(payload, "organism"),
            namespace=_raw_context_value(payload, "protein_namespace"),
            redistribution_status=_raw_context_value(
                payload,
                "redistribution_status",
            ),
        )


def _enum_value(value: object) -> object:
    if isinstance(value, (RedistributionStatus, RedistributionEvidenceType)):
        return value.value
    return value


def _raw_context_value(payload: dict[str, object], key: str) -> object:
    if key not in payload:
        return _MISSING
    return payload[key]


def _format_diagnostic_value(value: object) -> str:
    if value is _MISSING:
        return "<missing>"
    return repr(_enum_value(value))


def _format_digest_value(value: object) -> str:
    if value is _MISSING:
        return "<missing>"
    if isinstance(value, str):
        return value
    return repr(value)


def _format_release_validation_failure(
    diagnostic_context: _ReferenceDiagnosticContext,
    *,
    field: str,
    actual_value: object = _MISSING,
    reason: str,
    file_path: str | None = None,
    expected_digest: object = _MISSING,
    actual_digest: object = _MISSING,
) -> str:
    parts = [
        f"reference_id={_format_diagnostic_value(diagnostic_context.reference_id)}",
        f"display_name={_format_diagnostic_value(diagnostic_context.display_name)}",
        f"organism={_format_diagnostic_value(diagnostic_context.organism)}",
        f"namespace={_format_diagnostic_value(diagnostic_context.namespace)}",
        f"field={_format_diagnostic_value(field)}",
        "redistribution_status="
        f"{_format_diagnostic_value(diagnostic_context.redistribution_status)}",
        f"actual_value={_format_diagnostic_value(actual_value)}",
    ]
    if file_path is not None:
        parts.extend(
            [
                f"file={_format_diagnostic_value(file_path)}",
                f"expected_digest={_format_digest_value(expected_digest)}",
                f"actual_digest={_format_digest_value(actual_digest)}",
            ]
        )
    parts.append(f"reason={reason}")
    prefix = (
        "Reference manifest file validation failed: "
        if file_path is not None
        else "Reference release validation failed: "
    )
    return prefix + ", ".join(parts)


def _manifest_parse_error(
    diagnostic_context: _ReferenceDiagnosticContext,
    *,
    field: str,
    actual_value: object = _MISSING,
    reason: str,
) -> ReferenceManifestError:
    return ReferenceManifestError(
        _format_release_validation_failure(
            diagnostic_context,
            field=field,
            actual_value=actual_value,
            reason=reason,
        )
    )


def _format_release_gate_failure(
    manifest: ReferenceManifest,
    *,
    field: str,
    actual_value: object = _MISSING,
    reason: str,
) -> str:
    return _format_release_validation_failure(
        _ReferenceDiagnosticContext.from_manifest(manifest),
        field=field,
        actual_value=actual_value,
        reason=reason,
    )
