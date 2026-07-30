"""Reference manifest JSON schema parsing."""

from __future__ import annotations

import json
from datetime import date
from pathlib import Path
from typing import cast

from phospy.science.references.errors import ReferenceManifestError
from phospy.science.references.manifest import (
    RedistributionAttribution,
    RedistributionEvidence,
    RedistributionEvidenceType,
    RedistributionScope,
    RedistributionStatus,
    ReferenceFileManifest,
    ReferenceManifest,
    UpstreamPackageLicenseEvidence,
)
from phospy.science.references.validation._diagnostics import (
    _MISSING,
    _manifest_parse_error,
    _ReferenceDiagnosticContext,
)

_REQUIRED_MANIFEST_FIELDS = frozenset(
    {
        "reference_id",
        "display_name",
        "organism",
        "taxonomy_id",
        "protein_namespace",
        "reference_version",
        "source_name",
        "source_url",
        "source_version",
        "retrieved_at",
        "table_sha256",
        "license_name",
        "license_url",
        "redistribution_status",
        "redistribution_notes",
        "derived_from",
        "generated_by",
        "generated_at_utc",
        "manifest_schema_version",
        "files",
    }
)

_ALLOWED_MANIFEST_FIELDS = frozenset(
    {
        "reference_id",
        "display_name",
        "organism",
        "taxonomy_id",
        "protein_namespace",
        "reference_version",
        "source_name",
        "source_url",
        "source_version",
        "retrieved_at",
        "table_sha256",
        "source_publication",
        "license_name",
        "license_url",
        "redistribution_status",
        "redistribution_allowed",
        "redistribution_notes",
        "redistribution_evidence",
        "derived_from",
        "generated_by",
        "generated_at_utc",
        "manifest_schema_version",
        "files",
        "sequence_context_policy",
        "sequence_window_length",
        "sequence_center_index",
        "allowed_sequence_alphabet",
        "organism_common_name",
        "supports",
        "limitations",
    }
)

_REQUIRED_FILE_FIELDS = frozenset(
    {
        "relative_path",
        "role",
        "format",
        "sha256",
        "row_count",
        "column_names",
    }
)

_ALLOWED_FILE_FIELDS = frozenset(
    {
        "relative_path",
        "role",
        "format",
        "sha256",
        "row_count",
        "column_names",
    }
)

_REQUIRED_REDISTRIBUTION_EVIDENCE_FIELDS = frozenset(
    {
        "evidence_type",
        "upstream_package",
        "scope",
        "attribution",
        "independent_database_permission_claimed",
    }
)

_ALLOWED_REDISTRIBUTION_EVIDENCE_FIELDS = frozenset(
    {
        "evidence_type",
        "upstream_package",
        "scope",
        "attribution",
        "independent_database_permission_claimed",
        "evidence_url",
        "verified_at",
        "notes",
    }
)

_REQUIRED_UPSTREAM_PACKAGE_FIELDS = frozenset(
    {
        "package_name",
        "package_version",
        "license_name",
    }
)

_ALLOWED_UPSTREAM_PACKAGE_FIELDS = frozenset(
    {
        "package_name",
        "package_version",
        "license_name",
        "license_url",
    }
)

_REQUIRED_REDISTRIBUTION_SCOPE_FIELDS = frozenset(
    {
        "reference_id",
        "reference_version",
        "applies_to_exact_packaged_files",
        "packaged_files",
        "applies_to_future_bundles",
    }
)

_ALLOWED_REDISTRIBUTION_SCOPE_FIELDS = _REQUIRED_REDISTRIBUTION_SCOPE_FIELDS

_REQUIRED_REDISTRIBUTION_ATTRIBUTION_FIELDS = frozenset(
    {
        "repository_notice_path",
        "bundle_attribution_path",
    }
)

_ALLOWED_REDISTRIBUTION_ATTRIBUTION_FIELDS = _REQUIRED_REDISTRIBUTION_ATTRIBUTION_FIELDS


def load_reference_manifest(
    manifest_path: str | Path,
    *,
    bundle_root: str | Path | None = None,
    bundled: bool = False,
    require_redistribution_allowed: bool = False,
    require_all_files_listed: bool = False,
) -> ReferenceManifest:
    """Load, parse, and validate one reference manifest JSON file."""

    resolved_manifest_path = Path(manifest_path)
    try:
        payload = json.loads(resolved_manifest_path.read_text(encoding="utf-8"))
    except FileNotFoundError as exc:
        raise ReferenceManifestError(
            f"reference manifest does not exist: {resolved_manifest_path}"
        ) from exc
    except json.JSONDecodeError as exc:
        raise ReferenceManifestError(
            f"reference manifest is not valid JSON: {resolved_manifest_path}"
        ) from exc
    if not isinstance(payload, dict):
        raise ReferenceManifestError(
            f"reference manifest must decode to an object: {resolved_manifest_path}"
        )
    root = (
        Path(bundle_root) if bundle_root is not None else resolved_manifest_path.parent
    )
    manifest = parse_reference_manifest_payload(
        cast(dict[str, object], payload),
        context=str(resolved_manifest_path),
    )
    from phospy.science.references.validation.bundle_semantics import (
        validate_reference_manifest,
    )

    return validate_reference_manifest(
        manifest,
        bundle_root=root,
        bundled=bundled,
        require_redistribution_allowed=require_redistribution_allowed,
        require_all_files_listed=require_all_files_listed,
    )


def parse_reference_manifest_payload(
    payload: dict[str, object],
    *,
    context: str,
) -> ReferenceManifest:
    """Parse a JSON object into the typed manifest model."""

    diagnostic_context = _ReferenceDiagnosticContext.from_raw_payload(payload)
    _require_fields(
        payload,
        required_fields=_REQUIRED_MANIFEST_FIELDS,
        context=context,
        diagnostic_context=diagnostic_context,
    )
    _reject_unrecognized_fields(
        payload,
        allowed_fields=_ALLOWED_MANIFEST_FIELDS,
        context=context,
        diagnostic_context=diagnostic_context,
    )
    files = _parse_file_manifests(
        payload.get("files"),
        context=f"{context}.files",
        diagnostic_context=diagnostic_context,
    )
    redistribution_status = _require_redistribution_status(
        payload,
        key="redistribution_status",
        context=context,
        diagnostic_context=diagnostic_context,
    )
    return ReferenceManifest(
        reference_id=_require_string(
            payload,
            key="reference_id",
            context=context,
            diagnostic_context=diagnostic_context,
        ),
        display_name=_require_string(
            payload,
            key="display_name",
            context=context,
            diagnostic_context=diagnostic_context,
        ),
        organism=_require_string(
            payload,
            key="organism",
            context=context,
            diagnostic_context=diagnostic_context,
        ),
        taxonomy_id=_optional_int(
            payload,
            key="taxonomy_id",
            context=context,
            diagnostic_context=diagnostic_context,
        ),
        protein_namespace=_require_string(
            payload,
            key="protein_namespace",
            context=context,
            diagnostic_context=diagnostic_context,
        ),
        reference_version=_require_string(
            payload,
            key="reference_version",
            context=context,
            diagnostic_context=diagnostic_context,
        ),
        source_name=_require_string(
            payload,
            key="source_name",
            context=context,
            diagnostic_context=diagnostic_context,
        ),
        source_version=_optional_string(
            payload,
            key="source_version",
            context=context,
            diagnostic_context=diagnostic_context,
        ),
        source_url=_optional_string(
            payload,
            key="source_url",
            context=context,
            diagnostic_context=diagnostic_context,
        ),
        retrieved_at=_require_date(
            payload,
            key="retrieved_at",
            context=context,
            diagnostic_context=diagnostic_context,
        ),
        table_sha256=_require_string(
            payload,
            key="table_sha256",
            context=context,
            diagnostic_context=diagnostic_context,
        ),
        source_publication=_optional_string(
            payload,
            key="source_publication",
            context=context,
            diagnostic_context=diagnostic_context,
        ),
        license_name=_optional_string(
            payload,
            key="license_name",
            context=context,
            diagnostic_context=diagnostic_context,
        ),
        license_url=_optional_string(
            payload,
            key="license_url",
            context=context,
            diagnostic_context=diagnostic_context,
        ),
        redistribution_status=redistribution_status,
        redistribution_notes=_require_string(
            payload,
            key="redistribution_notes",
            context=context,
            diagnostic_context=diagnostic_context,
        ),
        redistribution_evidence=_optional_redistribution_evidence(
            payload,
            key="redistribution_evidence",
            context=context,
            diagnostic_context=diagnostic_context,
        ),
        derived_from=_require_string_tuple(
            payload,
            key="derived_from",
            context=context,
            diagnostic_context=diagnostic_context,
        ),
        generated_by=_require_string(
            payload,
            key="generated_by",
            context=context,
            diagnostic_context=diagnostic_context,
        ),
        generated_at_utc=_require_string(
            payload,
            key="generated_at_utc",
            context=context,
            diagnostic_context=diagnostic_context,
        ),
        manifest_schema_version=_require_string(
            payload,
            key="manifest_schema_version",
            context=context,
            diagnostic_context=diagnostic_context,
        ),
        files=files,
        sequence_context_policy=_optional_string(
            payload,
            key="sequence_context_policy",
            context=context,
            diagnostic_context=diagnostic_context,
        ),
        sequence_window_length=_optional_int(
            payload,
            key="sequence_window_length",
            context=context,
            diagnostic_context=diagnostic_context,
        ),
        sequence_center_index=_optional_int(
            payload,
            key="sequence_center_index",
            context=context,
            diagnostic_context=diagnostic_context,
        ),
        allowed_sequence_alphabet=_optional_string(
            payload,
            key="allowed_sequence_alphabet",
            context=context,
            diagnostic_context=diagnostic_context,
        ),
        organism_common_name=_optional_string(
            payload,
            key="organism_common_name",
            context=context,
            diagnostic_context=diagnostic_context,
        ),
        supports=_optional_string_tuple(
            payload,
            key="supports",
            context=context,
            diagnostic_context=diagnostic_context,
        ),
        limitations=_optional_string_tuple(
            payload,
            key="limitations",
            context=context,
            diagnostic_context=diagnostic_context,
        ),
        raw_redistribution_allowed=_optional_bool(
            payload,
            key="redistribution_allowed",
            context=context,
            diagnostic_context=diagnostic_context,
        ),
    )


def _parse_file_manifests(
    value: object,
    *,
    context: str,
    diagnostic_context: _ReferenceDiagnosticContext,
) -> tuple[ReferenceFileManifest, ...]:
    if not isinstance(value, list):
        raise _manifest_parse_error(
            diagnostic_context,
            field="files",
            actual_value=value,
            reason=f"{context} must be an array",
        )
    if not value:
        raise _manifest_parse_error(
            diagnostic_context,
            field="files",
            actual_value=value,
            reason=f"{context} must not be empty",
        )
    files: list[ReferenceFileManifest] = []
    for index, item in enumerate(value):
        file_context = f"{context}[{index}]"
        if not isinstance(item, dict):
            raise _manifest_parse_error(
                diagnostic_context,
                field=f"files[{index}]",
                actual_value=item,
                reason=f"{file_context} must be an object",
            )
        file_payload = cast(dict[str, object], item)
        _require_fields(
            file_payload,
            required_fields=_REQUIRED_FILE_FIELDS,
            context=file_context,
            field_path_prefix=f"files[{index}]",
            diagnostic_context=diagnostic_context,
        )
        _reject_unrecognized_fields(
            file_payload,
            allowed_fields=_ALLOWED_FILE_FIELDS,
            context=file_context,
            field_path_prefix=f"files[{index}]",
            diagnostic_context=diagnostic_context,
        )
        files.append(
            ReferenceFileManifest(
                relative_path=_require_string(
                    file_payload,
                    key="relative_path",
                    context=file_context,
                    field_path=f"files[{index}].relative_path",
                    diagnostic_context=diagnostic_context,
                ),
                role=_require_string(
                    file_payload,
                    key="role",
                    context=file_context,
                    field_path=f"files[{index}].role",
                    diagnostic_context=diagnostic_context,
                ),
                format=_require_string(
                    file_payload,
                    key="format",
                    context=file_context,
                    field_path=f"files[{index}].format",
                    diagnostic_context=diagnostic_context,
                ),
                sha256=_require_string(
                    file_payload,
                    key="sha256",
                    context=file_context,
                    field_path=f"files[{index}].sha256",
                    diagnostic_context=diagnostic_context,
                ),
                row_count=_optional_int(
                    file_payload,
                    key="row_count",
                    context=file_context,
                    field_path=f"files[{index}].row_count",
                    diagnostic_context=diagnostic_context,
                ),
                column_names=_optional_column_names(
                    file_payload,
                    key="column_names",
                    context=file_context,
                    field_path=f"files[{index}].column_names",
                    diagnostic_context=diagnostic_context,
                ),
            )
        )
    return tuple(files)


def _require_fields(
    payload: dict[str, object],
    *,
    required_fields: frozenset[str],
    context: str,
    diagnostic_context: _ReferenceDiagnosticContext,
    field_path_prefix: str = "",
) -> None:
    missing = sorted(field for field in required_fields if field not in payload)
    if missing:
        paths = (
            missing
            if not field_path_prefix
            else [f"{field_path_prefix}.{field}" for field in missing]
        )
        raise _manifest_parse_error(
            diagnostic_context,
            field=", ".join(paths),
            reason=(
                f"reference manifest is missing required field(s) for {context}: "
                f"{', '.join(missing)}"
            ),
        )


def _require_string(
    payload: dict[str, object],
    *,
    key: str,
    context: str,
    diagnostic_context: _ReferenceDiagnosticContext,
    field_path: str | None = None,
) -> str:
    value = payload.get(key)
    if not isinstance(value, str) or not value.strip():
        raise _manifest_parse_error(
            diagnostic_context,
            field=field_path or key,
            actual_value=value if key in payload else _MISSING,
            reason=(
                f"reference manifest {key} must be a non-empty string for {context}"
            ),
        )
    return value.strip()


def _optional_string(
    payload: dict[str, object],
    *,
    key: str,
    context: str,
    diagnostic_context: _ReferenceDiagnosticContext,
    field_path: str | None = None,
) -> str | None:
    if key not in payload or payload.get(key) is None:
        return None
    value = payload.get(key)
    if not isinstance(value, str):
        raise _manifest_parse_error(
            diagnostic_context,
            field=field_path or key,
            actual_value=value,
            reason=f"reference manifest {key} must be a string or null for {context}",
        )
    text = value.strip()
    return text if text else None


def _require_date(
    payload: dict[str, object],
    *,
    key: str,
    context: str,
    diagnostic_context: _ReferenceDiagnosticContext,
    field_path: str | None = None,
) -> date:
    value = payload.get(key)
    if not isinstance(value, str) or not value.strip():
        raise _manifest_parse_error(
            diagnostic_context,
            field=field_path or key,
            actual_value=value if key in payload else _MISSING,
            reason=f"reference manifest {key} must be YYYY-MM-DD for {context}",
        )
    try:
        return date.fromisoformat(value.strip())
    except ValueError as exc:
        raise _manifest_parse_error(
            diagnostic_context,
            field=field_path or key,
            actual_value=value,
            reason=f"reference manifest {key} must be YYYY-MM-DD for {context}",
        ) from exc


def _optional_date(
    payload: dict[str, object],
    *,
    key: str,
    context: str,
    diagnostic_context: _ReferenceDiagnosticContext,
    field_path: str | None = None,
) -> date | None:
    if key not in payload or payload.get(key) is None:
        return None
    return _require_date(
        payload,
        key=key,
        context=context,
        field_path=field_path,
        diagnostic_context=diagnostic_context,
    )


def _require_redistribution_status(
    payload: dict[str, object],
    *,
    key: str,
    context: str,
    diagnostic_context: _ReferenceDiagnosticContext,
) -> RedistributionStatus:
    value = payload.get(key)
    if not isinstance(value, str) or not value.strip():
        allowed = ", ".join(item.value for item in RedistributionStatus)
        raise _manifest_parse_error(
            diagnostic_context,
            field=key,
            actual_value=value if key in payload else _MISSING,
            reason=f"reference manifest {key} must be one of {allowed} for {context}",
        )
    try:
        return RedistributionStatus(value.strip())
    except ValueError as exc:
        allowed = ", ".join(item.value for item in RedistributionStatus)
        raise _manifest_parse_error(
            diagnostic_context,
            field=key,
            actual_value=value,
            reason=f"reference manifest {key} must be one of {allowed} for {context}",
        ) from exc


def _optional_redistribution_evidence(
    payload: dict[str, object],
    *,
    key: str,
    context: str,
    diagnostic_context: _ReferenceDiagnosticContext,
) -> RedistributionEvidence | None:
    if key not in payload or payload.get(key) is None:
        return None
    value = payload.get(key)
    if not isinstance(value, dict):
        raise _manifest_parse_error(
            diagnostic_context,
            field=key,
            actual_value=value,
            reason=(
                f"reference manifest {key} must be an object or null for {context}"
            ),
        )
    evidence_payload = cast(dict[str, object], value)
    evidence_context = f"{context}.{key}"
    _require_fields(
        evidence_payload,
        required_fields=_REQUIRED_REDISTRIBUTION_EVIDENCE_FIELDS,
        context=evidence_context,
        field_path_prefix=key,
        diagnostic_context=diagnostic_context,
    )
    _reject_unrecognized_fields(
        evidence_payload,
        allowed_fields=_ALLOWED_REDISTRIBUTION_EVIDENCE_FIELDS,
        context=evidence_context,
        field_path_prefix=key,
        diagnostic_context=diagnostic_context,
    )
    return RedistributionEvidence(
        evidence_type=_require_redistribution_evidence_type(
            evidence_payload,
            key="evidence_type",
            context=evidence_context,
            field_path="redistribution_evidence.evidence_type",
            diagnostic_context=diagnostic_context,
        ),
        upstream_package=_require_upstream_package_license_evidence(
            evidence_payload,
            context=evidence_context,
            diagnostic_context=diagnostic_context,
        ),
        scope=_require_redistribution_scope(
            evidence_payload,
            context=evidence_context,
            diagnostic_context=diagnostic_context,
        ),
        attribution=_require_redistribution_attribution(
            evidence_payload,
            context=evidence_context,
            diagnostic_context=diagnostic_context,
        ),
        independent_database_permission_claimed=_require_bool(
            evidence_payload,
            key="independent_database_permission_claimed",
            context=evidence_context,
            field_path=(
                "redistribution_evidence.independent_database_permission_claimed"
            ),
            diagnostic_context=diagnostic_context,
        ),
        evidence_url=_optional_string(
            evidence_payload,
            key="evidence_url",
            context=evidence_context,
            field_path="redistribution_evidence.evidence_url",
            diagnostic_context=diagnostic_context,
        ),
        verified_at=_optional_redistribution_evidence_verified_at(
            evidence_payload,
            context=evidence_context,
            diagnostic_context=diagnostic_context,
        ),
        notes=_optional_string(
            evidence_payload,
            key="notes",
            context=evidence_context,
            field_path="redistribution_evidence.notes",
            diagnostic_context=diagnostic_context,
        ),
    )


def _optional_redistribution_evidence_verified_at(
    payload: dict[str, object],
    *,
    context: str,
    diagnostic_context: _ReferenceDiagnosticContext,
) -> date | None:
    if "verified_at" not in payload or payload.get("verified_at") is None:
        return None
    return _require_date(
        payload,
        key="verified_at",
        context=context,
        field_path="redistribution_evidence.verified_at",
        diagnostic_context=diagnostic_context,
    )


def _require_upstream_package_license_evidence(
    payload: dict[str, object],
    *,
    context: str,
    diagnostic_context: _ReferenceDiagnosticContext,
) -> UpstreamPackageLicenseEvidence:
    upstream_context = f"{context}.upstream_package"
    upstream_payload = _require_object(
        payload,
        key="upstream_package",
        context=context,
        field_path="redistribution_evidence.upstream_package",
        diagnostic_context=diagnostic_context,
    )
    _require_fields(
        upstream_payload,
        required_fields=_REQUIRED_UPSTREAM_PACKAGE_FIELDS,
        context=upstream_context,
        field_path_prefix="redistribution_evidence.upstream_package",
        diagnostic_context=diagnostic_context,
    )
    _reject_unrecognized_fields(
        upstream_payload,
        allowed_fields=_ALLOWED_UPSTREAM_PACKAGE_FIELDS,
        context=upstream_context,
        field_path_prefix="redistribution_evidence.upstream_package",
        diagnostic_context=diagnostic_context,
    )
    return UpstreamPackageLicenseEvidence(
        package_name=_require_string(
            upstream_payload,
            key="package_name",
            context=upstream_context,
            field_path="redistribution_evidence.upstream_package.package_name",
            diagnostic_context=diagnostic_context,
        ),
        package_version=_require_string(
            upstream_payload,
            key="package_version",
            context=upstream_context,
            field_path="redistribution_evidence.upstream_package.package_version",
            diagnostic_context=diagnostic_context,
        ),
        license_name=_require_string(
            upstream_payload,
            key="license_name",
            context=upstream_context,
            field_path="redistribution_evidence.upstream_package.license_name",
            diagnostic_context=diagnostic_context,
        ),
        license_url=_optional_string(
            upstream_payload,
            key="license_url",
            context=upstream_context,
            field_path="redistribution_evidence.upstream_package.license_url",
            diagnostic_context=diagnostic_context,
        ),
    )


def _require_redistribution_scope(
    payload: dict[str, object],
    *,
    context: str,
    diagnostic_context: _ReferenceDiagnosticContext,
) -> RedistributionScope:
    scope_context = f"{context}.scope"
    scope_payload = _require_object(
        payload,
        key="scope",
        context=context,
        field_path="redistribution_evidence.scope",
        diagnostic_context=diagnostic_context,
    )
    _require_fields(
        scope_payload,
        required_fields=_REQUIRED_REDISTRIBUTION_SCOPE_FIELDS,
        context=scope_context,
        field_path_prefix="redistribution_evidence.scope",
        diagnostic_context=diagnostic_context,
    )
    _reject_unrecognized_fields(
        scope_payload,
        allowed_fields=_ALLOWED_REDISTRIBUTION_SCOPE_FIELDS,
        context=scope_context,
        field_path_prefix="redistribution_evidence.scope",
        diagnostic_context=diagnostic_context,
    )
    return RedistributionScope(
        reference_id=_require_string(
            scope_payload,
            key="reference_id",
            context=scope_context,
            field_path="redistribution_evidence.scope.reference_id",
            diagnostic_context=diagnostic_context,
        ),
        reference_version=_require_string(
            scope_payload,
            key="reference_version",
            context=scope_context,
            field_path="redistribution_evidence.scope.reference_version",
            diagnostic_context=diagnostic_context,
        ),
        applies_to_exact_packaged_files=_require_bool(
            scope_payload,
            key="applies_to_exact_packaged_files",
            context=scope_context,
            field_path=(
                "redistribution_evidence.scope.applies_to_exact_packaged_files"
            ),
            diagnostic_context=diagnostic_context,
        ),
        packaged_files=_require_string_tuple(
            scope_payload,
            key="packaged_files",
            context=scope_context,
            field_path="redistribution_evidence.scope.packaged_files",
            diagnostic_context=diagnostic_context,
        ),
        applies_to_future_bundles=_require_bool(
            scope_payload,
            key="applies_to_future_bundles",
            context=scope_context,
            field_path="redistribution_evidence.scope.applies_to_future_bundles",
            diagnostic_context=diagnostic_context,
        ),
    )


def _require_redistribution_attribution(
    payload: dict[str, object],
    *,
    context: str,
    diagnostic_context: _ReferenceDiagnosticContext,
) -> RedistributionAttribution:
    attribution_context = f"{context}.attribution"
    attribution_payload = _require_object(
        payload,
        key="attribution",
        context=context,
        field_path="redistribution_evidence.attribution",
        diagnostic_context=diagnostic_context,
    )
    _require_fields(
        attribution_payload,
        required_fields=_REQUIRED_REDISTRIBUTION_ATTRIBUTION_FIELDS,
        context=attribution_context,
        field_path_prefix="redistribution_evidence.attribution",
        diagnostic_context=diagnostic_context,
    )
    _reject_unrecognized_fields(
        attribution_payload,
        allowed_fields=_ALLOWED_REDISTRIBUTION_ATTRIBUTION_FIELDS,
        context=attribution_context,
        field_path_prefix="redistribution_evidence.attribution",
        diagnostic_context=diagnostic_context,
    )
    return RedistributionAttribution(
        repository_notice_path=_require_string(
            attribution_payload,
            key="repository_notice_path",
            context=attribution_context,
            field_path="redistribution_evidence.attribution.repository_notice_path",
            diagnostic_context=diagnostic_context,
        ),
        bundle_attribution_path=_require_string(
            attribution_payload,
            key="bundle_attribution_path",
            context=attribution_context,
            field_path="redistribution_evidence.attribution.bundle_attribution_path",
            diagnostic_context=diagnostic_context,
        ),
    )


def _require_object(
    payload: dict[str, object],
    *,
    key: str,
    context: str,
    diagnostic_context: _ReferenceDiagnosticContext,
    field_path: str | None = None,
) -> dict[str, object]:
    value = payload.get(key)
    if not isinstance(value, dict):
        raise _manifest_parse_error(
            diagnostic_context,
            field=field_path or key,
            actual_value=value if key in payload else _MISSING,
            reason=f"reference manifest {key} must be an object for {context}",
        )
    return cast(dict[str, object], value)


def _reject_unrecognized_fields(
    payload: dict[str, object],
    *,
    allowed_fields: frozenset[str],
    context: str,
    diagnostic_context: _ReferenceDiagnosticContext,
    field_path_prefix: str = "",
) -> None:
    extra = sorted(field for field in payload if field not in allowed_fields)
    if extra:
        paths = (
            extra
            if not field_path_prefix
            else [f"{field_path_prefix}.{field}" for field in extra]
        )
        actual_value = (
            payload[extra[0]]
            if len(extra) == 1
            else {
                path: payload[field] for path, field in zip(paths, extra, strict=True)
            }
        )
        raise _manifest_parse_error(
            diagnostic_context,
            field=", ".join(paths),
            actual_value=actual_value,
            reason=(
                f"reference manifest has unrecognized field(s) for {context}: "
                f"{', '.join(paths)}"
            ),
        )


def _require_redistribution_evidence_type(
    payload: dict[str, object],
    *,
    key: str,
    context: str,
    diagnostic_context: _ReferenceDiagnosticContext,
    field_path: str | None = None,
) -> RedistributionEvidenceType:
    value = payload.get(key)
    if not isinstance(value, str) or not value.strip():
        allowed = ", ".join(item.value for item in RedistributionEvidenceType)
        raise _manifest_parse_error(
            diagnostic_context,
            field=field_path or key,
            actual_value=value if key in payload else _MISSING,
            reason=f"reference manifest {key} must be one of {allowed} for {context}",
        )
    try:
        return RedistributionEvidenceType(value.strip())
    except ValueError as exc:
        allowed = ", ".join(item.value for item in RedistributionEvidenceType)
        raise _manifest_parse_error(
            diagnostic_context,
            field=field_path or key,
            actual_value=value,
            reason=f"reference manifest {key} must be one of {allowed} for {context}",
        ) from exc


def _require_bool(
    payload: dict[str, object],
    *,
    key: str,
    context: str,
    diagnostic_context: _ReferenceDiagnosticContext,
    field_path: str | None = None,
) -> bool:
    value = payload.get(key)
    if not isinstance(value, bool):
        raise _manifest_parse_error(
            diagnostic_context,
            field=field_path or key,
            actual_value=value if key in payload else _MISSING,
            reason=(
                f"reference manifest {key} must be a JSON Boolean "
                f"(true or false) for {context}; got {value!r} "
                f"(type {type(value).__name__})"
            ),
        )
    return value


def _optional_bool(
    payload: dict[str, object],
    *,
    key: str,
    context: str,
    diagnostic_context: _ReferenceDiagnosticContext,
    field_path: str | None = None,
) -> bool | None:
    if key not in payload:
        return None
    return _require_bool(
        payload,
        key=key,
        context=context,
        field_path=field_path,
        diagnostic_context=diagnostic_context,
    )


def _optional_int(
    payload: dict[str, object],
    *,
    key: str,
    context: str,
    diagnostic_context: _ReferenceDiagnosticContext,
    field_path: str | None = None,
) -> int | None:
    if key not in payload or payload.get(key) is None:
        return None
    value = payload.get(key)
    if isinstance(value, bool) or not isinstance(value, int):
        raise _manifest_parse_error(
            diagnostic_context,
            field=field_path or key,
            actual_value=value,
            reason=(
                f"reference manifest {key} must be an integer or null for {context}"
            ),
        )
    return int(value)


def _require_string_tuple(
    payload: dict[str, object],
    *,
    key: str,
    context: str,
    diagnostic_context: _ReferenceDiagnosticContext,
    field_path: str | None = None,
) -> tuple[str, ...]:
    if key not in payload:
        raise _manifest_parse_error(
            diagnostic_context,
            field=field_path or key,
            reason=f"reference manifest {key} is required for {context}",
        )
    return _string_tuple(
        payload.get(key),
        key=key,
        context=context,
        allow_empty=False,
        field_path=field_path,
        diagnostic_context=diagnostic_context,
    )


def _optional_string_tuple(
    payload: dict[str, object],
    *,
    key: str,
    context: str,
    diagnostic_context: _ReferenceDiagnosticContext,
    field_path: str | None = None,
) -> tuple[str, ...]:
    if key not in payload or payload.get(key) is None:
        return ()
    return _string_tuple(
        payload.get(key),
        key=key,
        context=context,
        allow_empty=True,
        field_path=field_path,
        diagnostic_context=diagnostic_context,
    )


def _string_tuple(
    value: object,
    *,
    key: str,
    context: str,
    allow_empty: bool,
    diagnostic_context: _ReferenceDiagnosticContext,
    field_path: str | None = None,
) -> tuple[str, ...]:
    resolved: list[str] = []
    base_field_path = field_path or key
    if not isinstance(value, list):
        raise _manifest_parse_error(
            diagnostic_context,
            field=base_field_path,
            actual_value=value,
            reason=(
                f"reference manifest {key} must be an array of strings for {context}"
            ),
        )
    for index, item in enumerate(value):
        if not isinstance(item, str) or not item.strip():
            raise _manifest_parse_error(
                diagnostic_context,
                field=f"{base_field_path}[{index}]",
                actual_value=item,
                reason=(
                    f"reference manifest {key}[{index}] must be a non-empty "
                    f"string for {context}"
                ),
            )
        resolved.append(item.strip())
    if not allow_empty and not resolved:
        raise _manifest_parse_error(
            diagnostic_context,
            field=base_field_path,
            actual_value=value,
            reason=f"reference manifest {key} must not be empty for {context}",
        )
    return tuple(resolved)


def _optional_column_names(
    payload: dict[str, object],
    *,
    key: str,
    context: str,
    diagnostic_context: _ReferenceDiagnosticContext,
    field_path: str | None = None,
) -> tuple[str, ...] | None:
    if key not in payload or payload.get(key) is None:
        return None
    value = payload.get(key)
    base_field_path = field_path or key
    if not isinstance(value, list):
        raise _manifest_parse_error(
            diagnostic_context,
            field=base_field_path,
            actual_value=value,
            reason=(
                f"reference manifest {key} must be an array of strings or null "
                f"for {context}"
            ),
        )
    resolved: list[str] = []
    for index, item in enumerate(value):
        if not isinstance(item, str) or not item.strip():
            raise _manifest_parse_error(
                diagnostic_context,
                field=f"{base_field_path}[{index}]",
                actual_value=item,
                reason=(
                    f"reference manifest {key}[{index}] must be a non-empty "
                    f"string for {context}"
                ),
            )
        resolved.append(item.strip())
    return tuple(resolved)
