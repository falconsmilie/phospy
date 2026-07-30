"""Environment provenance payload serialization."""

from __future__ import annotations

from collections.abc import Mapping

from phospy.provenance.models import (
    ENVIRONMENT_PROVENANCE_SCHEMA_VERSION_V2,
    EnvironmentProvenance,
)
from phospy.provenance.serialization._payload import (
    optional_str,
    raise_legacy_provenance_schema,
    require_int,
    require_mapping,
    require_str,
    to_json_safe,
    to_json_value,
)


def environment_to_payload(environment: EnvironmentProvenance) -> dict[str, object]:
    return {
        "schema_version": int(environment.schema_version),
        "package_name": environment.package_name,
        "package_version": environment.package_version,
        "python_version": environment.python_version,
        "dependency_versions": to_json_safe(environment.dependency_versions),
        "platform": to_json_safe(environment.platform),
        "blas_lapack": to_json_safe(environment.blas_lapack),
        "thread_environment": to_json_safe(environment.thread_environment),
        "timezone": environment.timezone,
        "locale": to_json_safe(environment.locale),
        "constraints_fingerprint": to_json_safe(environment.constraints_fingerprint),
    }


def environment_from_payload(payload: Mapping[str, object]) -> EnvironmentProvenance:
    dependency_versions = require_mapping(
        payload.get("dependency_versions"),
        field_name="provenance.environment.dependency_versions",
    )
    platform_payload = require_mapping(
        payload.get("platform", {}),
        field_name="provenance.environment.platform",
    )
    blas_lapack_payload = require_mapping(
        payload.get("blas_lapack", {}),
        field_name="provenance.environment.blas_lapack",
    )
    thread_environment_payload = require_mapping(
        payload.get("thread_environment", {}),
        field_name="provenance.environment.thread_environment",
    )
    locale_payload = require_mapping(
        payload.get("locale", {}),
        field_name="provenance.environment.locale",
    )
    constraints_fingerprint_payload = require_mapping(
        payload.get("constraints_fingerprint", {}),
        field_name="provenance.environment.constraints_fingerprint",
    )
    return EnvironmentProvenance(
        schema_version=_require_current_environment_schema_version(payload),
        package_name=require_str(
            payload.get("package_name"),
            field_name="provenance.environment.package_name",
        ),
        package_version=require_str(
            payload.get("package_version"),
            field_name="provenance.environment.package_version",
        ),
        python_version=require_str(
            payload.get("python_version"),
            field_name="provenance.environment.python_version",
        ),
        dependency_versions={
            key: (
                None
                if value is None
                else require_str(
                    value,
                    field_name=(f"provenance.environment.dependency_versions['{key}']"),
                )
            )
            for key, value in dependency_versions.items()
        },
        platform={
            key: require_str(
                value,
                field_name=f"provenance.environment.platform['{key}']",
            )
            for key, value in platform_payload.items()
        },
        blas_lapack={
            key: to_json_value(value) for key, value in blas_lapack_payload.items()
        },
        thread_environment={
            key: (
                None
                if value is None
                else require_str(
                    value,
                    field_name=f"provenance.environment.thread_environment['{key}']",
                )
            )
            for key, value in thread_environment_payload.items()
        },
        timezone=optional_str(
            payload.get("timezone"),
            field_name="provenance.environment.timezone",
        ),
        locale={
            key: (
                None
                if value is None
                else require_str(
                    value,
                    field_name=f"provenance.environment.locale['{key}']",
                )
            )
            for key, value in locale_payload.items()
        },
        constraints_fingerprint={
            key: (
                None
                if value is None
                else require_str(
                    value,
                    field_name=(
                        f"provenance.environment.constraints_fingerprint['{key}']"
                    ),
                )
            )
            for key, value in constraints_fingerprint_payload.items()
        },
    )


def _require_current_environment_schema_version(
    payload: Mapping[str, object],
) -> int:
    if "schema_version" not in payload:
        raise_legacy_provenance_schema()
    schema_version = require_int(
        payload.get("schema_version"),
        field_name="provenance.environment.schema_version",
    )
    if schema_version != ENVIRONMENT_PROVENANCE_SCHEMA_VERSION_V2:
        raise_legacy_provenance_schema()
    return schema_version
