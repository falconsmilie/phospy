"""Runtime environment provenance models."""

from __future__ import annotations

from dataclasses import dataclass, field

from phospy.provenance.immutability import freeze_json_mapping
from phospy.provenance.models._shared import JsonValue

ENVIRONMENT_PROVENANCE_SCHEMA_VERSION_V1 = 1

ENVIRONMENT_PROVENANCE_SCHEMA_VERSION_V2 = 2


def _empty_platform_provenance() -> dict[str, str]:
    return {}


@dataclass(frozen=True, slots=True)
class EnvironmentProvenance:
    """Runtime environment fingerprint for reproducibility."""

    package_name: str
    package_version: str
    python_version: str
    dependency_versions: dict[str, str | None]
    schema_version: int = ENVIRONMENT_PROVENANCE_SCHEMA_VERSION_V2
    platform: dict[str, str] = field(default_factory=_empty_platform_provenance)
    blas_lapack: dict[str, JsonValue] = field(default_factory=dict)
    thread_environment: dict[str, str | None] = field(default_factory=dict)
    timezone: str | None = None
    locale: dict[str, str | None] = field(default_factory=dict)
    constraints_fingerprint: dict[str, str | None] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "dependency_versions",
            freeze_json_mapping(
                self.dependency_versions,
                field_name="environment_provenance.dependency_versions",
            ),
        )
        object.__setattr__(
            self,
            "platform",
            freeze_json_mapping(
                self.platform,
                field_name="environment_provenance.platform",
            ),
        )
        object.__setattr__(
            self,
            "blas_lapack",
            freeze_json_mapping(
                self.blas_lapack,
                field_name="environment_provenance.blas_lapack",
            ),
        )
        object.__setattr__(
            self,
            "thread_environment",
            freeze_json_mapping(
                self.thread_environment,
                field_name="environment_provenance.thread_environment",
            ),
        )
        object.__setattr__(
            self,
            "locale",
            freeze_json_mapping(
                self.locale,
                field_name="environment_provenance.locale",
            ),
        )
        object.__setattr__(
            self,
            "constraints_fingerprint",
            freeze_json_mapping(
                self.constraints_fingerprint,
                field_name="environment_provenance.constraints_fingerprint",
            ),
        )
