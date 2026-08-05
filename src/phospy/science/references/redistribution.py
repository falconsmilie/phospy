"""Reference manifest redistribution-evidence value models."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import date

from phospy.provenance.models import JsonValue
from phospy.science.references.manifest_common import (
    _coerce_date,
    _optional_string,
    _required_string,
)
from phospy.science.references.manifest_policy import RedistributionEvidenceType


@dataclass(frozen=True, slots=True)
class UpstreamPackageLicenseEvidence:
    """Upstream package metadata that supplies the redistribution basis."""

    package_name: str
    package_version: str
    license_name: str
    license_url: str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "package_name",
            _required_string(
                self.package_name,
                "redistribution_evidence.upstream_package.package_name",
            ),
        )
        object.__setattr__(
            self,
            "package_version",
            _required_string(
                self.package_version,
                "redistribution_evidence.upstream_package.package_version",
            ),
        )
        object.__setattr__(
            self,
            "license_name",
            _required_string(
                self.license_name,
                "redistribution_evidence.upstream_package.license_name",
            ),
        )
        object.__setattr__(
            self,
            "license_url",
            _optional_string(self.license_url),
        )

    def to_payload(self) -> dict[str, JsonValue]:
        return {
            "package_name": self.package_name,
            "package_version": self.package_version,
            "license_name": self.license_name,
            "license_url": self.license_url,
        }


@dataclass(frozen=True, slots=True)
class RedistributionScope:
    """Exact PhosPy bundle snapshot and file scope covered by the evidence."""

    reference_id: str
    reference_version: str
    applies_to_exact_packaged_files: bool
    packaged_files: tuple[str, ...]
    applies_to_future_bundles: bool

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "reference_id",
            _required_string(
                self.reference_id,
                "redistribution_evidence.scope.reference_id",
            ),
        )
        object.__setattr__(
            self,
            "reference_version",
            _required_string(
                self.reference_version,
                "redistribution_evidence.scope.reference_version",
            ),
        )
        if not isinstance(self.applies_to_exact_packaged_files, bool):
            raise ValueError(
                "reference manifest "
                "redistribution_evidence.scope.applies_to_exact_packaged_files "
                "must be bool"
            )
        object.__setattr__(
            self,
            "applies_to_exact_packaged_files",
            self.applies_to_exact_packaged_files,
        )
        packaged_files = tuple(
            _required_string(
                item,
                "redistribution_evidence.scope.packaged_files",
            )
            for item in self.packaged_files
        )
        object.__setattr__(
            self,
            "packaged_files",
            packaged_files,
        )
        if not isinstance(self.applies_to_future_bundles, bool):
            raise ValueError(
                "reference manifest "
                "redistribution_evidence.scope.applies_to_future_bundles "
                "must be bool"
            )
        object.__setattr__(
            self,
            "applies_to_future_bundles",
            self.applies_to_future_bundles,
        )

    def to_payload(self) -> dict[str, JsonValue]:
        return {
            "reference_id": self.reference_id,
            "reference_version": self.reference_version,
            "applies_to_exact_packaged_files": bool(
                self.applies_to_exact_packaged_files
            ),
            "packaged_files": list(self.packaged_files),
            "applies_to_future_bundles": bool(self.applies_to_future_bundles),
        }


@dataclass(frozen=True, slots=True)
class RedistributionAttribution:
    """Repository and bundle-local attribution locations for a bundled reference."""

    repository_notice_path: str
    bundle_attribution_path: str

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "repository_notice_path",
            _required_string(
                self.repository_notice_path,
                "redistribution_evidence.attribution.repository_notice_path",
            ),
        )
        object.__setattr__(
            self,
            "bundle_attribution_path",
            _required_string(
                self.bundle_attribution_path,
                "redistribution_evidence.attribution.bundle_attribution_path",
            ),
        )

    def to_payload(self) -> dict[str, JsonValue]:
        return {
            "repository_notice_path": self.repository_notice_path,
            "bundle_attribution_path": self.bundle_attribution_path,
        }


@dataclass(frozen=True, slots=True)
class RedistributionEvidence:
    """Typed exact-snapshot redistribution evidence for release validation."""

    evidence_type: RedistributionEvidenceType | str
    upstream_package: UpstreamPackageLicenseEvidence
    scope: RedistributionScope
    attribution: RedistributionAttribution
    independent_database_permission_claimed: bool
    evidence_url: str | None = None
    verified_at: date | None = None
    notes: str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "evidence_type",
            _coerce_redistribution_evidence_type(self.evidence_type),
        )
        if not isinstance(self.upstream_package, UpstreamPackageLicenseEvidence):
            raise ValueError(
                "reference manifest redistribution_evidence.upstream_package "
                "must be UpstreamPackageLicenseEvidence"
            )
        if not isinstance(self.scope, RedistributionScope):
            raise ValueError(
                "reference manifest redistribution_evidence.scope "
                "must be RedistributionScope"
            )
        if not isinstance(self.attribution, RedistributionAttribution):
            raise ValueError(
                "reference manifest redistribution_evidence.attribution "
                "must be RedistributionAttribution"
            )
        if not isinstance(self.independent_database_permission_claimed, bool):
            raise ValueError(
                "reference manifest "
                "redistribution_evidence.independent_database_permission_claimed "
                "must be bool"
            )
        object.__setattr__(
            self,
            "independent_database_permission_claimed",
            self.independent_database_permission_claimed,
        )
        object.__setattr__(self, "evidence_url", _optional_string(self.evidence_url))
        object.__setattr__(
            self,
            "verified_at",
            None if self.verified_at is None else _coerce_date(self.verified_at),
        )
        object.__setattr__(self, "notes", _optional_string(self.notes))

    def to_payload(self) -> dict[str, JsonValue]:
        evidence_type = _coerce_redistribution_evidence_type(self.evidence_type)
        return {
            "evidence_type": evidence_type.value,
            "upstream_package": self.upstream_package.to_payload(),
            "scope": self.scope.to_payload(),
            "attribution": self.attribution.to_payload(),
            "independent_database_permission_claimed": bool(
                self.independent_database_permission_claimed
            ),
            "evidence_url": self.evidence_url,
            "verified_at": (
                None if self.verified_at is None else self.verified_at.isoformat()
            ),
            "notes": self.notes,
        }


def _coerce_redistribution_evidence_type(
    value: RedistributionEvidenceType | str,
) -> RedistributionEvidenceType:
    if isinstance(value, RedistributionEvidenceType):
        return value
    try:
        return RedistributionEvidenceType(str(value).strip())
    except ValueError as exc:
        allowed = ", ".join(item.value for item in RedistributionEvidenceType)
        raise ValueError(
            "reference manifest redistribution_evidence.evidence_type must be "
            f"one of: {allowed}"
        ) from exc


__all__ = [
    "RedistributionAttribution",
    "RedistributionEvidence",
    "RedistributionScope",
    "UpstreamPackageLicenseEvidence",
]
