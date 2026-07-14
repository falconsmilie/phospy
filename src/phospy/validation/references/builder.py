"""Validation for local reference-bundle build requests."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from datetime import date
from os import PathLike
from pathlib import Path
from typing import cast

from phospy.errors.references import ReferenceResolutionError
from phospy.science.references.models import (
    Organism,
    ReferenceBundleBuildRequest,
    SequenceWindowDefinition,
)


@dataclass(frozen=True, slots=True)
class ValidatedReferenceBundleBuildRequest:
    """Validated local-source reference build request."""

    organism: Organism
    kinase_substrate_path: Path
    site_sequence_path: Path
    source_name: str
    source_version: str
    retrieved_at: date
    license: str
    redistribution_status: str
    identifier_namespace: str
    sequence_window: SequenceWindowDefinition | None
    bundle_id: str | None
    organism_common_name: str | None
    supports: tuple[str, ...]
    limitations: tuple[str, ...]
    reference_version: str | None


class ReferenceBundleBuildRequestValidator:
    """Validate local-source reference builder request fields."""

    def run(
        self,
        request: ReferenceBundleBuildRequest,
    ) -> ValidatedReferenceBundleBuildRequest:
        if not isinstance(request, ReferenceBundleBuildRequest):
            raise ReferenceResolutionError(
                "reference bundle build request must be a ReferenceBundleBuildRequest"
            )
        organism = request.organism
        if not isinstance(cast(object, organism), Organism):
            raise ReferenceResolutionError(
                "reference bundle build request organism must be an Organism enum value"
            )
        sequence_window = request.sequence_window
        if sequence_window is not None:
            self._validate_sequence_window(sequence_window)
        return ValidatedReferenceBundleBuildRequest(
            organism=organism,
            kinase_substrate_path=self._require_local_path(
                request.kinase_substrate_path,
                field_name="reference bundle build request kinase_substrate_path",
            ),
            site_sequence_path=self._require_local_path(
                request.site_sequence_path,
                field_name="reference bundle build request site_sequence_path",
            ),
            source_name=self._require_non_empty_string(
                request.source_name,
                field_name="reference bundle build request source_name",
            ),
            source_version=self._require_non_empty_string(
                request.source_version,
                field_name="reference bundle build request source_version",
            ),
            retrieved_at=self._require_date(
                request.retrieved_at,
                field_name="reference bundle build request retrieved_at",
            ),
            license=self._require_non_empty_string(
                request.license,
                field_name="reference bundle build request license",
            ),
            redistribution_status=self._require_non_empty_string(
                request.redistribution_status,
                field_name="reference bundle build request redistribution_status",
            ),
            identifier_namespace=self._require_non_empty_string(
                request.identifier_namespace,
                field_name="reference bundle build request identifier_namespace",
            ),
            sequence_window=sequence_window,
            bundle_id=self._optional_non_empty_string(
                request.bundle_id,
                field_name="reference bundle build request bundle_id",
            ),
            organism_common_name=self._optional_non_empty_string(
                request.organism_common_name,
                field_name="reference bundle build request organism_common_name",
            ),
            supports=self._require_non_empty_string_sequence(
                request.supports,
                field_name="reference bundle build request supports",
            ),
            limitations=self._require_non_empty_string_sequence(
                request.limitations,
                field_name="reference bundle build request limitations",
            ),
            reference_version=self._optional_non_empty_string(
                request.reference_version,
                field_name="reference bundle build request reference_version",
            ),
        )

    @staticmethod
    def _require_local_path(value: object, *, field_name: str) -> Path:
        if isinstance(value, str):
            normalized = value.strip()
            if not normalized:
                raise ReferenceResolutionError(f"{field_name} must be non-empty")
            if "://" in normalized.lower():
                raise ReferenceResolutionError(
                    f"{field_name} must be a local filesystem path; remote URLs "
                    "are not supported"
                )
            return Path(normalized)
        if isinstance(value, Path):
            return value
        if isinstance(value, PathLike):
            return Path(cast(PathLike[str], value))
        raise ReferenceResolutionError(f"{field_name} must be a local filesystem path")

    @staticmethod
    def _require_non_empty_string(value: object, *, field_name: str) -> str:
        if not isinstance(value, str) or not value.strip():
            raise ReferenceResolutionError(f"{field_name} must be a non-empty string")
        return value.strip()

    @classmethod
    def _optional_non_empty_string(
        cls,
        value: object,
        *,
        field_name: str,
    ) -> str | None:
        if value is None:
            return None
        return cls._require_non_empty_string(value, field_name=field_name)

    @classmethod
    def _require_non_empty_string_sequence(
        cls,
        value: object,
        *,
        field_name: str,
    ) -> tuple[str, ...]:
        if not isinstance(value, Sequence) or isinstance(
            value,
            (str, bytes, bytearray),
        ):
            raise ReferenceResolutionError(
                f"{field_name} must be a non-empty sequence of strings"
            )
        resolved = tuple(
            cls._require_non_empty_string(
                item,
                field_name=f"{field_name}[{index}]",
            )
            for index, item in enumerate(value)
        )
        if not resolved:
            raise ReferenceResolutionError(
                f"{field_name} must be a non-empty sequence of strings"
            )
        return resolved

    @staticmethod
    def _require_date(value: object, *, field_name: str) -> date:
        if isinstance(value, date):
            return value
        if not isinstance(value, str) or not value.strip():
            raise ReferenceResolutionError(f"{field_name} must be a YYYY-MM-DD date")
        try:
            return date.fromisoformat(value.strip())
        except ValueError as exc:
            raise ReferenceResolutionError(
                f"{field_name} must be a YYYY-MM-DD date"
            ) from exc

    @staticmethod
    def _validate_sequence_window(value: object) -> None:
        if not isinstance(value, SequenceWindowDefinition):
            raise ReferenceResolutionError(
                "reference bundle build request sequence_window must be a "
                "SequenceWindowDefinition or None"
            )
        if value.upstream_residues < 0 or value.downstream_residues < 0:
            raise ReferenceResolutionError(
                "reference bundle build request sequence_window residue counts "
                "must be >= 0"
            )
