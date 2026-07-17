"""Build validated reference bundles from local source files."""

from __future__ import annotations

import re
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from datetime import date
from hashlib import sha256
from os import PathLike
from pathlib import Path
from typing import Protocol, cast

import pandas as pd

from phospy.errors.references import ReferenceResolutionError
from phospy.provenance.hashing import fingerprint_table
from phospy.provenance.models import ReferenceProvenance
from phospy.provenance.references import fingerprint_local_reference_source_file
from phospy.science.references.identifiers import (
    merge_reference_identifier_normalisation_reports,
)
from phospy.science.references.manifest import REFERENCE_MANIFEST_SCHEMA_VERSION
from phospy.science.references.models import (
    Organism,
    RedistributionStatus,
    ReferenceBuildPath,
    ReferenceBundle,
    ReferenceBundleBuildRequest,
    ReferenceContext,
    ReferenceFileManifest,
    ReferenceManifest,
    SequenceWindowDefinition,
)
from phospy.science.sites.identifiers import canonicalize_site_series
from phospy.science.tables.references import (
    KinaseSubstrateReference,
    SiteSequenceReference,
)

_COLUMN_TOKEN_PATTERN = re.compile(r"[^a-z0-9]+")
_CENTRAL_PHOSPHO_RESIDUES = frozenset({"S", "T", "Y"})
_KINASE_ALIASES = (
    "kinase",
    "kinase_gene",
    "kinase_gene_symbol",
    "kinase_symbol",
    "enzyme",
)
_SUBSTRATE_SITE_ALIASES = (
    "substrate_site",
    "site_id",
    "display_id",
    "phosphosite",
    "site_identifier",
)
_SITE_ID_ALIASES = (
    "site_id",
    "substrate_site",
    "display_id",
    "phosphosite",
    "site_identifier",
    "index",
    "unnamed_0",
)
_SITE_SEQUENCE_ALIASES = (
    "site_sequence",
    "centralized_sequence",
    "centralised_sequence",
    "sequence_window",
    "motif_sequence",
    "sequence",
)
_DISPLAY_ID_ALIASES = (
    "display_id",
    "display_site",
    "display_site_id",
)
_ORGANISM_ALIASES = ("organism", "species")
_GENE_SYMBOL_ALIASES = (
    "gene_symbol",
    "gene",
    "substrate_gene",
    "substrate_gene_symbol",
)
_PROTEIN_ACCESSION_ALIASES = (
    "protein_accession",
    "accession",
    "substrate_accession",
    "uniprot",
    "uniprot_id",
)
_PROTEIN_ID_ALIASES = (
    "protein_id",
    "substrate_protein",
    "substrate_protein_id",
)
_ORGANISM_TOKENS = {
    Organism.HUMAN: frozenset({"human", "homo sapiens", "homo_sapiens", "9606"}),
    Organism.MOUSE: frozenset({"mouse", "mus musculus", "mus_musculus", "10090"}),
    Organism.RAT: frozenset({"rat", "rattus norvegicus", "rattus_norvegicus", "10116"}),
}


class ValidatedReferenceBundleBuildRequestProtocol(Protocol):
    """Validated local-source reference build request shape."""

    @property
    def organism(self) -> Organism: ...

    @property
    def kinase_substrate_path(self) -> Path: ...

    @property
    def site_sequence_path(self) -> Path: ...

    @property
    def source_name(self) -> str: ...

    @property
    def source_version(self) -> str: ...

    @property
    def retrieved_at(self) -> date: ...

    @property
    def license(self) -> str: ...

    @property
    def redistribution_status(self) -> str: ...

    @property
    def identifier_namespace(self) -> str: ...

    @property
    def sequence_window(self) -> SequenceWindowDefinition | None: ...

    @property
    def bundle_id(self) -> str | None: ...

    @property
    def organism_common_name(self) -> str | None: ...

    @property
    def supports(self) -> tuple[str, ...]: ...

    @property
    def limitations(self) -> tuple[str, ...]: ...

    @property
    def reference_version(self) -> str | None: ...


class ReferenceBundleBuildRequestValidatorProtocol(Protocol):
    """Validator adapter consumed by reference bundle building."""

    def run(
        self,
        request: ReferenceBundleBuildRequest,
    ) -> ValidatedReferenceBundleBuildRequestProtocol: ...


@dataclass(frozen=True, slots=True)
class _ValidatedReferenceBundleBuildRequest:
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
    ) -> ValidatedReferenceBundleBuildRequestProtocol:
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
        return _ValidatedReferenceBundleBuildRequest(
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
                f"{field_name} must be a sequence of non-empty strings"
            )
        resolved: list[str] = []
        for item in value:
            resolved.append(cls._require_non_empty_string(item, field_name=field_name))
        return tuple(resolved)

    @staticmethod
    def _require_date(value: object, *, field_name: str) -> date:
        if isinstance(value, str):
            try:
                return date.fromisoformat(value.strip())
            except ValueError as exc:
                raise ReferenceResolutionError(
                    f"{field_name} must be an ISO date string"
                ) from exc
        if not isinstance(value, date):
            raise ReferenceResolutionError(f"{field_name} must be a datetime.date")
        return value

    @staticmethod
    def _validate_sequence_window(value: SequenceWindowDefinition) -> None:
        if not isinstance(value, SequenceWindowDefinition):
            raise ReferenceResolutionError(
                "reference bundle build request sequence_window must be a "
                "SequenceWindowDefinition"
            )


class ReferenceSourceTableReader(Protocol):
    """Reader protocol for local reference source tables."""

    def run(self, path: Path, *, field_name: str) -> pd.DataFrame: ...


class ReferenceBundleBuilder:
    """Build a validated ``ReferenceBundle`` from local kinase and sequence files."""

    def __init__(
        self,
        *,
        source_reader: ReferenceSourceTableReader,
        request_validator: ReferenceBundleBuildRequestValidatorProtocol | None = None,
    ) -> None:
        self._request_validator = (
            request_validator or ReferenceBundleBuildRequestValidator()
        )
        self._source_reader = source_reader

    def run(self, request: ReferenceBundleBuildRequest) -> ReferenceBundle:
        """Read local source files and return a normal validated reference bundle."""

        validated = self._request_validator.run(request)
        kinase_source = self._source_reader.run(
            validated.kinase_substrate_path,
            field_name="reference bundle kinase_substrate_path",
        )
        sequence_source = self._source_reader.run(
            validated.site_sequence_path,
            field_name="reference bundle site_sequence_path",
        )
        self._validate_source_organism(
            kinase_source,
            organism=validated.organism,
            field_name="reference bundle kinase-substrate source",
        )
        self._validate_source_organism(
            sequence_source,
            organism=validated.organism,
            field_name="reference bundle site-sequence source",
        )

        kinase_substrate_map = self._build_kinase_substrate_map(kinase_source)
        site_sequences = self._build_site_sequences(sequence_source)

        kinase_reference = KinaseSubstrateReference(frame=kinase_substrate_map)
        site_sequence_reference = SiteSequenceReference(frame=site_sequences)
        site_sequences = self._ensure_display_ids(site_sequence_reference.frame)
        sequence_window = self._resolve_sequence_window(
            validated,
            site_sequences=site_sequences,
        )
        manifest = self._build_manifest(
            validated,
            sequence_window=sequence_window,
            kinase_source=kinase_source,
            sequence_source=sequence_source,
        )
        identifier_normalisation = merge_reference_identifier_normalisation_reports(
            report
            for report in (
                kinase_reference.identifier_normalisation,
                site_sequence_reference.identifier_normalisation,
            )
            if report is not None
        )
        provenance = ReferenceProvenance(
            source_type="local",
            organism=validated.organism.value,
            bundle_id=manifest.bundle_id,
            source_name=manifest.source_name,
            source_version=manifest.source_version,
            retrieved_at=manifest.retrieved_at.isoformat(),
            identifier_namespace=manifest.identifier_namespace,
            sequence_window=manifest.sequence_window.to_payload(),
            manifest=manifest.to_payload(),
            table_fingerprints=(
                fingerprint_table(
                    kinase_reference.frame,
                    name="references.kinase_substrate_map",
                ),
                fingerprint_table(
                    site_sequences,
                    name="references.site_sequences",
                ),
            ),
            identifier_normalisation=identifier_normalisation,
            reference_context=ReferenceContext.from_manifest(manifest),
        )
        return ReferenceBundle._from_owned(  # pyright: ignore[reportPrivateUsage] - builder owns trusted transformation before public validation
            organism=validated.organism,
            kinase_substrate_map=kinase_reference.frame,
            site_sequences=site_sequences,
            provenance=provenance,
            manifest=manifest,
        )

    def _build_kinase_substrate_map(self, frame: pd.DataFrame) -> pd.DataFrame:
        field_name = "reference bundle kinase-substrate source"
        kinase_column = self._require_column(
            frame,
            aliases=_KINASE_ALIASES,
            semantic_name="kinase",
            field_name=field_name,
        )
        site_column = self._require_column(
            frame,
            aliases=_SUBSTRATE_SITE_ALIASES,
            semantic_name="substrate_site/site_id",
            field_name=field_name,
        )
        built = pd.DataFrame(
            {
                "kinase": self._normalise_text_series(
                    frame.loc[:, kinase_column],
                    uppercase=True,
                ),
                "substrate_site": self._normalise_text_series(
                    frame.loc[:, site_column],
                    uppercase=False,
                ),
            }
        )
        display_column = self._find_column(frame, _DISPLAY_ID_ALIASES)
        if display_column is not None:
            built.loc[:, "display_id"] = self._canonical_site_series(
                frame.loc[:, display_column],
                field_name=f"{field_name}.{display_column}",
            )
            self._require_matching_site_columns(
                left=built.loc[:, "substrate_site"],
                right=built.loc[:, "display_id"],
                left_name="substrate_site",
                right_name="display_id",
                field_name=field_name,
            )
        self._copy_optional_metadata_columns(frame, built)
        return built

    def _build_site_sequences(self, frame: pd.DataFrame) -> pd.DataFrame:
        field_name = "reference bundle site-sequence source"
        site_id_column = self._require_column(
            frame,
            aliases=_SITE_ID_ALIASES,
            semantic_name="site_id/substrate_site",
            field_name=field_name,
        )
        sequence_column = self._require_column(
            frame,
            aliases=_SITE_SEQUENCE_ALIASES,
            semantic_name="site_sequence",
            field_name=field_name,
        )
        built = pd.DataFrame(
            {
                "site_id": self._normalise_text_series(
                    frame.loc[:, site_id_column],
                    uppercase=False,
                ),
                "site_sequence": self._normalise_text_series(
                    frame.loc[:, sequence_column],
                    uppercase=True,
                ),
            }
        )
        display_column = self._find_column(frame, _DISPLAY_ID_ALIASES)
        if display_column is not None:
            built.loc[:, "display_id"] = self._canonical_site_series(
                frame.loc[:, display_column],
                field_name=f"{field_name}.{display_column}",
            )
            self._require_matching_site_columns(
                left=built.loc[:, "site_id"],
                right=built.loc[:, "display_id"],
                left_name="site_id",
                right_name="display_id",
                field_name=field_name,
            )
        self._copy_optional_metadata_columns(frame, built)
        built = built.set_index("site_id")
        built.index.name = "site_id"
        return built

    def _build_manifest(
        self,
        request: ValidatedReferenceBundleBuildRequestProtocol,
        *,
        sequence_window: SequenceWindowDefinition,
        kinase_source: pd.DataFrame,
        sequence_source: pd.DataFrame,
    ) -> ReferenceManifest:
        kinase_source_file = fingerprint_local_reference_source_file(
            request.kinase_substrate_path,
            role="kinase_substrate",
        )
        sequence_source_file = fingerprint_local_reference_source_file(
            request.site_sequence_path,
            role="site_sequences",
        )
        kinase_source_sha256 = str(kinase_source_file["sha256"])
        sequence_source_sha256 = str(sequence_source_file["sha256"])
        reference_version = request.reference_version
        if reference_version is None:
            reference_version = _generated_local_reference_version(
                kinase_sha256=kinase_source_sha256,
                sequence_sha256=sequence_source_sha256,
            )
        bundle_id = (
            request.bundle_id
            if request.bundle_id is not None
            else self._default_bundle_id(request)
        )
        files = (
            ReferenceFileManifest(
                relative_path=str(request.kinase_substrate_path),
                role="kinase_substrate",
                format=_reference_file_format(request.kinase_substrate_path),
                sha256=kinase_source_sha256,
                row_count=int(kinase_source.shape[0]),
                column_names=tuple(str(item) for item in kinase_source.columns),
            ),
            ReferenceFileManifest(
                relative_path=str(request.site_sequence_path),
                role="site_sequences",
                format=_reference_file_format(request.site_sequence_path),
                sha256=sequence_source_sha256,
                row_count=int(sequence_source.shape[0]),
                column_names=tuple(str(item) for item in sequence_source.columns),
            ),
        )
        return ReferenceManifest(
            reference_id=bundle_id,
            display_name=request.source_name,
            organism=(
                request.organism_common_name
                if request.organism_common_name is not None
                else request.organism.value
            ),
            taxonomy_id=_taxonomy_id(request.organism),
            organism_common_name=request.organism.value,
            protein_namespace=request.identifier_namespace,
            reference_version=reference_version,
            source_name=request.source_name,
            source_url=None,
            source_version=request.source_version,
            retrieved_at=request.retrieved_at,
            table_sha256=kinase_source_sha256,
            license_name=request.license,
            license_url=None,
            redistribution_status=_structured_redistribution_status(
                request.redistribution_status
            ),
            redistribution_notes=request.redistribution_status,
            derived_from=(
                str(request.kinase_substrate_path),
                str(request.site_sequence_path),
            ),
            generated_by="ReferenceBundleBuilder",
            generated_at_utc=f"{request.retrieved_at.isoformat()}T00:00:00Z",
            manifest_schema_version=REFERENCE_MANIFEST_SCHEMA_VERSION,
            files=files,
            sequence_context_policy=(
                "centered phosphosite sequence window"
                if sequence_window.central_residue_required
                else "sequence window"
            ),
            sequence_window_length=(
                int(sequence_window.upstream_residues)
                + 1
                + int(sequence_window.downstream_residues)
            ),
            sequence_center_index=int(sequence_window.upstream_residues),
            allowed_sequence_alphabet="ACDEFGHIKLMNPQRSTVWY",
            supports=request.supports,
            limitations=request.limitations,
        )

    @staticmethod
    def _resolve_sequence_window(
        request: ValidatedReferenceBundleBuildRequestProtocol,
        *,
        site_sequences: pd.DataFrame,
    ) -> SequenceWindowDefinition:
        if request.sequence_window is not None:
            return request.sequence_window
        sequence_values = [
            str(value).strip()
            for value in site_sequences.loc[:, "site_sequence"].tolist()
        ]
        lengths = {len(value) for value in sequence_values}
        if len(lengths) != 1:
            raise ReferenceResolutionError(
                "reference bundle build request sequence_window is required when "
                "site_sequence values do not all have the same length"
            )
        length = next(iter(lengths))
        if length <= 0 or length % 2 == 0:
            raise ReferenceResolutionError(
                "reference bundle build request sequence_window is required when "
                "site_sequence values are not non-empty odd-length centered windows"
            )
        center = length // 2
        central_residues = {value[center].upper() for value in sequence_values}
        if not central_residues.issubset(_CENTRAL_PHOSPHO_RESIDUES):
            residues = ", ".join(sorted(central_residues))
            raise ReferenceResolutionError(
                "reference bundle build request sequence_window is required when "
                "site_sequence central residues are not all S/T/Y; "
                f"central_residues={residues}"
            )
        return SequenceWindowDefinition(
            upstream_residues=center,
            downstream_residues=center,
            central_residue_required=True,
        )

    def _validate_source_organism(
        self,
        frame: pd.DataFrame,
        *,
        organism: Organism,
        field_name: str,
    ) -> None:
        organism_column = self._find_column(frame, _ORGANISM_ALIASES)
        if organism_column is None:
            return
        values = self._normalise_text_series(
            frame.loc[:, organism_column],
            uppercase=False,
        )
        blank_count = int((values == "").sum())
        if blank_count > 0:
            raise ReferenceResolutionError(
                f"{field_name}.{organism_column} contains blank organism values; "
                f"blank_count={blank_count}"
            )
        mismatches = [
            value
            for value in dict.fromkeys(values.tolist())
            if self._resolve_organism_token(value) is not organism
        ]
        if mismatches:
            preview = ", ".join(repr(value) for value in mismatches[:5])
            suffix = "" if len(mismatches) <= 5 else " ..."
            raise ReferenceResolutionError(
                f"{field_name}.{organism_column} does not match requested organism "
                f"{organism.value!r}; mismatched_values={preview}{suffix}"
            )

    @staticmethod
    def _ensure_display_ids(frame: pd.DataFrame) -> pd.DataFrame:
        if "display_id" in frame.columns:
            return frame
        enriched = frame.copy(deep=True)
        enriched.loc[:, "display_id"] = pd.Series(
            frame.index.astype(str).tolist(),
            index=frame.index.copy(),
            dtype="string",
        )
        return enriched

    def _copy_optional_metadata_columns(
        self,
        source: pd.DataFrame,
        target: pd.DataFrame,
    ) -> None:
        optional_columns: tuple[
            tuple[str, tuple[str, ...], Callable[[pd.Series], pd.Series]],
            ...,
        ] = (
            (
                "gene_symbol",
                _GENE_SYMBOL_ALIASES,
                lambda series: self._normalise_text_series(series, uppercase=True),
            ),
            (
                "protein_accession",
                _PROTEIN_ACCESSION_ALIASES,
                lambda series: self._normalise_text_series(series, uppercase=False),
            ),
            (
                "protein_id",
                _PROTEIN_ID_ALIASES,
                lambda series: self._normalise_text_series(series, uppercase=False),
            ),
        )
        for output_column, aliases, normaliser in optional_columns:
            if output_column in target.columns:
                continue
            source_column = self._find_column(source, aliases)
            if source_column is None:
                continue
            target.loc[:, output_column] = normaliser(source.loc[:, source_column])

    @staticmethod
    def _normalise_text_series(series: pd.Series, *, uppercase: bool) -> pd.Series:
        values = series.fillna("").astype(str).str.strip()
        if uppercase:
            values = values.str.upper()
        return pd.Series(values.tolist(), index=series.index.copy(), dtype="string")

    @staticmethod
    def _canonical_site_series(series: pd.Series, *, field_name: str) -> pd.Series:
        return canonicalize_site_series(
            series,
            field_name=field_name,
            error_type=ReferenceResolutionError,
        )

    def _require_matching_site_columns(
        self,
        *,
        left: pd.Series,
        right: pd.Series,
        left_name: str,
        right_name: str,
        field_name: str,
    ) -> None:
        left_canonical = self._canonical_site_series(
            left,
            field_name=f"{field_name}.{left_name}",
        )
        right_canonical = self._canonical_site_series(
            right,
            field_name=f"{field_name}.{right_name}",
        )
        mismatch_positions = [
            position
            for position, (left_value, right_value) in enumerate(
                zip(left_canonical.tolist(), right_canonical.tolist(), strict=False)
            )
            if left_value != right_value
        ]
        if not mismatch_positions:
            return
        preview = ", ".join(str(position) for position in mismatch_positions[:5])
        suffix = "" if len(mismatch_positions) <= 5 else " ..."
        raise ReferenceResolutionError(
            f"{field_name} has conflicting {left_name} and {right_name} values "
            f"after normalisation at row positions: {preview}{suffix}"
        )

    def _require_column(
        self,
        frame: pd.DataFrame,
        *,
        aliases: tuple[str, ...],
        semantic_name: str,
        field_name: str,
    ) -> str:
        column = self._find_column(frame, aliases)
        if column is not None:
            return column
        accepted = ", ".join(aliases)
        available = ", ".join(str(item) for item in frame.columns.tolist())
        if not available:
            available = "(none)"
        raise ReferenceResolutionError(
            f"{field_name} is missing required {semantic_name} column; "
            f"accepted aliases: {accepted}; available columns: {available}"
        )

    @staticmethod
    def _find_column(frame: pd.DataFrame, aliases: tuple[str, ...]) -> str | None:
        keyed_columns = {
            _normalise_column_token(column): str(column)
            for column in frame.columns.tolist()
        }
        for alias in aliases:
            column = keyed_columns.get(alias)
            if column is not None:
                return column
        return None

    @staticmethod
    def _resolve_organism_token(value: str) -> Organism | None:
        normalized = value.strip().lower()
        collapsed = _COLUMN_TOKEN_PATTERN.sub(" ", normalized).strip()
        underscored = collapsed.replace(" ", "_")
        for organism, tokens in _ORGANISM_TOKENS.items():
            if normalized in tokens or collapsed in tokens or underscored in tokens:
                return organism
        return None

    @staticmethod
    def _default_bundle_id(
        request: ValidatedReferenceBundleBuildRequestProtocol,
    ) -> str:
        raw = (
            f"local-{request.organism.value}-{request.source_name}-"
            f"{request.source_version}-{request.retrieved_at.isoformat()}"
        )
        normalized = _COLUMN_TOKEN_PATTERN.sub("-", raw.lower()).strip("-")
        return normalized or f"local-{request.organism.value}"


def _normalise_column_token(column: object) -> str:
    raw = str(column).strip().lower()
    return _COLUMN_TOKEN_PATTERN.sub("_", raw).strip("_")


def _taxonomy_id(organism: Organism) -> int:
    return {
        Organism.HUMAN: 9606,
        Organism.MOUSE: 10090,
        Organism.RAT: 10116,
    }[organism]


def _structured_redistribution_status(status: str) -> RedistributionStatus:
    try:
        resolved = RedistributionStatus(status.strip())
    except ValueError:
        return RedistributionStatus.UNRESOLVED
    if resolved is RedistributionStatus.APPROVED:
        return RedistributionStatus.UNRESOLVED
    return resolved


def _generated_local_reference_version(
    *,
    kinase_sha256: str,
    sequence_sha256: str,
) -> str:
    canonical = f"kinase_substrate:{kinase_sha256}\nsite_sequences:{sequence_sha256}\n"
    return f"local-snapshot-sha256-{sha256(canonical.encode('ascii')).hexdigest()}"


def _reference_file_format(path: ReferenceBuildPath) -> str:
    suffix = Path(cast(str | PathLike[str], path)).suffix.lower().lstrip(".")
    return suffix or "table"


__all__ = [
    "ReferenceBundleBuilder",
    "ReferenceBundleBuildRequestValidator",
    "ReferenceBundleBuildRequestValidatorProtocol",
    "ReferenceSourceTableReader",
    "ValidatedReferenceBundleBuildRequestProtocol",
]
