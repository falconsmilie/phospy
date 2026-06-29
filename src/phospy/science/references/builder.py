"""Build validated reference bundles from local source files."""

from __future__ import annotations

import re
from collections.abc import Callable
from pathlib import Path

import pandas as pd

from phospy.errors.references import ReferenceResolutionError
from phospy.io.bundles.reference_sources import ReferenceSourceTableReader
from phospy.provenance.hashing import fingerprint_table
from phospy.provenance.models import ReferenceProvenance
from phospy.provenance.references import fingerprint_local_reference_source_file
from phospy.science.references.identifiers import (
    merge_reference_identifier_normalisation_reports,
)
from phospy.science.references.models import (
    Organism,
    ReferenceBuildPath,
    ReferenceBundle,
    ReferenceBundleBuildRequest,
    ReferenceFileManifest,
    ReferenceManifest,
    SequenceWindowDefinition,
)
from phospy.science.sites.identifiers import canonicalize_site_series
from phospy.tables.references import KinaseSubstrateReference, SiteSequenceReference
from phospy.validation.references.builder import (
    ReferenceBundleBuildRequestValidator,
    ValidatedReferenceBundleBuildRequest,
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


class ReferenceBundleBuilder:
    """Build a validated ``ReferenceBundle`` from local kinase and sequence files."""

    def __init__(
        self,
        *,
        request_validator: ReferenceBundleBuildRequestValidator | None = None,
        source_reader: ReferenceSourceTableReader | None = None,
    ) -> None:
        self._request_validator = (
            request_validator or ReferenceBundleBuildRequestValidator()
        )
        self._source_reader = source_reader or ReferenceSourceTableReader()

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
            source_version=manifest.reference_version,
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
        request: ValidatedReferenceBundleBuildRequest,
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
        bundle_id = (
            request.bundle_id
            if request.bundle_id is not None
            else self._default_bundle_id(request)
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
            reference_version=request.source_version,
            source_name=request.source_name,
            source_version=request.source_version,
            source_license=request.license,
            redistribution_allowed=_redistribution_allowed_from_status(
                request.redistribution_status
            ),
            redistribution_notes=request.redistribution_status,
            derived_from=(
                str(request.kinase_substrate_path),
                str(request.site_sequence_path),
            ),
            generated_by="ReferenceBundleBuilder",
            generated_at_utc=f"{request.retrieved_at.isoformat()}T00:00:00Z",
            manifest_schema_version="1.0",
            files=(
                ReferenceFileManifest(
                    relative_path=str(request.kinase_substrate_path),
                    role="kinase_substrate",
                    format=_reference_file_format(request.kinase_substrate_path),
                    sha256=str(kinase_source_file["sha256"]),
                    row_count=int(kinase_source.shape[0]),
                    column_names=tuple(str(item) for item in kinase_source.columns),
                ),
                ReferenceFileManifest(
                    relative_path=str(request.site_sequence_path),
                    role="site_sequences",
                    format=_reference_file_format(request.site_sequence_path),
                    sha256=str(sequence_source_file["sha256"]),
                    row_count=int(sequence_source.shape[0]),
                    column_names=tuple(str(item) for item in sequence_source.columns),
                ),
            ),
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
        request: ValidatedReferenceBundleBuildRequest,
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
    def _default_bundle_id(request: ValidatedReferenceBundleBuildRequest) -> str:
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


def _redistribution_allowed_from_status(status: str) -> bool:
    normalized = status.strip().lower()
    disallowed_tokens = ("not", "unknown", "unclear", "restricted", "prohibited")
    allowed_tokens = ("redistributable", "redistribution allowed", "approved", "cc0")
    return any(token in normalized for token in allowed_tokens) and not any(
        token in normalized for token in disallowed_tokens
    )


def _reference_file_format(path: ReferenceBuildPath) -> str:
    suffix = Path(path).suffix.lower().lstrip(".")
    return suffix or "table"


__all__ = ["ReferenceBundleBuilder"]
