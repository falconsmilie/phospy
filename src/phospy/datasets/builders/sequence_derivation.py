"""Internal site-sequence derivation for dataset builder inputs."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

import pandas as pd

from phospy.errors.input import UnsupportedInputFormatError
from phospy.references.models import Organism
from phospy.references.resources import (
    bundled_reference_name_for_organism,
    load_bundled_reference_manifest,
    load_bundled_site_sequences,
)
from phospy.sites.identifiers import (
    canonicalize_site_components_series,
    canonicalize_site_index,
    canonicalize_site_series,
)

_SITE_SEQUENCE_DERIVATION_SCHEMA_VERSION = 1


@dataclass(frozen=True, slots=True)
class SiteSequenceDerivationReport:
    """Structured summary for builder-owned site-sequence derivation."""

    schema_version: int
    input_site_sequence_column_present: bool
    provided_sequence_count: int
    derived_sequence_count: int
    unresolved_sequence_count: int
    derivation_attempted: bool
    reference_source: str | None
    reference_bundle_id: str | None
    reference_manifest: Mapping[str, object] | None
    reference_support: str

    def to_payload(self) -> dict[str, object]:
        return {
            "schema_version": int(self.schema_version),
            "input_site_sequence_column_present": bool(
                self.input_site_sequence_column_present
            ),
            "provided_sequence_count": int(self.provided_sequence_count),
            "derived_sequence_count": int(self.derived_sequence_count),
            "unresolved_sequence_count": int(self.unresolved_sequence_count),
            "derivation_attempted": bool(self.derivation_attempted),
            "reference_source": self.reference_source,
            "reference_bundle_id": self.reference_bundle_id,
            "reference_manifest": (
                None
                if self.reference_manifest is None
                else dict(self.reference_manifest)
            ),
            "reference_support": str(self.reference_support),
        }


@dataclass(frozen=True, slots=True)
class _BundledDerivationAttempt:
    derived: pd.Series | None
    reference_source: str | None
    reference_bundle_id: str | None
    reference_manifest: Mapping[str, object] | None
    reference_support: str


class SiteSequenceDeriver:
    """Validate/provision `site_sequence` as optional builder enrichment."""

    def __init__(self) -> None:
        self._last_report: SiteSequenceDerivationReport | None = None

    @property
    def last_report(self) -> SiteSequenceDerivationReport | None:
        return self._last_report

    def run(
        self,
        site_metadata: pd.DataFrame,
        *,
        organism: Organism | None,
        allow_partial: bool = False,
        derive_missing_from_reference: bool = True,
    ) -> pd.DataFrame:
        normalized = site_metadata
        row_count = int(len(normalized.index))
        input_has_site_sequence = "site_sequence" in normalized.columns
        provided_sequence_count = 0
        derived_sequence_count = 0
        derivation_attempted = bool(
            derive_missing_from_reference and organism is not None
        )
        reference_source: str | None = None
        reference_bundle_id: str | None = None
        reference_manifest: Mapping[str, object] | None = None
        reference_support = "not_attempted"

        if "site_sequence" in normalized.columns and not allow_partial:
            existing = self._validated_existing_site_sequence(normalized)
            normalized.loc[:, "site_sequence"] = existing.astype(str)
            self._last_report = SiteSequenceDerivationReport(
                schema_version=_SITE_SEQUENCE_DERIVATION_SCHEMA_VERSION,
                input_site_sequence_column_present=input_has_site_sequence,
                provided_sequence_count=int(existing.shape[0]),
                derived_sequence_count=0,
                unresolved_sequence_count=0,
                derivation_attempted=derivation_attempted,
                reference_source=None,
                reference_bundle_id=None,
                reference_manifest=None,
                reference_support=(
                    "not_attempted"
                    if not derivation_attempted
                    else "available_without_missing_rows"
                ),
            )
            return normalized

        if "site_sequence" in normalized.columns:
            existing = self._normalized_optional_site_sequence(
                normalized.loc[:, "site_sequence"]
            )
            provided_sequence_count = int(existing.notna().sum())
            final_sequences = existing.copy(deep=True)
            if (
                derive_missing_from_reference
                and organism is not None
                and bool(existing.isna().any())
            ):
                attempt = self._derive_from_bundled_sequences_if_available(
                    site_metadata=normalized,
                    organism=organism,
                )
                reference_source = attempt.reference_source
                reference_bundle_id = attempt.reference_bundle_id
                reference_manifest = attempt.reference_manifest
                reference_support = attempt.reference_support
                if attempt.derived is not None:
                    derived_optional = self._normalized_optional_site_sequence(
                        attempt.derived
                    )
                    derived_mask = existing.isna() & derived_optional.notna()
                    derived_sequence_count = int(derived_mask.sum())
                    final_sequences = existing.where(
                        existing.notna(),
                        other=derived_optional,
                    )
            normalized.loc[:, "site_sequence"] = final_sequences.astype("string")
            unresolved_count = int(final_sequences.isna().sum())
            self._last_report = SiteSequenceDerivationReport(
                schema_version=_SITE_SEQUENCE_DERIVATION_SCHEMA_VERSION,
                input_site_sequence_column_present=input_has_site_sequence,
                provided_sequence_count=provided_sequence_count,
                derived_sequence_count=derived_sequence_count,
                unresolved_sequence_count=unresolved_count,
                derivation_attempted=derivation_attempted,
                reference_source=reference_source,
                reference_bundle_id=reference_bundle_id,
                reference_manifest=reference_manifest,
                reference_support=reference_support,
            )
            return normalized

        if not derive_missing_from_reference or organism is None:
            self._last_report = SiteSequenceDerivationReport(
                schema_version=_SITE_SEQUENCE_DERIVATION_SCHEMA_VERSION,
                input_site_sequence_column_present=input_has_site_sequence,
                provided_sequence_count=0,
                derived_sequence_count=0,
                unresolved_sequence_count=row_count,
                derivation_attempted=derivation_attempted,
                reference_source=None,
                reference_bundle_id=None,
                reference_manifest=None,
                reference_support="not_attempted",
            )
            return normalized

        attempt = self._derive_from_bundled_sequences_if_available(
            site_metadata=normalized,
            organism=organism,
        )
        reference_source = attempt.reference_source
        reference_bundle_id = attempt.reference_bundle_id
        reference_manifest = attempt.reference_manifest
        reference_support = attempt.reference_support
        if attempt.derived is None:
            self._last_report = SiteSequenceDerivationReport(
                schema_version=_SITE_SEQUENCE_DERIVATION_SCHEMA_VERSION,
                input_site_sequence_column_present=input_has_site_sequence,
                provided_sequence_count=0,
                derived_sequence_count=0,
                unresolved_sequence_count=row_count,
                derivation_attempted=derivation_attempted,
                reference_source=reference_source,
                reference_bundle_id=reference_bundle_id,
                reference_manifest=reference_manifest,
                reference_support=reference_support,
            )
            return normalized
        if allow_partial:
            derived_optional = self._normalized_optional_site_sequence(attempt.derived)
            derived_sequence_count = int(derived_optional.notna().sum())
            unresolved_count = int(derived_optional.isna().sum())
            normalized.loc[:, "site_sequence"] = derived_optional.astype("string")
            self._last_report = SiteSequenceDerivationReport(
                schema_version=_SITE_SEQUENCE_DERIVATION_SCHEMA_VERSION,
                input_site_sequence_column_present=input_has_site_sequence,
                provided_sequence_count=0,
                derived_sequence_count=derived_sequence_count,
                unresolved_sequence_count=unresolved_count,
                derivation_attempted=derivation_attempted,
                reference_source=reference_source,
                reference_bundle_id=reference_bundle_id,
                reference_manifest=reference_manifest,
                reference_support=reference_support,
            )
            return normalized
        unresolved = attempt.derived.isna() | (attempt.derived == "")
        if unresolved.any():
            self._last_report = SiteSequenceDerivationReport(
                schema_version=_SITE_SEQUENCE_DERIVATION_SCHEMA_VERSION,
                input_site_sequence_column_present=input_has_site_sequence,
                provided_sequence_count=0,
                derived_sequence_count=int((~unresolved).sum()),
                unresolved_sequence_count=int(unresolved.sum()),
                derivation_attempted=derivation_attempted,
                reference_source=reference_source,
                reference_bundle_id=reference_bundle_id,
                reference_manifest=reference_manifest,
                reference_support=reference_support,
            )
            return normalized
        normalized.loc[:, "site_sequence"] = attempt.derived.astype(str)
        self._last_report = SiteSequenceDerivationReport(
            schema_version=_SITE_SEQUENCE_DERIVATION_SCHEMA_VERSION,
            input_site_sequence_column_present=input_has_site_sequence,
            provided_sequence_count=0,
            derived_sequence_count=int(attempt.derived.shape[0]),
            unresolved_sequence_count=0,
            derivation_attempted=derivation_attempted,
            reference_source=reference_source,
            reference_bundle_id=reference_bundle_id,
            reference_manifest=reference_manifest,
            reference_support=reference_support,
        )
        return normalized

    @staticmethod
    def _validated_existing_site_sequence(site_metadata: pd.DataFrame) -> pd.Series:
        column = site_metadata.loc[:, "site_sequence"]
        as_string = column.astype("string").str.strip()
        invalid = as_string.isna() | (as_string == "")
        if invalid.any():
            preview = site_metadata.index[invalid.to_numpy()].astype(str).tolist()[:5]
            joined_preview = ", ".join(preview)
            raise UnsupportedInputFormatError(
                "dataset build request site_metadata.site_sequence must contain "
                "non-empty string values. remove blank values or drop the "
                f"'site_sequence' column to enable derivation; invalid sites: "
                f"{joined_preview}"
            )
        return as_string

    @staticmethod
    def _normalized_optional_site_sequence(column: pd.Series) -> pd.Series:
        as_string = column.astype("string").str.strip()
        missing = column.isna() | as_string.isna() | (as_string == "")
        return as_string.where(~missing, other=pd.NA)

    @staticmethod
    def _derive_from_bundled_sequences_if_available(
        *,
        site_metadata: pd.DataFrame,
        organism: Organism,
    ) -> _BundledDerivationAttempt:
        reference_name = bundled_reference_name_for_organism(organism)
        bundled_sequences = load_bundled_site_sequences(organism)
        manifest = load_bundled_reference_manifest(organism)
        sequence_map = bundled_sequences.loc[:, "site_sequence"]
        sequence_map.index = canonicalize_site_index(
            sequence_map.index,
            field_name="bundled site sequence index",
            error_type=UnsupportedInputFormatError,
        )
        site_lookup_index = _resolve_row_level_site_lookup_index(site_metadata)
        derived = sequence_map.reindex(
            pd.Index(site_lookup_index.astype("object").tolist())
        )
        derived.index = site_metadata.index.copy()
        derived.name = "site_sequence"
        return _BundledDerivationAttempt(
            derived=derived,
            reference_source=(
                f"bundled_reference:{organism.value}/{reference_name}"
                "/site_sequences.csv"
            ),
            reference_bundle_id=manifest.bundle_id,
            reference_manifest=manifest.to_payload(),
            reference_support="available",
        )


def _resolve_row_level_site_lookup_index(site_metadata: pd.DataFrame) -> pd.Series:
    metadata_site_ids = _resolve_site_ids_from_metadata_columns(site_metadata)
    index_site_ids = _resolve_site_ids_from_index(site_metadata.index)
    return metadata_site_ids.where(metadata_site_ids.notna(), other=index_site_ids)


def _resolve_site_ids_from_metadata_columns(site_metadata: pd.DataFrame) -> pd.Series:
    if (
        "gene_symbol" not in site_metadata.columns
        or "site" not in site_metadata.columns
    ):
        return pd.Series(pd.NA, index=site_metadata.index.copy(), dtype="string")
    gene_symbol = site_metadata.loc[:, "gene_symbol"].astype("string")
    site = site_metadata.loc[:, "site"].astype("string")
    normalized_gene_symbol = gene_symbol.str.strip()
    normalized_site = site.str.strip()
    has_tokens = (
        normalized_gene_symbol.notna()
        & normalized_site.notna()
        & (normalized_gene_symbol != "")
        & (normalized_site != "")
    )
    normalized = pd.Series(pd.NA, index=site_metadata.index.copy(), dtype="string")
    if not bool(has_tokens.any()):
        return normalized
    canonical = canonicalize_site_components_series(
        gene_symbol=gene_symbol.loc[has_tokens],
        site=site.loc[has_tokens],
        field_name=(
            "dataset build request site_metadata.gene_symbol/site for sequence lookup"
        ),
        error_type=UnsupportedInputFormatError,
        output_name="site_id",
    ).astype("string")
    normalized.loc[has_tokens] = canonical.to_numpy(dtype="object", copy=False)
    return normalized


def _resolve_site_ids_from_index(index: pd.Index) -> pd.Series:
    index_series = pd.Series(index.tolist(), index=index.copy(), dtype="string")
    normalized = index_series.str.strip()
    has_tokens = normalized.notna() & (normalized != "")
    site_like_tokens = has_tokens & normalized.str.contains(";", regex=False)
    resolved = pd.Series(pd.NA, index=index.copy(), dtype="string")
    if not bool(site_like_tokens.any()):
        return resolved
    canonical = canonicalize_site_series(
        normalized.loc[site_like_tokens],
        field_name="dataset build request site_metadata.index for sequence lookup",
        error_type=UnsupportedInputFormatError,
    ).astype("string")
    resolved.loc[site_like_tokens] = canonical.to_numpy(dtype="object", copy=False)
    return resolved
