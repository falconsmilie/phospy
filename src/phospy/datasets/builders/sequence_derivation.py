"""Internal site-sequence derivation for dataset builder inputs."""

from __future__ import annotations

import pandas as pd

from phospy.errors.input import UnsupportedInputFormatError
from phospy.errors.references import ReferenceResolutionError, UnsupportedOrganismError
from phospy.references.models import Organism
from phospy.references.resources import load_bundled_site_sequences
from phospy.site_ids import (
    canonicalize_site_components_series,
    canonicalize_site_index,
    canonicalize_site_series,
)


class SiteSequenceDeriver:
    """Validate/provision `site_sequence` as optional builder enrichment.

    In partial mode (`allow_partial=True`), enrichment is row-wise:
    rows with supplied or derivable sequence support keep their values, while
    unresolved rows remain missing for downstream policy-specific exclusion.
    """

    def run(
        self,
        site_metadata: pd.DataFrame,
        *,
        organism: Organism | None,
        allow_partial: bool = False,
    ) -> pd.DataFrame:
        normalized = site_metadata
        if "site_sequence" in normalized.columns:
            if allow_partial:
                existing = self._normalized_optional_site_sequence(
                    normalized.loc[:, "site_sequence"]
                )
                if organism is not None and bool(existing.isna().any()):
                    derived = self._derive_from_bundled_sequences_if_available(
                        site_metadata=normalized,
                        organism=organism,
                    )
                    if derived is not None:
                        derived_optional = self._normalized_optional_site_sequence(
                            derived
                        )
                        existing = existing.where(
                            existing.notna(), other=derived_optional
                        )
                normalized.loc[:, "site_sequence"] = existing.astype("string")
                return normalized
            existing = self._validated_existing_site_sequence(normalized)
            normalized.loc[:, "site_sequence"] = existing.astype(str)
            return normalized
        if organism is None:
            return normalized
        derived = self._derive_from_bundled_sequences_if_available(
            site_metadata=normalized,
            organism=organism,
        )
        if derived is None:
            return normalized
        if allow_partial:
            normalized.loc[:, "site_sequence"] = (
                self._normalized_optional_site_sequence(derived).astype("string")
            )
            return normalized
        unresolved = derived.isna() | (derived == "")
        if unresolved.any():
            return normalized
        normalized.loc[:, "site_sequence"] = derived.astype(str)
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
    ) -> pd.Series | None:
        try:
            bundled_sequences = load_bundled_site_sequences(organism)
        except (UnsupportedOrganismError, ReferenceResolutionError):
            return None
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
        return derived


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
    has_tokens = gene_symbol.notna() & site.notna() & (gene_symbol != "") & (site != "")
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
