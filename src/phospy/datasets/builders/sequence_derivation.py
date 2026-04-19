"""Internal site-sequence derivation for dataset builder inputs."""

from __future__ import annotations

import pandas as pd

from phospy.errors.input import UnsupportedInputFormatError
from phospy.errors.references import ReferenceResolutionError, UnsupportedOrganismError
from phospy.references.models import Organism
from phospy.references.resources import load_bundled_site_sequences
from phospy.site_ids import canonicalize_site_index


class SiteSequenceDeriver:
    """Ensure site metadata contains populated `site_sequence` values."""

    def run(
        self,
        site_metadata: pd.DataFrame,
        *,
        organism: Organism | None,
    ) -> pd.DataFrame:
        normalized = site_metadata
        if "site_sequence" in normalized.columns:
            existing = self._validated_existing_site_sequence(normalized)
            normalized.loc[:, "site_sequence"] = existing.astype(str)
            return normalized
        if organism is None:
            raise UnsupportedInputFormatError(
                "dataset build request site_metadata is missing site_sequence values. "
                "provide site_sequence explicitly or set organism to enable bundled "
                "site-sequence derivation"
            )
        derived = self._derive_from_bundled_sequences(
            site_index=normalized.index,
            organism=organism,
        )
        unresolved = derived.isna() | (derived == "")
        if unresolved.any():
            missing_sites = normalized.index[unresolved].astype(str).tolist()[:5]
            preview = ", ".join(missing_sites)
            raise UnsupportedInputFormatError(
                "dataset build request site_metadata is missing site_sequence for "
                f"{int(unresolved.sum())} sites after derivation: {preview}. provide "
                "site_sequence values or use supported identifiers for the selected "
                "organism"
            )
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
    def _derive_from_bundled_sequences(
        *,
        site_index: pd.Index,
        organism: Organism,
    ) -> pd.Series:
        try:
            bundled_sequences = load_bundled_site_sequences(organism)
        except (UnsupportedOrganismError, ReferenceResolutionError) as exc:
            raise UnsupportedInputFormatError(
                f"site_sequence derivation is unavailable for organism "
                f"'{organism.value}'. provide site_metadata.site_sequence explicitly"
            ) from exc
        sequence_map = bundled_sequences.loc[:, "site_sequence"]
        sequence_map.index = canonicalize_site_index(
            sequence_map.index,
            field_name="bundled site sequence index",
            error_type=UnsupportedInputFormatError,
        )
        normalized_site_index = canonicalize_site_index(
            site_index,
            field_name="dataset site_metadata.index",
            error_type=UnsupportedInputFormatError,
        )
        return sequence_map.reindex(normalized_site_index)
