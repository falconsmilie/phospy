"""Internal site-sequence derivation for dataset builder inputs."""

from __future__ import annotations

import pandas as pd

from phospy.errors.input import UnsupportedInputFormatError
from phospy.errors.references import ReferenceResolutionError, UnsupportedOrganismError
from phospy.references.models import Organism
from phospy.references.resources import load_bundled_site_sequences


class SiteSequenceDeriver:
    """Ensure site metadata contains populated `site_sequence` values."""

    def run(
        self,
        site_metadata: pd.DataFrame,
        *,
        organism: Organism | None,
    ) -> pd.DataFrame:
        normalized = site_metadata.copy(deep=True)
        existing = self._existing_site_sequence(normalized)
        missing = existing.isna() | (existing == "")
        if not missing.any():
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
        resolved = existing.copy()
        resolved.loc[missing] = derived.loc[missing]
        unresolved = resolved.isna() | (resolved == "")
        if unresolved.any():
            missing_sites = normalized.index[unresolved].astype(str).tolist()[:5]
            preview = ", ".join(missing_sites)
            raise UnsupportedInputFormatError(
                "dataset build request site_metadata is missing site_sequence for "
                f"{int(unresolved.sum())} sites after derivation: {preview}. provide "
                "site_sequence values or use supported identifiers for the selected "
                "organism"
            )
        normalized.loc[:, "site_sequence"] = resolved.astype(str)
        return normalized

    @staticmethod
    def _existing_site_sequence(site_metadata: pd.DataFrame) -> pd.Series:
        if "site_sequence" not in site_metadata.columns:
            return pd.Series(index=site_metadata.index, dtype="object")
        column = site_metadata.loc[:, "site_sequence"]
        as_string = column.astype("string").str.strip()
        return as_string.fillna("")

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
        sequence_map = bundled_sequences.loc[:, "site_sequence"].copy()
        sequence_map.index = pd.Index(
            sequence_map.index.astype(str).str.strip(),
            name=sequence_map.index.name,
        )
        normalized_site_index = pd.Index(
            site_index.astype(str).str.strip(), name=site_index.name
        )
        return sequence_map.reindex(normalized_site_index)
