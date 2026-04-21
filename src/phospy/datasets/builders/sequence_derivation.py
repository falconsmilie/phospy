"""Internal site-sequence derivation for dataset builder inputs."""

from __future__ import annotations

import pandas as pd

from phospy.errors.input import UnsupportedInputFormatError
from phospy.errors.references import ReferenceResolutionError, UnsupportedOrganismError
from phospy.references.models import Organism
from phospy.references.resources import load_bundled_site_sequences
from phospy.site_ids import canonicalize_site_index


class SiteSequenceDeriver:
    """Validate/provision `site_sequence` as optional builder enrichment."""

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
                        site_index=normalized.index,
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
            site_index=normalized.index,
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
        site_index: pd.Index,
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
        normalized_site_index = canonicalize_site_index(
            site_index,
            field_name="dataset site_metadata.index",
            error_type=UnsupportedInputFormatError,
        )
        return sequence_map.reindex(normalized_site_index)
