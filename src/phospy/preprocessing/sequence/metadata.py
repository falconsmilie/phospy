"""Metadata mutation collaborator for site-sequence stage."""

from __future__ import annotations

import pandas as pd


class SiteSequenceMetadataUpdater:
    """Own stateful updates to site_metadata.site_sequence for the stage."""

    def __init__(
        self,
        *,
        site_metadata: pd.DataFrame,
        existing_site_sequence: pd.Series,
    ) -> None:
        self._site_metadata = site_metadata.copy(deep=True)
        self._updated_site_sequence = existing_site_sequence.copy(deep=True)

    def assign(self, *, row_id: object, site_sequence: str) -> None:
        self._updated_site_sequence.at[row_id] = site_sequence

    def build(self) -> pd.DataFrame:
        result = self._site_metadata.copy(deep=True)
        result.loc[:, "site_sequence"] = self._updated_site_sequence.astype("string")
        return result


__all__ = ["SiteSequenceMetadataUpdater"]
