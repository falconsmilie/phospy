"""Protein identifier resolution for interpreted signalome sites."""

from __future__ import annotations

import pandas as pd

from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.workflows.signalome.boundary_errors import raise_signalome_boundary_error
from phospy.workflows.signalome.constants import (
    SIGNALOME_INTERPRETER_PROTEIN_MAPPING_SEAM,
    SIGNALOME_PROTEIN_RESOLUTION_SOURCE_SITE_METADATA,
)


class SignalomeProteinResolver:
    """Resolve retained sites to signalome protein grouping labels."""

    _PROTEIN_COLUMN = "protein_id"

    def run(
        self,
        *,
        dataset: AnalysisReadyPhosphoDataset,
        site_index: pd.Index,
        removed_by_score_preconditioning_count: int,
    ) -> pd.Series:
        metadata = dataset._borrow_site_metadata_frame()
        if self._PROTEIN_COLUMN not in metadata.columns:
            raise_signalome_boundary_error(
                seam=SIGNALOME_INTERPRETER_PROTEIN_MAPPING_SEAM,
                next_action=(
                    "populate signalome protein grouping metadata in "
                    "dataset.site_metadata.protein_id for retained signalome sites "
                    "after score preconditioning"
                ),
                protein_resolution_source=SIGNALOME_PROTEIN_RESOLUTION_SOURCE_SITE_METADATA,
                interpreted_sites=int(site_index.size),
                resolved_protein_sites=0,
                unresolved_protein_sites=int(site_index.size),
                removed_by_score_preconditioning=(
                    int(removed_by_score_preconditioning_count)
                ),
                retained_with_missing_protein_id=int(site_index.size),
                retained_and_valid=0,
            )
        aligned_metadata = metadata
        if not site_index.isin(metadata.index).all():
            raise_signalome_boundary_error(
                seam=SIGNALOME_INTERPRETER_PROTEIN_MAPPING_SEAM,
                next_action=(
                    "ensure retained signalome site_key labels align with "
                    "dataset.site_metadata.index; signalome assumes dataset "
                    "site_key identity is already valid and does not repair it"
                ),
                protein_resolution_source=SIGNALOME_PROTEIN_RESOLUTION_SOURCE_SITE_METADATA,
                interpreted_sites=int(site_index.size),
                resolved_protein_sites=0,
                unresolved_protein_sites=int(site_index.size),
                removed_by_score_preconditioning=(
                    int(removed_by_score_preconditioning_count)
                ),
                retained_with_missing_protein_id=int(site_index.size),
                retained_and_valid=0,
            )
        resolved = (
            aligned_metadata.reindex(site_index)
            .loc[:, self._PROTEIN_COLUMN]
            .fillna("")
            .astype(str)
            .str.strip()
        )
        unresolved_mask = resolved.astype(str).str.strip() == ""
        if unresolved_mask.any():
            resolved_sites = int((~unresolved_mask).sum())
            unresolved_sites = resolved.loc[unresolved_mask].index.astype(str).tolist()
            raise_signalome_boundary_error(
                seam=SIGNALOME_INTERPRETER_PROTEIN_MAPPING_SEAM,
                next_action=(
                    "populate signalome protein grouping metadata in "
                    "dataset.site_metadata.protein_id for retained signalome sites "
                    "after score preconditioning"
                ),
                protein_resolution_source=SIGNALOME_PROTEIN_RESOLUTION_SOURCE_SITE_METADATA,
                interpreted_sites=int(site_index.size),
                resolved_protein_sites=resolved_sites,
                unresolved_protein_sites=int(unresolved_mask.sum()),
                retained_with_missing_protein_id=int(unresolved_mask.sum()),
                retained_and_valid=resolved_sites,
                removed_by_score_preconditioning=(
                    int(removed_by_score_preconditioning_count)
                ),
                missing_protein_id_sites=unresolved_sites,
            )
        resolved.index = site_index.copy()
        resolved.name = self._PROTEIN_COLUMN
        return resolved.astype(str)


__all__ = ["SignalomeProteinResolver"]
