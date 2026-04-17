"""Reference bundle validator."""

from __future__ import annotations

import pandas as pd

from phospy.errors.validation import ReferenceValidationError
from phospy.references.models import Organism
from phospy.validation.common.dataframes import (
    require_dataframe,
    require_unique_index,
)


class ReferenceBundleValidator:
    """Validate the stable `ReferenceBundle` contract."""

    def run(
        self,
        *,
        organism: Organism,
        kinase_substrate_map: pd.DataFrame,
        site_sequences: pd.DataFrame,
        dataset_organism: Organism | None = None,
    ) -> None:
        if not isinstance(organism, Organism):
            raise ReferenceValidationError(
                "references.organism must be an Organism enum value"
            )

        kinase_substrate_map_frame = require_dataframe(
            kinase_substrate_map,
            field_name="references.kinase_substrate_map",
            allow_empty=False,
            error_type=ReferenceValidationError,
        )
        site_sequences_frame = require_dataframe(
            site_sequences,
            field_name="references.site_sequences",
            allow_empty=False,
            error_type=ReferenceValidationError,
        )
        require_unique_index(
            site_sequences_frame,
            field_name="references.site_sequences",
            error_type=ReferenceValidationError,
        )

        if kinase_substrate_map_frame.empty:
            raise ReferenceValidationError(
                "references.kinase_substrate_map must be non-empty"
            )
        if site_sequences_frame.empty:
            raise ReferenceValidationError(
                "references.site_sequences must be non-empty"
            )
        if dataset_organism is not None and dataset_organism is not organism:
            raise ReferenceValidationError(
                "references.organism must match dataset.organism when both are present"
            )
