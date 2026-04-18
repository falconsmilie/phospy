"""Reference bundle validator."""

from __future__ import annotations

import pandas as pd

from phospy.errors.validation import ReferenceValidationError
from phospy.references.models import Organism
from phospy.validation.common.dataframes import (
    require_canonical_site_index,
    require_canonical_site_series,
    require_canonical_string_column,
    require_columns,
    require_dataframe,
    require_non_empty_string_column,
    require_unique_index,
    require_unique_row_pairs,
)


class ReferenceBundleValidator:
    """Validate the stable `ReferenceBundle` contract."""

    def run(
        self,
        *,
        organism: Organism,
        kinase_substrate_map: pd.DataFrame,
        site_sequences: pd.DataFrame,
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
        require_columns(
            kinase_substrate_map_frame,
            field_name="references.kinase_substrate_map",
            required_columns=("kinase", "substrate_site"),
            error_type=ReferenceValidationError,
        )
        require_non_empty_string_column(
            kinase_substrate_map_frame,
            field_name="references.kinase_substrate_map",
            column_name="kinase",
            error_type=ReferenceValidationError,
        )
        require_canonical_string_column(
            kinase_substrate_map_frame,
            field_name="references.kinase_substrate_map",
            column_name="kinase",
            error_type=ReferenceValidationError,
        )
        require_non_empty_string_column(
            kinase_substrate_map_frame,
            field_name="references.kinase_substrate_map",
            column_name="substrate_site",
            error_type=ReferenceValidationError,
        )
        require_canonical_site_series(
            kinase_substrate_map_frame.loc[:, "substrate_site"],
            field_name="references.kinase_substrate_map.substrate_site",
            error_type=ReferenceValidationError,
        )
        require_unique_row_pairs(
            kinase_substrate_map_frame,
            field_name="references.kinase_substrate_map",
            column_names=("kinase", "substrate_site"),
            error_type=ReferenceValidationError,
        )
        require_columns(
            site_sequences_frame,
            field_name="references.site_sequences",
            required_columns=("site_sequence",),
            error_type=ReferenceValidationError,
        )
        require_non_empty_string_column(
            site_sequences_frame,
            field_name="references.site_sequences",
            column_name="site_sequence",
            error_type=ReferenceValidationError,
        )
        require_canonical_string_column(
            site_sequences_frame,
            field_name="references.site_sequences",
            column_name="site_sequence",
            error_type=ReferenceValidationError,
        )
        require_canonical_site_index(
            site_sequences_frame.index,
            field_name="references.site_sequences.index",
            error_type=ReferenceValidationError,
        )
        substrate_sites = set(
            kinase_substrate_map_frame.loc[:, "substrate_site"].tolist()
        )
        known_sites = set(site_sequences_frame.index.tolist())
        missing_sequences = sorted(substrate_sites.difference(known_sites))
        if missing_sequences:
            missing_sample = ", ".join(missing_sequences[:10])
            raise ReferenceValidationError(
                "references.site_sequences is missing sequence entries for "
                f"substrate sites in references.kinase_substrate_map: {missing_sample}"
            )
