"""Dataset validators."""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from phospy.validation.datasets.analysis_ready import (
        AnalysisReadyDatasetModelBoundaryValidator,
    )
    from phospy.validation.datasets.display_site_identity import (
        DISPLAY_SITE_CONTEXT_COLUMNS,
        enforce_unique_display_site_identity_rows,
    )
    from phospy.validation.datasets.inputs import DatasetInputSourceValidator
    from phospy.validation.datasets.preprocessing import (
        DatasetPreprocessingConfigValidator,
    )
    from phospy.validation.datasets.protein_scoped_site_identity import (
        enforce_display_id_column,
        enforce_site_key_column,
        enforce_site_key_index,
        enforce_site_key_matches_metadata,
        enforce_unique_site_key_identity,
    )

__all__ = [
    "AnalysisReadyDatasetModelBoundaryValidator",
    "DISPLAY_SITE_CONTEXT_COLUMNS",
    "DatasetInputSourceValidator",
    "DatasetPreprocessingConfigValidator",
    "enforce_display_id_column",
    "enforce_site_key_column",
    "enforce_site_key_index",
    "enforce_site_key_matches_metadata",
    "enforce_unique_display_site_identity_rows",
    "enforce_unique_site_key_identity",
]


def __getattr__(name: str) -> object:
    if name == "AnalysisReadyDatasetModelBoundaryValidator":
        from phospy.validation.datasets.analysis_ready import (
            AnalysisReadyDatasetModelBoundaryValidator,
        )

        return AnalysisReadyDatasetModelBoundaryValidator
    if name == "DatasetInputSourceValidator":
        from phospy.validation.datasets.inputs import DatasetInputSourceValidator

        return DatasetInputSourceValidator
    if name == "DISPLAY_SITE_CONTEXT_COLUMNS":
        from phospy.validation.datasets.display_site_identity import (
            DISPLAY_SITE_CONTEXT_COLUMNS,
        )

        return DISPLAY_SITE_CONTEXT_COLUMNS
    if name == "enforce_unique_display_site_identity_rows":
        from phospy.validation.datasets.display_site_identity import (
            enforce_unique_display_site_identity_rows,
        )

        return enforce_unique_display_site_identity_rows
    if name == "DatasetPreprocessingConfigValidator":
        from phospy.validation.datasets.preprocessing import (
            DatasetPreprocessingConfigValidator,
        )

        return DatasetPreprocessingConfigValidator
    if name == "enforce_site_key_column":
        from phospy.validation.datasets.protein_scoped_site_identity import (
            enforce_site_key_column,
        )

        return enforce_site_key_column
    if name == "enforce_display_id_column":
        from phospy.validation.datasets.protein_scoped_site_identity import (
            enforce_display_id_column,
        )

        return enforce_display_id_column
    if name == "enforce_unique_site_key_identity":
        from phospy.validation.datasets.protein_scoped_site_identity import (
            enforce_unique_site_key_identity,
        )

        return enforce_unique_site_key_identity
    if name == "enforce_site_key_matches_metadata":
        from phospy.validation.datasets.protein_scoped_site_identity import (
            enforce_site_key_matches_metadata,
        )

        return enforce_site_key_matches_metadata
    if name == "enforce_site_key_index":
        from phospy.validation.datasets.protein_scoped_site_identity import (
            enforce_site_key_index,
        )

        return enforce_site_key_index
    raise AttributeError(name)
