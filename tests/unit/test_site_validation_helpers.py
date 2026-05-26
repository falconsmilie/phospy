from __future__ import annotations

import pandas as pd
import pytest

from phospy.errors.validation import DatasetValidationError
from phospy.science.sites.validation import (
    require_canonical_site_index,
    require_canonical_site_series,
    require_site_identity_coherence,
)


def test_require_canonical_site_index_accepts_strict_canonical_ids() -> None:
    index = pd.Index(["MAPK14;Y182;", "AKT1;T308;"], name="site_id")
    validated = require_canonical_site_index(
        index,
        field_name="dataset.phospho.index",
        error_type=DatasetValidationError,
    )
    assert validated.equals(index)


def test_require_canonical_site_index_rejects_non_canonical_ids() -> None:
    index = pd.Index([" mapk14 ; y182 ; "], name="site_id")
    with pytest.raises(
        DatasetValidationError,
        match="must contain canonical site identifiers in 'GENE;SITE;' format",
    ):
        require_canonical_site_index(
            index,
            field_name="dataset.phospho.index",
            error_type=DatasetValidationError,
        )


def test_require_canonical_site_index_non_strict_rejects_stripped_collisions() -> None:
    index = pd.Index(["MAPK14;Y182;", " MAPK14;Y182; "], name="site_id")
    with pytest.raises(
        DatasetValidationError,
        match="contains colliding site identifiers when stripped",
    ):
        require_canonical_site_index(
            index,
            field_name="dataset.phospho.index",
            error_type=DatasetValidationError,
            strict_supported_format=False,
        )


def test_require_canonical_site_series_rejects_non_canonical_ids() -> None:
    series = pd.Series(["mapk14;y182;"], dtype="object")
    with pytest.raises(
        DatasetValidationError,
        match="must contain canonical site identifiers in 'GENE;SITE;' format",
    ):
        require_canonical_site_series(
            series,
            field_name="activity_result.target_table.site_id",
            error_type=DatasetValidationError,
        )


def test_require_site_identity_coherence_reports_unparseable_and_mismatched_rows() -> (
    None
):
    site_index = pd.Index(["MAPK14;Y182;", 7], dtype="object")
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14"],
            "site": ["Y181"],
        },
        index=pd.Index(["MAPK14;Y182;"], name="site_id"),
    )
    with pytest.raises(
        DatasetValidationError,
        match="dataset site-identity coherence failed",
    ) as exc_info:
        require_site_identity_coherence(
            site_index=site_index,
            site_metadata=site_metadata,
            site_index_field_name="dataset.site_metadata.index",
            site_metadata_field_name="dataset.site_metadata",
            error_type=DatasetValidationError,
        )

    message = str(exc_info.value)
    assert "unparseable site IDs" in message
    assert "mismatched rows:" in message
    assert "MAPK14;Y182;" in message
