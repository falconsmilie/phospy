from __future__ import annotations

import pytest

from phospy.errors import PhosPyInputError
from phospy.validation.datasets.batch_correction_controls import (
    normalize_applied_selected_site_key_rows,
    require_control_site_source_metadata,
    require_unique_selected_control_site_rows,
)


def test_selected_control_site_rows_normalize_unique_site_keys() -> None:
    rows = normalize_applied_selected_site_key_rows(
        (" AKT1_S473 ", "GSK3B_S9"),
    )

    require_unique_selected_control_site_rows(rows)

    assert rows == ("AKT1_S473", "GSK3B_S9")


def test_selected_control_site_rows_reject_duplicate_site_keys() -> None:
    rows = normalize_applied_selected_site_key_rows(("AKT1_S473", "AKT1_S473"))

    with pytest.raises(
        PhosPyInputError,
        match="duplicate selected control row identifiers.*AKT1_S473",
    ):
        require_unique_selected_control_site_rows(rows)


def test_selected_control_site_rows_reject_sentinels() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="selected_site_key_rows contains sentinel.*\\[1\\]",
    ):
        normalize_applied_selected_site_key_rows(("AKT1_S473", "not provided"))


def test_control_site_source_metadata_accepts_supported_caller_source() -> None:
    require_control_site_source_metadata(
        {
            "source_type": "caller_supplied",
            "organism": "rat",
            "identifier_namespace": "site_key",
            "source_name": "manual-curated-controls",
            "source_version": "manual-v1",
            "license": "caller local use",
            "redistribution": "not redistributed",
        },
        selected_site_key_rows=("AKT1_S473", "GSK3B_S9"),
    )


def test_control_site_source_metadata_accepts_supported_strict_source() -> None:
    require_control_site_source_metadata(
        {
            "source_type": "packaged_reference",
            "organism": "rat",
            "identifier_namespace": "site_key",
            "source_name": "packaged-control-reference",
            "source_version": "2026-08",
            "license": "GPL-3.0",
            "redistribution": "redistributable with package",
            "selection_method": "packaged curated reference",
        },
        selected_site_key_rows=("AKT1_S473", "GSK3B_S9"),
    )


def test_control_site_source_metadata_accepts_field_missing_reason() -> None:
    require_control_site_source_metadata(
        {
            "source_type": "caller_supplied",
            "organism": "rat",
            "identifier_namespace": "site_key",
            "source_name": "local notebook",
            "redistribution": "not redistributed",
            "metadata_missing_reason": {
                "source_version": "local notebook was not versioned",
                "license": "caller-only local material",
            },
        },
        selected_site_key_rows=("AKT1_S473", "GSK3B_S9"),
    )


def test_control_site_source_metadata_rejects_missing_reason_metadata() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="control source is missing 'license' without explicit rationale",
    ):
        require_control_site_source_metadata(
            {
                "source_type": "caller_supplied",
                "organism": "rat",
                "identifier_namespace": "site_key",
                "source_name": "manual-curated-controls",
                "source_version": "manual-v1",
                "redistribution": "not redistributed",
            },
            selected_site_key_rows=("AKT1_S473", "GSK3B_S9"),
        )


def test_control_site_source_metadata_rejects_conflicting_source_declarations() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="conflicting source declarations.*caller_supplied.*dataset_metadata",
    ):
        require_control_site_source_metadata(
            {
                "source_type": "caller_supplied",
                "source": "dataset_metadata",
                "organism": "rat",
                "identifier_namespace": "site_key",
                "source_name": "manual-curated-controls",
                "source_version": "manual-v1",
                "license": "caller local use",
                "redistribution": "not redistributed",
            },
            selected_site_key_rows=("AKT1_S473", "GSK3B_S9"),
        )
