from __future__ import annotations

import pandas as pd
import pytest

from phospy.errors import PhosPyInputError
from phospy.io.readers._table_parsing import (
    build_row_ids,
    build_unique_feature_ids,
    first_list_token,
    is_missing,
    optional_text,
    parse_flag,
    raise_for_forbidden_flags,
    required_text,
    resolve_column,
    resolve_flag_series,
    resolve_intensity_columns,
    resolve_required_column,
    split_multi_value,
)


def _validate_optional_column_name(
    value: object,
    *,
    field_name: str,
) -> str | None:
    if value is None:
        return None
    if not isinstance(value, str) or value.strip() == "":
        raise PhosPyInputError(f"{field_name} must be a non-empty string or None")
    return value.strip()


def test_resolve_required_column_matches_normalised_labels() -> None:
    columns = pd.Index(["  Gene   Name  ", "Intensity A"])

    resolved = resolve_required_column(
        columns,
        explicit=None,
        candidates=("Gene Name",),
        field_name="demo column_mapping.gene_symbol",
        importer_label="Demo",
        validate_column_name=_validate_optional_column_name,
    )

    assert resolved == "  Gene   Name  "


def test_resolve_column_reports_ambiguous_normalised_labels() -> None:
    columns = pd.Index(["Gene   Name", "gene name"])

    with pytest.raises(
        PhosPyInputError,
        match="Demo importer found ambiguous source columns for 'Gene Name'",
    ):
        resolve_column(
            columns,
            explicit=None,
            candidates=("Gene Name",),
            field_name="demo column_mapping.gene_symbol",
            importer_label="Demo",
            required=True,
            validate_column_name=_validate_optional_column_name,
        )


def test_text_and_multi_value_helpers_share_missing_handling() -> None:
    assert is_missing(pd.NA) is True
    assert split_multi_value(" P1; P2, P3 ") == ["P1", "P2", "P3"]
    assert (
        first_list_token(" P1; P2 ", field_name="Demo column", row_position=0) == "P1"
    )
    assert (
        required_text("  value  ", field_name="Demo column", row_position=0) == "value"
    )
    assert optional_text("  ") is None

    with pytest.raises(
        PhosPyInputError,
        match="Demo column must not contain missing values; row_position=1",
    ):
        required_text(pd.NA, field_name="Demo column", row_position=1)


def test_flag_helpers_parse_supported_values_and_report_bad_tokens() -> None:
    source = pd.DataFrame({"flag": ["+", "", 0, 1, "no", "yes"]})

    flags = resolve_flag_series(source, column="flag", field_name="Demo flag")

    assert flags is not None
    assert flags.tolist() == [True, False, False, True, False, True]
    assert parse_flag(None, field_name="Demo flag", row_position=0) is False

    with pytest.raises(
        PhosPyInputError,
        match="Demo flag row_position=2 contains unsupported flag value",
    ):
        parse_flag("maybe", field_name="Demo flag", row_position=2)

    with pytest.raises(
        PhosPyInputError,
        match="Demo importer encountered 1 decoy row",
    ):
        raise_for_forbidden_flags(
            pd.Series([False, True]),
            policy="error",
            error_policy="error",
            importer_label="Demo",
            label="decoy",
        )


def test_intensity_resolution_supports_duplicate_sample_policy() -> None:
    source = pd.DataFrame(
        {
            "Intensity Control": [1.0],
            "LFQ intensity Control": [2.0],
            "Intensity Stim": [3.0],
        }
    )

    with pytest.raises(
        PhosPyInputError,
        match="Demo importer inferred multiple intensity columns",
    ):
        resolve_intensity_columns(
            source,
            None,
            intensity_column_prefixes=("Intensity ", "LFQ intensity "),
            importer_label="Demo",
            request_label="demo",
            mapping_class_name="DemoColumnMapping",
            reject_duplicate_inferred_sample_ids=True,
        )

    inferred = resolve_intensity_columns(
        source.loc[:, ["Intensity Control", "Intensity Stim"]],
        None,
        intensity_column_prefixes=("Intensity ",),
        importer_label="Demo",
        request_label="demo",
        mapping_class_name="DemoColumnMapping",
        reject_duplicate_inferred_sample_ids=True,
    )

    assert inferred == {
        "Intensity Control": "Control",
        "Intensity Stim": "Stim",
    }


def test_feature_id_builders_preserve_explicit_value_errors() -> None:
    source = pd.DataFrame({"row_id": ["r1", ""], "feature_id": ["f1", "f2"]})

    assert build_row_ids(
        source=source,
        explicit_column=None,
        protein_values=["P1", "P2"],
        site_values=["S1", "T2"],
        source_row_numbers=[1, 2],
        importer_label="Demo",
        generated_prefix="demo",
    ) == ["demo:P1:S1:row1", "demo:P2:T2:row2"]
    assert build_unique_feature_ids(
        source=source,
        explicit_column=None,
        source_row_numbers=[1, 2],
        importer_label="Demo",
        generated_prefix="demo",
    ) == ["demo_feature_1", "demo_feature_2"]

    with pytest.raises(
        PhosPyInputError,
        match="Demo row_id must contain non-empty values; row_position=1",
    ):
        build_row_ids(
            source=source,
            explicit_column="row_id",
            protein_values=["P1", "P2"],
            site_values=["S1", "T2"],
            source_row_numbers=[1, 2],
            importer_label="Demo",
            generated_prefix="demo",
        )
