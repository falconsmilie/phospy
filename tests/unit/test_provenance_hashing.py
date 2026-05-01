from __future__ import annotations

import numpy as np
import pandas as pd

from phospy.provenance.hashing import fingerprint_table, hash_table


def _base_table() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample_a": [1.0, np.nan, 3.0],
            "sample_b": [4.0, 5.0, 6.0],
        },
        index=pd.Index(["A;S1;", "B;S2;", "C;S3;"], name="site_id"),
    )


def test_hash_is_deterministic_for_identical_table() -> None:
    table = _base_table()
    assert hash_table(table, name="dataset.phospho") == hash_table(
        table.copy(deep=True),
        name="dataset.phospho",
    )


def test_hash_changes_when_row_order_changes() -> None:
    table = _base_table()
    reordered = table.iloc[[2, 1, 0], :]
    assert hash_table(table, name="dataset.phospho") != hash_table(
        reordered,
        name="dataset.phospho",
    )


def test_hash_changes_when_column_order_changes() -> None:
    table = _base_table()
    reordered = table.loc[:, ["sample_b", "sample_a"]]
    assert hash_table(table, name="dataset.phospho") != hash_table(
        reordered,
        name="dataset.phospho",
    )


def test_hash_changes_when_value_changes() -> None:
    table = _base_table()
    changed = table.copy(deep=True)
    changed.loc["A;S1;", "sample_a"] = 9.0
    assert hash_table(table, name="dataset.phospho") != hash_table(
        changed,
        name="dataset.phospho",
    )


def test_hash_changes_when_dtype_changes() -> None:
    table = _base_table()
    changed = table.astype({"sample_b": "int64"})
    assert hash_table(table, name="dataset.phospho") != hash_table(
        changed,
        name="dataset.phospho",
    )


def test_hash_distinguishes_numeric_and_string_column_labels() -> None:
    numeric_columns = pd.DataFrame(
        {1: [1.0, 2.0]},
        index=pd.Index(["A", "B"], name="row_id"),
    )
    string_columns = pd.DataFrame(
        {"1": [1.0, 2.0]},
        index=pd.Index(["A", "B"], name="row_id"),
    )
    assert hash_table(numeric_columns, name="table") != hash_table(
        string_columns, name="table"
    )


def test_hash_distinguishes_numeric_and_string_index_labels() -> None:
    numeric_index = pd.DataFrame(
        {"x": [1.0, 2.0]},
        index=pd.Index([1, 2], name="row_id"),
    )
    string_index = pd.DataFrame(
        {"x": [1.0, 2.0]},
        index=pd.Index(["1", "2"], name="row_id"),
    )
    assert hash_table(numeric_index, name="table") != hash_table(
        string_index, name="table"
    )


def test_hash_changes_when_row_index_name_changes() -> None:
    first = pd.DataFrame(
        {"x": [1.0, 2.0]},
        index=pd.Index(["A", "B"], name="row_id"),
    )
    second = first.copy(deep=True)
    second.index = second.index.rename("site_id")
    assert hash_table(first, name="table") != hash_table(second, name="table")


def test_hash_changes_when_column_index_name_changes() -> None:
    first = pd.DataFrame(
        {"sample_a": [1.0], "sample_b": [2.0]},
        index=pd.Index(["A"], name="row_id"),
    )
    first.columns = first.columns.rename("sample_id")
    second = first.copy(deep=True)
    second.columns = second.columns.rename("run_id")
    assert hash_table(first, name="table") != hash_table(second, name="table")


def test_hash_distinguishes_range_index_from_equivalent_integer_index() -> None:
    range_index_table = pd.DataFrame({"x": [1.0, 2.0, 3.0]})
    integer_index_table = range_index_table.copy(deep=True)
    integer_index_table.index = pd.Index([0, 1, 2], dtype="int64")
    assert hash_table(range_index_table, name="table") != hash_table(
        integer_index_table, name="table"
    )


def test_hash_supports_multiindex_rows() -> None:
    first = pd.DataFrame(
        {"x": [1.0, 2.0]},
        index=pd.MultiIndex.from_tuples(
            [(1, "A"), (2, "B")], names=["batch", "replicate"]
        ),
    )
    second = first.copy(deep=True)
    second.index = pd.MultiIndex.from_tuples(
        [("1", "A"), (2, "B")], names=["batch", "replicate"]
    )
    assert hash_table(first, name="table") != hash_table(second, name="table")


def test_hash_supports_multiindex_columns() -> None:
    columns_first = pd.MultiIndex.from_tuples(
        [(1, "treated"), ("1", "control")], names=["sample_id", "group"]
    )
    columns_second = pd.MultiIndex.from_tuples(
        [("1", "treated"), ("1", "control")], names=["sample_id", "group"]
    )
    first = pd.DataFrame(
        [[1.0, 2.0], [3.0, 4.0]],
        index=pd.Index(["A", "B"], name="row_id"),
        columns=columns_first,
    )
    second = pd.DataFrame(
        [[1.0, 2.0], [3.0, 4.0]],
        index=pd.Index(["A", "B"], name="row_id"),
        columns=columns_second,
    )
    assert hash_table(first, name="table") != hash_table(second, name="table")


def test_hash_changes_when_display_is_identical_but_dtype_differs() -> None:
    int64_table = pd.DataFrame(
        {"x": pd.Series([1, 2], dtype="int64")},
        index=pd.Index(["A", "B"], name="row_id"),
    )
    nullable_int_table = pd.DataFrame(
        {"x": pd.Series([1, 2], dtype="Int64")},
        index=pd.Index(["A", "B"], name="row_id"),
    )
    assert hash_table(int64_table, name="table") != hash_table(
        nullable_int_table, name="table"
    )


def test_missing_value_representation_is_stable() -> None:
    first = pd.DataFrame(
        {"x": [1.0, np.nan]},
        index=pd.Index(["A", "B"], name="row_id"),
    )
    second = pd.DataFrame(
        {"x": [1.0, None]},
        index=pd.Index(["A", "B"], name="row_id"),
    )
    assert hash_table(first, name="table") == hash_table(second, name="table")


def test_fingerprint_captures_structural_metadata() -> None:
    table = _base_table()
    fingerprint = fingerprint_table(table, name="dataset.phospho")
    assert fingerprint.name == "dataset.phospho"
    assert fingerprint.rows == 3
    assert fingerprint.columns == 2
    assert fingerprint.index_name == "site_id"
    assert fingerprint.column_names == ("sample_a", "sample_b")
    assert fingerprint.index_structure is not None
    assert fingerprint.index_structure["type"] == "index"
    assert fingerprint.column_index_structure is not None
    assert fingerprint.column_index_structure["type"] == "index"
    assert fingerprint.hash_algorithm == "sha256"
    assert isinstance(fingerprint.hash_value, str)
    assert len(fingerprint.hash_value) == 64
