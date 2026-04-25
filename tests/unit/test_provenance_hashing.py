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
    assert fingerprint.hash_algorithm == "sha256"
    assert isinstance(fingerprint.hash_value, str)
    assert len(fingerprint.hash_value) == 64
