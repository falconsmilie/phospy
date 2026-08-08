from __future__ import annotations

import pandas as pd
import pandas.testing as pdt
import pytest

from phospy.advanced import (
    filter_differential_results,
    rank_differential_results,
)
from phospy.errors import PhosPyInputError


def _result_table() -> pd.DataFrame:
    index = pd.Index(
        ["site_a", "site_b", "site_c", "site_d", "site_e"],
        name="site_key",
    )
    return pd.DataFrame(
        {
            "site_key": index.tolist(),
            "display_id": [
                "GENE1;S1;",
                "GENE2;S2;",
                "GENE3;S3;",
                "GENE4;S4;",
                "GENE5;S5;",
            ],
            "logFC": [1.2, -0.7, 2.5, -1.8, -1.2],
            "P.Value": [0.02, 0.20, 0.03, 0.03, 0.10],
            "adj.P.Val": [0.04, 0.20, 0.08, 0.03, 0.10],
        },
        index=index,
    )


def test_filter_differential_results_by_adjusted_p_value() -> None:
    filtered = filter_differential_results(
        _result_table(),
        adjusted_p_value_max=0.05,
    )

    assert filtered.index.tolist() == ["site_a", "site_d"]


def test_filter_differential_results_by_raw_p_value() -> None:
    filtered = filter_differential_results(
        _result_table(),
        p_value_max=0.05,
    )

    assert filtered.index.tolist() == ["site_a", "site_c", "site_d"]


def test_filter_differential_results_by_effect_size_threshold() -> None:
    filtered = filter_differential_results(
        _result_table(),
        min_abs_effect_size=1.5,
    )

    assert filtered.index.tolist() == ["site_c", "site_d"]


def test_filter_differential_results_combines_p_value_and_effect_size() -> None:
    filtered = filter_differential_results(
        _result_table(),
        adjusted_p_value_max=0.05,
        min_abs_effect_size=1.5,
    )

    assert filtered.index.tolist() == ["site_d"]


def test_filter_differential_results_does_not_mutate_input() -> None:
    table = _result_table()
    expected = table.copy(deep=True)

    filtered = filter_differential_results(
        table,
        adjusted_p_value_max=0.05,
    )
    filtered.loc[:, "logFC"] = 99.0

    pdt.assert_frame_equal(table, expected)


def test_rank_differential_results_ascending_and_descending() -> None:
    table = _result_table()

    ascending = rank_differential_results(table, by="P.Value")
    descending = rank_differential_results(table, by="logFC", ascending=False)

    assert ascending.index.tolist() == [
        "site_a",
        "site_c",
        "site_d",
        "site_e",
        "site_b",
    ]
    assert descending.index.tolist() == [
        "site_c",
        "site_a",
        "site_b",
        "site_e",
        "site_d",
    ]


def test_rank_differential_results_by_absolute_value() -> None:
    ranked = rank_differential_results(
        _result_table(),
        by="logFC",
        ascending=False,
        absolute=True,
    )

    assert ranked.index.tolist() == [
        "site_c",
        "site_d",
        "site_a",
        "site_e",
        "site_b",
    ]


def test_differential_result_filter_missing_column_error() -> None:
    table = _result_table().drop(columns=["adj.P.Val"])

    with pytest.raises(PhosPyInputError, match="missing required columns: adj.P.Val"):
        filter_differential_results(table, adjusted_p_value_max=0.05)


def test_differential_result_rank_missing_column_error() -> None:
    table = _result_table().drop(columns=["logFC"])

    with pytest.raises(PhosPyInputError, match="missing required columns: logFC"):
        rank_differential_results(table, by="logFC")


def test_filter_differential_results_can_return_empty_table() -> None:
    filtered = filter_differential_results(
        _result_table(),
        adjusted_p_value_max=0.001,
    )

    assert filtered.empty
    assert filtered.columns.tolist() == _result_table().columns.tolist()
