from __future__ import annotations

import pandas as pd
import pandas.testing as pdt
import pytest

from phospy.api.configs import (
    DATASET_SITE_MATRIX_DUPLICATE_STRATEGY_AGGREGATE_MEAN,
    DATASET_SITE_MATRIX_DUPLICATE_STRATEGY_AGGREGATE_MEDIAN,
    DATASET_SITE_MATRIX_DUPLICATE_STRATEGY_FIRST,
    DATASET_SITE_MATRIX_DUPLICATE_STRATEGY_MAX_MEAN_SIGNAL,
)
from phospy.datasets.preprocessing.models import DuplicateSiteResolutionResult
from phospy.datasets.preprocessing.stages.site_matrix import (
    _apply_duplicate_site_policy,
)


def _duplicate_policy_inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.Series]:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 3.0, 5.0],
            "sample_b": [2.0, 4.0, 6.0],
        },
        index=pd.Index(["row_a", "row_b", "row_c"], name="source_row"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "MAPK14", "AKT1"],
            "site": ["Y182", "Y182", "T308"],
            "site_sequence": ["SEQ_A", "SEQ_B", "SEQ_C"],
            "protein_id": ["PROT_A", "PROT_B", "PROT_C"],
            "uid": ["A", "B", "C"],
        },
        index=phospho.index.copy(),
    )
    constructed_site_id = pd.Series(
        ["MAPK14;Y182;", "MAPK14;Y182;", "AKT1;T308;"],
        index=phospho.index.copy(),
        name="site_id",
    )
    return phospho, site_metadata, constructed_site_id


def test_apply_duplicate_site_policy_returns_structured_result() -> None:
    phospho, site_metadata, constructed_site_id = _duplicate_policy_inputs()
    result = _apply_duplicate_site_policy(
        phospho=phospho,
        site_metadata=site_metadata,
        constructed_site_id=constructed_site_id,
        duplicate_site_strategy=DATASET_SITE_MATRIX_DUPLICATE_STRATEGY_FIRST,
    )
    assert isinstance(result, DuplicateSiteResolutionResult)


def test_duplicate_site_policy_first_reports_retained_and_dropped_rows() -> None:
    phospho, site_metadata, constructed_site_id = _duplicate_policy_inputs()
    result = _apply_duplicate_site_policy(
        phospho=phospho,
        site_metadata=site_metadata,
        constructed_site_id=constructed_site_id,
        duplicate_site_strategy=DATASET_SITE_MATRIX_DUPLICATE_STRATEGY_FIRST,
    )

    table = result.duplicate_site_resolution.sort_values("source_row_id")
    assert table.shape[0] == 2
    assert {
        "site_id",
        "source_row_id",
        "retained",
        "resolution_strategy",
        "retained_reason",
        "dropped_reason",
        "observed_values",
        "mean_signal",
    }.issubset(set(table.columns))
    assert table.loc[table["source_row_id"] == "row_a", "retained"].item()
    assert not table.loc[table["source_row_id"] == "row_b", "retained"].item()
    assert "input order" in str(
        table.loc[table["source_row_id"] == "row_a", "retained_reason"].item()
    )


def test_duplicate_site_policy_max_mean_signal_reports_mean_and_selection() -> None:
    phospho, site_metadata, constructed_site_id = _duplicate_policy_inputs()
    result = _apply_duplicate_site_policy(
        phospho=phospho,
        site_metadata=site_metadata,
        constructed_site_id=constructed_site_id,
        duplicate_site_strategy=DATASET_SITE_MATRIX_DUPLICATE_STRATEGY_MAX_MEAN_SIGNAL,
    )

    table = result.duplicate_site_resolution.sort_values("source_row_id")
    retained_row_id = table.loc[table["retained"], "source_row_id"].astype(str).tolist()
    assert retained_row_id == ["row_b"]
    assert table.loc[
        table["source_row_id"] == "row_a", "mean_signal"
    ].item() == pytest.approx(1.5)
    assert table.loc[
        table["source_row_id"] == "row_b", "mean_signal"
    ].item() == pytest.approx(3.5)


def test_duplicate_site_policy_aggregate_mean_reports_all_contributing_rows() -> None:
    phospho, site_metadata, constructed_site_id = _duplicate_policy_inputs()
    result = _apply_duplicate_site_policy(
        phospho=phospho,
        site_metadata=site_metadata,
        constructed_site_id=constructed_site_id,
        duplicate_site_strategy=DATASET_SITE_MATRIX_DUPLICATE_STRATEGY_AGGREGATE_MEAN,
    )

    table = result.duplicate_site_resolution.sort_values("source_row_id")
    assert table["source_row_id"].tolist() == ["row_a", "row_b"]
    assert table["retained"].tolist() == [True, True]
    assert table["n_aggregated_rows"].tolist() == [2, 2]
    assert (table["resolution_strategy"] == "aggregate_mean").all()


def test_duplicate_site_policy_aggregate_median_reports_all_contributing_rows() -> None:
    phospho, site_metadata, constructed_site_id = _duplicate_policy_inputs()
    result = _apply_duplicate_site_policy(
        phospho=phospho,
        site_metadata=site_metadata,
        constructed_site_id=constructed_site_id,
        duplicate_site_strategy=DATASET_SITE_MATRIX_DUPLICATE_STRATEGY_AGGREGATE_MEDIAN,
    )

    table = result.duplicate_site_resolution.sort_values("source_row_id")
    assert table["source_row_id"].tolist() == ["row_a", "row_b"]
    assert table["retained"].tolist() == [True, True]
    assert table["n_aggregated_rows"].tolist() == [2, 2]
    assert (table["resolution_strategy"] == "aggregate_median").all()


def test_duplicate_site_policy_reports_metadata_conflicts() -> None:
    phospho, site_metadata, constructed_site_id = _duplicate_policy_inputs()
    result = _apply_duplicate_site_policy(
        phospho=phospho,
        site_metadata=site_metadata,
        constructed_site_id=constructed_site_id,
        duplicate_site_strategy=DATASET_SITE_MATRIX_DUPLICATE_STRATEGY_AGGREGATE_MEAN,
    )

    assert not result.metadata_conflicts.empty
    protein_conflict = result.metadata_conflicts.loc[
        result.metadata_conflicts["field"] == "protein_id"
    ]
    assert protein_conflict.shape[0] == 1
    assert protein_conflict.iloc[0]["site_id"] == "MAPK14;Y182;"
    assert protein_conflict.iloc[0]["n_distinct_values"] == 2
    assert protein_conflict.iloc[0]["source_row_ids"] == ("row_a", "row_b")
    assert (
        result.duplicate_site_resolution.loc[:, "metadata_conflict_detected"]
        .astype(bool)
        .all()
    )


@pytest.mark.parametrize(
    ("strategy", "expected_phospho", "expected_site_metadata"),
    [
        (
            DATASET_SITE_MATRIX_DUPLICATE_STRATEGY_FIRST,
            pd.DataFrame(
                {
                    "sample_a": [1.0, 5.0],
                    "sample_b": [2.0, 6.0],
                },
                index=pd.Index(["MAPK14;Y182;", "AKT1;T308;"], name="site_id"),
            ),
            pd.DataFrame(
                {
                    "gene_symbol": ["MAPK14", "AKT1"],
                    "site": ["Y182", "T308"],
                    "site_sequence": ["SEQ_A", "SEQ_C"],
                    "protein_id": ["PROT_A", "PROT_C"],
                    "uid": ["A", "C"],
                },
                index=pd.Index(["MAPK14;Y182;", "AKT1;T308;"], name="site_id"),
            ),
        ),
        (
            DATASET_SITE_MATRIX_DUPLICATE_STRATEGY_MAX_MEAN_SIGNAL,
            pd.DataFrame(
                {
                    "sample_a": [5.0, 3.0],
                    "sample_b": [6.0, 4.0],
                },
                index=pd.Index(["AKT1;T308;", "MAPK14;Y182;"], name="site_id"),
            ),
            pd.DataFrame(
                {
                    "gene_symbol": ["AKT1", "MAPK14"],
                    "site": ["T308", "Y182"],
                    "site_sequence": ["SEQ_C", "SEQ_B"],
                    "protein_id": ["PROT_C", "PROT_B"],
                    "uid": ["C", "B"],
                },
                index=pd.Index(["AKT1;T308;", "MAPK14;Y182;"], name="site_id"),
            ),
        ),
        (
            DATASET_SITE_MATRIX_DUPLICATE_STRATEGY_AGGREGATE_MEAN,
            pd.DataFrame(
                {
                    "sample_a": [2.0, 5.0],
                    "sample_b": [3.0, 6.0],
                },
                index=pd.Index(["MAPK14;Y182;", "AKT1;T308;"], name="site_id"),
            ),
            pd.DataFrame(
                {
                    "gene_symbol": ["MAPK14", "AKT1"],
                    "site": ["Y182", "T308"],
                    "site_sequence": ["SEQ_A", "SEQ_C"],
                    "protein_id": ["PROT_A", "PROT_C"],
                    "uid": ["A", "C"],
                },
                index=pd.Index(["MAPK14;Y182;", "AKT1;T308;"], name="site_id"),
            ),
        ),
        (
            DATASET_SITE_MATRIX_DUPLICATE_STRATEGY_AGGREGATE_MEDIAN,
            pd.DataFrame(
                {
                    "sample_a": [2.0, 5.0],
                    "sample_b": [3.0, 6.0],
                },
                index=pd.Index(["MAPK14;Y182;", "AKT1;T308;"], name="site_id"),
            ),
            pd.DataFrame(
                {
                    "gene_symbol": ["MAPK14", "AKT1"],
                    "site": ["Y182", "T308"],
                    "site_sequence": ["SEQ_A", "SEQ_C"],
                    "protein_id": ["PROT_A", "PROT_C"],
                    "uid": ["A", "C"],
                },
                index=pd.Index(["MAPK14;Y182;", "AKT1;T308;"], name="site_id"),
            ),
        ),
    ],
)
def test_duplicate_site_policy_preserves_existing_outputs_for_current_strategies(
    strategy: str,
    expected_phospho: pd.DataFrame,
    expected_site_metadata: pd.DataFrame,
) -> None:
    phospho, site_metadata, constructed_site_id = _duplicate_policy_inputs()
    result = _apply_duplicate_site_policy(
        phospho=phospho,
        site_metadata=site_metadata,
        constructed_site_id=constructed_site_id,
        duplicate_site_strategy=strategy,
    )

    pdt.assert_frame_equal(result.phospho, expected_phospho)
    pdt.assert_frame_equal(result.site_metadata, expected_site_metadata)
