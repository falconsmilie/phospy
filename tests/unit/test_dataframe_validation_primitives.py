from __future__ import annotations

import pandas as pd
import pytest

from phospy.errors.validation import DatasetValidationError
from phospy.validation.common.dataframes import (
    require_aligned_dataframe_shape,
    require_columns,
    require_dataframe,
    require_exact_index_match,
    require_finite_numeric_dataframe,
    require_no_duplicate_labels,
    require_non_empty_dataframe,
    require_non_empty_index_intersection,
    require_numeric_dataframe,
    require_string_index,
    require_unique_columns,
    require_unique_index,
)


def test_require_dataframe_rejects_non_dataframe_input() -> None:
    with pytest.raises(DatasetValidationError, match="dataset.phospho"):
        require_dataframe(
            {"sample_a": [1.0]},
            field_name="dataset.phospho",
            allow_empty=False,
            error_type=DatasetValidationError,
        )


def test_require_non_empty_dataframe_rejects_empty_frame() -> None:
    empty = pd.DataFrame(columns=["sample_a"])
    with pytest.raises(DatasetValidationError, match="must be non-empty"):
        require_non_empty_dataframe(
            empty,
            field_name="dataset.phospho",
            error_type=DatasetValidationError,
        )


def test_require_unique_index_rejects_duplicate_index_labels() -> None:
    frame = pd.DataFrame(
        {"sample_a": [1.0, 2.0]},
        index=pd.Index(["MAPK14;Y182;", "MAPK14;Y182;"], name="site_id"),
    )
    with pytest.raises(DatasetValidationError, match="index must be unique"):
        require_unique_index(
            frame,
            field_name="dataset.phospho",
            error_type=DatasetValidationError,
        )


def test_require_unique_columns_rejects_duplicate_column_labels() -> None:
    frame = pd.DataFrame([[1.0, 2.0]], columns=["sample_a", "sample_a"])
    with pytest.raises(DatasetValidationError, match="columns must be unique"):
        require_unique_columns(
            frame,
            field_name="dataset.phospho",
            error_type=DatasetValidationError,
        )


def test_require_columns_rejects_missing_required_columns() -> None:
    frame = pd.DataFrame({"gene_symbol": ["MAPK14"]})
    with pytest.raises(DatasetValidationError, match="missing required columns: site"):
        require_columns(
            frame,
            field_name="dataset.site_metadata",
            required_columns=("gene_symbol", "site"),
            error_type=DatasetValidationError,
        )


def test_require_numeric_dataframe_rejects_non_numeric_columns() -> None:
    frame = pd.DataFrame({"sample_a": [1.0], "sample_b": ["x"]})
    with pytest.raises(DatasetValidationError, match="non-numeric columns: sample_b"):
        require_numeric_dataframe(
            frame,
            field_name="dataset.phospho",
            error_type=DatasetValidationError,
        )


def test_require_numeric_dataframe_rejects_boolean_columns_for_scientific_data() -> (
    None
):
    frame = pd.DataFrame({"sample_a": [1.0], "sample_b": [True]})
    with pytest.raises(DatasetValidationError, match="boolean columns are invalid"):
        require_numeric_dataframe(
            frame,
            field_name="dataset.phospho",
            error_type=DatasetValidationError,
        )


def test_require_finite_numeric_dataframe_rejects_infinite_values() -> None:
    frame = pd.DataFrame({"sample_a": [1.0, float("inf")]})
    with pytest.raises(DatasetValidationError, match="finite numeric values"):
        require_finite_numeric_dataframe(
            frame,
            field_name="dataset.phospho",
            error_type=DatasetValidationError,
            allow_missing=True,
        )


def test_require_non_empty_index_intersection_rejects_disjoint_indexes() -> None:
    with pytest.raises(DatasetValidationError, match="shared_count=0"):
        require_non_empty_index_intersection(
            left=pd.Index(["A;S1;", "B;S2;"]),
            right=pd.Index(["C;S3;", "D;S4;"]),
            left_name="prediction_result.pred_mat.index",
            right_name="dataset.phospho.index",
            error_type=DatasetValidationError,
        )


def test_require_exact_index_match_rejects_mismatched_index_order() -> None:
    with pytest.raises(DatasetValidationError, match="must exactly match"):
        require_exact_index_match(
            left=pd.Index(["A;S1;", "B;S2;"]),
            right=pd.Index(["B;S2;", "A;S1;"]),
            left_name="prediction_result.pred_mat.index",
            right_name="dataset.phospho.index",
            error_type=DatasetValidationError,
        )


def test_require_no_duplicate_labels_rejects_duplicate_label_sequences() -> None:
    with pytest.raises(DatasetValidationError, match="duplicate_count=2"):
        require_no_duplicate_labels(
            pd.Index(["sample_a", "sample_a"]),
            field_name="dataset.phospho.columns",
            error_type=DatasetValidationError,
        )


def test_require_string_index_rejects_non_string_labels() -> None:
    with pytest.raises(DatasetValidationError, match="non_string_label_count=1"):
        require_string_index(
            pd.Index(["MAPK14;Y182;", 3]),
            field_name="dataset.phospho.index",
            error_type=DatasetValidationError,
        )


def test_require_aligned_dataframe_shape_rejects_shape_mismatch() -> None:
    left = pd.DataFrame({"a": [1.0], "b": [2.0]})
    right = pd.DataFrame({"a": [1.0]})
    with pytest.raises(DatasetValidationError, match="shape must align"):
        require_aligned_dataframe_shape(
            left=left,
            right=right,
            left_name="prediction_result.pred_mat",
            right_name="scoring_result.profile_scores",
            error_type=DatasetValidationError,
        )
