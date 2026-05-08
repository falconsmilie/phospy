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


@pytest.mark.parametrize(
    ("factory", "pattern"),
    [
        pytest.param(
            lambda: require_non_empty_dataframe(
                pd.DataFrame(columns=["sample_a"]),
                field_name="dataset.phospho",
                error_type=DatasetValidationError,
            ),
            "must be non-empty",
            id="non-empty-frame",
        ),
        pytest.param(
            lambda: require_unique_index(
                pd.DataFrame(
                    {"sample_a": [1.0, 2.0]},
                    index=pd.Index(["MAPK14;Y182;", "MAPK14;Y182;"], name="site_id"),
                ),
                field_name="dataset.phospho",
                error_type=DatasetValidationError,
            ),
            "index must be unique",
            id="unique-index",
        ),
        pytest.param(
            lambda: require_unique_columns(
                pd.DataFrame([[1.0, 2.0]], columns=["sample_a", "sample_a"]),
                field_name="dataset.phospho",
                error_type=DatasetValidationError,
            ),
            "columns must be unique",
            id="unique-columns",
        ),
        pytest.param(
            lambda: require_columns(
                pd.DataFrame({"gene_symbol": ["MAPK14"]}),
                field_name="dataset.site_metadata",
                required_columns=("gene_symbol", "site"),
                error_type=DatasetValidationError,
            ),
            "missing required columns: site",
            id="required-columns",
        ),
    ],
)
def test_dataframe_schema_primitives_reject_invalid_shapes(
    factory: object, pattern: str
) -> None:
    assert callable(factory)
    with pytest.raises(DatasetValidationError, match=pattern):
        factory()


@pytest.mark.parametrize(
    ("factory", "pattern"),
    [
        pytest.param(
            lambda: require_numeric_dataframe(
                pd.DataFrame({"sample_a": [1.0], "sample_b": ["x"]}),
                field_name="dataset.phospho",
                error_type=DatasetValidationError,
            ),
            "non-numeric columns: sample_b",
            id="numeric-only-non-numeric-column",
        ),
        pytest.param(
            lambda: require_numeric_dataframe(
                pd.DataFrame({"sample_a": [1.0], "sample_b": [True]}),
                field_name="dataset.phospho",
                error_type=DatasetValidationError,
            ),
            "boolean columns are invalid",
            id="numeric-only-boolean-column",
        ),
        pytest.param(
            lambda: require_finite_numeric_dataframe(
                pd.DataFrame({"sample_a": [1.0, float("inf")]}),
                field_name="dataset.phospho",
                error_type=DatasetValidationError,
                allow_missing=True,
            ),
            "finite numeric values",
            id="finite-values-infinite-rejected",
        ),
    ],
)
def test_numeric_dataframe_primitives_reject_invalid_content(
    factory: object, pattern: str
) -> None:
    assert callable(factory)
    with pytest.raises(DatasetValidationError, match=pattern):
        factory()


@pytest.mark.parametrize(
    ("factory", "pattern"),
    [
        pytest.param(
            lambda: require_non_empty_index_intersection(
                left=pd.Index(["A;S1;", "B;S2;"]),
                right=pd.Index(["C;S3;", "D;S4;"]),
                left_name="prediction_result.pred_mat.index",
                right_name="dataset.phospho.index",
                error_type=DatasetValidationError,
            ),
            "shared_count=0",
            id="non-empty-index-intersection",
        ),
        pytest.param(
            lambda: require_exact_index_match(
                left=pd.Index(["A;S1;", "B;S2;"]),
                right=pd.Index(["B;S2;", "A;S1;"]),
                left_name="prediction_result.pred_mat.index",
                right_name="dataset.phospho.index",
                error_type=DatasetValidationError,
            ),
            "must exactly match",
            id="exact-index-match-order",
        ),
        pytest.param(
            lambda: require_no_duplicate_labels(
                pd.Index(["sample_a", "sample_a"]),
                field_name="dataset.phospho.columns",
                error_type=DatasetValidationError,
            ),
            "duplicate_count=2",
            id="no-duplicate-labels",
        ),
        pytest.param(
            lambda: require_string_index(
                pd.Index(["MAPK14;Y182;", 3]),
                field_name="dataset.phospho.index",
                error_type=DatasetValidationError,
            ),
            "non_string_label_count=1",
            id="string-index",
        ),
        pytest.param(
            lambda: require_aligned_dataframe_shape(
                left=pd.DataFrame({"a": [1.0], "b": [2.0]}),
                right=pd.DataFrame({"a": [1.0]}),
                left_name="prediction_result.pred_mat",
                right_name="scoring_result.profile_scores",
                error_type=DatasetValidationError,
            ),
            "shape must align",
            id="aligned-dataframe-shape",
        ),
    ],
)
def test_dataframe_alignment_primitives_reject_mismatch_conditions(
    factory: object, pattern: str
) -> None:
    assert callable(factory)
    with pytest.raises(DatasetValidationError, match=pattern):
        factory()
