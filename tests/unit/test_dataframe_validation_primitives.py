from __future__ import annotations

import pandas as pd
import pytest

from phospy.errors.validation import DatasetValidationError
from phospy.validation.common.dataframes import (
    format_label_examples,
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
    summarise_column_mismatch,
    summarise_index_mismatch,
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


def test_format_label_examples_caps_large_sequences() -> None:
    labels = [f"S{idx:03d}" for idx in range(1, 8)]
    formatted = format_label_examples(labels, max_examples=3)
    assert formatted == "'S001', 'S002', 'S003', +4 more"


def test_summarise_index_mismatch_reports_order_only() -> None:
    summary = summarise_index_mismatch(
        left=pd.Index(["S001", "S003", "S002"]),
        right=pd.Index(["S001", "S002", "S003"]),
        left_name="matrix.columns",
        right_name="metadata.index",
    )
    assert "Only in matrix.columns: (none)" in summary
    assert "Only in metadata.index: (none)" in summary
    assert "First positional mismatch: position 1" in summary
    assert "Labels match as a set but order differs: true" in summary


@pytest.mark.parametrize(
    ("left", "right", "left_fragment", "right_fragment"),
    [
        pytest.param(
            pd.Index(["S001"]),
            pd.Index(["S001", "S002"]),
            "Only in matrix.columns: (none)",
            "Only in metadata.index: 'S002'",
            id="missing-label",
        ),
        pytest.param(
            pd.Index(["S001", "S002"]),
            pd.Index(["S001"]),
            "Only in matrix.columns: 'S002'",
            "Only in metadata.index: (none)",
            id="extra-label",
        ),
        pytest.param(
            pd.Index(["S001", "S003"]),
            pd.Index(["S001", "S004"]),
            "Only in matrix.columns: 'S003'",
            "Only in metadata.index: 'S004'",
            id="missing-and-extra",
        ),
    ],
)
def test_summarise_index_mismatch_reports_missing_and_extra_labels(
    left: pd.Index,
    right: pd.Index,
    left_fragment: str,
    right_fragment: str,
) -> None:
    summary = summarise_index_mismatch(
        left=left,
        right=right,
        left_name="matrix.columns",
        right_name="metadata.index",
    )
    assert left_fragment in summary
    assert right_fragment in summary


def test_require_exact_index_match_reports_positional_and_set_diagnostics() -> None:
    with pytest.raises(DatasetValidationError) as exc_info:
        require_exact_index_match(
            left=pd.Index(["S001", "S003", "S004"]),
            right=pd.Index(["S001", "Sample_3", "Sample_4"]),
            left_name="matrix.columns",
            right_name="metadata.index",
            error_type=DatasetValidationError,
        )
    message = str(exc_info.value)
    assert "matrix.columns must exactly match metadata.index" in message
    assert "expected_length=3, actual_length=3" in message
    assert "Only in matrix.columns: 'S003', 'S004'" in message
    assert "Only in metadata.index: 'Sample_3', 'Sample_4'" in message
    assert (
        "First positional mismatch: position 1, matrix.columns='S003', "
        "metadata.index='Sample_3'"
    ) in message
    assert "Labels match as a set but order differs: false" in message


def test_require_exact_index_match_caps_large_examples() -> None:
    left = pd.Index([f"L{idx:03d}" for idx in range(1, 11)])
    right = pd.Index([f"R{idx:03d}" for idx in range(1, 11)])
    with pytest.raises(DatasetValidationError) as exc_info:
        require_exact_index_match(
            left=left,
            right=right,
            left_name="dataset.comparisons.index",
            right_name="dataset.phospho.index",
            error_type=DatasetValidationError,
        )
    message = str(exc_info.value)
    assert (
        "Only in dataset.comparisons.index: 'L001', 'L002', 'L003', 'L004', 'L005', +5 more"
        in message
    )
    assert (
        "Only in dataset.phospho.index: 'R001', 'R002', 'R003', 'R004', 'R005', +5 more"
        in message
    )


def test_summarise_column_mismatch_uses_column_names() -> None:
    summary = summarise_column_mismatch(
        left=pd.Index(["sample_a", "sample_x"]),
        right=pd.Index(["sample_a", "sample_b"]),
        left_name="dataset.total.columns",
        right_name="dataset.phospho.columns",
    )
    assert "Only in dataset.total.columns: 'sample_x'" in summary
    assert "Only in dataset.phospho.columns: 'sample_b'" in summary
