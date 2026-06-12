from __future__ import annotations

from types import MappingProxyType

import pandas as pd
import pytest

from phospy.errors import WorkflowValidationError
from phospy.science.design import (
    PAIRED_DESIGN_POLICY_FIXED_BLOCK,
    PAIRED_DESIGN_POLICY_REJECT,
    BatchCovariate,
    CategoricalCovariate,
    ContinuousCovariate,
    DesignMatrixBuilder,
    ExperimentalDesign,
    SampleDesignRecord,
)


def _expected_frame(
    data: dict[str, list[float]],
    *,
    samples: tuple[str, ...],
) -> pd.DataFrame:
    frame = pd.DataFrame(data, index=pd.Index(samples, name="sample"), dtype=float)
    frame.columns = pd.Index(tuple(data), name="coefficient")
    return frame


def _condition_only_design() -> ExperimentalDesign:
    return ExperimentalDesign(
        samples=(
            SampleDesignRecord(sample_id="A_1", condition="A"),
            SampleDesignRecord(sample_id="A_2", condition="A"),
            SampleDesignRecord(sample_id="B_1", condition="B"),
            SampleDesignRecord(sample_id="B_2", condition="B"),
        )
    )


def _continuous_design_from_raw_values(
    values: tuple[object, ...],
) -> ExperimentalDesign:
    samples = tuple(
        SampleDesignRecord(
            sample_id=sample_id,
            condition=condition,
            covariates={"dose": 0.0},
        )
        for sample_id, condition in (
            ("A_1", "A"),
            ("A_2", "A"),
            ("B_1", "B"),
            ("B_2", "B"),
        )
    )
    for record, value in zip(samples, values, strict=True):
        object.__setattr__(
            record,
            "covariates",
            MappingProxyType({"dose": value}),
        )
    return ExperimentalDesign(
        samples=samples,
        fixed_effects=(ContinuousCovariate("dose"),),
    )


def test_condition_only_design_matrix_matches_existing_expected_shape() -> None:
    result = DesignMatrixBuilder().run(design=_condition_only_design())

    expected = _expected_frame(
        {
            "A": [1.0, 1.0, 0.0, 0.0],
            "B": [0.0, 0.0, 1.0, 1.0],
        },
        samples=("A_1", "A_2", "B_1", "B_2"),
    )
    pd.testing.assert_frame_equal(result.frame, expected)
    assert result.sample_labels == ("A_1", "A_2", "B_1", "B_2")
    assert result.coefficient_labels == ("A", "B")


def test_categorical_covariate_creates_stable_columns() -> None:
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(
                sample_id="A_1",
                condition="A",
                covariates={"sex": "F"},
            ),
            SampleDesignRecord(
                sample_id="A_2",
                condition="A",
                covariates={"sex": "M"},
            ),
            SampleDesignRecord(
                sample_id="B_1",
                condition="B",
                covariates={"sex": "F"},
            ),
            SampleDesignRecord(
                sample_id="B_2",
                condition="B",
                covariates={"sex": "M"},
            ),
        ),
        fixed_effects=(CategoricalCovariate("sex"),),
    )

    result = DesignMatrixBuilder().run(design=design)

    expected = _expected_frame(
        {
            "A": [1.0, 1.0, 0.0, 0.0],
            "B": [0.0, 0.0, 1.0, 1.0],
            "sex[M]": [0.0, 1.0, 0.0, 1.0],
        },
        samples=("A_1", "A_2", "B_1", "B_2"),
    )
    pd.testing.assert_frame_equal(result.frame, expected)
    assert dict(result.categorical_levels) == {"sex": ("F", "M")}
    assert dict(result.reference_levels) == {"sex": "F"}
    assert dict(result.unused_levels) == {"sex": ()}


def test_continuous_covariate_creates_raw_numeric_column() -> None:
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(
                sample_id="A_1",
                condition="A",
                covariates={"dose": 0.0},
            ),
            SampleDesignRecord(
                sample_id="A_2",
                condition="A",
                covariates={"dose": 2.5},
            ),
            SampleDesignRecord(
                sample_id="B_1",
                condition="B",
                covariates={"dose": -1.0},
            ),
            SampleDesignRecord(
                sample_id="B_2",
                condition="B",
                covariates={"dose": 10.0},
            ),
        ),
        fixed_effects=(ContinuousCovariate("dose"),),
    )

    result = DesignMatrixBuilder().run(design=design)

    expected = _expected_frame(
        {
            "A": [1.0, 1.0, 0.0, 0.0],
            "B": [0.0, 0.0, 1.0, 1.0],
            "dose": [0.0, 2.5, -1.0, 10.0],
        },
        samples=("A_1", "A_2", "B_1", "B_2"),
    )
    pd.testing.assert_frame_equal(result.frame, expected)
    assert result.encoded_covariates == ("dose",)
    assert dict(result.continuous_columns) == {"dose": "dose"}
    assert dict(result.categorical_levels) == {}


def test_batch_covariate_creates_stable_columns() -> None:
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(sample_id="A_1", condition="A", batch="batch_1"),
            SampleDesignRecord(sample_id="A_2", condition="A", batch="batch_2"),
            SampleDesignRecord(sample_id="B_1", condition="B", batch="batch_1"),
            SampleDesignRecord(sample_id="B_2", condition="B", batch="batch_2"),
        ),
        fixed_effects=(BatchCovariate(),),
    )

    result = DesignMatrixBuilder().run(design=design)

    expected = _expected_frame(
        {
            "A": [1.0, 1.0, 0.0, 0.0],
            "B": [0.0, 0.0, 1.0, 1.0],
            "batch[batch_2]": [0.0, 1.0, 0.0, 1.0],
        },
        samples=("A_1", "A_2", "B_1", "B_2"),
    )
    pd.testing.assert_frame_equal(result.frame, expected)
    assert dict(result.categorical_levels) == {"batch": ("batch_1", "batch_2")}
    assert dict(result.reference_levels) == {"batch": "batch_1"}


def test_fixed_block_two_condition_design_matrix_uses_reference_block() -> None:
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(sample_id="pair_1_A", condition="A", block_id="pair_1"),
            SampleDesignRecord(sample_id="pair_1_B", condition="B", block_id="pair_1"),
            SampleDesignRecord(sample_id="pair_2_A", condition="A", block_id="pair_2"),
            SampleDesignRecord(sample_id="pair_2_B", condition="B", block_id="pair_2"),
        )
    )

    result = DesignMatrixBuilder().run(
        design=design,
        paired_design_policy=PAIRED_DESIGN_POLICY_FIXED_BLOCK,
    )

    expected = _expected_frame(
        {
            "A": [1.0, 0.0, 1.0, 0.0],
            "B": [0.0, 1.0, 0.0, 1.0],
            "block[pair_2]": [0.0, 0.0, 1.0, 1.0],
        },
        samples=("pair_1_A", "pair_1_B", "pair_2_A", "pair_2_B"),
    )
    pd.testing.assert_frame_equal(result.frame, expected)
    assert result.frame.shape == (4, 3)
    assert result.coefficient_labels == ("A", "B", "block[pair_2]")
    assert result.block_levels == ("pair_1", "pair_2")
    assert result.block_reference_level == "pair_1"
    assert dict(result.block_columns) == {"pair_2": "block[pair_2]"}
    assert result.formula == "~0 + condition + block"


def test_fixed_block_design_matrix_supports_multiple_blocks() -> None:
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(sample_id="pair_1_A", condition="A", block_id="pair_1"),
            SampleDesignRecord(sample_id="pair_1_B", condition="B", block_id="pair_1"),
            SampleDesignRecord(sample_id="pair_2_A", condition="A", block_id="pair_2"),
            SampleDesignRecord(sample_id="pair_2_B", condition="B", block_id="pair_2"),
            SampleDesignRecord(sample_id="pair_3_A", condition="A", block_id="pair_3"),
            SampleDesignRecord(sample_id="pair_3_B", condition="B", block_id="pair_3"),
        )
    )

    result = DesignMatrixBuilder().run(
        design=design,
        paired_design_policy=PAIRED_DESIGN_POLICY_FIXED_BLOCK,
    )

    assert result.frame.shape == (6, 4)
    assert result.frame.columns.tolist() == [
        "A",
        "B",
        "block[pair_2]",
        "block[pair_3]",
    ]
    assert result.block_levels == ("pair_1", "pair_2", "pair_3")
    assert result.block_reference_level == "pair_1"
    assert dict(result.block_columns) == {
        "pair_2": "block[pair_2]",
        "pair_3": "block[pair_3]",
    }


def test_fixed_block_column_names_are_deterministic_from_block_levels() -> None:
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(sample_id="s_c_A", condition="A", block_id="pair_c"),
            SampleDesignRecord(sample_id="s_a_A", condition="A", block_id="pair_a"),
            SampleDesignRecord(sample_id="s_b_B", condition="B", block_id="pair_b"),
            SampleDesignRecord(sample_id="s_c_B", condition="B", block_id="pair_c"),
        )
    )

    result = DesignMatrixBuilder().run(
        design=design,
        paired_design_policy=PAIRED_DESIGN_POLICY_FIXED_BLOCK,
    )

    assert result.frame.columns.tolist() == [
        "A",
        "B",
        "block[pair_b]",
        "block[pair_c]",
    ]
    assert result.block_levels == ("pair_a", "pair_b", "pair_c")
    assert result.block_reference_level == "pair_a"


def test_fixed_block_design_matrix_preserves_sample_order() -> None:
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(sample_id="B_pair_2", condition="B", block_id="pair_2"),
            SampleDesignRecord(sample_id="A_pair_1", condition="A", block_id="pair_1"),
            SampleDesignRecord(sample_id="B_pair_1", condition="B", block_id="pair_1"),
            SampleDesignRecord(sample_id="A_pair_2", condition="A", block_id="pair_2"),
        )
    )

    result = DesignMatrixBuilder().run(
        design=design,
        paired_design_policy=PAIRED_DESIGN_POLICY_FIXED_BLOCK,
    )

    assert result.frame.index.tolist() == [
        "B_pair_2",
        "A_pair_1",
        "B_pair_1",
        "A_pair_2",
    ]
    assert result.sample_labels == (
        "B_pair_2",
        "A_pair_1",
        "B_pair_1",
        "A_pair_2",
    )
    assert result.frame.loc["B_pair_2", "block[pair_2]"] == 1.0
    assert result.frame.loc["A_pair_1", "block[pair_2]"] == 0.0


def test_categorical_covariate_and_fixed_block_columns_are_both_encoded() -> None:
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(
                sample_id="pair_1_A",
                condition="A",
                block_id="pair_1",
                covariates={"sex": "F"},
            ),
            SampleDesignRecord(
                sample_id="pair_1_B",
                condition="B",
                block_id="pair_1",
                covariates={"sex": "M"},
            ),
            SampleDesignRecord(
                sample_id="pair_2_A",
                condition="A",
                block_id="pair_2",
                covariates={"sex": "M"},
            ),
            SampleDesignRecord(
                sample_id="pair_2_B",
                condition="B",
                block_id="pair_2",
                covariates={"sex": "F"},
            ),
            SampleDesignRecord(
                sample_id="pair_3_A",
                condition="A",
                block_id="pair_3",
                covariates={"sex": "F"},
            ),
            SampleDesignRecord(
                sample_id="pair_3_B",
                condition="B",
                block_id="pair_3",
                covariates={"sex": "M"},
            ),
        ),
        fixed_effects=(CategoricalCovariate("sex"),),
    )

    result = DesignMatrixBuilder().run(
        design=design,
        paired_design_policy=PAIRED_DESIGN_POLICY_FIXED_BLOCK,
    )

    expected = _expected_frame(
        {
            "A": [1.0, 0.0, 1.0, 0.0, 1.0, 0.0],
            "B": [0.0, 1.0, 0.0, 1.0, 0.0, 1.0],
            "sex[M]": [0.0, 1.0, 1.0, 0.0, 0.0, 1.0],
            "block[pair_2]": [0.0, 0.0, 1.0, 1.0, 0.0, 0.0],
            "block[pair_3]": [0.0, 0.0, 0.0, 0.0, 1.0, 1.0],
        },
        samples=(
            "pair_1_A",
            "pair_1_B",
            "pair_2_A",
            "pair_2_B",
            "pair_3_A",
            "pair_3_B",
        ),
    )
    pd.testing.assert_frame_equal(result.frame, expected)
    assert result.encoded_covariates == ("sex",)
    assert dict(result.covariate_columns) == {"sex": ("sex[M]",)}
    assert dict(result.block_columns) == {
        "pair_2": "block[pair_2]",
        "pair_3": "block[pair_3]",
    }
    assert result.formula == "~0 + condition + sex + block"


def test_block_values_do_not_construct_matrix_when_policy_is_reject() -> None:
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(sample_id="pair_1_A", condition="A", block_id="pair_1"),
            SampleDesignRecord(sample_id="pair_1_B", condition="B", block_id="pair_1"),
        )
    )

    with pytest.raises(
        WorkflowValidationError,
        match="paired_design_policy='reject'",
    ):
        DesignMatrixBuilder().run(
            design=design,
            paired_design_policy=PAIRED_DESIGN_POLICY_REJECT,
        )


def test_fixed_block_policy_rejects_missing_block_ids_without_dropping_samples() -> (
    None
):
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(sample_id="pair_1_A", condition="A", block_id="pair_1"),
            SampleDesignRecord(sample_id="pair_1_B", condition="B"),
            SampleDesignRecord(sample_id="pair_2_A", condition="A", block_id="pair_2"),
        )
    )

    with pytest.raises(
        WorkflowValidationError,
        match="fixed_block.*missing block_id.*pair_1_B",
    ):
        DesignMatrixBuilder().run(
            design=design,
            paired_design_policy=PAIRED_DESIGN_POLICY_FIXED_BLOCK,
        )


def test_sample_ordering_follows_design_sample_order() -> None:
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(sample_id="B_2", condition="B"),
            SampleDesignRecord(sample_id="A_1", condition="A"),
            SampleDesignRecord(sample_id="B_1", condition="B"),
            SampleDesignRecord(sample_id="A_2", condition="A"),
        )
    )

    result = DesignMatrixBuilder().run(design=design)

    assert result.frame.index.tolist() == ["B_2", "A_1", "B_1", "A_2"]
    assert result.frame.columns.tolist() == ["B", "A"]
    assert result.frame.loc["B_2", "B"] == 1.0
    assert result.frame.loc["A_1", "A"] == 1.0


def test_continuous_covariate_numeric_string_is_rejected_without_parsing() -> None:
    design = _continuous_design_from_raw_values((0.0, "2.5", 5.0, 10.0))

    with pytest.raises(
        WorkflowValidationError,
        match="continuous fixed-effect covariate 'dose' must be numeric.*A_2",
    ):
        DesignMatrixBuilder().run(design=design)


def test_continuous_covariate_non_numeric_value_is_rejected() -> None:
    design = _continuous_design_from_raw_values((0.0, "unknown", 5.0, 10.0))

    with pytest.raises(
        WorkflowValidationError,
        match="continuous fixed-effect covariate 'dose' must be numeric.*A_2",
    ):
        DesignMatrixBuilder().run(design=design)


def test_continuous_covariate_missing_value_is_rejected() -> None:
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(
                sample_id="A_1",
                condition="A",
                covariates={"dose": 0.0},
            ),
            SampleDesignRecord(sample_id="A_2", condition="A"),
            SampleDesignRecord(
                sample_id="B_1",
                condition="B",
                covariates={"dose": 5.0},
            ),
            SampleDesignRecord(
                sample_id="B_2",
                condition="B",
                covariates={"dose": 10.0},
            ),
        ),
        fixed_effects=(ContinuousCovariate("dose"),),
    )

    with pytest.raises(
        WorkflowValidationError,
        match="continuous fixed-effect covariate 'dose' has missing values.*A_2",
    ):
        DesignMatrixBuilder().run(design=design)


def test_continuous_covariate_nan_is_rejected() -> None:
    design = _continuous_design_from_raw_values((0.0, float("nan"), 5.0, 10.0))

    with pytest.raises(
        WorkflowValidationError,
        match="continuous fixed-effect covariate 'dose' must be finite.*A_2",
    ):
        DesignMatrixBuilder().run(design=design)


def test_continuous_covariate_infinite_value_is_rejected() -> None:
    design = _continuous_design_from_raw_values((0.0, float("inf"), 5.0, 10.0))

    with pytest.raises(
        WorkflowValidationError,
        match="continuous fixed-effect covariate 'dose' must be finite.*A_2",
    ):
        DesignMatrixBuilder().run(design=design)


def test_mixed_covariates_have_deterministic_column_order() -> None:
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(
                sample_id="A_1",
                condition="A",
                covariates={"dose": 0.0, "sex": "F"},
            ),
            SampleDesignRecord(
                sample_id="A_2",
                condition="A",
                covariates={"dose": 2.5, "sex": "M"},
            ),
            SampleDesignRecord(
                sample_id="B_1",
                condition="B",
                covariates={"dose": 5.0, "sex": "F"},
            ),
            SampleDesignRecord(
                sample_id="B_2",
                condition="B",
                covariates={"dose": 10.0, "sex": "M"},
            ),
        ),
        fixed_effects=(
            ContinuousCovariate("dose"),
            CategoricalCovariate("sex"),
        ),
    )

    result = DesignMatrixBuilder().run(design=design)

    assert result.frame.columns.tolist() == ["A", "B", "dose", "sex[M]"]
    assert result.encoded_covariates == ("dose", "sex")


def test_unused_explicit_categories_are_reported_and_encoded() -> None:
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(
                sample_id="A_1",
                condition="A",
                covariates={"sex": "F"},
            ),
            SampleDesignRecord(
                sample_id="A_2",
                condition="A",
                covariates={"sex": "M"},
            ),
            SampleDesignRecord(
                sample_id="B_1",
                condition="B",
                covariates={"sex": "F"},
            ),
            SampleDesignRecord(
                sample_id="B_2",
                condition="B",
                covariates={"sex": "M"},
            ),
        ),
        fixed_effects=(CategoricalCovariate("sex"),),
    )

    result = DesignMatrixBuilder().run(
        design=design,
        categorical_levels={"sex": ("F", "M", "not_observed")},
    )

    assert result.frame.columns.tolist() == ["A", "B", "sex[M]", "sex[not_observed]"]
    assert result.frame.loc[:, "sex[not_observed]"].tolist() == [
        0.0,
        0.0,
        0.0,
        0.0,
    ]
    assert dict(result.unused_levels) == {"sex": ("not_observed",)}


def test_missing_covariate_value_fails_explicitly() -> None:
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(
                sample_id="A_1",
                condition="A",
                covariates={"sex": "F"},
            ),
            SampleDesignRecord(sample_id="A_2", condition="A"),
            SampleDesignRecord(
                sample_id="B_1",
                condition="B",
                covariates={"sex": "M"},
            ),
            SampleDesignRecord(
                sample_id="B_2",
                condition="B",
                covariates={"sex": "F"},
            ),
        ),
        fixed_effects=(CategoricalCovariate("sex", required=False),),
    )

    with pytest.raises(WorkflowValidationError, match="has missing values"):
        DesignMatrixBuilder().run(design=design)
