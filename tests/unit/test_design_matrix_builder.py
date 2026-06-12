from __future__ import annotations

import pandas as pd
import pytest

from phospy.errors import WorkflowValidationError
from phospy.science.design import (
    BatchCovariate,
    CategoricalCovariate,
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
