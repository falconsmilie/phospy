from __future__ import annotations

import pytest

from phospy.api.configs import (
    DatasetGroupCoverageFilterConfig,
    DatasetPreprocessingConfig,
)
from phospy.errors.input import PhosPyInputError
from phospy.science.datasets.preprocessing.models import PreprocessingPlan


def test_group_coverage_filter_accepts_count_threshold_config() -> None:
    config = DatasetGroupCoverageFilterConfig(
        enabled=True,
        group_column="condition",
        min_finite_observations_per_group=2,
        min_groups_passing_threshold=1,
    )

    assert config.enabled is True
    assert config.group_column == "condition"
    assert config.min_finite_observations_per_group == 2
    assert config.min_finite_fraction_per_group is None
    assert config.min_groups_passing_threshold == 1


def test_group_coverage_filter_accepts_fraction_threshold_config() -> None:
    config = DatasetGroupCoverageFilterConfig(
        enabled=True,
        group_column="condition",
        min_finite_fraction_per_group=0.75,
        min_groups_passing_threshold=2,
    )

    assert config.min_finite_observations_per_group is None
    assert config.min_finite_fraction_per_group == pytest.approx(0.75)
    assert config.min_groups_passing_threshold == 2


@pytest.mark.parametrize(
    ("kwargs", "pattern"),
    [
        pytest.param(
            {"min_finite_observations_per_group": 0},
            "min_finite_observations_per_group must be greater than or equal to 1",
            id="count-below-one",
        ),
        pytest.param(
            {"min_finite_observations_per_group": True},
            "min_finite_observations_per_group must be an int",
            id="count-bool",
        ),
        pytest.param(
            {"min_finite_fraction_per_group": 0.0},
            "min_finite_fraction_per_group must satisfy 0 < value <= 1",
            id="fraction-zero",
        ),
        pytest.param(
            {"min_finite_fraction_per_group": 1.2},
            "min_finite_fraction_per_group must satisfy 0 < value <= 1",
            id="fraction-above-one",
        ),
        pytest.param(
            {"min_finite_fraction_per_group": float("nan")},
            "min_finite_fraction_per_group must satisfy 0 < value <= 1",
            id="fraction-nan",
        ),
        pytest.param(
            {"min_finite_fraction_per_group": True},
            "min_finite_fraction_per_group must be a float",
            id="fraction-bool",
        ),
    ],
)
def test_group_coverage_filter_rejects_invalid_threshold_values(
    kwargs: dict[str, object],
    pattern: str,
) -> None:
    with pytest.raises(PhosPyInputError, match=pattern):
        DatasetGroupCoverageFilterConfig(
            enabled=True,
            group_column="condition",
            min_groups_passing_threshold=1,
            **kwargs,
        )


@pytest.mark.parametrize(
    "group_column",
    [None, "", "   "],
)
def test_group_coverage_filter_rejects_missing_group_column_when_enabled(
    group_column: object,
) -> None:
    with pytest.raises(PhosPyInputError, match="group_column must be a non-empty"):
        DatasetGroupCoverageFilterConfig(
            enabled=True,
            group_column=group_column,  # type: ignore[arg-type]
            min_finite_observations_per_group=2,
            min_groups_passing_threshold=1,
        )


@pytest.mark.parametrize(
    ("value", "pattern"),
    [
        pytest.param(
            0,
            "min_groups_passing_threshold must be greater than or equal to 1",
            id="zero",
        ),
        pytest.param(
            -1,
            "min_groups_passing_threshold must be greater than or equal to 1",
            id="negative",
        ),
        pytest.param(True, "min_groups_passing_threshold must be an int", id="bool"),
    ],
)
def test_group_coverage_filter_rejects_invalid_minimum_number_of_groups(
    value: object,
    pattern: str,
) -> None:
    with pytest.raises(PhosPyInputError, match=pattern):
        DatasetGroupCoverageFilterConfig(
            enabled=True,
            group_column="condition",
            min_finite_observations_per_group=2,
            min_groups_passing_threshold=value,  # type: ignore[arg-type]
        )


def test_group_coverage_filter_default_is_disabled_and_unscheduled() -> None:
    config = DatasetPreprocessingConfig()
    plan = PreprocessingPlan.from_config(config)

    assert config.group_coverage_filter == DatasetGroupCoverageFilterConfig()
    assert config.group_coverage_filter.enabled is False
    assert config.group_coverage_filter.group_column is None
    assert config.group_coverage_filter.min_finite_observations_per_group is None
    assert config.group_coverage_filter.min_finite_fraction_per_group is None
    assert plan.stage_order == ("localisation_confidence", "missing_data")


def test_group_coverage_filter_enabled_config_schedules_pre_missing_stage() -> None:
    config = DatasetPreprocessingConfig(
        group_coverage_filter=DatasetGroupCoverageFilterConfig(
            enabled=True,
            group_column="condition",
            min_finite_observations_per_group=2,
            min_groups_passing_threshold=1,
        )
    )

    assert PreprocessingPlan.from_config(config).stage_order == (
        "localisation_confidence",
        "group_coverage_filter",
        "missing_data",
    )


def test_group_coverage_filter_rejects_count_and_fraction_together() -> None:
    with pytest.raises(PhosPyInputError, match="mutually exclusive"):
        DatasetGroupCoverageFilterConfig(
            enabled=True,
            group_column="condition",
            min_finite_observations_per_group=2,
            min_finite_fraction_per_group=0.5,
            min_groups_passing_threshold=1,
        )


def test_group_coverage_filter_rejects_enabled_config_without_threshold() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="must set exactly one finite-observation threshold when enabled",
    ):
        DatasetGroupCoverageFilterConfig(
            enabled=True,
            group_column="condition",
        )


def test_dataset_preprocessing_config_rejects_wrong_group_coverage_filter_type() -> (
    None
):
    with pytest.raises(
        PhosPyInputError,
        match=(
            "preprocessing_config.group_coverage_filter must be a "
            "DatasetGroupCoverageFilterConfig"
        ),
    ):
        DatasetPreprocessingConfig(group_coverage_filter=object())  # type: ignore[arg-type]
