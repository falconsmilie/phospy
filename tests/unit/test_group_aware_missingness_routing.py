from __future__ import annotations

from dataclasses import dataclass
from typing import cast

import numpy as np
import pandas as pd
import pandas.testing as pdt
import pytest

from phospy.errors.input import PhosPyInputError
from phospy.science.datasets.preprocessing.stages.missing_data import (
    GroupAwareMissingnessRouter,
    GroupAwareRoutingOutcome,
    GroupMissingnessClassification,
    GroupMissingnessRoute,
    GroupRoutingAssumption,
)


@dataclass(frozen=True)
class RoutingCase:
    name: str
    groups: tuple[str, ...]
    values: tuple[float, ...]
    expected: tuple[GroupMissingnessClassification, ...]
    retained: bool
    knn_cells: int
    minprob_cells: int


NAN = float("nan")


CASES = (
    RoutingCase(
        "case_a_complete_groups",
        ("A", "A", "A", "B", "B", "B"),
        (10.0, 10.2, 10.1, 11.0, 11.1, 10.9),
        (
            GroupMissingnessClassification.COMPLETE,
            GroupMissingnessClassification.COMPLETE,
        ),
        True,
        0,
        0,
    ),
    RoutingCase(
        "case_b_supported_partial_observation",
        ("A", "A", "A", "B", "B", "B"),
        (10.0, NAN, 10.2, 11.0, 11.1, 10.9),
        (
            GroupMissingnessClassification.SUPPORTED_PARTIAL,
            GroupMissingnessClassification.COMPLETE,
        ),
        True,
        1,
        0,
    ),
    RoutingCase(
        "case_c_unsupported_partial_observation",
        ("A", "A", "A", "B", "B", "B"),
        (10.0, NAN, NAN, 11.0, 11.1, 10.9),
        (
            GroupMissingnessClassification.UNSUPPORTED_PARTIAL,
            GroupMissingnessClassification.COMPLETE,
        ),
        False,
        0,
        0,
    ),
    RoutingCase(
        "case_d_asymmetric_complete_absence",
        ("A", "A", "A", "B", "B", "B"),
        (NAN, NAN, NAN, 11.0, 11.1, 10.9),
        (
            GroupMissingnessClassification.SUPPORTED_FULLY_MISSING,
            GroupMissingnessClassification.COMPLETE,
        ),
        True,
        0,
        3,
    ),
    RoutingCase(
        "case_e_unsupported_absence",
        ("A", "A", "A", "B", "B", "B"),
        (NAN, NAN, NAN, 11.0, NAN, NAN),
        (
            GroupMissingnessClassification.UNSUPPORTED_FULLY_MISSING,
            GroupMissingnessClassification.UNSUPPORTED_PARTIAL,
        ),
        False,
        0,
        0,
    ),
    RoutingCase(
        "fully_missing_every_group",
        ("A", "A", "B", "B"),
        (NAN, NAN, NAN, NAN),
        (
            GroupMissingnessClassification.UNSUPPORTED_FULLY_MISSING,
            GroupMissingnessClassification.UNSUPPORTED_FULLY_MISSING,
        ),
        False,
        0,
        0,
    ),
    RoutingCase(
        "multiple_fully_missing_one_strong_reference",
        ("A", "A", "B", "B", "C", "C"),
        (NAN, NAN, NAN, NAN, 1.0, 2.0),
        (
            GroupMissingnessClassification.SUPPORTED_FULLY_MISSING,
            GroupMissingnessClassification.SUPPORTED_FULLY_MISSING,
            GroupMissingnessClassification.COMPLETE,
        ),
        True,
        0,
        4,
    ),
    RoutingCase(
        "multiple_partially_observed_groups",
        ("A", "A", "B", "B"),
        (1.0, NAN, NAN, 2.0),
        (
            GroupMissingnessClassification.SUPPORTED_PARTIAL,
            GroupMissingnessClassification.SUPPORTED_PARTIAL,
        ),
        True,
        2,
        0,
    ),
    RoutingCase(
        "case_f_mixed_knn_and_minprob",
        ("A", "A", "A", "B", "B", "B", "C", "C", "C"),
        (10.0, NAN, 10.2, NAN, NAN, NAN, 11.0, 11.2, 11.1),
        (
            GroupMissingnessClassification.SUPPORTED_PARTIAL,
            GroupMissingnessClassification.SUPPORTED_FULLY_MISSING,
            GroupMissingnessClassification.COMPLETE,
        ),
        True,
        1,
        3,
    ),
    RoutingCase(
        "unequal_group_sizes",
        ("A", "A", "A", "A", "B", "B"),
        (1.0, 2.0, NAN, NAN, 3.0, 4.0),
        (
            GroupMissingnessClassification.SUPPORTED_PARTIAL,
            GroupMissingnessClassification.COMPLETE,
        ),
        True,
        2,
        0,
    ),
    RoutingCase(
        "threshold_boundary_equality",
        ("A", "A", "B", "B", "C", "C", "C", "C"),
        (1.0, NAN, NAN, NAN, 2.0, 3.0, 4.0, NAN),
        (
            GroupMissingnessClassification.SUPPORTED_PARTIAL,
            GroupMissingnessClassification.SUPPORTED_FULLY_MISSING,
            GroupMissingnessClassification.SUPPORTED_PARTIAL,
        ),
        True,
        2,
        2,
    ),
)


def _route(
    values: tuple[float, ...],
    groups: tuple[str, ...],
    *,
    columns: tuple[str, ...] | None = None,
    min_partial_observed_fraction: float = 0.5,
    min_reference_observed_fraction: float = 0.75,
) -> tuple[pd.DataFrame, GroupAwareRoutingOutcome]:
    sample_ids = columns or tuple(
        f"sample_{position}" for position in range(len(values))
    )
    phospho = pd.DataFrame(
        [values],
        index=pd.Index(["site_1"], name="site_key"),
        columns=pd.Index(sample_ids, name="sample"),
    )
    metadata = pd.DataFrame(
        {"condition": groups, "irrelevant_group": tuple(reversed(groups))},
        index=phospho.columns.copy(),
    )
    outcome = GroupAwareMissingnessRouter().run(
        phospho=phospho,
        sample_metadata=metadata,
        group_column="condition",
        min_partial_observed_fraction=min_partial_observed_fraction,
        min_reference_observed_fraction=min_reference_observed_fraction,
    )
    return phospho, outcome


@pytest.mark.parametrize("case", CASES, ids=lambda case: case.name)
def test_group_aware_classification_table(case: RoutingCase) -> None:
    phospho, outcome = _route(case.values, case.groups)

    assert tuple(fact.classification for fact in outcome.group_facts) == case.expected
    assert bool(outcome.retained_row_mask.iat[0]) is case.retained
    assert outcome.knn_target_cell_count == case.knn_cells
    assert outcome.minprob_target_cell_count == case.minprob_cells

    original_missing = pd.DataFrame(
        ~np.isfinite(phospho.to_numpy(dtype=float)),
        index=phospho.index,
        columns=phospho.columns,
    )
    assigned = outcome.knn_target_mask | outcome.minprob_target_mask
    assert not bool(
        (outcome.knn_target_mask & outcome.minprob_target_mask).to_numpy().any()
    )
    assert not bool((assigned & ~original_missing).to_numpy().any())
    if case.retained:
        pdt.assert_frame_equal(assigned, original_missing)
    else:
        assert not bool(assigned.to_numpy().any())
        assert outcome.dropped_row_ids == ("site_1",)
        assert outcome.dropped_row_reasons[0].reasons_by_group


def test_group_facts_record_counts_fractions_and_routes() -> None:
    _, outcome = _route(
        (1.0, NAN, NAN, NAN, 2.0, 3.0),
        ("A", "A", "B", "B", "C", "C"),
    )
    facts = {fact.group_label: fact for fact in outcome.group_facts}

    assert (
        facts["A"].group_sample_count,
        facts["A"].observed_finite_count,
        facts["A"].missing_count,
        facts["A"].observed_fraction,
        facts["A"].route,
    ) == (2, 1, 1, 0.5, GroupMissingnessRoute.KNN)
    assert facts["B"].route is GroupMissingnessRoute.MINPROB
    assert facts["C"].route is GroupMissingnessRoute.NONE
    assert tuple(outcome.knn_target_mask.loc["site_1"]) == (
        False,
        True,
        False,
        False,
        False,
        False,
    )
    assert tuple(outcome.minprob_target_mask.loc["site_1"]) == (
        False,
        False,
        True,
        True,
        False,
        False,
    )
    assert outcome.min_partial_observed_fraction == 0.5
    assert outcome.min_reference_observed_fraction == 0.75
    assert outcome.original_missingness_mask_hash


def test_dropped_row_clears_supported_block_routes_and_target_assignments() -> None:
    _, outcome = _route(
        (NAN, NAN, 1.0, 2.0, NAN, NAN),
        ("A", "A", "B", "B", "B", "B"),
    )
    facts = {fact.group_label: fact for fact in outcome.group_facts}

    assert not bool(outcome.retained_row_mask.loc["site_1"])
    assert (
        facts["A"].classification
        is GroupMissingnessClassification.UNSUPPORTED_FULLY_MISSING
    )
    assert facts["B"].classification is GroupMissingnessClassification.SUPPORTED_PARTIAL
    assert all(fact.route is GroupMissingnessRoute.NONE for fact in facts.values())
    assert all(
        fact.routing_assumption is GroupRoutingAssumption.NONE
        for fact in facts.values()
    )
    assert not bool(outcome.knn_target_mask.loc["site_1"].any())
    assert not bool(outcome.minprob_target_mask.loc["site_1"].any())
    assert outcome.knn_target_cell_count == 0
    assert outcome.minprob_target_cell_count == 0


def test_sample_order_does_not_change_group_classification_or_routes() -> None:
    values = (1.0, NAN, NAN, NAN, 2.0, 3.0)
    groups = ("A", "A", "B", "B", "C", "C")
    columns = tuple(f"sample_{position}" for position in range(6))
    phospho, baseline = _route(values, groups, columns=columns)
    order = (5, 2, 0, 4, 1, 3)
    _, reordered = _route(
        tuple(values[position] for position in order),
        tuple(groups[position] for position in order),
        columns=tuple(columns[position] for position in order),
    )

    baseline_classes = {
        fact.group_label: fact.classification for fact in baseline.group_facts
    }
    reordered_classes = {
        fact.group_label: fact.classification for fact in reordered.group_facts
    }
    assert reordered_classes == baseline_classes
    pdt.assert_frame_equal(
        reordered.knn_target_mask.loc[:, phospho.columns], baseline.knn_target_mask
    )
    pdt.assert_frame_equal(
        reordered.minprob_target_mask.loc[:, phospho.columns],
        baseline.minprob_target_mask,
    )


def test_sample_metadata_row_order_does_not_change_routing() -> None:
    columns = pd.Index(
        ["a1", "a2", "a3", "b1", "b2", "b3", "c1", "c2", "c3"],
        name="sample",
    )
    phospho = pd.DataFrame(
        [[10.0, NAN, 10.2, NAN, NAN, NAN, 11.0, 11.2, 11.1]],
        index=pd.Index(["mixed"], name="site_key"),
        columns=columns,
    )
    metadata = pd.DataFrame(
        {"condition": ["A", "A", "A", "B", "B", "B", "C", "C", "C"]},
        index=columns.copy(),
    )
    router = GroupAwareMissingnessRouter()
    baseline = router.run(
        phospho=phospho,
        sample_metadata=metadata,
        group_column="condition",
        min_partial_observed_fraction=0.5,
        min_reference_observed_fraction=0.75,
    )
    reordered = router.run(
        phospho=phospho,
        sample_metadata=metadata.loc[
            ["c3", "a2", "b1", "c1", "a1", "b3", "c2", "b2", "a3"]
        ],
        group_column="condition",
        min_partial_observed_fraction=0.5,
        min_reference_observed_fraction=0.75,
    )

    assert reordered.resolved_groups == baseline.resolved_groups
    assert reordered.group_facts == baseline.group_facts
    assert (
        reordered.original_missingness_mask_hash
        == baseline.original_missingness_mask_hash
    )
    pdt.assert_frame_equal(reordered.knn_target_mask, baseline.knn_target_mask)
    pdt.assert_frame_equal(reordered.minprob_target_mask, baseline.minprob_target_mask)


def test_group_labels_come_only_from_sample_metadata() -> None:
    misleading_names = ("A_rep1", "B_rep1", "A_rep2", "B_rep2")
    _, outcome = _route(
        (1.0, NAN, 2.0, NAN),
        ("left", "left", "right", "right"),
        columns=misleading_names,
    )

    assert outcome.resolved_groups.group_labels == ("left", "right")
    assert {fact.group_label for fact in outcome.group_facts} == {"left", "right"}


def test_row_without_missing_values_stays_untargeted_when_other_rows_are_routed() -> (
    None
):
    phospho = pd.DataFrame(
        {
            "a1": [1.0, 1.0],
            "a2": [2.0, NAN],
            "b1": [3.0, 3.0],
            "b2": [4.0, 4.0],
        },
        index=pd.Index(["complete", "partial"], name="site_key"),
    )
    original = phospho.copy(deep=True)
    outcome = GroupAwareMissingnessRouter().run(
        phospho=phospho,
        sample_metadata=pd.DataFrame(
            {"condition": ("A", "A", "B", "B")}, index=phospho.columns
        ),
        group_column="condition",
        min_partial_observed_fraction=0.5,
        min_reference_observed_fraction=0.75,
    )

    pdt.assert_frame_equal(phospho, original)
    assert bool(outcome.retained_row_mask.loc["complete"])
    assert not bool(outcome.knn_target_mask.loc["complete"].any())
    assert not bool(outcome.minprob_target_mask.loc["complete"].any())


def test_unsupported_rows_are_removed_before_numerical_inputs_are_formed() -> None:
    phospho = pd.DataFrame(
        {
            "a1": [1.0, 1.0],
            "a2": [NAN, NAN],
            "a3": [NAN, 2.0],
            "b1": [3.0, 3.0],
            "b2": [4.0, 4.0],
        },
        index=pd.Index(["drop", "retain"], name="site_key"),
    )
    metadata = pd.DataFrame(
        {"condition": ("A", "A", "A", "B", "B")},
        index=phospho.columns,
    )
    outcome = GroupAwareMissingnessRouter().run(
        phospho=phospho,
        sample_metadata=metadata,
        group_column="condition",
        min_partial_observed_fraction=0.5,
        min_reference_observed_fraction=0.75,
    )

    numerical_input = outcome.retain_rows(phospho)

    assert numerical_input.index.tolist() == ["retain"]
    assert not bool(outcome.knn_target_mask.loc["drop"].any())
    assert not bool(outcome.minprob_target_mask.loc["drop"].any())
    assert outcome.unsupported_group_count == 1
    assert outcome.unsupported_partial_group_count == 1
    assert outcome.unsupported_fully_missing_group_count == 0


def test_fully_missing_is_described_as_pattern_and_left_censored_route() -> None:
    _, outcome = _route((NAN, NAN, 1.0, 2.0), ("A", "A", "B", "B"))
    fact = outcome.group_facts[0]

    assert fact.classification.value == "supported_fully_missing"
    assert fact.route.value == "minprob"
    assert fact.routing_assumption is GroupRoutingAssumption.ASYMMETRIC_LEFT_CENSORED
    assert fact.qualifying_reference_group_labels == ("B",)
    assert "mar" not in fact.classification.value
    assert "mnar" not in fact.classification.value


def test_infinite_values_are_not_reinterpreted_as_missing_cells() -> None:
    with pytest.raises(PhosPyInputError, match="finite or missing"):
        _route((1.0, float("inf"), 2.0, 3.0), ("A", "A", "B", "B"))


@pytest.mark.parametrize(
    "threshold_name",
    ("min_partial_observed_fraction", "min_reference_observed_fraction"),
)
@pytest.mark.parametrize(
    "invalid_value",
    (
        True,
        False,
        "0.5",
        0.0,
        -0.1,
        1.1,
        float("nan"),
        float("inf"),
        float("-inf"),
    ),
    ids=(
        "true",
        "false",
        "numeric-string",
        "zero",
        "negative",
        "above-one",
        "nan",
        "infinity",
        "negative-infinity",
    ),
)
def test_router_rejects_invalid_threshold_types_and_ranges(
    threshold_name: str,
    invalid_value: object,
) -> None:
    partial_threshold = 0.5
    reference_threshold = 0.75
    if threshold_name == "min_partial_observed_fraction":
        partial_threshold = cast(float, invalid_value)
    else:
        reference_threshold = cast(float, invalid_value)

    with pytest.raises(PhosPyInputError, match=threshold_name):
        _route(
            (1.0, NAN, 2.0, 3.0),
            ("A", "A", "B", "B"),
            min_partial_observed_fraction=partial_threshold,
            min_reference_observed_fraction=reference_threshold,
        )


@pytest.mark.parametrize("threshold", (0.5, 1))
def test_router_accepts_valid_numeric_thresholds(threshold: float) -> None:
    _, outcome = _route(
        (1.0, NAN, 2.0, 3.0),
        ("A", "A", "B", "B"),
        min_partial_observed_fraction=threshold,
        min_reference_observed_fraction=threshold,
    )

    assert outcome.min_partial_observed_fraction == float(threshold)
    assert outcome.min_reference_observed_fraction == float(threshold)
