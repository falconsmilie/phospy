"""Classify and route group-aware missingness before numerical imputation."""

from __future__ import annotations

import math
from dataclasses import replace

import numpy as np
import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.science.datasets.preprocessing.sample_group_metadata import (
    SampleGroupMetadataResolver,
)

from .diagnostics import hash_missingness_mask
from .models import (
    DroppedRowRoutingRecord,
    GroupAwareRoutingOutcome,
    GroupMissingnessClassification,
    GroupMissingnessRoute,
    GroupRoutingAssumption,
    GroupRoutingFact,
)


class GroupAwareMissingnessRouter:
    """Route cells using only the input matrix's observed/missing pattern."""

    def __init__(
        self,
        *,
        sample_group_resolver: SampleGroupMetadataResolver | None = None,
    ) -> None:
        self._sample_group_resolver = (
            sample_group_resolver or SampleGroupMetadataResolver()
        )

    def run(
        self,
        *,
        phospho: pd.DataFrame,
        sample_metadata: pd.DataFrame | None,
        group_column: str | None,
        min_partial_observed_fraction: float,
        min_reference_observed_fraction: float,
    ) -> GroupAwareRoutingOutcome:
        """Classify all row/group blocks once and assign eligible missing cells."""

        partial_threshold = _validate_fraction_threshold(
            min_partial_observed_fraction,
            field_name="min_partial_observed_fraction",
        )
        reference_threshold = _validate_fraction_threshold(
            min_reference_observed_fraction,
            field_name="min_reference_observed_fraction",
        )
        _require_numeric_matrix(phospho)
        resolved_groups = self._sample_group_resolver.run(
            phospho=phospho,
            sample_metadata=sample_metadata,
            group_column=group_column,
        )

        values = phospho.to_numpy(dtype="float64", copy=True, na_value=np.nan)
        finite_values = np.isfinite(values)
        missing_values = np.isnan(values)
        infinite_values = (~finite_values) & (~missing_values)
        if bool(infinite_values.any()):
            raise PhosPyInputError(
                "group-aware missingness routing requires every phospho cell to "
                "be finite or missing; infinite cell count="
                f"{int(infinite_values.sum())}"
            )
        original_missing_mask = pd.DataFrame(
            missing_values,
            index=phospho.index.copy(),
            columns=phospho.columns.copy(),
        )
        knn_target_mask = _false_mask_like(phospho)
        minprob_target_mask = _false_mask_like(phospho)
        retained_values = np.ones(len(phospho.index), dtype=bool)
        facts: list[GroupRoutingFact] = []
        dropped_records: list[DroppedRowRoutingRecord] = []

        groups = resolved_groups.group_labels
        group_positions = tuple(
            (
                group_label,
                tuple(
                    position
                    for position, assigned_group in enumerate(
                        resolved_groups.group_by_column_position
                    )
                    if assigned_group == group_label
                ),
            )
            for group_label in groups
        )

        for row_position, row_label in enumerate(phospho.index):
            observed_counts = tuple(
                int(finite_values[row_position, positions].sum())
                for _, positions in group_positions
            )
            observed_fractions = tuple(
                observed_count / float(len(positions))
                for observed_count, (_, positions) in zip(
                    observed_counts, group_positions, strict=True
                )
            )
            row_facts: list[GroupRoutingFact] = []
            unsupported: dict[str, GroupMissingnessClassification] = {}

            for group_index, (group_label, positions) in enumerate(group_positions):
                observed_count = observed_counts[group_index]
                sample_count = len(positions)
                qualifying_references = tuple(
                    other_group_label
                    for other_index, (other_group_label, _) in enumerate(
                        group_positions
                    )
                    if other_index != group_index
                    and observed_fractions[other_index] >= reference_threshold
                )
                classification, route = _classify_group(
                    observed_count=observed_count,
                    sample_count=sample_count,
                    observed_fraction=observed_fractions[group_index],
                    has_qualifying_reference=bool(qualifying_references),
                    min_partial_observed_fraction=partial_threshold,
                )
                row_facts.append(
                    GroupRoutingFact(
                        row_id=str(row_label),
                        group_label=group_label,
                        group_sample_count=sample_count,
                        observed_finite_count=observed_count,
                        missing_count=sample_count - observed_count,
                        observed_fraction=observed_fractions[group_index],
                        classification=classification,
                        route=route,
                        routing_assumption=_routing_assumption(classification),
                        qualifying_reference_group_labels=(
                            qualifying_references
                            if classification
                            is GroupMissingnessClassification.SUPPORTED_FULLY_MISSING
                            else ()
                        ),
                    )
                )
                if classification in {
                    GroupMissingnessClassification.UNSUPPORTED_PARTIAL,
                    GroupMissingnessClassification.UNSUPPORTED_FULLY_MISSING,
                }:
                    unsupported[group_label] = classification

            if unsupported:
                row_facts = [
                    replace(
                        fact,
                        route=GroupMissingnessRoute.NONE,
                        routing_assumption=GroupRoutingAssumption.NONE,
                    )
                    for fact in row_facts
                ]
            facts.extend(row_facts)
            if unsupported:
                retained_values[row_position] = False
                dropped_records.append(
                    DroppedRowRoutingRecord(
                        row_id=str(row_label),
                        reasons_by_group=unsupported,
                    )
                )
                continue

            for fact, (_, positions) in zip(row_facts, group_positions, strict=True):
                missing_positions = [
                    position
                    for position in positions
                    if not finite_values[row_position, position]
                ]
                if fact.route is GroupMissingnessRoute.KNN:
                    knn_target_mask.iloc[row_position, missing_positions] = True
                elif fact.route is GroupMissingnessRoute.MINPROB:
                    minprob_target_mask.iloc[row_position, missing_positions] = True

        retained_row_mask = pd.Series(
            retained_values,
            index=phospho.index.copy(),
            name="group_aware_routing_retained",
            dtype=bool,
        )
        _validate_routing_invariants(
            original_missing_mask=original_missing_mask,
            retained_row_mask=retained_row_mask,
            knn_target_mask=knn_target_mask,
            minprob_target_mask=minprob_target_mask,
        )
        unsupported_partial_count = sum(
            fact.classification is GroupMissingnessClassification.UNSUPPORTED_PARTIAL
            for fact in facts
        )
        unsupported_fully_missing_count = sum(
            fact.classification
            is GroupMissingnessClassification.UNSUPPORTED_FULLY_MISSING
            for fact in facts
        )
        return GroupAwareRoutingOutcome(
            retained_row_mask=retained_row_mask,
            retained_row_ids=tuple(
                str(row_id) for row_id in phospho.index[retained_values].tolist()
            ),
            dropped_row_ids=tuple(record.row_id for record in dropped_records),
            dropped_row_reasons=tuple(dropped_records),
            resolved_groups=resolved_groups,
            group_facts=tuple(facts),
            knn_target_mask=knn_target_mask,
            minprob_target_mask=minprob_target_mask,
            knn_target_cell_count=int(knn_target_mask.to_numpy().sum()),
            minprob_target_cell_count=int(minprob_target_mask.to_numpy().sum()),
            unsupported_group_count=(
                unsupported_partial_count + unsupported_fully_missing_count
            ),
            unsupported_partial_group_count=unsupported_partial_count,
            unsupported_fully_missing_group_count=unsupported_fully_missing_count,
            min_partial_observed_fraction=partial_threshold,
            min_reference_observed_fraction=reference_threshold,
            original_missingness_mask_hash=hash_missingness_mask(original_missing_mask),
        )


def route_group_aware_missingness(
    *,
    phospho: pd.DataFrame,
    sample_metadata: pd.DataFrame | None,
    group_column: str | None,
    min_partial_observed_fraction: float,
    min_reference_observed_fraction: float,
) -> GroupAwareRoutingOutcome:
    """Convenience entry point for the independently testable router."""

    return GroupAwareMissingnessRouter().run(
        phospho=phospho,
        sample_metadata=sample_metadata,
        group_column=group_column,
        min_partial_observed_fraction=min_partial_observed_fraction,
        min_reference_observed_fraction=min_reference_observed_fraction,
    )


def _classify_group(
    *,
    observed_count: int,
    sample_count: int,
    observed_fraction: float,
    has_qualifying_reference: bool,
    min_partial_observed_fraction: float,
) -> tuple[GroupMissingnessClassification, GroupMissingnessRoute]:
    if observed_count == sample_count:
        return (
            GroupMissingnessClassification.COMPLETE,
            GroupMissingnessRoute.NONE,
        )
    if observed_count > 0:
        if observed_fraction >= min_partial_observed_fraction:
            return (
                GroupMissingnessClassification.SUPPORTED_PARTIAL,
                GroupMissingnessRoute.KNN,
            )
        return (
            GroupMissingnessClassification.UNSUPPORTED_PARTIAL,
            GroupMissingnessRoute.NONE,
        )
    if has_qualifying_reference:
        return (
            GroupMissingnessClassification.SUPPORTED_FULLY_MISSING,
            GroupMissingnessRoute.MINPROB,
        )
    return (
        GroupMissingnessClassification.UNSUPPORTED_FULLY_MISSING,
        GroupMissingnessRoute.NONE,
    )


def _routing_assumption(
    classification: GroupMissingnessClassification,
) -> GroupRoutingAssumption:
    if classification is GroupMissingnessClassification.SUPPORTED_PARTIAL:
        return GroupRoutingAssumption.SIMILARITY_BASED_ELIGIBILITY
    if classification is GroupMissingnessClassification.SUPPORTED_FULLY_MISSING:
        return GroupRoutingAssumption.ASYMMETRIC_LEFT_CENSORED
    return GroupRoutingAssumption.NONE


def _validate_fraction_threshold(value: float, *, field_name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise PhosPyInputError(
            f"group-aware missingness routing {field_name} must be a finite "
            "number satisfying 0 < value <= 1"
        )
    resolved = float(value)
    if not math.isfinite(resolved) or not (0.0 < resolved <= 1.0):
        raise PhosPyInputError(
            f"group-aware missingness routing {field_name} must be a finite "
            "number satisfying 0 < value <= 1"
        )
    return resolved


def _require_numeric_matrix(phospho: pd.DataFrame) -> None:
    invalid_columns = [
        str(column)
        for column in phospho.columns
        if not pd.api.types.is_numeric_dtype(phospho[column])
        or pd.api.types.is_bool_dtype(phospho[column])
    ]
    if invalid_columns:
        preview = ", ".join(repr(column) for column in invalid_columns[:5])
        raise PhosPyInputError(
            "group-aware missingness routing requires numeric, non-boolean "
            f"phospho sample columns; invalid columns: {preview}"
        )


def _false_mask_like(phospho: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame(
        False, index=phospho.index.copy(), columns=phospho.columns.copy()
    )


def _validate_routing_invariants(
    *,
    original_missing_mask: pd.DataFrame,
    retained_row_mask: pd.Series,
    knn_target_mask: pd.DataFrame,
    minprob_target_mask: pd.DataFrame,
) -> None:
    if bool((knn_target_mask & minprob_target_mask).to_numpy().any()):
        raise RuntimeError("group-aware KNN and MinProb target masks overlap")
    assigned = knn_target_mask | minprob_target_mask
    if bool((assigned & ~original_missing_mask).to_numpy().any()):
        raise RuntimeError("group-aware routing targeted an originally observed cell")
    retained_missing = original_missing_mask.loc[retained_row_mask]
    retained_assigned = assigned.loc[retained_row_mask]
    if not retained_assigned.equals(retained_missing):
        raise RuntimeError(
            "group-aware routing left missing cells unresolved on retained rows"
        )
    if bool(assigned.loc[~retained_row_mask].to_numpy().any()):
        raise RuntimeError("group-aware routing targeted cells on dropped rows")


__all__ = [
    "GroupAwareMissingnessRouter",
    "route_group_aware_missingness",
]
