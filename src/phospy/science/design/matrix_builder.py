"""Deterministic experimental-design matrix construction."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from types import MappingProxyType
from typing import cast

import pandas as pd

from phospy.errors.validation import WorkflowValidationError
from phospy.science.design.models import (
    FIXED_EFFECT_COVARIATE_KIND_BATCH,
    FIXED_EFFECT_COVARIATE_KIND_CATEGORICAL,
    FIXED_EFFECT_COVARIATE_KIND_CONTINUOUS,
    ExperimentalDesign,
    FixedEffectCovariate,
    SampleDesignRecord,
)

CategoricalLevelSpec = Mapping[str, Sequence[str | int | float]]


def _empty_level_mapping() -> Mapping[str, tuple[str, ...]]:
    return MappingProxyType({})


def _empty_reference_mapping() -> Mapping[str, str]:
    return MappingProxyType({})


@dataclass(frozen=True, slots=True)
class DesignMatrixBuildResult:
    """Owned design matrix plus encoding metadata."""

    frame: pd.DataFrame
    condition_labels: tuple[str, ...]
    sample_labels: tuple[str, ...]
    coefficient_labels: tuple[str, ...]
    encoded_covariates: tuple[str, ...] = ()
    categorical_levels: Mapping[str, tuple[str, ...]] = field(
        default_factory=_empty_level_mapping
    )
    reference_levels: Mapping[str, str] = field(
        default_factory=_empty_reference_mapping
    )
    unused_levels: Mapping[str, tuple[str, ...]] = field(
        default_factory=_empty_level_mapping
    )

    def __post_init__(self) -> None:
        object.__setattr__(self, "frame", self.frame.copy(deep=True))
        object.__setattr__(self, "condition_labels", tuple(self.condition_labels))
        object.__setattr__(self, "sample_labels", tuple(self.sample_labels))
        object.__setattr__(self, "coefficient_labels", tuple(self.coefficient_labels))
        object.__setattr__(self, "encoded_covariates", tuple(self.encoded_covariates))
        object.__setattr__(
            self,
            "categorical_levels",
            MappingProxyType(
                {
                    str(name): tuple(levels)
                    for name, levels in self.categorical_levels.items()
                }
            ),
        )
        object.__setattr__(
            self,
            "reference_levels",
            MappingProxyType(
                {str(name): str(level) for name, level in self.reference_levels.items()}
            ),
        )
        object.__setattr__(
            self,
            "unused_levels",
            MappingProxyType(
                {
                    str(name): tuple(levels)
                    for name, levels in self.unused_levels.items()
                }
            ),
        )


class DesignMatrixBuilder:
    """Build deterministic condition and categorical fixed-effect matrices."""

    def run(
        self,
        *,
        design: ExperimentalDesign,
        condition_labels: Sequence[str] | None = None,
        categorical_levels: CategoricalLevelSpec | None = None,
    ) -> DesignMatrixBuildResult:
        if not isinstance(cast(object, design), ExperimentalDesign):
            raise WorkflowValidationError(
                "design matrix builder requires an ExperimentalDesign"
            )

        resolved_condition_labels = _resolve_condition_labels(
            design=design,
            condition_labels=condition_labels,
        )
        sample_labels = design.sample_ids()
        records = design.samples
        explicit_levels = _normalise_explicit_levels(categorical_levels)
        data: dict[str, list[float]] = {
            condition: [
                1.0 if record.condition == condition else 0.0 for record in records
            ]
            for condition in resolved_condition_labels
        }

        observed_conditions = {record.condition for record in records}
        unknown_conditions = sorted(
            observed_conditions - set(resolved_condition_labels)
        )
        if unknown_conditions:
            raise WorkflowValidationError(
                "design matrix builder condition labels are missing observed "
                "conditions: " + ", ".join(unknown_conditions)
            )

        categorical_levels_by_name: dict[str, tuple[str, ...]] = {}
        reference_levels: dict[str, str] = {}
        unused_levels: dict[str, tuple[str, ...]] = {}
        encoded_covariates: list[str] = []
        modelled_effect_names: set[str] = set()
        for covariate in design.fixed_effects:
            if not covariate.include_in_model:
                continue
            modelled_effect_names.add(covariate.name)
            if covariate.kind == FIXED_EFFECT_COVARIATE_KIND_CONTINUOUS:
                raise WorkflowValidationError(
                    "design matrix builder does not support continuous fixed effects"
                )
            if covariate.kind not in {
                FIXED_EFFECT_COVARIATE_KIND_BATCH,
                FIXED_EFFECT_COVARIATE_KIND_CATEGORICAL,
            }:
                raise WorkflowValidationError(
                    "design matrix builder does not support fixed-effect covariate "
                    f"kind {covariate.kind!r}"
                )

            values = _collect_categorical_values(
                records=records,
                covariate=covariate,
            )
            observed_levels = _unique_in_order(values)
            levels = explicit_levels.get(covariate.name, observed_levels)
            _validate_observed_levels_are_declared(
                covariate_name=covariate.name,
                observed_levels=observed_levels,
                declared_levels=levels,
            )
            reference_level = levels[0]
            categorical_levels_by_name[covariate.name] = levels
            reference_levels[covariate.name] = reference_level
            unused_levels[covariate.name] = tuple(
                level for level in levels if level not in set(observed_levels)
            )
            encoded_covariates.append(covariate.name)
            for level in levels[1:]:
                column_name = _categorical_column_name(covariate.name, level)
                data[column_name] = [1.0 if value == level else 0.0 for value in values]

        unknown_level_specs = sorted(set(explicit_levels) - modelled_effect_names)
        if unknown_level_specs:
            raise WorkflowValidationError(
                "design matrix builder received categorical level order for "
                "non-modelled covariates: " + ", ".join(unknown_level_specs)
            )

        frame = pd.DataFrame(
            data,
            index=pd.Index(sample_labels, name="sample"),
            dtype=float,
        )
        frame.columns = pd.Index(tuple(data), name="coefficient")
        return DesignMatrixBuildResult(
            frame=frame,
            condition_labels=resolved_condition_labels,
            sample_labels=sample_labels,
            coefficient_labels=tuple(data),
            encoded_covariates=tuple(encoded_covariates),
            categorical_levels=categorical_levels_by_name,
            reference_levels=reference_levels,
            unused_levels=unused_levels,
        )


def _resolve_condition_labels(
    *,
    design: ExperimentalDesign,
    condition_labels: Sequence[str] | None,
) -> tuple[str, ...]:
    if condition_labels is None:
        return design.condition_labels()
    if isinstance(condition_labels, str):
        raise WorkflowValidationError(
            "design matrix builder condition_labels must be a sequence of strings"
        )
    resolved = tuple(
        _normalise_level_label(
            value,
            field_name="design matrix builder condition_labels",
        )
        for value in condition_labels
    )
    if not resolved:
        raise WorkflowValidationError(
            "design matrix builder condition_labels must not be empty"
        )
    duplicates = sorted({label for label in resolved if resolved.count(label) > 1})
    if duplicates:
        raise WorkflowValidationError(
            "design matrix builder condition_labels contain duplicate values: "
            + ", ".join(duplicates)
        )
    return resolved


def _normalise_explicit_levels(
    categorical_levels: CategoricalLevelSpec | None,
) -> dict[str, tuple[str, ...]]:
    if categorical_levels is None:
        return {}
    if not isinstance(categorical_levels, Mapping):
        raise WorkflowValidationError(
            "design matrix builder categorical_levels must be a mapping"
        )
    normalised: dict[str, tuple[str, ...]] = {}
    for raw_name, raw_levels in categorical_levels.items():
        name = _normalise_level_label(
            raw_name,
            field_name="design matrix builder categorical_levels keys",
        )
        if name in normalised:
            raise WorkflowValidationError(
                "design matrix builder categorical_levels contain duplicate "
                f"covariate name: {name}"
            )
        if isinstance(raw_levels, str) or not isinstance(raw_levels, Sequence):
            raise WorkflowValidationError(
                "design matrix builder categorical_levels values must be sequences"
            )
        levels = tuple(
            _normalise_level_label(
                level,
                field_name=(f"design matrix builder categorical_levels[{name!r}]"),
            )
            for level in raw_levels
        )
        if not levels:
            raise WorkflowValidationError(
                "design matrix builder categorical_levels values must not be empty"
            )
        duplicates = sorted({level for level in levels if levels.count(level) > 1})
        if duplicates:
            raise WorkflowValidationError(
                "design matrix builder categorical_levels contain duplicate levels "
                f"for {name!r}: " + ", ".join(duplicates)
            )
        normalised[name] = levels
    return normalised


def _collect_categorical_values(
    *,
    records: tuple[SampleDesignRecord, ...],
    covariate: FixedEffectCovariate,
) -> tuple[str, ...]:
    values: list[str] = []
    missing_samples: list[str] = []
    for record in records:
        if covariate.kind == FIXED_EFFECT_COVARIATE_KIND_BATCH:
            value = record.batch
            missing = value is None
        else:
            missing = covariate.name not in record.covariates
            value = None if missing else record.covariates[covariate.name]
        if missing:
            missing_samples.append(record.sample_id)
            continue
        values.append(
            _normalise_level_label(
                value,
                field_name=(
                    f"design matrix builder fixed-effect covariate {covariate.name!r}"
                ),
            )
        )
    if missing_samples:
        raise WorkflowValidationError(
            "design matrix builder fixed-effect covariate "
            f"{covariate.name!r} has missing values for samples: "
            + ", ".join(missing_samples)
        )
    if not values:
        raise WorkflowValidationError(
            "design matrix builder fixed-effect covariate "
            f"{covariate.name!r} has no observed levels"
        )
    return tuple(values)


def _unique_in_order(values: Sequence[str]) -> tuple[str, ...]:
    seen: set[str] = set()
    ordered: list[str] = []
    for value in values:
        if value in seen:
            continue
        seen.add(value)
        ordered.append(value)
    return tuple(ordered)


def _validate_observed_levels_are_declared(
    *,
    covariate_name: str,
    observed_levels: tuple[str, ...],
    declared_levels: tuple[str, ...],
) -> None:
    missing = [level for level in observed_levels if level not in set(declared_levels)]
    if missing:
        raise WorkflowValidationError(
            "design matrix builder observed levels are missing from explicit level "
            f"order for {covariate_name!r}: " + ", ".join(missing)
        )


def _categorical_column_name(covariate_name: str, level: str) -> str:
    return f"{covariate_name}[{level}]"


def _normalise_level_label(value: object, *, field_name: str) -> str:
    if isinstance(value, bool) or not isinstance(value, str | int | float):
        raise WorkflowValidationError(f"{field_name} values must be strings or numbers")
    label = str(value).strip()
    if label == "":
        raise WorkflowValidationError(f"{field_name} values must be non-empty")
    return label


__all__ = ["DesignMatrixBuildResult", "DesignMatrixBuilder"]
