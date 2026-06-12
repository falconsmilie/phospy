"""Typed experimental-design and contrast contracts."""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Literal, cast

from phospy.errors.validation import WorkflowValidationError

FIXED_EFFECT_COVARIATE_KIND_CATEGORICAL = "categorical"
FIXED_EFFECT_COVARIATE_KIND_CONTINUOUS = "continuous"
FIXED_EFFECT_COVARIATE_KIND_BATCH = "batch"
FixedEffectCovariateKind = Literal["categorical", "continuous", "batch"]
SUPPORTED_FIXED_EFFECT_COVARIATE_KINDS = frozenset(
    {
        FIXED_EFFECT_COVARIATE_KIND_CATEGORICAL,
        FIXED_EFFECT_COVARIATE_KIND_CONTINUOUS,
        FIXED_EFFECT_COVARIATE_KIND_BATCH,
    }
)

_RESERVED_SAMPLE_FIELD_NAMES = {
    "biological_replicate_id",
    "block",
    "condition",
    "sample",
    "sample_id",
    "technical_replicate_id",
}


def _require_non_empty_text(value: str, *, field_name: str) -> str:
    if not isinstance(value, str):
        raise WorkflowValidationError(f"{field_name} must be a string")
    normalized = value.strip()
    if normalized == "":
        raise WorkflowValidationError(f"{field_name} must be a non-empty string")
    return normalized


def _normalize_optional_text(value: str | None, *, field_name: str) -> str | None:
    if value is None:
        return None
    if not isinstance(value, str):
        raise WorkflowValidationError(f"{field_name} must be a string or None")
    normalized = value.strip()
    if normalized == "":
        raise WorkflowValidationError(
            f"{field_name} must be a non-empty string when provided"
        )
    return normalized


def _require_bool(value: bool, *, field_name: str) -> bool:
    if not isinstance(value, bool):
        raise WorkflowValidationError(f"{field_name} must be a bool")
    return value


def _normalize_fixed_effect_kind(
    value: str,
    *,
    field_name: str,
) -> FixedEffectCovariateKind:
    if not isinstance(value, str):
        raise WorkflowValidationError(
            f"{field_name} has unsupported covariate kind: {value!r}"
        )
    normalized = value.strip()
    if normalized not in SUPPORTED_FIXED_EFFECT_COVARIATE_KINDS:
        supported = ", ".join(sorted(SUPPORTED_FIXED_EFFECT_COVARIATE_KINDS))
        raise WorkflowValidationError(
            f"{field_name} has unsupported covariate kind: {value!r}; "
            f"supported kinds: {supported}"
        )
    return cast(FixedEffectCovariateKind, normalized)


def _normalize_covariate_value(
    value: str | int | float,
    *,
    field_name: str,
) -> str | float:
    if isinstance(value, bool) or not isinstance(value, str | int | float):
        raise WorkflowValidationError(
            f"{field_name} values must be strings or finite numeric values"
        )
    if isinstance(value, str):
        normalized = value.strip()
        if normalized == "":
            raise WorkflowValidationError(
                f"{field_name} values must be non-empty strings"
            )
        return normalized
    numeric = float(value)
    if not math.isfinite(numeric):
        raise WorkflowValidationError(f"{field_name} values must be finite")
    return numeric


def _normalize_covariate_mapping(
    value: Mapping[str, str | int | float] | None,
    *,
    field_name: str,
) -> Mapping[str, str | float]:
    if value is None:
        return MappingProxyType({})
    if not isinstance(value, Mapping):
        raise WorkflowValidationError(f"{field_name} must be a mapping")
    normalized: dict[str, str | float] = {}
    for raw_key, raw_value in value.items():
        key = _require_non_empty_text(str(raw_key), field_name=f"{field_name}.keys")
        if key in normalized:
            raise WorkflowValidationError(
                f"{field_name} contains duplicate covariate name: {key}"
            )
        normalized[key] = _normalize_covariate_value(
            raw_value,
            field_name=f"{field_name}[{key!r}]",
        )
    return MappingProxyType(normalized)


@dataclass(frozen=True, slots=True)
class FixedEffectCovariate:
    """Explicit fixed-effect covariate declaration.

    The declaration records user intent only. It does not build a design matrix
    and does not make adjusted differential models executable.
    """

    name: str
    kind: FixedEffectCovariateKind
    required: bool = True
    include_in_model: bool = True

    def __post_init__(self) -> None:
        normalized_name = _require_non_empty_text(
            self.name,
            field_name="experimental_design.fixed_effects[].name",
        )
        normalized_kind = _normalize_fixed_effect_kind(
            self.kind,
            field_name="experimental_design.fixed_effects[].kind",
        )
        if (
            normalized_kind == FIXED_EFFECT_COVARIATE_KIND_BATCH
            and normalized_name != "batch"
        ):
            raise WorkflowValidationError(
                "experimental_design batch fixed-effect covariate must use name='batch'"
            )
        if (
            normalized_name == "batch"
            and normalized_kind != FIXED_EFFECT_COVARIATE_KIND_BATCH
        ):
            raise WorkflowValidationError(
                "experimental_design covariate name 'batch' is reserved for "
                "batch fixed effects"
            )
        object.__setattr__(self, "name", normalized_name)
        object.__setattr__(self, "kind", normalized_kind)
        object.__setattr__(
            self,
            "required",
            _require_bool(
                self.required,
                field_name="experimental_design.fixed_effects[].required",
            ),
        )
        object.__setattr__(
            self,
            "include_in_model",
            _require_bool(
                self.include_in_model,
                field_name="experimental_design.fixed_effects[].include_in_model",
            ),
        )

    @classmethod
    def categorical(
        cls,
        name: str,
        *,
        required: bool = True,
        include_in_model: bool = True,
    ) -> FixedEffectCovariate:
        return cls(
            name=name,
            kind=FIXED_EFFECT_COVARIATE_KIND_CATEGORICAL,
            required=required,
            include_in_model=include_in_model,
        )

    @classmethod
    def continuous(
        cls,
        name: str,
        *,
        required: bool = True,
        include_in_model: bool = True,
    ) -> FixedEffectCovariate:
        return cls(
            name=name,
            kind=FIXED_EFFECT_COVARIATE_KIND_CONTINUOUS,
            required=required,
            include_in_model=include_in_model,
        )

    @classmethod
    def batch(
        cls,
        *,
        required: bool = True,
        include_in_model: bool = True,
    ) -> FixedEffectCovariate:
        return cls(
            name="batch",
            kind=FIXED_EFFECT_COVARIATE_KIND_BATCH,
            required=required,
            include_in_model=include_in_model,
        )


@dataclass(frozen=True, slots=True)
class CategoricalCovariate(FixedEffectCovariate):
    """Explicit categorical fixed-effect covariate declaration."""

    kind: FixedEffectCovariateKind = field(
        default=FIXED_EFFECT_COVARIATE_KIND_CATEGORICAL,
        init=False,
    )


@dataclass(frozen=True, slots=True)
class ContinuousCovariate(FixedEffectCovariate):
    """Explicit continuous fixed-effect covariate declaration."""

    kind: FixedEffectCovariateKind = field(
        default=FIXED_EFFECT_COVARIATE_KIND_CONTINUOUS,
        init=False,
    )


@dataclass(frozen=True, slots=True)
class BatchCovariate(FixedEffectCovariate):
    """Explicit batch-as-fixed-effect covariate declaration."""

    name: str = "batch"
    kind: FixedEffectCovariateKind = field(
        default=FIXED_EFFECT_COVARIATE_KIND_BATCH,
        init=False,
    )


def _normalize_fixed_effects(
    values: Iterable[FixedEffectCovariate],
) -> tuple[FixedEffectCovariate, ...]:
    if isinstance(values, str) or not isinstance(values, Iterable):
        raise WorkflowValidationError(
            "experimental_design.fixed_effects must be an iterable of "
            "FixedEffectCovariate values"
        )
    normalized: list[FixedEffectCovariate] = []
    for value in values:
        if not isinstance(value, FixedEffectCovariate):
            raise WorkflowValidationError(
                "experimental_design.fixed_effects must contain "
                "FixedEffectCovariate values"
            )
        normalized.append(value)
    return tuple(normalized)


@dataclass(frozen=True, slots=True)
class SampleDesignRecord:
    """Design metadata for a single sample."""

    sample_id: str
    condition: str
    biological_replicate_id: str | None = None
    technical_replicate_id: str | None = None
    batch: str | None = None
    block: str | None = None
    covariates: Mapping[str, str | int | float] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "sample_id",
            _require_non_empty_text(
                self.sample_id,
                field_name="experimental_design.samples[].sample_id",
            ),
        )
        object.__setattr__(
            self,
            "condition",
            _require_non_empty_text(
                self.condition,
                field_name="experimental_design.samples[].condition",
            ),
        )
        object.__setattr__(
            self,
            "biological_replicate_id",
            _normalize_optional_text(
                self.biological_replicate_id,
                field_name="experimental_design.samples[].biological_replicate_id",
            ),
        )
        object.__setattr__(
            self,
            "technical_replicate_id",
            _normalize_optional_text(
                self.technical_replicate_id,
                field_name="experimental_design.samples[].technical_replicate_id",
            ),
        )
        object.__setattr__(
            self,
            "batch",
            _normalize_optional_text(
                self.batch,
                field_name="experimental_design.samples[].batch",
            ),
        )
        object.__setattr__(
            self,
            "block",
            _normalize_optional_text(
                self.block,
                field_name="experimental_design.samples[].block",
            ),
        )
        object.__setattr__(
            self,
            "covariates",
            _normalize_covariate_mapping(
                self.covariates,
                field_name="experimental_design.samples[].covariates",
            ),
        )


@dataclass(frozen=True, slots=True)
class ExperimentalDesign:
    """Typed contract for workflow-consumed sample-to-condition design.

    Experimental design is independent from passive dataset sample metadata.
    """

    samples: tuple[SampleDesignRecord, ...]
    fixed_effects: tuple[FixedEffectCovariate, ...] = ()

    def __post_init__(self) -> None:
        samples = tuple(self.samples)
        if not samples:
            raise WorkflowValidationError(
                "experimental_design.samples must contain at least one sample"
            )
        for record in samples:
            if not isinstance(record, SampleDesignRecord):
                raise WorkflowValidationError(
                    "experimental_design.samples must contain SampleDesignRecord values"
                )
        sample_ids = [record.sample_id for record in samples]
        duplicate_sample_ids = sorted(
            {sample_id for sample_id in sample_ids if sample_ids.count(sample_id) > 1}
        )
        if duplicate_sample_ids:
            joined = ", ".join(duplicate_sample_ids)
            raise WorkflowValidationError(
                f"experimental_design contains duplicate sample IDs: {joined}"
            )

        fixed_effects = _normalize_fixed_effects(self.fixed_effects)
        fixed_effect_names = [covariate.name for covariate in fixed_effects]
        duplicate_fixed_effects = sorted(
            {name for name in fixed_effect_names if fixed_effect_names.count(name) > 1}
        )
        if duplicate_fixed_effects:
            raise WorkflowValidationError(
                "experimental_design contains duplicate covariate names: "
                + ", ".join(duplicate_fixed_effects)
            )
        reserved = sorted(
            covariate.name
            for covariate in fixed_effects
            if covariate.name in _RESERVED_SAMPLE_FIELD_NAMES
        )
        if reserved:
            raise WorkflowValidationError(
                "experimental_design covariate names are reserved: "
                + ", ".join(reserved)
            )

        object.__setattr__(self, "samples", samples)
        object.__setattr__(self, "fixed_effects", fixed_effects)

    def sample_ids(self) -> tuple[str, ...]:
        return tuple(record.sample_id for record in self.samples)

    def condition_labels(self) -> tuple[str, ...]:
        seen: set[str] = set()
        ordered: list[str] = []
        for record in self.samples:
            if record.condition in seen:
                continue
            seen.add(record.condition)
            ordered.append(record.condition)
        return tuple(ordered)

    @classmethod
    def from_records(cls, records: Iterable[SampleDesignRecord]) -> ExperimentalDesign:
        return cls(samples=tuple(records))


@dataclass(frozen=True, slots=True)
class Contrast:
    """Typed condition-vs-condition contrast definition."""

    name: str
    numerator_condition: str
    denominator_condition: str

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "name",
            _require_non_empty_text(
                self.name,
                field_name="contrast.name",
            ),
        )
        object.__setattr__(
            self,
            "numerator_condition",
            _require_non_empty_text(
                self.numerator_condition,
                field_name="contrast.numerator_condition",
            ),
        )
        object.__setattr__(
            self,
            "denominator_condition",
            _require_non_empty_text(
                self.denominator_condition,
                field_name="contrast.denominator_condition",
            ),
        )
        if self.numerator_condition == self.denominator_condition:
            raise WorkflowValidationError(
                "contrast numerator and denominator conditions must differ"
            )


__all__ = [
    "BatchCovariate",
    "CategoricalCovariate",
    "Contrast",
    "ContinuousCovariate",
    "FIXED_EFFECT_COVARIATE_KIND_BATCH",
    "FIXED_EFFECT_COVARIATE_KIND_CATEGORICAL",
    "FIXED_EFFECT_COVARIATE_KIND_CONTINUOUS",
    "FixedEffectCovariate",
    "FixedEffectCovariateKind",
    "ExperimentalDesign",
    "SUPPORTED_FIXED_EFFECT_COVARIATE_KINDS",
    "SampleDesignRecord",
]
