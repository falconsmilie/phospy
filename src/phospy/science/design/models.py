"""Typed experimental-design and contrast contracts."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass

from phospy.errors.validation import WorkflowValidationError


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


@dataclass(frozen=True, slots=True)
class SampleDesignRecord:
    """Design metadata for a single sample."""

    sample_id: str
    condition: str
    biological_replicate_id: str | None = None
    technical_replicate_id: str | None = None
    batch: str | None = None
    block: str | None = None

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


@dataclass(frozen=True, slots=True)
class ExperimentalDesign:
    """Typed contract for sample-to-condition experimental design."""

    samples: tuple[SampleDesignRecord, ...]

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
        object.__setattr__(self, "samples", samples)

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
    "Contrast",
    "ExperimentalDesign",
    "SampleDesignRecord",
]
