"""Structured kinase workflow site attrition metrics."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass

import pandas as pd

from phospy.errors.workflows import PhosPyWorkflowError


@dataclass(frozen=True, slots=True)
class KinaseAttritionMetrics:
    """Pre-scoring kinase workflow site attrition counts and fractions."""

    total_dataset_sites: int
    reference_overlap_sites: int
    sequence_supported_sites: int
    scored_sites: int
    reference_overlap_fraction: float
    sequence_supported_fraction: float
    scored_fraction: float

    @classmethod
    def from_counts(
        cls,
        *,
        total_dataset_sites: int,
        reference_overlap_sites: int,
        sequence_supported_sites: int,
        scored_sites: int,
    ) -> KinaseAttritionMetrics:
        total = int(total_dataset_sites)
        if total <= 0:
            raise PhosPyWorkflowError(
                "kinase workflow internal attrition validation failed at "
                "seam=kinase.attrition_metrics.total_dataset_sites; "
                "total_dataset_sites must be greater than zero before "
                "fractional attrition metrics can be computed"
            )
        reference_overlap = _validate_count(
            value=reference_overlap_sites,
            field_name="reference_overlap_sites",
            total=total,
        )
        sequence_supported = _validate_count(
            value=sequence_supported_sites,
            field_name="sequence_supported_sites",
            total=total,
        )
        scored = _validate_count(
            value=scored_sites,
            field_name="scored_sites",
            total=total,
        )
        return cls(
            total_dataset_sites=total,
            reference_overlap_sites=reference_overlap,
            sequence_supported_sites=sequence_supported,
            scored_sites=scored,
            reference_overlap_fraction=reference_overlap / total,
            sequence_supported_fraction=sequence_supported / total,
            scored_fraction=scored / total,
        )

    def to_payload(self) -> dict[str, object]:
        return {
            "total_dataset_sites": int(self.total_dataset_sites),
            "reference_overlap_sites": int(self.reference_overlap_sites),
            "sequence_supported_sites": int(self.sequence_supported_sites),
            "scored_sites": int(self.scored_sites),
            "reference_overlap_fraction": float(self.reference_overlap_fraction),
            "sequence_supported_fraction": float(self.sequence_supported_fraction),
            "scored_fraction": float(self.scored_fraction),
        }


def build_kinase_attrition_metrics(
    *,
    dataset_site_index: pd.Index,
    reference_site_index: Iterable[object],
    sequence_supported_site_index: pd.Index,
) -> KinaseAttritionMetrics:
    """Build attrition metrics from resolved site identity sets."""

    dataset_sites = set(_string_values(dataset_site_index))
    reference_sites = set(_string_values(reference_site_index))
    sequence_supported_sites = set(_string_values(sequence_supported_site_index))
    reference_overlap_sites = dataset_sites.intersection(reference_sites)
    scored_sites = reference_overlap_sites.intersection(sequence_supported_sites)
    return KinaseAttritionMetrics.from_counts(
        total_dataset_sites=len(dataset_sites),
        reference_overlap_sites=len(reference_overlap_sites),
        sequence_supported_sites=len(
            dataset_sites.intersection(sequence_supported_sites)
        ),
        scored_sites=len(scored_sites),
    )


def build_kinase_attrition_metrics_from_overlap(
    *,
    total_dataset_sites: int,
    reference_overlap_site_ids: Iterable[object],
    sequence_supported_site_index: pd.Index,
) -> KinaseAttritionMetrics:
    """Build attrition metrics from pre-resolved reference-overlap identities."""

    reference_overlap_sites = set(_string_values(reference_overlap_site_ids))
    sequence_supported_sites = set(_string_values(sequence_supported_site_index))
    scored_sites = reference_overlap_sites.intersection(sequence_supported_sites)
    return KinaseAttritionMetrics.from_counts(
        total_dataset_sites=int(total_dataset_sites),
        reference_overlap_sites=len(reference_overlap_sites),
        sequence_supported_sites=len(sequence_supported_sites),
        scored_sites=len(scored_sites),
    )


def _validate_count(*, value: int, field_name: str, total: int) -> int:
    count = int(value)
    if count < 0:
        raise PhosPyWorkflowError(
            "kinase workflow internal attrition validation failed at "
            f"seam=kinase.attrition_metrics.{field_name}; "
            f"{field_name} must be greater than or equal to zero"
        )
    if count > total:
        raise PhosPyWorkflowError(
            "kinase workflow internal attrition validation failed at "
            f"seam=kinase.attrition_metrics.{field_name}; "
            f"{field_name} must not exceed total_dataset_sites; "
            f"{field_name}={count}; total_dataset_sites={total}"
        )
    return count


def _string_values(values: Iterable[object]) -> list[str]:
    return [str(value) for value in values]


__all__ = [
    "KinaseAttritionMetrics",
    "build_kinase_attrition_metrics",
    "build_kinase_attrition_metrics_from_overlap",
]
