"""Allocation of peptide signals to resolved sites and site-level summarisation."""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.science.evidence.dataset_resolution.mapping import ResolvedMappingFractions
from phospy.science.evidence.dataset_resolution.models import (
    MAPPING_FRACTION_COLUMN,
    PeptideToSiteAggregationPolicy,
)


@dataclass(frozen=True, slots=True)
class AllocatedEvidence:
    """Resolved mapping rows with sample signals allocated by mapping fraction."""

    rows: pd.DataFrame
    sample_columns: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class SiteSignalSummary:
    """Site-level phospho matrix produced from allocated evidence rows."""

    phospho: pd.DataFrame


def allocate_peptide_signals_to_resolved_sites(
    *,
    resolved_mapping: ResolvedMappingFractions,
    sample_columns: tuple[str, ...],
    aggregation_policy: PeptideToSiteAggregationPolicy,
) -> AllocatedEvidence:
    """Allocate peptide signals under a typed peptide-to-site estimand policy."""

    if not isinstance(aggregation_policy, PeptideToSiteAggregationPolicy):
        raise PhosPyInputError(
            "peptide-to-site signal allocation requires a typed "
            "PeptideToSiteAggregationPolicy; mapping fractions cannot be applied "
            "from untyped scale or column metadata"
        )

    allocated = resolved_mapping.rows.copy(deep=True)
    aggregation_policy.validate_allocation_rows(allocated)
    mapping_fractions = allocated.loc[:, MAPPING_FRACTION_COLUMN].to_numpy(dtype=float)
    for sample_column in sample_columns:
        allocated.loc[:, sample_column] = (
            pd.to_numeric(allocated.loc[:, sample_column], errors="coerce")
            * mapping_fractions
        )
    return AllocatedEvidence(rows=allocated, sample_columns=sample_columns)


def summarise_allocated_site_signals(
    *,
    allocated_evidence: AllocatedEvidence,
) -> SiteSignalSummary:
    """Summarise allocated peptide rows to site-level arithmetic means."""

    matrix = (
        allocated_evidence.rows.groupby("site_id", sort=True)[
            list(allocated_evidence.sample_columns)
        ]
        .mean(numeric_only=True)
        .astype(float)
    )
    matrix.index = pd.Index(matrix.index.astype(str), name="site_id")
    return SiteSignalSummary(phospho=matrix)


__all__ = [
    "AllocatedEvidence",
    "SiteSignalSummary",
    "allocate_peptide_signals_to_resolved_sites",
    "summarise_allocated_site_signals",
]
