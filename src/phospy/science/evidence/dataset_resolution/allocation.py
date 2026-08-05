"""Allocation of peptide signals to resolved sites and site-level summarisation."""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from phospy.science.evidence.dataset_resolution.contracts import (
    MAPPING_FRACTION_COLUMN,
)
from phospy.science.evidence.dataset_resolution.mapping import ResolvedMappingFractions


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
) -> AllocatedEvidence:
    """Multiply each peptide signal by its resolved mapping fraction."""

    allocated = resolved_mapping.rows.copy(deep=True)
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
