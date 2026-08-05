"""Allocation of peptide signals to resolved sites and site-level summarisation."""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.science.evidence.dataset_resolution.contracts import (
    MAPPING_FRACTION_COLUMN,
    MAPPING_WEIGHT_SUM_TOLERANCE,
)
from phospy.science.evidence.dataset_resolution.mapping import ResolvedMappingFractions
from phospy.science.transformations.models import IntensityScaleKind


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
    input_intensity_scale: IntensityScaleKind,
) -> AllocatedEvidence:
    """Multiply each peptide signal by its resolved mapping fraction."""

    allocated = resolved_mapping.rows.copy(deep=True)
    mapping_fractions = allocated.loc[:, MAPPING_FRACTION_COLUMN].to_numpy(dtype=float)
    _reject_fractional_allocation_for_non_linear_input(
        rows=allocated,
        input_intensity_scale=input_intensity_scale,
    )
    for sample_column in sample_columns:
        allocated.loc[:, sample_column] = (
            pd.to_numeric(allocated.loc[:, sample_column], errors="coerce")
            * mapping_fractions
        )
    return AllocatedEvidence(rows=allocated, sample_columns=sample_columns)


def _reject_fractional_allocation_for_non_linear_input(
    *,
    rows: pd.DataFrame,
    input_intensity_scale: IntensityScaleKind,
) -> None:
    if input_intensity_scale is IntensityScaleKind.LINEAR:
        return
    mapping_fractions = pd.to_numeric(
        rows.loc[:, MAPPING_FRACTION_COLUMN],
        errors="coerce",
    )
    fractional_mask = mapping_fractions.sub(1.0).abs() > MAPPING_WEIGHT_SUM_TOLERANCE
    if not bool(fractional_mask.any()):
        return

    fractional_rows = rows.loc[fractional_mask, ["peptide_row_id", "site_id"]].head(5)
    preview = ", ".join(
        f"{str(row.peptide_row_id)!r}->{str(row.site_id)!r}"
        for row in fractional_rows.itertuples(index=False)
    )
    suffix = "" if int(fractional_mask.sum()) <= 5 else " ..."
    raise PhosPyInputError(
        "dataset peptide-evidence mode declared input_intensity_scale="
        f"{input_intensity_scale.value!r} cannot use fractional allocation: "
        "mapping_fraction contains non-unit values for peptide-to-site mappings "
        f"({preview}{suffix}). Supported corrective action: provide linear "
        "peptide evidence before applying split/fractional allocation, or use only "
        "unit/unambiguous peptide-to-site mappings until a validated scale-aware "
        "peptide allocation estimator is implemented."
    )


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
