"""Peptide-row to resolved-site mapping and mapping-fraction validation."""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.science.evidence.dataset_resolution.contracts import (
    DATASET_PEPTIDE_MAPPING_WEIGHT_SOURCE_DERIVED_EQUAL,
    DATASET_PEPTIDE_MAPPING_WEIGHT_SOURCE_EXPLICIT,
    MAPPING_FRACTION_COLUMN,
    MAPPING_WEIGHT_SUM_TOLERANCE,
)


@dataclass(frozen=True, slots=True)
class JoinedMappingRows:
    """Peptide evidence rows joined to resolved site-mapping rows."""

    rows: pd.DataFrame


@dataclass(frozen=True, slots=True)
class ResolvedMappingFractions:
    """Resolved mapping rows with validated allocation fractions."""

    rows: pd.DataFrame
    mapping_weight_source: str


def join_peptide_rows_to_site_mapping(
    *,
    evidence_frame: pd.DataFrame,
    mapping: pd.DataFrame,
    sample_columns: tuple[str, ...],
) -> JoinedMappingRows:
    """Join validated peptide evidence rows onto validated site mapping rows."""

    peptide_fields = [
        "peptide_row_id",
        "protein_accession",
        "site_string",
        "peptide_sequence",
        "multi_site",
    ]
    if "site_sequence" in evidence_frame.columns:
        peptide_fields.append("site_sequence")
    if "localisation_confidence" in evidence_frame.columns:
        peptide_fields.append("localisation_confidence")
    peptide_rows = evidence_frame.loc[:, peptide_fields + list(sample_columns)].copy(
        deep=True
    )
    return JoinedMappingRows(
        rows=mapping.merge(peptide_rows, how="inner", on="peptide_row_id")
    )


def resolve_and_validate_mapping_fractions(
    *,
    joined_mapping: JoinedMappingRows,
) -> ResolvedMappingFractions:
    """Resolve explicit or derived mapping fractions and validate per-row totals."""

    resolved = joined_mapping.rows.copy(deep=True)
    mapping_weight_source = DATASET_PEPTIDE_MAPPING_WEIGHT_SOURCE_EXPLICIT
    if "mapping_weight" not in resolved.columns:
        counts = resolved.groupby("peptide_row_id", sort=False).size().astype(float)
        resolved.loc[:, "mapping_weight"] = resolved.loc[:, "peptide_row_id"].map(
            lambda peptide_row_id: float(1.0 / counts.loc[peptide_row_id])
        )
        mapping_weight_source = DATASET_PEPTIDE_MAPPING_WEIGHT_SOURCE_DERIVED_EQUAL
    mapping_fractions = pd.to_numeric(
        resolved.loc[:, "mapping_weight"],
        errors="coerce",
    )
    if mapping_fractions.isna().any() or (mapping_fractions <= 0.0).any():
        raise PhosPyInputError(
            "dataset build request peptide_evidence site mapping contains "
            "non-positive or non-numeric mapping_weight values"
        )
    per_peptide_weight_sum = mapping_fractions.groupby(
        resolved.loc[:, "peptide_row_id"]
    ).sum()
    invalid_weight_rows = per_peptide_weight_sum.loc[
        (per_peptide_weight_sum - 1.0).abs() > MAPPING_WEIGHT_SUM_TOLERANCE
    ]
    if not invalid_weight_rows.empty:
        preview = ", ".join(
            f"{str(peptide_row_id)!r}={float(total_weight):.6f}"
            for peptide_row_id, total_weight in invalid_weight_rows.iloc[:5].items()
        )
        suffix = "" if int(invalid_weight_rows.shape[0]) <= 5 else " ..."
        raise PhosPyInputError(
            "dataset build request peptide_evidence mapping_weight values must sum "
            "to 1.0 per peptide_row_id; invalid totals: "
            f"{preview}{suffix}"
        )
    resolved.loc[:, MAPPING_FRACTION_COLUMN] = mapping_fractions.to_numpy(dtype=float)
    return ResolvedMappingFractions(
        rows=resolved,
        mapping_weight_source=mapping_weight_source,
    )


__all__ = [
    "JoinedMappingRows",
    "ResolvedMappingFractions",
    "join_peptide_rows_to_site_mapping",
    "resolve_and_validate_mapping_fractions",
]
