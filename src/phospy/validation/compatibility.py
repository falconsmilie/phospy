from __future__ import annotations

from collections.abc import Mapping, Sequence

import pandas as pd

from .errors import InputCompatibilityError
from .tables import SiteMatrixSchema


def validate_pred_mat_overlap(
    pred_mat: pd.DataFrame,
    phospho_matrix: pd.DataFrame,
    *,
    pred_context: str = "pred_mat",
    matrix_context: str = "phospho_matrix",
) -> None:
    overlap = pred_mat.index.intersection(phospho_matrix.index)
    if overlap.empty:
        msg = f"{pred_context} and {matrix_context} have no overlapping phosphosite IDs"
        raise InputCompatibilityError(msg)


def validate_workflow_inputs(
    phospho_matrix: pd.DataFrame,
    substrate_map: Mapping[str, Sequence[str]],
    site_sequences: Mapping[str, str] | Sequence[str] | pd.Series | None,
    motif_sequences: Mapping[str, Sequence[str]] | None,
    *,
    context: str = "Kinase workflow inputs",
) -> pd.DataFrame:
    validated_matrix = SiteMatrixSchema.validate(
        phospho_matrix,
        context="phospho_matrix",
    )

    overlapping_sites = {
        site
        for sites in substrate_map.values()
        for site in sites
        if site in validated_matrix.index
    }
    if not overlapping_sites:
        msg = f"{context} contain no overlap between substrate_map and phospho_matrix"
        raise InputCompatibilityError(msg)

    if site_sequences is not None:
        sequence_index = _extract_sequence_index(site_sequences, validated_matrix.index)
        missing = [
            site for site in validated_matrix.index if site not in sequence_index
        ]
        if missing:
            missing_preview = ", ".join(missing[:5])
            msg = (
                f"site_sequences is missing entries for phosphosites: {missing_preview}"
            )
            raise InputCompatibilityError(msg)

    if motif_sequences is not None:
        widths: set[int] = set()
        for kinase, sequences in motif_sequences.items():
            if not sequences:
                continue
            kinase_widths = {len(str(sequence)) for sequence in sequences}
            if len(kinase_widths) != 1:
                msg = (
                    f"motif_sequences for kinase {kinase} must use a consistent "
                    "sequence width"
                )
                raise InputCompatibilityError(msg)
            widths.update(kinase_widths)
        if len(widths) > 1:
            msg = "motif_sequences must use the same sequence width across kinases"
            raise InputCompatibilityError(msg)

    return validated_matrix


def _extract_sequence_index(
    site_sequences: Mapping[str, str] | Sequence[str] | pd.Series,
    site_index: pd.Index,
) -> set[str]:
    if isinstance(site_sequences, pd.Series):
        return {str(value) for value in site_sequences.index}
    if isinstance(site_sequences, Mapping):
        return {str(value) for value in site_sequences}
    if len(site_sequences) != len(site_index):
        msg = (
            "site_sequences must match phospho_matrix length when passed as a sequence"
        )
        raise InputCompatibilityError(msg)
    return {str(value) for value in site_index}
