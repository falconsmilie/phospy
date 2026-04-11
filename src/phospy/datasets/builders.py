from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import pandas as pd

from ..errors import InputCompatibilityError
from ..internal.constants import (
    CENTRALIZED_SEQUENCE_COLUMN,
    GENE_P_SITE_COLUMN,
    SITE_MATRIX_ID_COLUMN,
)
from ..preprocessing.site_matrix import SiteMatrixBuilder, SiteMatrixResult
from .schema import DatasetSchema

__all__ = ["DatasetSiteMatrix"]


@dataclass(frozen=True, slots=True)
class DatasetSiteMatrix:
    """Bound site-matrix builder for a validated phosphoproteomics dataset."""

    schema: DatasetSchema

    def _builder(self) -> SiteMatrixBuilder:
        return SiteMatrixBuilder(value_cols=self.schema.corrected_cols)

    def build(
        self,
        corrected_df: pd.DataFrame,
        *,
        gene_p_site_col: str = GENE_P_SITE_COLUMN,
        sequence_col: str = CENTRALIZED_SEQUENCE_COLUMN,
    ) -> SiteMatrixResult:
        return self._builder().build(
            corrected_df,
            gene_p_site_col=gene_p_site_col,
            sequence_col=sequence_col,
        )


def _build_site_metadata(
    *,
    phosr_input: pd.DataFrame,
    corrected_cols: Sequence[str],
) -> pd.DataFrame:
    if SITE_MATRIX_ID_COLUMN not in phosr_input.columns:
        msg = (
            "Analysis-ready site metadata requires the site-matrix source table to "
            f"include '{SITE_MATRIX_ID_COLUMN}'."
        )
        raise InputCompatibilityError(msg)

    metadata = phosr_input.drop(
        columns=[*tuple(corrected_cols), CENTRALIZED_SEQUENCE_COLUMN],
        errors="ignore",
    ).copy(deep=True)
    metadata = metadata.set_index(SITE_MATRIX_ID_COLUMN)
    metadata.index = pd.Index(
        metadata.index.astype("string"),
        name=SITE_MATRIX_ID_COLUMN,
    )
    return metadata


def _validate_analysis_ready_alignment(
    *,
    phospho_matrix: pd.DataFrame,
    site_metadata: pd.DataFrame,
    site_sequences: pd.Series,
) -> None:
    matrix_index = pd.Index(
        phospho_matrix.index.astype("string"),
        name=SITE_MATRIX_ID_COLUMN,
    )
    metadata_index = pd.Index(
        site_metadata.index.astype("string"),
        name=SITE_MATRIX_ID_COLUMN,
    )
    sequence_index = pd.Index(
        site_sequences.index.astype("string"),
        name=SITE_MATRIX_ID_COLUMN,
    )

    duplicate_index_contexts = {
        "phospho_matrix": matrix_index,
        "site_metadata": metadata_index,
        "site_sequences": sequence_index,
    }
    for name, index in duplicate_index_contexts.items():
        if not index.is_unique:
            msg = (
                "AnalysisReadyPhosphoDataset requires unique site identifiers; "
                f"{name} contains duplicate site IDs."
            )
            raise InputCompatibilityError(msg)

    if not matrix_index.equals(metadata_index) or not matrix_index.equals(
        sequence_index
    ):
        msg = (
            "AnalysisReadyPhosphoDataset requires phospho_matrix, site_metadata, "
            "and site_sequences to share the same aligned site_id index."
        )
        raise InputCompatibilityError(msg)

    phospho_matrix.index = matrix_index
    site_metadata.index = metadata_index
    site_sequences.index = sequence_index
