"""Total-protein correction stage for dataset preprocessing."""

from __future__ import annotations

from dataclasses import replace

import pandas as pd

from phospy.api.configs import (
    DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_NONE,
    DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_RATIO_TO_TOTAL,
)
from phospy.datasets.preprocessing.models import (
    DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION,
    PreprocessingState,
)
from phospy.errors.input import PhosPyInputError

_GENE_SYMBOL_COLUMN = "gene_symbol"


class TotalProteinCorrectionStage:
    """Apply historical-baseline phospho-to-total correction when requested."""

    stage_key = DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION

    def run(self, state: PreprocessingState) -> PreprocessingState:
        policy = state.plan.total_protein_correction_policy
        if policy == DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_NONE:
            return state
        if policy != DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_RATIO_TO_TOTAL:
            raise PhosPyInputError(
                "dataset build request preprocessing_config contains an unsupported "
                "total_protein_correction.policy"
            )

        if state.total is None:
            raise PhosPyInputError(
                "dataset build request preprocessing_config.total_protein_correction."
                "policy='ratio_to_total' requires total input data"
            )
        if _GENE_SYMBOL_COLUMN not in state.site_metadata.columns:
            raise PhosPyInputError(
                "dataset build request preprocessing total/protein correction requires "
                "site_metadata.gene_symbol"
            )
        if not state.total.columns.equals(state.phospho.columns):
            raise PhosPyInputError(
                "dataset build request preprocessing total/protein correction requires "
                "total columns to exactly match phospho columns"
            )
        non_numeric_phospho_columns = [
            str(column)
            for column in state.phospho.columns
            if not pd.api.types.is_numeric_dtype(state.phospho[column])
        ]
        if non_numeric_phospho_columns:
            raise PhosPyInputError(
                "dataset build request preprocessing total/protein correction requires "
                "numeric phospho columns; non-numeric columns: "
                + ", ".join(non_numeric_phospho_columns)
            )
        non_numeric_total_columns = [
            str(column)
            for column in state.total.columns
            if not pd.api.types.is_numeric_dtype(state.total[column])
        ]
        if non_numeric_total_columns:
            raise PhosPyInputError(
                "dataset build request preprocessing total/protein correction requires "
                "numeric total columns; non-numeric columns: "
                + ", ".join(non_numeric_total_columns)
            )

        phospho_gene_key = _normalize_identifier_series(
            state.site_metadata.loc[:, _GENE_SYMBOL_COLUMN]
        )
        total_gene_key = _normalize_identifier_series(
            pd.Series(state.total.index, index=state.total.index, dtype="string")
        )
        duplicate_total_gene_mask = total_gene_key.duplicated(keep=False)
        if bool(duplicate_total_gene_mask.any()):
            duplicate_preview = ", ".join(
                total_gene_key.loc[duplicate_total_gene_mask]
                .dropna()
                .astype(str)
                .unique()
                .tolist()[:5]
            )
            raise PhosPyInputError(
                "dataset build request preprocessing total/protein correction requires "
                "unique total protein identifiers after normalization; duplicate "
                f"identifiers: {duplicate_preview}"
            )

        total_lookup = state.total.set_index(total_gene_key)
        matched_mask = phospho_gene_key.isin(total_lookup.index)
        matched_rows = int(matched_mask.sum())
        if matched_rows == 0:
            raise PhosPyInputError(
                "dataset build request preprocessing total/protein correction found "
                "no overlapping identifiers between site_metadata.gene_symbol and "
                "total.index"
            )

        input_rows = int(len(phospho_gene_key))
        unmatched_rows = input_rows - matched_rows
        if unmatched_rows:
            unmatched_preview = ", ".join(
                phospho_gene_key.loc[~matched_mask]
                .dropna()
                .astype(str)
                .unique()
                .tolist()[:5]
            )
            percentage = (unmatched_rows / input_rows) * 100.0
            raise PhosPyInputError(
                "dataset build request preprocessing total/protein correction "
                f"requires complete phospho/total matching but would drop "
                f"{unmatched_rows} of {input_rows} rows ({percentage:.1f}%): "
                f"{unmatched_preview}"
            )

        matched_total = total_lookup.loc[phospho_gene_key.tolist()]
        corrected = pd.DataFrame(
            state.phospho.to_numpy(copy=False) - matched_total.to_numpy(copy=False),
            index=state.phospho.index,
            columns=state.phospho.columns,
        )
        return replace(state, phospho=corrected)


def _normalize_identifier_series(series: pd.Series) -> pd.Series:
    return series.astype("string").str.strip().str.upper()


__all__ = ["TotalProteinCorrectionStage"]
