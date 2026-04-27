"""Total-protein correction stage for dataset preprocessing."""

from __future__ import annotations

from dataclasses import replace

import pandas as pd

from phospy.api.configs import (
    DATASET_INTENSITY_TRANSFORM_POLICY_LOG2,
    DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_NONE,
    DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_SUBTRACT_LOG_TOTAL,
    resolve_dataset_total_protein_correction_policy,
)
from phospy.datasets.preprocessing.models import (
    DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION,
    PreprocessingStageResult,
    PreprocessingState,
)
from phospy.datasets.processing_state import (
    TOTAL_PROTEIN_CORRECTION_DIAGNOSTICS_SCHEMA_VERSION_V1,
)
from phospy.errors.input import PhosPyInputError
from phospy.provenance.hashing import hash_table

_GENE_SYMBOL_COLUMN = "gene_symbol"


class TotalProteinCorrectionStage:
    """Apply log-scale phospho-to-total subtraction when requested."""

    stage_key = DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION

    def run(self, state: PreprocessingState) -> PreprocessingStageResult:
        requested_policy = state.plan.total_protein_correction_policy
        if requested_policy == DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_NONE:
            return PreprocessingStageResult(
                state=state,
                diagnostics={
                    "dropped_row_ids": (),
                    "dropped_row_count": 0,
                    "imputed_cell_count": 0,
                    "imputed_row_ids": (),
                    "notes": "stage executed",
                    "diagnostics": {
                        "diagnostics_schema_version": (
                            TOTAL_PROTEIN_CORRECTION_DIAGNOSTICS_SCHEMA_VERSION_V1
                        ),
                        "policy": str(requested_policy),
                        "requested_policy": str(requested_policy),
                        "resolved_policy": str(requested_policy),
                    },
                },
            )
        resolved_policy = resolve_dataset_total_protein_correction_policy(
            requested_policy
        )
        if (
            resolved_policy
            != DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_SUBTRACT_LOG_TOTAL
        ):
            raise PhosPyInputError(
                "dataset build request preprocessing_config contains an unsupported "
                "total_protein_correction.policy"
            )
        if (
            state.plan.intensity_transform_policy
            != DATASET_INTENSITY_TRANSFORM_POLICY_LOG2
        ):
            raise PhosPyInputError(
                "dataset build request "
                f"preprocessing_config.total_protein_correction.policy={requested_policy!r} "
                "requires log2-scale phospho and total values. Configure "
                "preprocessing_config.intensity_transform.policy='log2', or disable "
                "total-protein correction."
            )

        if state.total is None:
            raise PhosPyInputError(
                "dataset build request preprocessing_config.total_protein_correction."
                f"policy={requested_policy!r} requires total input data"
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
        # Subtractive correction on log2-scale data yields a log2 ratio:
        # corrected = log2_phospho - log2_total.
        corrected = pd.DataFrame(
            state.phospho.to_numpy(copy=False) - matched_total.to_numpy(copy=False),
            index=state.phospho.index,
            columns=state.phospho.columns,
        )
        next_state = replace(state, phospho=corrected)
        return PreprocessingStageResult(
            state=next_state,
            diagnostics={
                "dropped_row_ids": (),
                "dropped_row_count": 0,
                "imputed_cell_count": 0,
                "imputed_row_ids": (),
                "notes": "stage executed",
                "diagnostics": {
                    "diagnostics_schema_version": (
                        TOTAL_PROTEIN_CORRECTION_DIAGNOSTICS_SCHEMA_VERSION_V1
                    ),
                    "policy": str(resolved_policy),
                    "requested_policy": str(requested_policy),
                    "resolved_policy": str(resolved_policy),
                    "formula": "log2_phospho - log2_total",
                    "requires_log_scale": True,
                    "input_scale": "log2",
                    "output_scale": "log2_ratio",
                    "quantitative_meaning": "phospho_total_log_ratio",
                    "matched_rows": matched_rows,
                    "total_table_hash": (
                        None
                        if state.total is None
                        else hash_table(
                            state.total,
                            name="total_protein_correction.total",
                        )
                    ),
                    "input_phospho_hash": hash_table(
                        state.phospho,
                        name="total_protein_correction.input.phospho",
                    ),
                    "output_phospho_hash": hash_table(
                        corrected,
                        name="total_protein_correction.output.phospho",
                    ),
                },
            },
        )


def _normalize_identifier_series(series: pd.Series) -> pd.Series:
    return series.astype("string").str.strip().str.upper()


__all__ = ["TotalProteinCorrectionStage"]
