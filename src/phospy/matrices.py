from __future__ import annotations

from collections.abc import Mapping, Sequence

import pandas as pd

from .validation.errors import TableSchemaError


def format_row_drop_diagnostics(row_drop_stats: Mapping[str, int]) -> str:
    """Format human-readable site-matrix row-drop diagnostics."""

    stats = {
        "input_rows": int(row_drop_stats.get("input_rows", 0)),
        "dropped_missing_sequence": int(
            row_drop_stats.get("dropped_missing_sequence", 0)
        ),
        "dropped_incomplete_values": int(
            row_drop_stats.get("dropped_incomplete_values", 0)
        ),
        "deduplicated_site_rows": int(row_drop_stats.get("deduplicated_site_rows", 0)),
        "retained_rows": int(row_drop_stats.get("retained_rows", 0)),
    }
    known_drops = (
        stats["dropped_missing_sequence"]
        + stats["dropped_incomplete_values"]
        + stats["deduplicated_site_rows"]
    )
    other_dropped_rows = max(
        stats["input_rows"] - stats["retained_rows"] - known_drops, 0
    )
    return (
        "row-drop diagnostics: "
        f"input_rows={stats['input_rows']}, "
        f"dropped_missing_sequence={stats['dropped_missing_sequence']}, "
        f"dropped_incomplete_values={stats['dropped_incomplete_values']}, "
        f"deduplicated_site_rows={stats['deduplicated_site_rows']}, "
        f"other_dropped_rows={other_dropped_rows}, "
        f"retained_rows={stats['retained_rows']}"
    )


def build_site_matrix(
    df: pd.DataFrame,
    gene_p_site_col: str,
    sequence_col: str,
    value_cols: Sequence[str],
    gene_col_name: str = "gene",
    p_site_col_name: str = "p_site",
) -> tuple[pd.DataFrame, pd.DataFrame, pd.Series]:
    work = df.copy()

    split_cols = work[gene_p_site_col].astype("string").str.split("_", n=1, expand=True)
    invalid_mask = (
        (split_cols.shape[1] < 2)
        or split_cols[0].isna().any()
        or split_cols[1].isna().any()
    )
    if invalid_mask:
        raise TableSchemaError(
            f"{gene_p_site_col} must contain values in the form GENE_SITE, for example PRKACA_S339"
        )
    work[gene_col_name] = split_cols[0].astype("string")
    work[p_site_col_name] = split_cols[1].astype("string")
    work["site_id"] = (
        work[gene_col_name].str.upper() + ";" + work[p_site_col_name].str.upper() + ";"
    )

    base_cols = [
        col
        for col in [
            "gene_names",
            gene_col_name,
            p_site_col_name,
            "uid",
            sequence_col,
            "site_id",
        ]
        if col in work.columns
    ]
    keep_cols: list[str] = [*base_cols, *list(value_cols)]
    phosr_input = work.loc[:, keep_cols].copy()

    total_rows = len(phosr_input)
    with_sequence = phosr_input[phosr_input[sequence_col].notna()].copy()
    dropped_missing_sequence = total_rows - len(with_sequence)

    complete_cases = with_sequence[
        with_sequence.loc[:, list(value_cols)].notna().all(axis=1)
    ].copy()
    dropped_incomplete_values = len(with_sequence) - len(complete_cases)

    complete_cases["__mean_signal"] = complete_cases.loc[:, list(value_cols)].mean(
        axis=1, skipna=True
    )
    idx = complete_cases.groupby("site_id")["__mean_signal"].idxmax()
    phosr_input = (
        complete_cases.loc[idx]
        .drop(columns="__mean_signal")
        .copy()
        .reset_index(drop=True)
    )
    deduplicated_site_rows = len(complete_cases) - len(phosr_input)

    row_drop_stats = {
        "input_rows": total_rows,
        "dropped_missing_sequence": dropped_missing_sequence,
        "dropped_incomplete_values": dropped_incomplete_values,
        "deduplicated_site_rows": deduplicated_site_rows,
        "retained_rows": len(phosr_input),
    }

    matrix = phosr_input.loc[:, ["site_id", *value_cols]].set_index("site_id")
    matrix.attrs["row_drop_stats"] = row_drop_stats.copy()
    sequences = phosr_input.set_index("site_id")[sequence_col].copy()
    sequences.attrs["row_drop_stats"] = row_drop_stats.copy()
    phosr_input.attrs["row_drop_stats"] = row_drop_stats.copy()
    return phosr_input, matrix, sequences
