from __future__ import annotations

from collections.abc import Sequence

import pandas as pd

from .validation.tables import split_gene_p_site


def build_site_matrix(
    df: pd.DataFrame,
    gene_p_site_col: str,
    sequence_col: str,
    value_cols: Sequence[str],
    gene_col_name: str = "gene",
    p_site_col_name: str = "p_site",
) -> tuple[pd.DataFrame, pd.DataFrame, pd.Series]:
    work = df.copy()

    split_cols = split_gene_p_site(work[gene_p_site_col], context=gene_p_site_col)
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
    complete_cases["__row_order"] = range(len(complete_cases))

    sort_cols = ["site_id", "__mean_signal"]
    ascending = [True, False]
    if "uid" in complete_cases.columns:
        complete_cases["__uid_sort"] = complete_cases["uid"].astype("string")
        sort_cols.append("__uid_sort")
        ascending.append(True)
    sort_cols.append("__row_order")
    ascending.append(True)

    phosr_input = (
        complete_cases.sort_values(sort_cols, ascending=ascending, kind="mergesort")
        .drop_duplicates(subset=["site_id"], keep="first")
        .drop(
            columns=[
                col
                for col in ["__mean_signal", "__uid_sort", "__row_order"]
                if col in complete_cases.columns
            ]
        )
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
