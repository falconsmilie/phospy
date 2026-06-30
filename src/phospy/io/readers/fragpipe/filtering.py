"""FragPipe contaminant and decoy filtering."""
# pyright: reportMissingTypeStubs=false

from __future__ import annotations

import pandas as pd

from phospy.io.readers._table_parsing import (
    raise_for_forbidden_flags,
    resolve_flag_series,
    split_multi_value,
)
from phospy.io.readers.fragpipe.constants import (
    _CONTAMINANT_PREFIXES,
    _DECOY_PREFIXES,
    _FRAGPIPE_CONTAMINANT_OUTPUT_COLUMN,
    _FRAGPIPE_DECOY_OUTPUT_COLUMN,
)
from phospy.io.readers.fragpipe.models import _ResolvedFragPipeColumns
from phospy.validation.datasets.fragpipe import (
    FRAGPIPE_FLAG_POLICY_ERROR,
    FRAGPIPE_FLAG_POLICY_REMOVE,
)


def _apply_flag_policies(
    source: pd.DataFrame,
    *,
    resolved: _ResolvedFragPipeColumns,
    contaminant_policy: str,
    decoy_policy: str,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object], tuple[str, ...]]:
    protein_values = source.loc[:, resolved.protein_accession]
    contaminant_flags = resolve_flag_series(
        source,
        column=resolved.contaminant,
        field_name="FragPipe contaminant flag",
    )
    decoy_flags = resolve_flag_series(
        source,
        column=resolved.decoy,
        field_name="FragPipe decoy flag",
    )
    contaminant_prefix_flags = protein_values.map(
        lambda value: _has_any_prefixed_token(value, prefixes=_CONTAMINANT_PREFIXES)
    )
    decoy_prefix_flags = protein_values.map(
        lambda value: _has_any_prefixed_token(value, prefixes=_DECOY_PREFIXES)
    )

    flags = pd.DataFrame(index=source.index.copy())
    flags[_FRAGPIPE_CONTAMINANT_OUTPUT_COLUMN] = (
        contaminant_prefix_flags.astype(bool)
        if contaminant_flags is None
        else (contaminant_flags.astype(bool) | contaminant_prefix_flags.astype(bool))
    ).tolist()
    flags[_FRAGPIPE_DECOY_OUTPUT_COLUMN] = (
        decoy_prefix_flags.astype(bool)
        if decoy_flags is None
        else (decoy_flags.astype(bool) | decoy_prefix_flags.astype(bool))
    ).tolist()

    _raise_for_forbidden_flags(
        flags[_FRAGPIPE_CONTAMINANT_OUTPUT_COLUMN],
        policy=contaminant_policy,
        label="contaminant",
    )
    _raise_for_forbidden_flags(
        flags[_FRAGPIPE_DECOY_OUTPUT_COLUMN],
        policy=decoy_policy,
        label="decoy",
    )
    keep_mask = pd.Series(True, index=source.index.copy(), dtype=bool)
    if contaminant_policy == FRAGPIPE_FLAG_POLICY_REMOVE:
        keep_mask &= ~flags[_FRAGPIPE_CONTAMINANT_OUTPUT_COLUMN].astype(bool)
    if decoy_policy == FRAGPIPE_FLAG_POLICY_REMOVE:
        keep_mask &= ~flags[_FRAGPIPE_DECOY_OUTPUT_COLUMN].astype(bool)

    diagnostics = {
        "input_row_count": int(source.shape[0]),
        "contaminant_column": resolved.contaminant,
        "decoy_column": resolved.decoy,
        "contaminant_policy": contaminant_policy,
        "decoy_policy": decoy_policy,
        "contaminant_rows": int(
            flags[_FRAGPIPE_CONTAMINANT_OUTPUT_COLUMN].astype(bool).sum()
        ),
        "decoy_rows": int(flags[_FRAGPIPE_DECOY_OUTPUT_COLUMN].astype(bool).sum()),
        "removed_rows": int((~keep_mask).sum()),
        "retained_row_count": int(keep_mask.sum()),
        "contaminant_prefix_rows": int(contaminant_prefix_flags.astype(bool).sum()),
        "decoy_prefix_rows": int(decoy_prefix_flags.astype(bool).sum()),
    }
    warnings: list[str] = []
    if resolved.contaminant is None:
        warnings.append(
            "FragPipe contaminant column was not found; contaminant filtering used "
            "protein accession prefixes only"
        )
    if resolved.decoy is None:
        warnings.append(
            "FragPipe decoy column was not found; decoy filtering used protein "
            "accession prefixes only"
        )
    filtered = source.loc[keep_mask, :].copy(deep=True)
    filtered_flags = flags.loc[keep_mask, :].copy(deep=True)
    return filtered, filtered_flags, diagnostics, tuple(warnings)


def _raise_for_forbidden_flags(
    values: pd.Series,
    *,
    policy: str,
    label: str,
) -> None:
    raise_for_forbidden_flags(
        values,
        policy=policy,
        error_policy=FRAGPIPE_FLAG_POLICY_ERROR,
        importer_label="FragPipe",
        label=label,
    )


def _has_any_prefixed_token(value: object, *, prefixes: tuple[str, ...]) -> bool:
    for token in split_multi_value(value):
        upper = token.upper()
        if any(upper.startswith(prefix) for prefix in prefixes):
            return True
    return False


__all__ = ["_apply_flag_policies"]
