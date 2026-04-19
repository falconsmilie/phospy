"""Signalome module table domain services."""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import pandas as pd

from phospy.api.configs import (
    SIGNALOME_ASSIGNMENT_POLICY_CUTOFF_BINARY,
    SIGNALOME_ASSIGNMENT_POLICY_WEIGHTED_TOP,
    SignalomeAssignmentPolicy,
)
from phospy.errors.workflows import WorkflowStageError
from phospy.signalomes.assignments import _normalize_top_kinase_weights
from phospy.signalomes.constants import (
    KINASE_COLUMN,
    MODULE_ID_COLUMN,
    PROTEIN_COLUMN,
    SITE_ID_COLUMN,
    SUPPORT_WEIGHT_COLUMN,
    TOP_KINASE_WEIGHTS_COLUMN,
)


def build_signalome_module_table(
    *,
    module_assignments: pd.DataFrame,
    kinase_substrates: Mapping[str, Sequence[str]],
    kinase_order: Sequence[str],
    assignment_policy: SignalomeAssignmentPolicy = (
        SIGNALOME_ASSIGNMENT_POLICY_CUTOFF_BINARY
    ),
) -> pd.DataFrame:
    """Build module-by-kinase signalome table as percent shares per module."""

    module_index = pd.Index(
        sorted(
            {
                int(value)
                for value in module_assignments.loc[:, MODULE_ID_COLUMN]
                if int(value) > 0
            }
        ),
        name=MODULE_ID_COLUMN,
    )
    kinase_index = pd.Index(
        [str(kinase) for kinase in kinase_order], name=KINASE_COLUMN
    )
    module_table = pd.DataFrame(
        0.0, index=module_index.copy(), columns=kinase_index.copy()
    )

    protein_to_module = (
        module_assignments.loc[:, [PROTEIN_COLUMN, MODULE_ID_COLUMN]]
        .drop_duplicates(subset=[PROTEIN_COLUMN])
        .set_index(PROTEIN_COLUMN)
        .loc[:, MODULE_ID_COLUMN]
        .astype("int64")
    )
    protein_to_module = protein_to_module.loc[protein_to_module > 0]

    if assignment_policy == SIGNALOME_ASSIGNMENT_POLICY_CUTOFF_BINARY:
        site_to_protein = module_assignments.loc[:, PROTEIN_COLUMN].astype(str)
        site_to_protein.index = pd.Index(
            site_to_protein.index.astype(str), name=SITE_ID_COLUMN
        )
        unique_kinases = tuple(dict.fromkeys(kinase_index.tolist()))
        kinase_site_map = pd.DataFrame(
            {
                KINASE_COLUMN: unique_kinases,
                SITE_ID_COLUMN: [
                    tuple(
                        str(site_id)
                        for site_id in kinase_substrates.get(str(kinase), ())
                    )
                    for kinase in unique_kinases
                ],
            }
        )
        if not kinase_site_map.empty:
            kinase_site_map = kinase_site_map.loc[
                kinase_site_map.loc[:, SITE_ID_COLUMN].map(len) > 0
            ]
            if not kinase_site_map.empty:
                kinase_site_map = kinase_site_map.explode(
                    SITE_ID_COLUMN, ignore_index=True
                )
                kinase_site_map.loc[:, SITE_ID_COLUMN] = kinase_site_map.loc[
                    :, SITE_ID_COLUMN
                ].astype(str)
                kinase_site_map.loc[:, PROTEIN_COLUMN] = kinase_site_map.loc[
                    :, SITE_ID_COLUMN
                ].map(site_to_protein)
                kinase_site_map = kinase_site_map.dropna(subset=[PROTEIN_COLUMN])
                if not kinase_site_map.empty:
                    kinase_site_map.loc[:, PROTEIN_COLUMN] = kinase_site_map.loc[
                        :, PROTEIN_COLUMN
                    ].astype(str)
                    kinase_site_map = kinase_site_map.drop_duplicates(
                        subset=[KINASE_COLUMN, PROTEIN_COLUMN],
                        keep="first",
                    )
                    kinase_site_map.loc[:, MODULE_ID_COLUMN] = kinase_site_map.loc[
                        :, PROTEIN_COLUMN
                    ].map(protein_to_module)
                    kinase_site_map = kinase_site_map.dropna(subset=[MODULE_ID_COLUMN])
                if not kinase_site_map.empty:
                    kinase_site_map.loc[:, MODULE_ID_COLUMN] = kinase_site_map.loc[
                        :, MODULE_ID_COLUMN
                    ].astype("int64")
                    module_hits = (
                        kinase_site_map.groupby(
                            [MODULE_ID_COLUMN, KINASE_COLUMN], sort=False
                        )
                        .size()
                        .unstack(KINASE_COLUMN, fill_value=0)
                        .astype(float)
                    )
                    module_table = module_hits.reindex(
                        index=module_index,
                        columns=kinase_index,
                        fill_value=0.0,
                    )
    elif assignment_policy == SIGNALOME_ASSIGNMENT_POLICY_WEIGHTED_TOP:
        module_table = _build_weighted_top_module_table(
            module_assignments=module_assignments,
            module_index=module_index,
            kinase_index=kinase_index,
            protein_to_module=protein_to_module,
        )
    else:
        allowed = ", ".join(
            sorted(
                (
                    SIGNALOME_ASSIGNMENT_POLICY_CUTOFF_BINARY,
                    SIGNALOME_ASSIGNMENT_POLICY_WEIGHTED_TOP,
                )
            )
        )
        raise WorkflowStageError(
            f"unsupported assignment_policy '{assignment_policy}'; expected one of: "
            f"{allowed}"
        )

    row_totals = module_table.sum(axis=1)
    non_zero_rows = row_totals > 0.0
    if non_zero_rows.any():
        module_table.loc[non_zero_rows] = (
            module_table.loc[non_zero_rows].div(row_totals.loc[non_zero_rows], axis=0)
            * 100.0
        )
    return module_table.astype(float).round(3)


def _build_weighted_top_module_table(
    *,
    module_assignments: pd.DataFrame,
    module_index: pd.Index,
    kinase_index: pd.Index,
    protein_to_module: pd.Series,
) -> pd.DataFrame:
    if TOP_KINASE_WEIGHTS_COLUMN not in module_assignments.columns:
        raise WorkflowStageError(
            "module assignments are missing top_kinase_weights required for "
            "assignment_policy='weighted_top'"
        )

    weighted_rows: list[dict[str, object]] = []
    site_payload = module_assignments.loc[
        :, [PROTEIN_COLUMN, TOP_KINASE_WEIGHTS_COLUMN]
    ].copy()
    site_payload.index = pd.Index(site_payload.index.astype(str), name=SITE_ID_COLUMN)
    for site_id, row in site_payload.iterrows():
        protein_id = str(row[PROTEIN_COLUMN])
        if protein_id not in protein_to_module.index:
            continue
        module_id = int(protein_to_module.loc[protein_id])
        for kinase, weight in _normalize_top_kinase_weights(
            row[TOP_KINASE_WEIGHTS_COLUMN],
            site_id=site_id,
        ):
            if kinase not in kinase_index:
                continue
            weighted_rows.append(
                {
                    MODULE_ID_COLUMN: module_id,
                    KINASE_COLUMN: kinase,
                    PROTEIN_COLUMN: protein_id,
                    SUPPORT_WEIGHT_COLUMN: float(weight),
                }
            )

    if not weighted_rows:
        return pd.DataFrame(0.0, index=module_index.copy(), columns=kinase_index.copy())

    weighted_hits = pd.DataFrame.from_records(weighted_rows).astype(
        {
            MODULE_ID_COLUMN: "int64",
            KINASE_COLUMN: str,
            PROTEIN_COLUMN: str,
            SUPPORT_WEIGHT_COLUMN: float,
        }
    )
    protein_level_weights = (
        weighted_hits.groupby(
            [MODULE_ID_COLUMN, KINASE_COLUMN, PROTEIN_COLUMN],
            sort=False,
        )[SUPPORT_WEIGHT_COLUMN]
        .max()
        .astype(float)
        .reset_index()
    )
    module_hits = (
        protein_level_weights.groupby([MODULE_ID_COLUMN, KINASE_COLUMN], sort=False)[
            SUPPORT_WEIGHT_COLUMN
        ]
        .sum()
        .astype(float)
        .unstack(KINASE_COLUMN, fill_value=0.0)
    )
    return module_hits.reindex(
        index=module_index.copy(),
        columns=kinase_index.copy(),
        fill_value=0.0,
    ).astype(float)


__all__ = ["build_signalome_module_table"]
