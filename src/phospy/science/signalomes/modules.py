"""Score-derived signalome module table domain services."""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import numpy as np
import pandas as pd

from phospy.errors.workflows import WorkflowStageError
from phospy.science.configs import (
    SIGNALOME_ASSIGNMENT_POLICY_CUTOFF_BINARY,
    SIGNALOME_ASSIGNMENT_POLICY_WEIGHTED_TOP,
    SignalomeAssignmentPolicy,
)
from phospy.science.signalomes.assignments import _normalize_top_kinase_weights
from phospy.science.signalomes.constants import (
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
    """Build candidate kinase-supported module shares as percentages."""

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
        module_table = _build_cutoff_binary_module_table(
            module_assignments=module_assignments,
            kinase_substrates=kinase_substrates,
            module_index=module_index,
            kinase_index=kinase_index,
            protein_to_module=protein_to_module,
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

    weighted_rows: list[tuple[int, str, str, float]] = []
    protein_to_module_lookup = {
        str(protein_id): int(module_id)
        for protein_id, module_id in protein_to_module.items()
        if int(module_id) > 0
    }
    kinase_membership = set(str(kinase) for kinase in kinase_index.tolist())
    site_payload = module_assignments.loc[
        :, [PROTEIN_COLUMN, TOP_KINASE_WEIGHTS_COLUMN]
    ].copy()
    site_payload.index = pd.Index(site_payload.index.astype(str), name=SITE_ID_COLUMN)
    site_ids = site_payload.index.to_numpy(dtype=object, copy=False)
    protein_values = site_payload.loc[:, PROTEIN_COLUMN].to_numpy(
        dtype=object, copy=False
    )
    top_weight_values = site_payload.loc[:, TOP_KINASE_WEIGHTS_COLUMN].to_numpy(
        dtype=object,
        copy=False,
    )
    for site_id, protein_value, top_weight_value in zip(
        site_ids,
        protein_values,
        top_weight_values,
        strict=True,
    ):
        protein_id = str(protein_value)
        module_id = protein_to_module_lookup.get(protein_id)
        if module_id is None:
            continue
        for kinase, weight in _normalize_top_kinase_weights(
            top_weight_value,
            site_id=str(site_id),
        ):
            if kinase not in kinase_membership:
                continue
            weighted_rows.append((module_id, kinase, protein_id, float(weight)))

    if not weighted_rows:
        return pd.DataFrame(0.0, index=module_index.copy(), columns=kinase_index.copy())

    weighted_hits = pd.DataFrame.from_records(
        weighted_rows,
        columns=[
            MODULE_ID_COLUMN,
            KINASE_COLUMN,
            PROTEIN_COLUMN,
            SUPPORT_WEIGHT_COLUMN,
        ],
    ).astype(
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
        .unstack(KINASE_COLUMN, fill_value=0)
    )
    return module_hits.reindex(
        index=module_index.copy(),
        columns=kinase_index.copy(),
        fill_value=0.0,
    ).astype(float)


def _build_cutoff_binary_module_table(
    *,
    module_assignments: pd.DataFrame,
    kinase_substrates: Mapping[str, Sequence[str]],
    module_index: pd.Index,
    kinase_index: pd.Index,
    protein_to_module: pd.Series,
) -> pd.DataFrame:
    unique_kinases = tuple(
        dict.fromkeys(str(kinase) for kinase in kinase_index.tolist())
    )
    unique_kinase_index = pd.Index(unique_kinases, name=KINASE_COLUMN)
    if module_index.empty or not unique_kinases:
        return pd.DataFrame(
            0.0,
            index=module_index.copy(),
            columns=kinase_index.copy(),
        ).astype(float)

    module_positions = {
        int(module_id): int(position)
        for position, module_id in enumerate(
            module_index.to_numpy(dtype=np.int64, copy=False)
        )
    }
    site_to_protein = module_assignments.loc[:, PROTEIN_COLUMN].astype(str)
    site_to_protein.index = pd.Index(
        site_to_protein.index.astype(str),
        name=SITE_ID_COLUMN,
    )
    site_to_protein_lookup = site_to_protein.to_dict()
    protein_to_module_lookup = {
        str(protein_id): int(module_id)
        for protein_id, module_id in protein_to_module.items()
        if int(module_id) > 0
    }

    counts_matrix = np.zeros(
        (int(module_index.size), len(unique_kinases)),
        dtype=float,
    )
    for kinase_position, kinase in enumerate(unique_kinases):
        substrate_sites = kinase_substrates.get(kinase, ())
        if not substrate_sites:
            continue
        seen_proteins: set[str] = set()
        for site_id in substrate_sites:
            protein_id = site_to_protein_lookup.get(str(site_id))
            if protein_id is None:
                continue
            protein_key = str(protein_id)
            if protein_key in seen_proteins:
                continue
            seen_proteins.add(protein_key)
            module_id = protein_to_module_lookup.get(protein_key)
            if module_id is None:
                continue
            module_position = module_positions.get(module_id)
            if module_position is None:
                continue
            counts_matrix[module_position, kinase_position] += 1.0

    module_hits = pd.DataFrame(
        counts_matrix,
        index=module_index.copy(),
        columns=unique_kinase_index,
    )
    return module_hits.reindex(
        index=module_index.copy(),
        columns=kinase_index.copy(),
        fill_value=0.0,
    ).astype(float)


__all__ = ["build_signalome_module_table"]
