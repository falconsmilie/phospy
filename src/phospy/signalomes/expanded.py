"""Expanded signalome domain services."""

from __future__ import annotations

import json
from collections.abc import Mapping, Sequence

import numpy as np
import pandas as pd

from phospy.api.configs import (
    SIGNALOME_ASSIGNMENT_POLICY_CUTOFF_BINARY,
    SIGNALOME_ASSIGNMENT_POLICY_WEIGHTED_TOP,
    SignalomeAssignmentPolicy,
)
from phospy.errors.workflows import WorkflowStageError
from phospy.signalomes.assignments import _normalize_top_kinase_weights
from phospy.signalomes.constants import (
    EXPANDED_SIGNALOME_ASSIGNMENT_POLICY_COLUMN,
    EXPANDED_SIGNALOME_KINASE_COLUMN,
    EXPANDED_SIGNALOME_LINKED_KINASES_COLUMN,
    EXPANDED_SIGNALOME_REGULATED_MODULE_IDS_COLUMN,
    EXPANDED_SIGNALOME_ROW_KIND_COLUMN,
    EXPANDED_SIGNALOME_SITE_ORDER_COLUMN,
    EXPANDED_SIGNALOME_SUPPORT_KINASES_COLUMN,
    JSON_EMPTY_ARRAY,
    MIN_EXPANDED_MODULE_SHARE_PERCENT,
    MODULE_ID_COLUMN,
    PROTEIN_COLUMN,
    SITE_ID_COLUMN,
    SOURCE_KINASE_COLUMN,
    SUPPORT_WEIGHT_COLUMN,
    TARGET_KINASE_COLUMN,
    TOP_KINASE_COLUMN,
    TOP_KINASE_WEIGHTS_COLUMN,
    TOP_SCORE_COLUMN,
)


def build_expanded_signalome_table(
    *,
    module_assignments: pd.DataFrame,
    signalome_modules: pd.DataFrame,
    kinase_network_edges: pd.DataFrame,
    kinase_substrates: Mapping[str, Sequence[str]],
    assignment_policy: SignalomeAssignmentPolicy = (
        SIGNALOME_ASSIGNMENT_POLICY_CUTOFF_BINARY
    ),
) -> pd.DataFrame:
    """Build a flattened expanded-signalome table for all supported kinases."""

    site_index = pd.Index(module_assignments.index.astype(str), name=SITE_ID_COLUMN)
    indexed_assignments = module_assignments.copy(deep=False)
    indexed_assignments.index = site_index

    module_id_values = indexed_assignments.loc[:, MODULE_ID_COLUMN].astype("int64")
    protein_ids = indexed_assignments.loc[:, PROTEIN_COLUMN].astype(str)
    top_kinases = (
        indexed_assignments.loc[:, TOP_KINASE_COLUMN].astype(str)
        if TOP_KINASE_COLUMN in indexed_assignments.columns
        else pd.Series("", index=site_index, dtype=str)
    )
    top_scores = (
        indexed_assignments.loc[:, TOP_SCORE_COLUMN].astype(float)
        if TOP_SCORE_COLUMN in indexed_assignments.columns
        else pd.Series(np.nan, index=site_index, dtype=float)
    )

    kinase_order = [str(kinase) for kinase in signalome_modules.columns.astype(str)]
    neighbor_map = _build_kinase_neighbor_map(
        kinase_network_edges=kinase_network_edges,
        kinase_order=kinase_order,
    )
    support_by_kinase = _build_site_support_by_kinase(
        module_assignments=indexed_assignments,
        site_index=site_index,
        kinase_substrates=kinase_substrates,
        assignment_policy=assignment_policy,
    )

    site_positions = np.arange(site_index.size, dtype=np.int64)
    site_module_ids = module_id_values.to_numpy(dtype=np.int64, copy=False)
    site_proteins = protein_ids.to_numpy(dtype=object, copy=False)
    site_top_kinases = top_kinases.to_numpy(dtype=object, copy=False)
    site_top_scores = top_scores.to_numpy(dtype=float, copy=False)
    site_ids = site_index.to_numpy(dtype=object, copy=False)

    expanded_rows: list[dict[str, object]] = []
    for focal_kinase in kinase_order:
        linked_kinases = tuple(
            dict.fromkeys((focal_kinase, *neighbor_map.get(focal_kinase, ())))
        )
        regulated_module_ids = tuple(
            int(module_id)
            for module_id, share in signalome_modules.loc[:, focal_kinase].items()
            if float(share) > MIN_EXPANDED_MODULE_SHARE_PERCENT
        )
        regulated_module_set = set(regulated_module_ids)

        linked_kinases_json = json.dumps(
            list(linked_kinases),
            separators=(",", ":"),
            ensure_ascii=True,
        )
        regulated_module_ids_json = json.dumps(
            list(regulated_module_ids),
            separators=(",", ":"),
            ensure_ascii=True,
        )

        matched_site_count = 0
        for position, site_id, module_id in zip(
            site_positions,
            site_ids,
            site_module_ids,
            strict=True,
        ):
            if int(module_id) not in regulated_module_set:
                continue
            support_kinases: list[str] = []
            support_weight = 0.0
            for linked_kinase in linked_kinases:
                kinase_support = support_by_kinase.get(linked_kinase)
                if kinase_support is None:
                    continue
                weight = float(kinase_support[int(position)])
                if weight <= 0.0:
                    continue
                support_kinases.append(linked_kinase)
                support_weight += weight
            if support_weight <= 0.0:
                continue
            matched_site_count += 1
            expanded_rows.append(
                {
                    EXPANDED_SIGNALOME_KINASE_COLUMN: focal_kinase,
                    EXPANDED_SIGNALOME_ROW_KIND_COLUMN: "site",
                    EXPANDED_SIGNALOME_ASSIGNMENT_POLICY_COLUMN: assignment_policy,
                    EXPANDED_SIGNALOME_LINKED_KINASES_COLUMN: linked_kinases_json,
                    EXPANDED_SIGNALOME_REGULATED_MODULE_IDS_COLUMN: regulated_module_ids_json,
                    SITE_ID_COLUMN: str(site_id),
                    EXPANDED_SIGNALOME_SITE_ORDER_COLUMN: int(position),
                    PROTEIN_COLUMN: str(site_proteins[int(position)]),
                    MODULE_ID_COLUMN: int(module_id),
                    EXPANDED_SIGNALOME_SUPPORT_KINASES_COLUMN: json.dumps(
                        support_kinases,
                        separators=(",", ":"),
                        ensure_ascii=True,
                    ),
                    SUPPORT_WEIGHT_COLUMN: float(support_weight),
                    TOP_KINASE_COLUMN: str(site_top_kinases[int(position)]),
                    TOP_SCORE_COLUMN: float(site_top_scores[int(position)]),
                }
            )
        if matched_site_count == 0:
            expanded_rows.append(
                {
                    EXPANDED_SIGNALOME_KINASE_COLUMN: focal_kinase,
                    EXPANDED_SIGNALOME_ROW_KIND_COLUMN: "summary",
                    EXPANDED_SIGNALOME_ASSIGNMENT_POLICY_COLUMN: assignment_policy,
                    EXPANDED_SIGNALOME_LINKED_KINASES_COLUMN: linked_kinases_json,
                    EXPANDED_SIGNALOME_REGULATED_MODULE_IDS_COLUMN: regulated_module_ids_json,
                    SITE_ID_COLUMN: "",
                    EXPANDED_SIGNALOME_SITE_ORDER_COLUMN: -1,
                    PROTEIN_COLUMN: "",
                    MODULE_ID_COLUMN: 0,
                    EXPANDED_SIGNALOME_SUPPORT_KINASES_COLUMN: JSON_EMPTY_ARRAY,
                    SUPPORT_WEIGHT_COLUMN: 0.0,
                    TOP_KINASE_COLUMN: "",
                    TOP_SCORE_COLUMN: np.nan,
                }
            )

    expanded = pd.DataFrame.from_records(
        expanded_rows,
        columns=[
            EXPANDED_SIGNALOME_KINASE_COLUMN,
            EXPANDED_SIGNALOME_ROW_KIND_COLUMN,
            EXPANDED_SIGNALOME_ASSIGNMENT_POLICY_COLUMN,
            EXPANDED_SIGNALOME_LINKED_KINASES_COLUMN,
            EXPANDED_SIGNALOME_REGULATED_MODULE_IDS_COLUMN,
            SITE_ID_COLUMN,
            EXPANDED_SIGNALOME_SITE_ORDER_COLUMN,
            PROTEIN_COLUMN,
            MODULE_ID_COLUMN,
            EXPANDED_SIGNALOME_SUPPORT_KINASES_COLUMN,
            SUPPORT_WEIGHT_COLUMN,
            TOP_KINASE_COLUMN,
            TOP_SCORE_COLUMN,
        ],
    )
    expanded = expanded.astype(
        {
            EXPANDED_SIGNALOME_KINASE_COLUMN: str,
            EXPANDED_SIGNALOME_ROW_KIND_COLUMN: str,
            EXPANDED_SIGNALOME_ASSIGNMENT_POLICY_COLUMN: str,
            EXPANDED_SIGNALOME_LINKED_KINASES_COLUMN: str,
            EXPANDED_SIGNALOME_REGULATED_MODULE_IDS_COLUMN: str,
            SITE_ID_COLUMN: str,
            EXPANDED_SIGNALOME_SITE_ORDER_COLUMN: "int64",
            PROTEIN_COLUMN: str,
            MODULE_ID_COLUMN: "int64",
            EXPANDED_SIGNALOME_SUPPORT_KINASES_COLUMN: str,
            SUPPORT_WEIGHT_COLUMN: float,
            TOP_KINASE_COLUMN: str,
            TOP_SCORE_COLUMN: float,
        }
    )
    return expanded.sort_values(
        [
            EXPANDED_SIGNALOME_KINASE_COLUMN,
            EXPANDED_SIGNALOME_ROW_KIND_COLUMN,
            EXPANDED_SIGNALOME_SITE_ORDER_COLUMN,
            SITE_ID_COLUMN,
        ],
        ascending=[True, True, True, True],
        kind="stable",
    ).reset_index(drop=True)


def _build_kinase_neighbor_map(
    *,
    kinase_network_edges: pd.DataFrame,
    kinase_order: Sequence[str],
) -> dict[str, tuple[str, ...]]:
    neighbor_sets: dict[str, set[str]] = {str(kinase): set() for kinase in kinase_order}
    if kinase_network_edges.empty:
        return {kinase: () for kinase in kinase_order}
    required = {SOURCE_KINASE_COLUMN, TARGET_KINASE_COLUMN}
    missing = sorted(
        required.difference(str(column) for column in kinase_network_edges.columns)
    )
    if missing:
        preview = ", ".join(missing)
        raise WorkflowStageError(
            "kinase network edges are missing required columns for expanded signalome: "
            f"{preview}"
        )
    for row in kinase_network_edges.loc[
        :, [SOURCE_KINASE_COLUMN, TARGET_KINASE_COLUMN]
    ].itertuples(index=False):
        source_kinase = str(row.source_kinase)
        target_kinase = str(row.target_kinase)
        neighbor_sets.setdefault(source_kinase, set()).add(target_kinase)
        neighbor_sets.setdefault(target_kinase, set()).add(source_kinase)
    return {
        kinase: tuple(sorted(neighbor_sets.get(kinase, set())))
        for kinase in kinase_order
    }


def _build_site_support_by_kinase(
    *,
    module_assignments: pd.DataFrame,
    site_index: pd.Index,
    kinase_substrates: Mapping[str, Sequence[str]],
    assignment_policy: SignalomeAssignmentPolicy,
) -> dict[str, np.ndarray]:
    site_size = int(site_index.size)
    site_positions = {
        str(site_id): int(position)
        for position, site_id in enumerate(
            site_index.to_numpy(dtype=object, copy=False)
        )
    }
    if assignment_policy == SIGNALOME_ASSIGNMENT_POLICY_CUTOFF_BINARY:
        support: dict[str, np.ndarray] = {}
        for kinase, substrates in kinase_substrates.items():
            weights = np.zeros(site_size, dtype=float)
            for site_id in substrates:
                position = site_positions.get(str(site_id))
                if position is None:
                    continue
                weights[position] = 1.0
            support[str(kinase)] = weights
        return support
    if assignment_policy == SIGNALOME_ASSIGNMENT_POLICY_WEIGHTED_TOP:
        return _build_weighted_top_site_support(
            module_assignments=module_assignments,
            site_index=site_index,
        )
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


def _build_weighted_top_site_support(
    *,
    module_assignments: pd.DataFrame,
    site_index: pd.Index,
) -> dict[str, np.ndarray]:
    if TOP_KINASE_WEIGHTS_COLUMN not in module_assignments.columns:
        raise WorkflowStageError(
            "module assignments are missing top_kinase_weights required for "
            "assignment_policy='weighted_top'"
        )
    support_by_kinase: dict[str, np.ndarray] = {}
    weight_values = module_assignments.loc[:, TOP_KINASE_WEIGHTS_COLUMN].to_numpy(
        copy=False
    )
    for row_position, value in enumerate(weight_values):
        site_id = str(site_index[int(row_position)])
        normalized_weights = _normalize_top_kinase_weights(value, site_id=site_id)
        for kinase, weight in normalized_weights:
            kinase_support = support_by_kinase.get(kinase)
            if kinase_support is None:
                kinase_support = np.zeros(int(site_index.size), dtype=float)
                support_by_kinase[kinase] = kinase_support
            kinase_support[int(row_position)] = float(weight)
    return support_by_kinase


__all__ = ["build_expanded_signalome_table"]
