"""Expanded signalome domain services."""

from __future__ import annotations

import json
from collections.abc import Mapping, Sequence

import numpy as np
import pandas as pd

from phospy.contracts.configs import (
    SIGNALOME_ASSIGNMENT_POLICY_CUTOFF_BINARY,
    SIGNALOME_ASSIGNMENT_POLICY_WEIGHTED_TOP,
    SignalomeAssignmentPolicy,
)
from phospy.errors.workflows import WorkflowStageError
from phospy.science.signalomes.assignments import _normalize_top_kinase_weights
from phospy.science.signalomes.constants import (
    DISPLAY_ID_COLUMN,
    EXPANDED_SIGNALOME_ASSIGNMENT_POLICY_COLUMN,
    EXPANDED_SIGNALOME_KINASE_COLUMN,
    EXPANDED_SIGNALOME_LINKED_KINASES_COLUMN,
    EXPANDED_SIGNALOME_REGULATED_MODULE_IDS_COLUMN,
    EXPANDED_SIGNALOME_ROW_KIND_COLUMN,
    EXPANDED_SIGNALOME_ROW_KIND_SITE,
    EXPANDED_SIGNALOME_ROW_KIND_SUMMARY,
    EXPANDED_SIGNALOME_SITE_ORDER_COLUMN,
    EXPANDED_SIGNALOME_SUPPORT_KINASES_COLUMN,
    GENE_SYMBOL_COLUMN,
    ISOFORM_ID_COLUMN,
    JSON_EMPTY_ARRAY,
    MIN_EXPANDED_MODULE_SHARE_PERCENT,
    MODULE_ID_COLUMN,
    PROTEIN_ACCESSION_COLUMN,
    PROTEIN_COLUMN,
    SITE_COLUMN,
    SITE_ID_COLUMN,
    SITE_KEY_COLUMN,
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
    """Build a flattened score-derived signalome table for supported kinases."""

    required_identity_columns = (
        SITE_KEY_COLUMN,
        DISPLAY_ID_COLUMN,
        GENE_SYMBOL_COLUMN,
        SITE_COLUMN,
        PROTEIN_COLUMN,
        PROTEIN_ACCESSION_COLUMN,
        ISOFORM_ID_COLUMN,
    )
    missing_identity_columns = [
        column
        for column in required_identity_columns
        if column not in module_assignments.columns
    ]
    if missing_identity_columns:
        joined = ", ".join(missing_identity_columns)
        raise WorkflowStageError(
            f"expanded signalome requires module assignment identity columns: {joined}"
        )
    site_key_values = (
        module_assignments.loc[:, SITE_KEY_COLUMN].fillna("").astype(str).str.strip()
    )
    if (site_key_values == "").any():
        raise WorkflowStageError(
            "expanded signalome requires non-empty module assignment site_key values"
        )
    if site_key_values.tolist() != module_assignments.index.astype(str).tolist():
        raise WorkflowStageError(
            "expanded signalome requires module assignment site_key values to "
            "match the assignment index"
        )
    display_id_values = (
        module_assignments.loc[:, DISPLAY_ID_COLUMN].fillna("").astype(str).str.strip()
    )
    if (display_id_values == "").any():
        raise WorkflowStageError(
            "expanded signalome requires non-empty module assignment display_id values"
        )
    site_index = pd.Index(site_key_values.tolist(), name=SITE_KEY_COLUMN)
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

    site_module_ids = module_id_values.to_numpy(dtype=np.int64, copy=False)
    site_proteins = protein_ids.to_numpy(dtype=object, copy=False)
    site_keys = site_key_values.to_numpy(dtype=object, copy=False)
    site_display_ids = display_id_values.to_numpy(dtype=object, copy=False)
    site_gene_symbols = (
        indexed_assignments.loc[:, GENE_SYMBOL_COLUMN]
        .fillna("")
        .astype(str)
        .str.strip()
        .to_numpy(dtype=object, copy=False)
    )
    site_tokens = (
        indexed_assignments.loc[:, SITE_COLUMN]
        .fillna("")
        .astype(str)
        .str.strip()
        .to_numpy(dtype=object, copy=False)
    )
    site_protein_accessions = (
        indexed_assignments.loc[:, PROTEIN_ACCESSION_COLUMN]
        .fillna("")
        .astype(str)
        .str.strip()
        .to_numpy(dtype=object, copy=False)
    )
    site_isoform_ids = (
        indexed_assignments.loc[:, ISOFORM_ID_COLUMN]
        .fillna("")
        .astype(str)
        .str.strip()
        .to_numpy(dtype=object, copy=False)
    )
    site_top_kinases = top_kinases.to_numpy(dtype=object, copy=False)
    site_top_scores = top_scores.to_numpy(dtype=float, copy=False)
    site_ids = site_display_ids
    module_site_positions = _build_module_site_positions(site_module_ids)
    supported_modules_by_kinase = _build_supported_modules_by_kinase(
        signalome_modules=signalome_modules,
        kinase_order=kinase_order,
    )
    site_size = int(site_index.size)

    row_columns = [
        EXPANDED_SIGNALOME_KINASE_COLUMN,
        EXPANDED_SIGNALOME_ROW_KIND_COLUMN,
        EXPANDED_SIGNALOME_ASSIGNMENT_POLICY_COLUMN,
        EXPANDED_SIGNALOME_LINKED_KINASES_COLUMN,
        EXPANDED_SIGNALOME_REGULATED_MODULE_IDS_COLUMN,
        SITE_KEY_COLUMN,
        DISPLAY_ID_COLUMN,
        SITE_ID_COLUMN,
        EXPANDED_SIGNALOME_SITE_ORDER_COLUMN,
        GENE_SYMBOL_COLUMN,
        SITE_COLUMN,
        PROTEIN_COLUMN,
        PROTEIN_ACCESSION_COLUMN,
        ISOFORM_ID_COLUMN,
        MODULE_ID_COLUMN,
        EXPANDED_SIGNALOME_SUPPORT_KINASES_COLUMN,
        SUPPORT_WEIGHT_COLUMN,
        TOP_KINASE_COLUMN,
        TOP_SCORE_COLUMN,
    ]
    expanded_rows: list[tuple[object, ...]] = []
    for focal_kinase in kinase_order:
        linked_kinases = tuple(
            dict.fromkeys((focal_kinase, *neighbor_map.get(focal_kinase, ())))
        )
        linked_support_arrays = tuple(
            (linked_kinase, support_array)
            for linked_kinase, support_array in (
                (linked_kinase, support_by_kinase.get(linked_kinase))
                for linked_kinase in linked_kinases
            )
            if support_array is not None
        )
        supported_module_ids = supported_modules_by_kinase[focal_kinase]

        linked_kinases_json = json.dumps(
            list(linked_kinases),
            separators=(",", ":"),
            ensure_ascii=True,
        )
        supported_module_ids_json = json.dumps(
            list(supported_module_ids),
            separators=(",", ":"),
            ensure_ascii=True,
        )

        candidate_positions = _collect_supported_site_positions(
            supported_module_ids=supported_module_ids,
            module_site_positions=module_site_positions,
            site_size=site_size,
        )
        if candidate_positions.size == 0 or not linked_support_arrays:
            expanded_rows.append(
                _summary_row(
                    focal_kinase=focal_kinase,
                    assignment_policy=assignment_policy,
                    linked_kinases_json=linked_kinases_json,
                    supported_module_ids_json=supported_module_ids_json,
                )
            )
            continue

        candidate_support_weights = np.zeros(candidate_positions.size, dtype=float)
        for _, support_array in linked_support_arrays:
            candidate_support_weights += support_array[candidate_positions]
        supported_mask = candidate_support_weights > 0.0
        if not supported_mask.any():
            expanded_rows.append(
                _summary_row(
                    focal_kinase=focal_kinase,
                    assignment_policy=assignment_policy,
                    linked_kinases_json=linked_kinases_json,
                    supported_module_ids_json=supported_module_ids_json,
                )
            )
            continue

        matched_positions = candidate_positions[supported_mask]
        matched_support_weights = candidate_support_weights[supported_mask]
        for position, support_weight in zip(
            matched_positions,
            matched_support_weights,
            strict=True,
        ):
            row_position = int(position)
            support_kinases = [
                linked_kinase
                for linked_kinase, support_array in linked_support_arrays
                if float(support_array[row_position]) > 0.0
            ]
            expanded_rows.append(
                (
                    focal_kinase,
                    EXPANDED_SIGNALOME_ROW_KIND_SITE,
                    assignment_policy,
                    linked_kinases_json,
                    supported_module_ids_json,
                    str(site_keys[row_position]),
                    str(site_display_ids[row_position]),
                    str(site_ids[row_position]),
                    row_position,
                    str(site_gene_symbols[row_position]),
                    str(site_tokens[row_position]),
                    str(site_proteins[row_position]),
                    str(site_protein_accessions[row_position]),
                    str(site_isoform_ids[row_position]),
                    int(site_module_ids[row_position]),
                    json.dumps(
                        support_kinases,
                        separators=(",", ":"),
                        ensure_ascii=True,
                    ),
                    float(support_weight),
                    str(site_top_kinases[row_position]),
                    float(site_top_scores[row_position]),
                )
            )

    expanded = pd.DataFrame.from_records(
        expanded_rows,
        columns=row_columns,
    )
    expanded = expanded.astype(
        {
            EXPANDED_SIGNALOME_KINASE_COLUMN: str,
            EXPANDED_SIGNALOME_ROW_KIND_COLUMN: str,
            EXPANDED_SIGNALOME_ASSIGNMENT_POLICY_COLUMN: str,
            EXPANDED_SIGNALOME_LINKED_KINASES_COLUMN: str,
            EXPANDED_SIGNALOME_REGULATED_MODULE_IDS_COLUMN: str,
            SITE_KEY_COLUMN: str,
            DISPLAY_ID_COLUMN: str,
            SITE_ID_COLUMN: str,
            EXPANDED_SIGNALOME_SITE_ORDER_COLUMN: "int64",
            GENE_SYMBOL_COLUMN: str,
            SITE_COLUMN: str,
            PROTEIN_COLUMN: str,
            PROTEIN_ACCESSION_COLUMN: str,
            ISOFORM_ID_COLUMN: str,
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
            SITE_KEY_COLUMN,
            SITE_ID_COLUMN,
        ],
        ascending=[True, True, True, True, True],
        kind="stable",
    ).reset_index(drop=True)


def _build_module_site_positions(site_module_ids: np.ndarray) -> dict[int, np.ndarray]:
    module_site_positions: dict[int, list[int]] = {}
    for position, module_id in enumerate(site_module_ids):
        module_site_positions.setdefault(int(module_id), []).append(int(position))
    return {
        int(module_id): np.asarray(site_positions, dtype=np.int64)
        for module_id, site_positions in module_site_positions.items()
    }


def _build_supported_modules_by_kinase(
    *,
    signalome_modules: pd.DataFrame,
    kinase_order: Sequence[str],
) -> dict[str, tuple[int, ...]]:
    module_values = signalome_modules.loc[:, kinase_order].astype(float)
    module_ids = module_values.index.to_numpy(dtype=np.int64, copy=False)
    module_matrix = module_values.to_numpy(dtype=float, copy=False)
    supported_modules: dict[str, tuple[int, ...]] = {}
    for column_position, kinase in enumerate(kinase_order):
        supported = (
            module_matrix[:, column_position] > MIN_EXPANDED_MODULE_SHARE_PERCENT
        )
        supported_modules[str(kinase)] = tuple(
            int(module_id) for module_id in module_ids[supported].tolist()
        )
    return supported_modules


def _collect_supported_site_positions(
    *,
    supported_module_ids: Sequence[int],
    module_site_positions: Mapping[int, np.ndarray],
    site_size: int,
) -> np.ndarray:
    if not supported_module_ids or site_size <= 0:
        return np.empty(0, dtype=np.int64)
    candidate_mask = np.zeros(site_size, dtype=bool)
    for module_id in supported_module_ids:
        site_positions = module_site_positions.get(int(module_id))
        if site_positions is None or site_positions.size == 0:
            continue
        candidate_mask[site_positions] = True
    return np.flatnonzero(candidate_mask).astype(np.int64, copy=False)


def _summary_row(
    *,
    focal_kinase: str,
    assignment_policy: SignalomeAssignmentPolicy,
    linked_kinases_json: str,
    supported_module_ids_json: str,
) -> tuple[object, ...]:
    return (
        focal_kinase,
        EXPANDED_SIGNALOME_ROW_KIND_SUMMARY,
        assignment_policy,
        linked_kinases_json,
        supported_module_ids_json,
        "",
        "",
        "",
        -1,
        "",
        "",
        "",
        "",
        "",
        0,
        JSON_EMPTY_ARRAY,
        0.0,
        "",
        np.nan,
    )


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
