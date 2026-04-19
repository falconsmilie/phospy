"""Minimal signalome science helpers for workflow execution."""

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

SITE_ID_COLUMN = "site_id"
KINASE_COLUMN = "kinase"
PROTEIN_COLUMN = "protein_id"
MODULE_ID_COLUMN = "module_id"
TOP_KINASE_COLUMN = "top_kinase"
TOP_SCORE_COLUMN = "top_score"
TOP_KINASE_CANDIDATES_COLUMN = "top_kinase_candidates"
TOP_KINASE_WEIGHTS_COLUMN = "top_kinase_weights"
TOP_KINASE_TIE_COUNT_COLUMN = "top_kinase_tie_count"
TOP_KINASE_IS_AMBIGUOUS_COLUMN = "top_kinase_is_ambiguous"
TOP_KINASE_SELECTION_POLICY_COLUMN = "top_kinase_selection_policy"
MODULE_TOP_KINASE_COLUMN = "module_top_kinase"
MODULE_TOP_KINASE_CANDIDATES_COLUMN = "module_top_kinase_candidates"
MODULE_TOP_KINASE_TIE_COUNT_COLUMN = "module_top_kinase_tie_count"
MODULE_TOP_KINASE_IS_AMBIGUOUS_COLUMN = "module_top_kinase_is_ambiguous"
MODULE_TOP_KINASE_SELECTION_POLICY_COLUMN = "module_top_kinase_selection_policy"
LEXICOGRAPHIC_TIE_BREAK_POLICY = "max_score_then_lexicographic_tiebreak"
NO_SUPPORT_SELECTION_POLICY = "no_support"
UNSUPPORTED_KINASE = "__UNSUPPORTED__"
DEGREE_COLUMN = "degree"
N_SUBSTRATES_COLUMN = "n_substrates"
SOURCE_KINASE_COLUMN = "source_kinase"
TARGET_KINASE_COLUMN = "target_kinase"
CORRELATION_COLUMN = "correlation"
SUPPORT_WEIGHT_COLUMN = "support_weight"

EXPANDED_SIGNALOME_KINASE_COLUMN = "kinase"
EXPANDED_SIGNALOME_LINKED_KINASES_COLUMN = "linked_kinases"
EXPANDED_SIGNALOME_REGULATED_MODULE_IDS_COLUMN = "regulated_module_ids"
EXPANDED_SIGNALOME_SITE_ORDER_COLUMN = "site_order"
EXPANDED_SIGNALOME_SUPPORT_KINASES_COLUMN = "support_kinases"
EXPANDED_SIGNALOME_ROW_KIND_COLUMN = "row_kind"
EXPANDED_SIGNALOME_ASSIGNMENT_POLICY_COLUMN = "assignment_policy"

_JSON_EMPTY_ARRAY = "[]"
_MIN_EXPANDED_MODULE_SHARE_PERCENT = 1.0


def build_module_assignments(
    *,
    prediction_matrix: pd.DataFrame,
    site_to_protein: pd.Series,
    protein_modules: pd.Series | None = None,
) -> pd.DataFrame:
    """Build site-level module assignments with explicit tie-handling metadata."""

    if prediction_matrix.shape[1] == 0:
        raise WorkflowStageError("prediction matrix must contain at least one kinase")
    site_index = _as_unique_string_index(prediction_matrix.index, context="pred_mat")
    resolved_site_to_protein = _resolve_site_to_protein(
        site_index=site_index,
        site_to_protein=site_to_protein,
    )

    sorted_kinase_columns = sorted(str(kinase) for kinase in prediction_matrix.columns)
    sorted_predictions = prediction_matrix.loc[:, sorted_kinase_columns].astype(float)
    top_scores = sorted_predictions.max(axis=1).astype(float)
    top_score_mask = sorted_predictions.eq(top_scores, axis=0).to_numpy(
        dtype=bool, copy=False
    )
    kinase_names = sorted_predictions.columns.to_numpy(dtype=object, copy=False)
    tie_counts = top_score_mask.sum(axis=1).astype("int64")

    top_kinase_candidates: list[tuple[str, ...]] = []
    top_kinase_weights: list[tuple[tuple[str, float], ...]] = []
    top_kinases: list[str] = []
    top_kinase_policies: list[str] = []
    for mask_row, tie_count in zip(top_score_mask, tie_counts, strict=True):
        candidates = tuple(str(kinase) for kinase in kinase_names[mask_row])
        top_kinase_candidates.append(candidates)
        if tie_count == 0:
            top_kinase_weights.append(())
            top_kinases.append(UNSUPPORTED_KINASE)
            top_kinase_policies.append(NO_SUPPORT_SELECTION_POLICY)
            continue
        weight = 1.0 / float(tie_count)
        top_kinase_weights.append(tuple((kinase, weight) for kinase in candidates))
        top_kinases.append(candidates[0])
        top_kinase_policies.append(LEXICOGRAPHIC_TIE_BREAK_POLICY)

    top_kinase_series = pd.Series(
        top_kinases,
        index=site_index.copy(),
        dtype=object,
        name=TOP_KINASE_COLUMN,
    )
    site_module_resolution = _resolve_site_module_resolution(
        top_kinases=top_kinase_series,
        top_kinase_weights=top_kinase_weights,
        site_to_protein=resolved_site_to_protein,
        protein_modules=protein_modules,
    )

    assignments = pd.DataFrame(
        {
            PROTEIN_COLUMN: resolved_site_to_protein.to_numpy(dtype=object, copy=False),
            MODULE_ID_COLUMN: site_module_resolution.loc[:, MODULE_ID_COLUMN].to_numpy(
                dtype=np.int64, copy=False
            ),
            TOP_KINASE_COLUMN: top_kinase_series.to_numpy(dtype=object, copy=False),
            TOP_SCORE_COLUMN: top_scores.to_numpy(dtype=float, copy=False),
            TOP_KINASE_CANDIDATES_COLUMN: top_kinase_candidates,
            TOP_KINASE_WEIGHTS_COLUMN: top_kinase_weights,
            TOP_KINASE_TIE_COUNT_COLUMN: tie_counts,
            TOP_KINASE_IS_AMBIGUOUS_COLUMN: tie_counts > 1,
            TOP_KINASE_SELECTION_POLICY_COLUMN: top_kinase_policies,
            MODULE_TOP_KINASE_COLUMN: site_module_resolution.loc[
                :, MODULE_TOP_KINASE_COLUMN
            ].to_numpy(dtype=object, copy=False),
            MODULE_TOP_KINASE_CANDIDATES_COLUMN: site_module_resolution.loc[
                :, MODULE_TOP_KINASE_CANDIDATES_COLUMN
            ].to_numpy(dtype=object, copy=False),
            MODULE_TOP_KINASE_TIE_COUNT_COLUMN: site_module_resolution.loc[
                :, MODULE_TOP_KINASE_TIE_COUNT_COLUMN
            ].to_numpy(dtype=np.int64, copy=False),
            MODULE_TOP_KINASE_IS_AMBIGUOUS_COLUMN: site_module_resolution.loc[
                :, MODULE_TOP_KINASE_IS_AMBIGUOUS_COLUMN
            ].to_numpy(dtype=bool, copy=False),
            MODULE_TOP_KINASE_SELECTION_POLICY_COLUMN: site_module_resolution.loc[
                :, MODULE_TOP_KINASE_SELECTION_POLICY_COLUMN
            ].to_numpy(dtype=object, copy=False),
        },
        index=site_index.copy(),
    )
    return assignments.astype(
        {
            PROTEIN_COLUMN: str,
            MODULE_ID_COLUMN: "int64",
            TOP_KINASE_COLUMN: str,
            TOP_SCORE_COLUMN: float,
            TOP_KINASE_TIE_COUNT_COLUMN: "int64",
            TOP_KINASE_IS_AMBIGUOUS_COLUMN: bool,
            TOP_KINASE_SELECTION_POLICY_COLUMN: str,
            MODULE_TOP_KINASE_COLUMN: str,
            MODULE_TOP_KINASE_TIE_COUNT_COLUMN: "int64",
            MODULE_TOP_KINASE_IS_AMBIGUOUS_COLUMN: bool,
            MODULE_TOP_KINASE_SELECTION_POLICY_COLUMN: str,
        }
    )


def select_kinase_substrates(
    *,
    prediction_matrix: pd.DataFrame,
    cutoff: float,
) -> dict[str, tuple[str, ...]]:
    """Select phosphosites supported per kinase above `cutoff`."""

    site_ids = prediction_matrix.index.astype(str).to_numpy(dtype=object, copy=False)
    kinase_names = prediction_matrix.columns.astype(str).to_numpy(
        dtype=object, copy=False
    )
    support_mask = prediction_matrix.to_numpy(dtype=float, copy=False) > float(cutoff)
    return {
        str(kinase): tuple(site_ids[support_mask[:, index]].tolist())
        for index, kinase in enumerate(kinase_names)
    }


def build_signalome_module_table(
    *,
    module_assignments: pd.DataFrame,
    kinase_substrates: Mapping[str, Sequence[str]],
    kinase_order: Sequence[str],
    assignment_policy: SignalomeAssignmentPolicy = (
        SIGNALOME_ASSIGNMENT_POLICY_CUTOFF_BINARY
    ),
) -> pd.DataFrame:
    """Build module-by-kinase signalome table as percent shares per module.

    `assignment_policy` controls support attribution:

    - `cutoff_binary`: binary kinase support from `kinase_substrates`.
    - `weighted_top`: fractional support propagated from
      `module_assignments.top_kinase_weights`.
    """

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


def _normalize_top_kinase_weights(
    value: object,
    *,
    site_id: str,
) -> tuple[tuple[str, float], ...]:
    if isinstance(value, dict):
        pairs = tuple((str(kinase), float(weight)) for kinase, weight in value.items())
    elif isinstance(value, (tuple, list)):
        pairs = _normalize_top_kinase_weight_pairs(value, site_id=site_id)
    elif value is None:
        pairs = ()
    else:
        raise WorkflowStageError(
            "top_kinase_weights entries must be dicts or (kinase, weight) sequences; "
            f"received {type(value).__name__} at site_id='{site_id}'"
        )
    if not pairs:
        return ()
    positive_pairs = tuple((kinase, weight) for kinase, weight in pairs if weight > 0.0)
    if not positive_pairs:
        return ()
    return positive_pairs


def _normalize_top_kinase_weight_pairs(
    values: tuple[object, ...] | list[object],
    *,
    site_id: str,
) -> tuple[tuple[str, float], ...]:
    normalized: list[tuple[str, float]] = []
    for value in values:
        if not isinstance(value, (tuple, list)) or len(value) != 2:
            raise WorkflowStageError(
                "top_kinase_weights entries must be (kinase, weight) pairs; "
                f"received invalid entry at site_id='{site_id}'"
            )
        kinase, weight = value
        try:
            normalized.append((str(kinase), float(weight)))
        except (TypeError, ValueError) as exc:
            raise WorkflowStageError(
                "top_kinase_weights entries must contain float-compatible weights at "
                f"site_id='{site_id}'"
            ) from exc
    return tuple(normalized)


def build_kinase_network(
    *,
    downstream_score_matrix: pd.DataFrame,
    kinase_order: Sequence[str],
    kinase_substrates: Mapping[str, Sequence[str]],
    threshold: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build deterministic edge and node tables for kinase network output.

    Missing score entries are consumed via pairwise-complete Pearson
    correlations after dropping rows with no finite score support.
    """

    kinase_index = pd.Index(
        [str(kinase) for kinase in kinase_order], name=KINASE_COLUMN
    )
    kinase_index = pd.Index(
        list(dict.fromkeys(kinase_index.tolist())),
        name=KINASE_COLUMN,
    )
    if kinase_index.empty:
        raise WorkflowStageError("kinase network requires at least one kinase")
    available_kinases = set(downstream_score_matrix.columns.astype(str).tolist())
    missing_kinases = [
        kinase for kinase in kinase_index if kinase not in available_kinases
    ]
    if missing_kinases:
        preview = ", ".join(missing_kinases[:3])
        suffix = "..." if len(missing_kinases) > 3 else ""
        raise WorkflowStageError(
            "downstream score matrix is missing kinases required for signalome network: "
            f"{preview}{suffix}"
        )

    aligned_scores = _precondition_network_scores(
        downstream_score_matrix=downstream_score_matrix,
        kinase_index=kinase_index,
    )
    correlation_matrix = aligned_scores.corr(method="pearson", min_periods=2).fillna(
        0.0
    )
    correlation_matrix = correlation_matrix.loc[kinase_index, kinase_index]
    correlation_matrix.index = kinase_index.copy()
    correlation_matrix.columns = kinase_index.copy()

    correlation_values = correlation_matrix.to_numpy(dtype=float, copy=True)
    np.fill_diagonal(correlation_values, 0.0)

    source_positions, target_positions = np.triu_indices(len(kinase_index), k=1)
    pair_correlations = correlation_values[source_positions, target_positions]
    edge_mask = np.abs(pair_correlations) >= float(threshold)

    selected_source = source_positions[edge_mask]
    selected_target = target_positions[edge_mask]
    selected_correlations = pair_correlations[edge_mask]
    edges = pd.DataFrame(
        {
            SOURCE_KINASE_COLUMN: kinase_index.to_numpy(dtype=object, copy=False)[
                selected_source
            ],
            TARGET_KINASE_COLUMN: kinase_index.to_numpy(dtype=object, copy=False)[
                selected_target
            ],
            CORRELATION_COLUMN: selected_correlations.astype(float, copy=False),
        }
    )
    if edges.empty:
        edges = pd.DataFrame(
            columns=[SOURCE_KINASE_COLUMN, TARGET_KINASE_COLUMN, CORRELATION_COLUMN]
        )
    edges = edges.astype(
        {
            SOURCE_KINASE_COLUMN: str,
            TARGET_KINASE_COLUMN: str,
            CORRELATION_COLUMN: float,
        }
    )
    edges = edges.sort_values(
        [SOURCE_KINASE_COLUMN, TARGET_KINASE_COLUMN],
        ascending=[True, True],
        kind="stable",
    ).reset_index(drop=True)

    degree_values = np.zeros(len(kinase_index), dtype=np.int64)
    np.add.at(degree_values, selected_source, 1)
    np.add.at(degree_values, selected_target, 1)
    node_substrates = np.asarray(
        [
            len(tuple(kinase_substrates.get(str(kinase), ())))
            for kinase in kinase_index.to_numpy(dtype=object, copy=False)
        ],
        dtype=np.int64,
    )
    nodes = pd.DataFrame(
        {
            DEGREE_COLUMN: degree_values,
            N_SUBSTRATES_COLUMN: node_substrates,
        },
        index=kinase_index.copy(),
    )
    nodes.index.name = KINASE_COLUMN
    nodes = nodes.astype({DEGREE_COLUMN: "int64", N_SUBSTRATES_COLUMN: "int64"})
    return edges, nodes


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
    """Build a flattened expanded-signalome table for all supported kinases.

    One site-level row is emitted for every focal kinase/site membership where:
    1) the site's module is regulated by the focal kinase, and
    2) the site is supported by at least one linked kinase under the selected
       assignment policy.

    If a focal kinase has no selected site memberships, one `row_kind=summary`
    row is emitted to preserve linked-kinase and regulated-module metadata.
    """

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
            if float(share) > _MIN_EXPANDED_MODULE_SHARE_PERCENT
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
                    EXPANDED_SIGNALOME_SUPPORT_KINASES_COLUMN: _JSON_EMPTY_ARRAY,
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


def _as_unique_string_index(index: pd.Index, *, context: str) -> pd.Index:
    resolved = pd.Index(index.astype(str), name=SITE_ID_COLUMN)
    if not resolved.has_duplicates:
        return resolved
    duplicates = sorted(
        {str(site_id) for site_id in resolved[resolved.duplicated(keep=False)]}
    )
    preview = ", ".join(duplicates[:3])
    suffix = "..." if len(duplicates) > 3 else ""
    raise WorkflowStageError(
        f"{context} contains duplicate site identifiers: {preview}{suffix}"
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


def _precondition_network_scores(
    *,
    downstream_score_matrix: pd.DataFrame,
    kinase_index: pd.Index,
) -> pd.DataFrame:
    aligned_scores = downstream_score_matrix.loc[:, kinase_index].astype(float)
    score_values = aligned_scores.to_numpy(dtype=float, copy=False)
    infinite_mask = np.isinf(score_values)
    if infinite_mask.any():
        raise WorkflowStageError(
            "downstream score matrix contains infinite values after interpreter "
            "preconditioning"
        )
    supported_row_mask = (
        aligned_scores.notna().any(axis=1).to_numpy(dtype=bool, copy=False)
    )
    if supported_row_mask.all():
        return aligned_scores
    return aligned_scores.iloc[supported_row_mask, :]


def _resolve_site_to_protein(
    *,
    site_index: pd.Index,
    site_to_protein: pd.Series,
) -> pd.Series:
    resolved = site_to_protein.copy()
    resolved.index = pd.Index(resolved.index.astype(str), name=SITE_ID_COLUMN)
    missing = [site_id for site_id in site_index if site_id not in resolved.index]
    if missing:
        preview = ", ".join(missing[:3])
        suffix = "..." if len(missing) > 3 else ""
        raise WorkflowStageError(
            f"site-to-protein mapping is missing prediction sites: {preview}{suffix}"
        )
    aligned = resolved.loc[site_index].astype(str).str.strip()
    if (aligned == "").any():
        raise WorkflowStageError(
            "site-to-protein mapping contains empty protein identifiers"
        )
    aligned.index = site_index.copy()
    aligned.name = PROTEIN_COLUMN
    return aligned


def _resolve_site_module_resolution(
    *,
    top_kinases: pd.Series,
    top_kinase_weights: Sequence[tuple[tuple[str, float], ...]],
    site_to_protein: pd.Series,
    protein_modules: pd.Series | None,
) -> pd.DataFrame:
    if protein_modules is None:
        protein_resolution = _derive_protein_modules_by_top_kinase(
            top_kinases=top_kinases,
            top_kinase_weights=top_kinase_weights,
            site_to_protein=site_to_protein,
        )
    else:
        protein_resolution = _build_protein_resolution_from_modules(
            top_kinases=top_kinases,
            top_kinase_weights=top_kinase_weights,
            site_to_protein=site_to_protein,
            protein_modules=protein_modules,
        )
    site_proteins = site_to_protein.to_numpy(dtype=object, copy=False)
    site_module_resolution = protein_resolution.loc[
        site_proteins,
        [
            MODULE_ID_COLUMN,
            MODULE_TOP_KINASE_COLUMN,
            MODULE_TOP_KINASE_CANDIDATES_COLUMN,
            MODULE_TOP_KINASE_TIE_COUNT_COLUMN,
            MODULE_TOP_KINASE_IS_AMBIGUOUS_COLUMN,
            MODULE_TOP_KINASE_SELECTION_POLICY_COLUMN,
        ],
    ]
    site_module_resolution.index = site_to_protein.index.copy()
    return site_module_resolution


def _build_protein_resolution_from_modules(
    *,
    top_kinases: pd.Series,
    top_kinase_weights: Sequence[tuple[tuple[str, float], ...]],
    site_to_protein: pd.Series,
    protein_modules: pd.Series,
) -> pd.DataFrame:
    if len(top_kinases) != len(top_kinase_weights):
        raise WorkflowStageError(
            "top kinase weights must align one-to-one with prediction sites"
        )
    normalized_modules = _normalize_protein_modules(protein_modules)
    protein_index = pd.Index(
        sorted(set(site_to_protein.astype(str).tolist())),
        name=PROTEIN_COLUMN,
        dtype=object,
    )
    module_ids = pd.Series(
        np.zeros(len(protein_index), dtype=np.int64),
        index=protein_index.copy(),
        dtype="int64",
    )
    shared_proteins = protein_index.intersection(normalized_modules.index)
    if not shared_proteins.empty:
        module_ids.loc[shared_proteins] = normalized_modules.loc[shared_proteins]

    top_table = pd.DataFrame(
        {
            PROTEIN_COLUMN: site_to_protein.to_numpy(dtype=object, copy=False),
            TOP_KINASE_COLUMN: top_kinases.to_numpy(dtype=object, copy=False),
            TOP_KINASE_WEIGHTS_COLUMN: list(top_kinase_weights),
        },
        index=site_to_protein.index.copy(),
    )
    top_table.loc[:, MODULE_ID_COLUMN] = (
        top_table.loc[:, PROTEIN_COLUMN]
        .astype(str)
        .map(module_ids)
        .fillna(0)
        .astype("int64")
    )
    module_resolution = _derive_module_top_kinase_resolution(top_table)
    protein_resolution = pd.DataFrame(
        {
            MODULE_ID_COLUMN: module_ids.to_numpy(dtype=np.int64, copy=False),
        },
        index=protein_index.copy(),
    )
    return protein_resolution.join(module_resolution, on=MODULE_ID_COLUMN).astype(
        {
            MODULE_ID_COLUMN: "int64",
            MODULE_TOP_KINASE_COLUMN: str,
            MODULE_TOP_KINASE_TIE_COUNT_COLUMN: "int64",
            MODULE_TOP_KINASE_IS_AMBIGUOUS_COLUMN: bool,
            MODULE_TOP_KINASE_SELECTION_POLICY_COLUMN: str,
        }
    )


def _normalize_protein_modules(protein_modules: pd.Series) -> pd.Series:
    resolved = protein_modules.copy()
    resolved.index = pd.Index(resolved.index.astype(str), name=PROTEIN_COLUMN)
    if resolved.index.has_duplicates:
        duplicated = sorted(
            {str(value) for value in resolved.index[resolved.index.duplicated()]}
        )
        preview = ", ".join(duplicated[:3])
        suffix = "..." if len(duplicated) > 3 else ""
        raise WorkflowStageError(
            f"protein_modules contains duplicate protein identifiers: {preview}{suffix}"
        )
    module_values = pd.to_numeric(resolved, errors="coerce")
    if module_values.isna().any():
        raise WorkflowStageError("protein_modules must contain integer module IDs")
    rounded = np.floor(module_values.to_numpy(dtype=float, copy=False))
    if not np.allclose(module_values.to_numpy(dtype=float, copy=False), rounded):
        raise WorkflowStageError("protein_modules must contain integer module IDs")
    integer_values = module_values.astype("int64")
    integer_values.loc[integer_values < 0] = 0
    return integer_values.astype("int64")


def _derive_protein_modules_by_top_kinase(
    *,
    top_kinases: pd.Series,
    top_kinase_weights: Sequence[tuple[tuple[str, float], ...]],
    site_to_protein: pd.Series,
) -> pd.DataFrame:
    if len(top_kinases) != len(top_kinase_weights):
        raise WorkflowStageError(
            "top kinase weights must align one-to-one with prediction sites"
        )
    top_table = pd.DataFrame(
        {
            PROTEIN_COLUMN: site_to_protein.to_numpy(dtype=object, copy=False),
            TOP_KINASE_COLUMN: top_kinases.to_numpy(dtype=object, copy=False),
            TOP_KINASE_WEIGHTS_COLUMN: list(top_kinase_weights),
        },
        index=site_to_protein.index.copy(),
    )

    protein_resolution_rows: list[dict[str, object]] = []
    for protein_id, group in top_table.groupby(PROTEIN_COLUMN, sort=True):
        supported_group = group.loc[
            group.loc[:, TOP_KINASE_COLUMN].astype(str) != UNSUPPORTED_KINASE
        ]
        counts = supported_group.loc[:, TOP_KINASE_COLUMN].astype(str).value_counts()
        if not counts.empty:
            max_count = int(counts.iloc[0])
            tied_kinases = tuple(
                sorted(kinase for kinase in counts[counts == max_count].index.to_list())
            )
            dominant_kinase = tied_kinases[0]
            selection_policy = LEXICOGRAPHIC_TIE_BREAK_POLICY
        else:
            tied_kinases = ()
            dominant_kinase = UNSUPPORTED_KINASE
            selection_policy = NO_SUPPORT_SELECTION_POLICY
        protein_resolution_rows.append(
            {
                PROTEIN_COLUMN: str(protein_id),
                MODULE_TOP_KINASE_COLUMN: dominant_kinase,
                MODULE_TOP_KINASE_CANDIDATES_COLUMN: tied_kinases,
                MODULE_TOP_KINASE_TIE_COUNT_COLUMN: len(tied_kinases),
                MODULE_TOP_KINASE_IS_AMBIGUOUS_COLUMN: len(tied_kinases) > 1,
                MODULE_TOP_KINASE_SELECTION_POLICY_COLUMN: selection_policy,
            }
        )

    protein_resolution = pd.DataFrame(protein_resolution_rows).set_index(PROTEIN_COLUMN)
    protein_resolution.index = pd.Index(
        protein_resolution.index.astype(str), name=PROTEIN_COLUMN
    )
    dominant_kinases = protein_resolution.loc[:, MODULE_TOP_KINASE_COLUMN].astype(str)
    supported_mask = dominant_kinases != UNSUPPORTED_KINASE
    module_by_kinase = {
        kinase: module_id
        for module_id, kinase in enumerate(
            sorted(set(dominant_kinases.loc[supported_mask].tolist())), start=1
        )
    }
    module_ids = pd.Series(
        np.zeros(len(protein_resolution), dtype=np.int64),
        index=protein_resolution.index.copy(),
        dtype="int64",
    )
    if supported_mask.any():
        module_ids.loc[supported_mask] = (
            dominant_kinases.loc[supported_mask].map(module_by_kinase).astype("int64")
        )
    protein_resolution.loc[:, MODULE_ID_COLUMN] = module_ids.to_numpy(
        dtype=np.int64, copy=False
    )
    return protein_resolution.astype(
        {
            MODULE_TOP_KINASE_COLUMN: str,
            MODULE_TOP_KINASE_TIE_COUNT_COLUMN: "int64",
            MODULE_TOP_KINASE_IS_AMBIGUOUS_COLUMN: bool,
            MODULE_TOP_KINASE_SELECTION_POLICY_COLUMN: str,
            MODULE_ID_COLUMN: "int64",
        }
    )


def _derive_module_top_kinase_resolution(top_table: pd.DataFrame) -> pd.DataFrame:
    module_resolution_rows: list[dict[str, object]] = []
    module_ids = sorted(
        {int(module_id) for module_id in top_table.loc[:, MODULE_ID_COLUMN].tolist()}
    )
    if 0 not in module_ids:
        module_ids = [0, *module_ids]
    for module_id in module_ids:
        module_group = top_table.loc[
            top_table.loc[:, MODULE_ID_COLUMN].astype("int64") == int(module_id)
        ]
        supported_group = module_group.loc[
            module_group.loc[:, TOP_KINASE_COLUMN].astype(str) != UNSUPPORTED_KINASE
        ]
        counts = supported_group.loc[:, TOP_KINASE_COLUMN].astype(str).value_counts()
        if not counts.empty and int(module_id) > 0:
            max_count = int(counts.iloc[0])
            tied_kinases = tuple(
                sorted(kinase for kinase in counts[counts == max_count].index.to_list())
            )
            dominant_kinase = tied_kinases[0]
            selection_policy = LEXICOGRAPHIC_TIE_BREAK_POLICY
        else:
            tied_kinases = ()
            dominant_kinase = UNSUPPORTED_KINASE
            selection_policy = NO_SUPPORT_SELECTION_POLICY
        module_resolution_rows.append(
            {
                MODULE_ID_COLUMN: int(module_id),
                MODULE_TOP_KINASE_COLUMN: dominant_kinase,
                MODULE_TOP_KINASE_CANDIDATES_COLUMN: tied_kinases,
                MODULE_TOP_KINASE_TIE_COUNT_COLUMN: len(tied_kinases),
                MODULE_TOP_KINASE_IS_AMBIGUOUS_COLUMN: len(tied_kinases) > 1,
                MODULE_TOP_KINASE_SELECTION_POLICY_COLUMN: selection_policy,
            }
        )
    module_resolution = pd.DataFrame(module_resolution_rows).set_index(MODULE_ID_COLUMN)
    module_resolution.index = pd.Index(
        module_resolution.index.astype("int64"),
        name=MODULE_ID_COLUMN,
    )
    return module_resolution
