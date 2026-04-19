"""Minimal signalome science helpers for workflow execution."""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import numpy as np
import pandas as pd

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


def build_module_assignments(
    *,
    prediction_matrix: pd.DataFrame,
    site_to_protein: pd.Series,
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
    protein_module_resolution = _derive_protein_modules(
        top_kinases=top_kinase_series,
        top_kinase_weights=top_kinase_weights,
        site_to_protein=resolved_site_to_protein,
    )
    site_proteins = resolved_site_to_protein.to_numpy(dtype=object, copy=False)
    site_module_resolution = protein_module_resolution.loc[
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
    site_module_resolution.index = site_index.copy()

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
    site_to_protein = module_assignments.loc[:, PROTEIN_COLUMN].astype(str)
    site_to_protein.index = pd.Index(
        site_to_protein.index.astype(str), name=SITE_ID_COLUMN
    )

    for kinase in kinase_index:
        substrate_sites = pd.Index(
            [str(site_id) for site_id in kinase_substrates.get(str(kinase), ())],
            name=SITE_ID_COLUMN,
        )
        if substrate_sites.empty:
            continue
        substrate_proteins = (
            site_to_protein.reindex(substrate_sites).dropna().astype(str)
        )
        if substrate_proteins.empty:
            continue
        unique_proteins = pd.Index(sorted(set(substrate_proteins.tolist())))
        module_hits = (
            protein_to_module.reindex(unique_proteins).dropna().astype("int64")
        )
        if module_hits.empty:
            continue
        counts = module_hits.value_counts().astype(float)
        module_table.loc[counts.index.astype(int), kinase] = counts.to_numpy(
            dtype=float, copy=False
        )

    row_totals = module_table.sum(axis=1)
    non_zero_rows = row_totals > 0.0
    if non_zero_rows.any():
        module_table.loc[non_zero_rows] = (
            module_table.loc[non_zero_rows].div(row_totals.loc[non_zero_rows], axis=0)
            * 100.0
        )
    return module_table.astype(float).round(3)


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


def _derive_protein_modules(
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
