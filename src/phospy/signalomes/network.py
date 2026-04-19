"""Signalome kinase-network domain services."""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import numpy as np
import pandas as pd

from phospy.api.configs import (
    SIGNALOME_KINASE_NETWORK_POLICIES,
    SIGNALOME_KINASE_NETWORK_POLICY_ABSOLUTE_THRESHOLD,
    SIGNALOME_KINASE_NETWORK_POLICY_POSITIVE_ONLY,
    SIGNALOME_KINASE_NETWORK_POLICY_SIGNED,
    SignalomeKinaseNetworkPolicy,
)
from phospy.errors.workflows import WorkflowStageError
from phospy.signalomes.constants import (
    CORRELATION_COLUMN,
    DEGREE_COLUMN,
    KINASE_COLUMN,
    N_SUBSTRATES_COLUMN,
    SOURCE_KINASE_COLUMN,
    TARGET_KINASE_COLUMN,
)


def build_kinase_network(
    *,
    downstream_score_matrix: pd.DataFrame,
    kinase_order: Sequence[str],
    kinase_substrates: Mapping[str, Sequence[str]],
    threshold: float,
    network_policy: SignalomeKinaseNetworkPolicy = (
        SIGNALOME_KINASE_NETWORK_POLICY_SIGNED
    ),
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build deterministic edge and node tables for kinase network output."""

    kinase_index = pd.Index(
        [str(kinase) for kinase in kinase_order], name=KINASE_COLUMN
    )
    kinase_index = pd.Index(
        list(dict.fromkeys(kinase_index.tolist())), name=KINASE_COLUMN
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
    edge_correlations, edge_mask = _resolve_network_edges_by_policy(
        pair_correlations=pair_correlations,
        threshold=float(threshold),
        network_policy=network_policy,
    )

    selected_source = source_positions[edge_mask]
    selected_target = target_positions[edge_mask]
    selected_correlations = edge_correlations[edge_mask]
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


def _resolve_network_edges_by_policy(
    *,
    pair_correlations: np.ndarray,
    threshold: float,
    network_policy: SignalomeKinaseNetworkPolicy,
) -> tuple[np.ndarray, np.ndarray]:
    if network_policy == SIGNALOME_KINASE_NETWORK_POLICY_POSITIVE_ONLY:
        edge_values = pair_correlations
        edge_mask = pair_correlations >= float(threshold)
        return edge_values, edge_mask
    if network_policy == SIGNALOME_KINASE_NETWORK_POLICY_ABSOLUTE_THRESHOLD:
        edge_values = np.abs(pair_correlations)
        edge_mask = edge_values >= float(threshold)
        return edge_values, edge_mask
    if network_policy == SIGNALOME_KINASE_NETWORK_POLICY_SIGNED:
        edge_values = pair_correlations
        edge_mask = np.abs(pair_correlations) >= float(threshold)
        return edge_values, edge_mask
    allowed = ", ".join(sorted(SIGNALOME_KINASE_NETWORK_POLICIES))
    raise WorkflowStageError(
        f"unsupported network_policy '{network_policy}'; expected one of: {allowed}"
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


__all__ = ["build_kinase_network"]
