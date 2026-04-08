from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass

import numpy as np
import pandas as pd
from sklearn.cluster import AgglomerativeClustering

from .validation.errors import InputCompatibilityError

__all__ = [
    "ExpandedSignalome",
    "SignalomeResult",
    "build_signalome_result",
]


@dataclass(frozen=True, slots=True)
class ExpandedSignalome:
    """Expanded signalome view for one kinase of interest.

    Each expanded signalome contains the subset of the expression matrix and the
    aligned site annotation rows that support one kinase-of-interest view.
    """

    kinase: str
    linked_kinases: tuple[str, ...]
    regulated_module_ids: tuple[int, ...]
    expression_matrix: pd.DataFrame
    site_annotations: pd.DataFrame


@dataclass(slots=True)
class SignalomeResult:
    """Structured signalome outputs derived from scoring and prediction tables.

    The result keeps the major downstream signalome seams explicit:

    ``signalome_modules``
        Module-by-kinase percentage table derived from predicted substrates.
    ``site_assignments``
        Site-level annotation table with protein, module, and top-kinase labels.
    ``protein_modules``
        Protein-to-module assignment series.
    ``kinase_substrates``
        Canonical kinase-to-site mapping after applying the signalome cutoff.
    ``expanded_signalomes``
        Kinase-of-interest specific views over the stable module assignments.
    """

    scoring_matrix: pd.DataFrame
    pred_mat: pd.DataFrame
    expression_matrix: pd.DataFrame
    signalome_modules: pd.DataFrame
    site_assignments: pd.DataFrame
    protein_modules: pd.Series
    kinase_substrates: dict[str, tuple[str, ...]]
    kinase_network: dict[str, tuple[str, ...]]
    kinase_correlation_matrix: pd.DataFrame
    expanded_signalomes: dict[str, ExpandedSignalome]

    @property
    def kinases_of_interest(self) -> tuple[str, ...]:
        """Return the kinases of interest included in ``expanded_signalomes``."""

        return tuple(self.expanded_signalomes)


@dataclass(frozen=True, slots=True)
class _SignalomePlan:
    scoring_matrix: pd.DataFrame
    pred_mat: pd.DataFrame
    expression_matrix: pd.DataFrame
    kinases_of_interest: tuple[str, ...]
    kinase_network_threshold: float
    signalome_cutoff: float
    module_count: int | None
    min_kinase_module_share_percent: float


class _SignalomeRunner:
    def execute(self, plan: _SignalomePlan) -> SignalomeResult:
        scoring_matrix = plan.scoring_matrix
        pred_mat = plan.pred_mat
        expression_matrix = plan.expression_matrix

        site_clusters = _cluster_sites(
            scoring_matrix=scoring_matrix,
            requested_module_count=plan.module_count,
        )
        protein_modules = _derive_protein_modules(
            site_clusters=site_clusters,
            site_ids=scoring_matrix.index,
        )
        site_assignments = _build_site_assignments(
            pred_mat=pred_mat,
            protein_modules=protein_modules,
        )
        kinase_substrates = _select_kinase_substrates(
            pred_mat=pred_mat,
            cutoff=plan.signalome_cutoff,
        )
        kinase_network, kinase_correlation_matrix = _build_kinase_network(
            scoring_matrix=scoring_matrix,
            threshold=plan.kinase_network_threshold,
        )
        signalome_modules = _build_signalome_module_table(
            site_assignments=site_assignments,
            kinase_substrates=kinase_substrates,
        )
        expanded_signalomes = _build_expanded_signalomes(
            kinases_of_interest=plan.kinases_of_interest,
            kinase_network=kinase_network,
            kinase_substrates=kinase_substrates,
            signalome_modules=signalome_modules,
            site_assignments=site_assignments,
            expression_matrix=expression_matrix,
            min_kinase_module_share_percent=plan.min_kinase_module_share_percent,
        )

        return SignalomeResult(
            scoring_matrix=scoring_matrix,
            pred_mat=pred_mat,
            expression_matrix=expression_matrix,
            signalome_modules=signalome_modules,
            site_assignments=site_assignments,
            protein_modules=protein_modules,
            kinase_substrates=kinase_substrates,
            kinase_network=kinase_network,
            kinase_correlation_matrix=kinase_correlation_matrix,
            expanded_signalomes=expanded_signalomes,
        )


def build_signalome_result(
    *,
    scoring_matrix: pd.DataFrame,
    pred_mat: pd.DataFrame,
    expression_matrix: pd.DataFrame,
    kinases_of_interest: Sequence[str],
    kinase_network_threshold: float = 0.9,
    signalome_cutoff: float = 0.5,
    module_count: int | None = None,
    min_kinase_module_share_percent: float = 1.0,
) -> SignalomeResult:
    """Build a structured signalome result from trusted aligned inputs."""

    plan = _SignalomePlan(
        scoring_matrix=scoring_matrix,
        pred_mat=pred_mat,
        expression_matrix=expression_matrix,
        kinases_of_interest=tuple(kinases_of_interest),
        kinase_network_threshold=kinase_network_threshold,
        signalome_cutoff=signalome_cutoff,
        module_count=module_count,
        min_kinase_module_share_percent=min_kinase_module_share_percent,
    )
    return _SignalomeRunner().execute(plan)


def _cluster_sites(
    *,
    scoring_matrix: pd.DataFrame,
    requested_module_count: int | None,
) -> pd.Series:
    n_sites = scoring_matrix.shape[0]
    if n_sites == 1:
        return pd.Series(
            [1], index=scoring_matrix.index, dtype=int, name="site_cluster"
        )

    module_count = (
        requested_module_count
        if requested_module_count is not None
        else _select_module_count(scoring_matrix)
    )
    module_count = max(1, min(module_count, n_sites))

    if module_count == 1:
        labels = np.ones(n_sites, dtype=int)
    else:
        labels = (
            AgglomerativeClustering(
                n_clusters=module_count,
                linkage="ward",
            )
            .fit_predict(scoring_matrix.to_numpy(dtype=float))
            .astype(int)
            + 1
        )

    return pd.Series(labels, index=scoring_matrix.index, dtype=int, name="site_cluster")


def _select_module_count(scoring_matrix: pd.DataFrame) -> int:
    n_sites = scoring_matrix.shape[0]
    if n_sites <= 1:
        return 1

    max_clusters = min(10, n_sites)
    if max_clusters < 2:
        return 1

    site_correlations = np.corrcoef(scoring_matrix.to_numpy(dtype=float))
    candidates = _score_cluster_candidates(
        scoring_matrix=scoring_matrix,
        site_correlations=site_correlations,
        threshold=0.5,
        cluster_range=range(2, max_clusters + 1),
    )
    if not candidates:
        candidates = _score_cluster_candidates(
            scoring_matrix=scoring_matrix,
            site_correlations=site_correlations,
            threshold=0.1,
            cluster_range=range(2, max_clusters + 1),
        )
    if not candidates:
        return 1

    return max(candidates.items(), key=lambda item: (item[1], -item[0]))[0]


def _score_cluster_candidates(
    *,
    scoring_matrix: pd.DataFrame,
    site_correlations: np.ndarray,
    threshold: float,
    cluster_range: Iterable[int],
) -> dict[int, float]:
    candidates: dict[int, float] = {}
    for cluster_count in cluster_range:
        labels = (
            AgglomerativeClustering(
                n_clusters=cluster_count,
                linkage="ward",
            )
            .fit_predict(scoring_matrix.to_numpy(dtype=float))
            .astype(int)
        )
        cluster_medians = [
            _cluster_median_correlation(site_correlations, labels, label)
            for label in np.unique(labels)
        ]
        if cluster_medians and all(median >= threshold for median in cluster_medians):
            candidates[cluster_count] = float(np.mean(cluster_medians))
    return candidates


def _cluster_median_correlation(
    site_correlations: np.ndarray,
    labels: np.ndarray,
    label: int,
) -> float:
    cluster_positions = np.flatnonzero(labels == label)
    if cluster_positions.size <= 1:
        return 0.0

    cluster_correlations = site_correlations[
        np.ix_(cluster_positions, cluster_positions)
    ]
    cluster_correlations = cluster_correlations.copy()
    np.fill_diagonal(cluster_correlations, np.nan)
    values = cluster_correlations[~np.isnan(cluster_correlations)]
    if values.size == 0:
        return 0.0
    return float(np.median(values))


def _derive_protein_modules(
    *,
    site_clusters: pd.Series,
    site_ids: pd.Index,
) -> pd.Series:
    site_labels = pd.Index(site_ids.astype(str), dtype=object)
    proteins = pd.Index(
        [_protein_id_from_site(site_id) for site_id in site_labels], dtype=object
    )
    membership = pd.crosstab(site_clusters, proteins)
    membership = (membership > 0).astype(int)

    assignments: dict[str, int] = {}
    pattern_to_module: dict[tuple[int, ...], int] = {}
    next_module_id = 1

    for protein in membership.columns:
        pattern = tuple(int(value) for value in membership.loc[:, protein].tolist())
        if pattern not in pattern_to_module:
            pattern_to_module[pattern] = next_module_id
            next_module_id += 1
        assignments[str(protein)] = pattern_to_module[pattern]

    protein_modules = pd.Series(assignments, dtype=int, name="module_id")
    protein_modules.index.name = "protein_id"
    return protein_modules


def _build_site_assignments(
    *,
    pred_mat: pd.DataFrame,
    protein_modules: pd.Series,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for site_id in pred_mat.index.astype(str):
        protein_id = _protein_id_from_site(site_id)
        scores = pred_mat.loc[site_id]
        top_score = float(scores.max())
        top_kinases = [str(kinase) for kinase in scores.index[scores == top_score]]
        rows.append(
            {
                "site_id": site_id,
                "protein_id": protein_id,
                "module_id": int(protein_modules.loc[protein_id]),
                "top_kinase": top_kinases[0],
                "top_score": top_score,
            }
        )

    site_assignments = pd.DataFrame.from_records(rows).set_index("site_id")
    site_assignments.index.name = "site_id"
    return site_assignments


def _select_kinase_substrates(
    *,
    pred_mat: pd.DataFrame,
    cutoff: float,
) -> dict[str, tuple[str, ...]]:
    selected: dict[str, tuple[str, ...]] = {}
    for kinase in pred_mat.columns.astype(str):
        mask = pred_mat.loc[:, kinase] > cutoff
        selected[kinase] = tuple(pred_mat.index[mask].astype(str).tolist())
    return selected


def _build_kinase_network(
    *,
    scoring_matrix: pd.DataFrame,
    threshold: float,
) -> tuple[dict[str, tuple[str, ...]], pd.DataFrame]:
    kinase_correlation_matrix = scoring_matrix.corr().fillna(0.0)
    if kinase_correlation_matrix.empty:
        msg = "scoring_matrix must contain at least one kinase column"
        raise InputCompatibilityError(msg)

    kinase_correlation_matrix = kinase_correlation_matrix.copy(deep=True)
    diagonal_positions = np.arange(len(kinase_correlation_matrix))
    kinase_correlation_matrix.iloc[diagonal_positions, diagonal_positions] = 0.0

    kinase_network: dict[str, tuple[str, ...]] = {}
    for kinase in kinase_correlation_matrix.columns.astype(str):
        correlated = kinase_correlation_matrix.index[
            kinase_correlation_matrix.loc[:, kinase] > threshold
        ].astype(str)
        kinase_network[kinase] = tuple(correlated.tolist())

    return kinase_network, kinase_correlation_matrix


def _build_signalome_module_table(
    *,
    site_assignments: pd.DataFrame,
    kinase_substrates: Mapping[str, Sequence[str]],
) -> pd.DataFrame:
    module_ids = sorted(
        {int(value) for value in site_assignments["module_id"].tolist()}
    )
    kinase_names = list(kinase_substrates)
    module_table = pd.DataFrame(
        0.0,
        index=pd.Index(module_ids, name="module_id"),
        columns=pd.Index(kinase_names, name="kinase"),
    )

    protein_to_module = (
        site_assignments.loc[:, ["protein_id", "module_id"]]
        .drop_duplicates(subset=["protein_id"])
        .set_index("protein_id")
        .loc[:, "module_id"]
    )
    site_to_protein = site_assignments.loc[:, "protein_id"]

    for kinase, substrates in kinase_substrates.items():
        proteins = (
            pd.Index(site_to_protein.loc[list(substrates)].drop_duplicates())
            if substrates
            else pd.Index([], dtype=object)
        )
        proteins = proteins.intersection(protein_to_module.index)
        if proteins.empty:
            continue

        counts = protein_to_module.loc[proteins].value_counts().sort_index()
        for module_id, count in counts.items():
            module_table.loc[int(module_id), kinase] = float(count)

    row_totals = module_table.sum(axis=1)
    non_zero = row_totals > 0
    if non_zero.any():
        module_table.loc[non_zero] = (
            module_table.loc[non_zero].div(
                row_totals.loc[non_zero],
                axis=0,
            )
            * 100.0
        )

    return module_table.round(3)


def _build_expanded_signalomes(
    *,
    kinases_of_interest: Sequence[str],
    kinase_network: Mapping[str, Sequence[str]],
    kinase_substrates: Mapping[str, Sequence[str]],
    signalome_modules: pd.DataFrame,
    site_assignments: pd.DataFrame,
    expression_matrix: pd.DataFrame,
    min_kinase_module_share_percent: float,
) -> dict[str, ExpandedSignalome]:
    expanded: dict[str, ExpandedSignalome] = {}
    available_sites = pd.Index(site_assignments.index.astype(str), dtype=object)

    for kinase in kinases_of_interest:
        linked_kinases = tuple(dict.fromkeys((kinase, *kinase_network.get(kinase, ()))))
        regulated_module_ids = tuple(
            int(module_id)
            for module_id in signalome_modules.index[
                signalome_modules.loc[:, kinase] > min_kinase_module_share_percent
            ].tolist()
        )
        substrate_site_ids = tuple(
            site_id
            for linked_kinase in linked_kinases
            for site_id in kinase_substrates.get(linked_kinase, ())
        )
        substrate_site_index = available_sites.intersection(
            pd.Index(substrate_site_ids, dtype=object)
        )
        site_mask = site_assignments.loc[substrate_site_index, "module_id"].isin(
            regulated_module_ids
        )
        selected_site_ids = site_mask.index[site_mask]

        expanded[kinase] = ExpandedSignalome(
            kinase=kinase,
            linked_kinases=linked_kinases,
            regulated_module_ids=regulated_module_ids,
            expression_matrix=expression_matrix.loc[selected_site_ids].copy(deep=True),
            site_annotations=site_assignments.loc[selected_site_ids].copy(deep=True),
        )

    return expanded


def _protein_id_from_site(site_id: str) -> str:
    return site_id.split(";", 1)[0] if ";" in site_id else site_id
