from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from sklearn.cluster import AgglomerativeClustering

from .validation.errors import InputCompatibilityError

if TYPE_CHECKING:
    from .signalome_maps import SignalomeMapData
    from .signalome_networks import SignalomeNetworkData

__all__ = [
    "ExpandedSignalome",
    "SignalomeAssignments",
    "SignalomeKinaseNetwork",
    "SignalomeModules",
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


@dataclass(frozen=True, slots=True)
class SignalomeModules:
    """Canonical module-centric signalome outputs.

    ``module_table`` is the wide module-by-kinase percentage matrix.
    ``kinase_module_relationships`` is the graph-friendly long table derived from
    the non-zero cells of ``module_table``.
    """

    module_table: pd.DataFrame
    kinase_module_relationships: pd.DataFrame

    def to_frame(self, *, copy: bool = True) -> pd.DataFrame:
        """Return the canonical wide signalome module table."""

        if copy:
            return self.module_table.copy(deep=True)
        return self.module_table

    def to_relationship_table(self, *, copy: bool = True) -> pd.DataFrame:
        """Return the canonical long kinase-to-module relationship table."""

        if copy:
            return self.kinase_module_relationships.copy(deep=True)
        return self.kinase_module_relationships


@dataclass(frozen=True, slots=True)
class SignalomeAssignments:
    """Canonical site- and protein-level signalome assignments."""

    site_assignments: pd.DataFrame
    protein_assignments: pd.DataFrame

    @property
    def protein_modules(self) -> pd.Series:
        """Return the protein-to-module assignment series for compatibility."""

        return self.protein_assignments.loc[:, "module_id"]

    def sites(self, *, copy: bool = True) -> pd.DataFrame:
        """Return the canonical site assignment table."""

        if copy:
            return self.site_assignments.copy(deep=True)
        return self.site_assignments

    def proteins(self, *, copy: bool = True) -> pd.DataFrame:
        """Return the canonical protein assignment table."""

        if copy:
            return self.protein_assignments.copy(deep=True)
        return self.protein_assignments


@dataclass(frozen=True, slots=True)
class SignalomeKinaseNetwork:
    """Canonical network-centric signalome outputs."""

    correlation_matrix: pd.DataFrame
    node_table: pd.DataFrame
    edge_table: pd.DataFrame
    neighbor_map: dict[str, tuple[str, ...]]

    def adjacency(self, *, copy: bool = True) -> pd.DataFrame:
        """Return the kinase correlation matrix used to derive network edges."""

        if copy:
            return self.correlation_matrix.copy(deep=True)
        return self.correlation_matrix

    def nodes(self, *, copy: bool = True) -> pd.DataFrame:
        """Return the canonical kinase network node table."""

        if copy:
            return self.node_table.copy(deep=True)
        return self.node_table

    def edges(self, *, copy: bool = True) -> pd.DataFrame:
        """Return the canonical kinase network edge table."""

        if copy:
            return self.edge_table.copy(deep=True)
        return self.edge_table


@dataclass(slots=True)
class SignalomeResult:
    """Structured signalome outputs with stable access and export contracts.

    Canonical access paths:

    ``modules``
        Wide module matrix plus long kinase-to-module relationships.
    ``assignments``
        Site-level and protein-level module assignments.
    ``network``
        Correlation matrix plus graph-friendly kinase node and edge tables.
    ``expanded_signalomes``
        Kinase-of-interest views derived from the canonical module assignments.

    The raw aligned scoring, prediction, and expression inputs remain available
    for read-oriented workflows, but they are not part of the default export set.
    """

    scoring_matrix: pd.DataFrame
    pred_mat: pd.DataFrame
    expression_matrix: pd.DataFrame
    modules: SignalomeModules
    assignments: SignalomeAssignments
    network: SignalomeKinaseNetwork
    kinase_substrate_map: dict[str, tuple[str, ...]]
    expanded_signalomes: dict[str, ExpandedSignalome]

    @property
    def kinases_of_interest(self) -> tuple[str, ...]:
        """Return the kinases of interest included in ``expanded_signalomes``."""

        return tuple(self.expanded_signalomes)

    @property
    def signalome_modules(self) -> pd.DataFrame:
        """Return the canonical module-by-kinase matrix."""

        return self.modules.module_table

    @property
    def kinase_module_relationships(self) -> pd.DataFrame:
        """Return the canonical long kinase-to-module relationship table."""

        return self.modules.kinase_module_relationships

    @property
    def site_assignments(self) -> pd.DataFrame:
        """Return the canonical site assignment table."""

        return self.assignments.site_assignments

    @property
    def protein_assignments(self) -> pd.DataFrame:
        """Return the canonical protein assignment table."""

        return self.assignments.protein_assignments

    @property
    def protein_modules(self) -> pd.Series:
        """Return the protein-to-module assignment series for compatibility."""

        return self.assignments.protein_modules

    @property
    def kinase_substrates(self) -> dict[str, tuple[str, ...]]:
        """Return kinase-to-site mappings derived from the relationship table."""

        return dict(self.kinase_substrate_map)

    @property
    def kinase_network(self) -> dict[str, tuple[str, ...]]:
        """Return the canonical kinase neighbor mapping for compatibility."""

        return dict(self.network.neighbor_map)

    @property
    def kinase_correlation_matrix(self) -> pd.DataFrame:
        """Return the canonical kinase correlation matrix."""

        return self.network.correlation_matrix

    @property
    def kinase_network_nodes(self) -> pd.DataFrame:
        """Return the canonical kinase network node table."""

        return self.network.node_table

    @property
    def kinase_network_edges(self) -> pd.DataFrame:
        """Return the canonical kinase network edge table."""

        return self.network.edge_table

    def to_map_data(self) -> SignalomeMapData:
        """Build serialisable map-ready plotting data from this result."""

        from .signalome_maps import build_signalome_map_data

        return build_signalome_map_data(self)

    def to_network_data(self) -> SignalomeNetworkData:
        """Build graph-friendly kinase-network data from this result."""

        from .signalome_networks import build_signalome_network_data

        return build_signalome_network_data(self)

    def to_frames(
        self,
        *,
        copy: bool = True,
        include_inputs: bool = False,
    ) -> dict[str, pd.DataFrame]:
        """Return the canonical signalome tables as named pandas objects.

        By default, this returns only the stable user-facing outputs. Pass
        ``include_inputs=True`` to also include the aligned input matrices that
        fed signalome construction.
        """

        frames = {
            "signalome_modules": self.modules.to_frame(copy=copy),
            "kinase_module_relationships": self.modules.to_relationship_table(
                copy=copy
            ),
            "site_assignments": self.assignments.sites(copy=copy),
            "protein_assignments": self.assignments.proteins(copy=copy),
            "kinase_network_nodes": self.network.nodes(copy=copy),
            "kinase_network_edges": self.network.edges(copy=copy),
            "kinase_correlation_matrix": self.network.adjacency(copy=copy),
        }
        if include_inputs:
            frames.update(
                {
                    "scoring_matrix": self.scoring_matrix.copy(deep=True)
                    if copy
                    else self.scoring_matrix,
                    "pred_mat": self.pred_mat.copy(deep=True)
                    if copy
                    else self.pred_mat,
                    "expression_matrix": self.expression_matrix.copy(deep=True)
                    if copy
                    else self.expression_matrix,
                }
            )
        return frames

    def to_csv(
        self,
        directory: str | Path,
        *,
        include_inputs: bool = False,
    ) -> dict[str, Path]:
        """Export the canonical signalome tables to a directory of CSV files."""

        target_dir = Path(directory)
        target_dir.mkdir(parents=True, exist_ok=True)

        written_paths: dict[str, Path] = {}
        for name, frame in self.to_frames(
            copy=True, include_inputs=include_inputs
        ).items():
            path = target_dir / f"{name}.csv"
            frame.to_csv(
                path,
                encoding="utf-8",
                float_format="%.17g",
                lineterminator="\n",
            )
            written_paths[name] = path
        return written_paths


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
        selected_kinase_substrates = _select_kinase_substrates(
            pred_mat=pred_mat,
            cutoff=plan.signalome_cutoff,
        )
        kinase_network, kinase_correlation_matrix = _build_kinase_network(
            scoring_matrix=scoring_matrix,
            threshold=plan.kinase_network_threshold,
        )
        signalome_modules = _build_signalome_module_table(
            site_assignments=site_assignments,
            kinase_substrates=selected_kinase_substrates,
        )
        protein_assignments = _build_protein_assignment_table(
            site_assignments=site_assignments,
        )
        kinase_module_relationships = _build_kinase_module_relationship_table(
            module_table=signalome_modules,
        )
        network = _build_kinase_network_view(
            kinase_network=kinase_network,
            kinase_correlation_matrix=kinase_correlation_matrix,
            kinase_substrates=selected_kinase_substrates,
        )
        expanded_signalomes = _build_expanded_signalomes(
            kinases_of_interest=plan.kinases_of_interest,
            kinase_network=kinase_network,
            kinase_substrates=selected_kinase_substrates,
            signalome_modules=signalome_modules,
            site_assignments=site_assignments,
            expression_matrix=expression_matrix,
            min_kinase_module_share_percent=plan.min_kinase_module_share_percent,
        )

        return SignalomeResult(
            scoring_matrix=scoring_matrix,
            pred_mat=pred_mat,
            expression_matrix=expression_matrix,
            modules=SignalomeModules(
                module_table=signalome_modules,
                kinase_module_relationships=kinase_module_relationships,
            ),
            assignments=SignalomeAssignments(
                site_assignments=site_assignments,
                protein_assignments=protein_assignments,
            ),
            network=network,
            kinase_substrate_map=selected_kinase_substrates,
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


def _build_protein_assignment_table(*, site_assignments: pd.DataFrame) -> pd.DataFrame:
    protein_assignments = (
        site_assignments.groupby("protein_id", sort=True)
        .agg(
            module_id=("module_id", "first"),
            site_count=("module_id", "size"),
        )
        .astype({"module_id": int, "site_count": int})
        .sort_index()
    )
    protein_assignments.index.name = "protein_id"
    return protein_assignments


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

    kinase_correlation_values = kinase_correlation_matrix.to_numpy(copy=True)
    np.fill_diagonal(kinase_correlation_values, 0.0)
    kinase_correlation_matrix = pd.DataFrame(
        kinase_correlation_values,
        index=kinase_correlation_matrix.index.copy(),
        columns=kinase_correlation_matrix.columns.copy(),
    )

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


def _build_kinase_module_relationship_table(
    *,
    module_table: pd.DataFrame,
) -> pd.DataFrame:
    try:
        relationships = (
            module_table.stack(future_stack=True).rename("share_percent").reset_index()
        )
    except TypeError:
        relationships = (
            module_table.stack(dropna=False).rename("share_percent").reset_index()
        )

    relationships = relationships.loc[relationships["share_percent"] > 0.0]
    relationships.columns = ["module_id", "kinase", "share_percent"]
    relationships = relationships.astype(
        {
            "module_id": int,
            "kinase": str,
            "share_percent": float,
        }
    )
    relationships = relationships.sort_values(
        ["module_id", "share_percent", "kinase"],
        ascending=[True, False, True],
        kind="stable",
    ).reset_index(drop=True)
    return relationships


def _build_kinase_network_view(
    *,
    kinase_network: Mapping[str, Sequence[str]],
    kinase_correlation_matrix: pd.DataFrame,
    kinase_substrates: Mapping[str, Sequence[str]],
) -> SignalomeKinaseNetwork:
    node_rows = [
        {
            "kinase": str(kinase),
            "degree": len(tuple(neighbors)),
            "n_substrates": len(tuple(kinase_substrates.get(str(kinase), ()))),
        }
        for kinase, neighbors in sorted(kinase_network.items())
    ]
    node_table = pd.DataFrame.from_records(node_rows).set_index("kinase")
    node_table = node_table.astype({"degree": int, "n_substrates": int})
    node_table.index.name = "kinase"

    edge_rows: list[dict[str, object]] = []
    columns = kinase_correlation_matrix.columns.astype(str)
    for left_position, source in enumerate(columns):
        for target in columns[left_position + 1 :]:
            target_name = str(target)
            if target_name not in set(kinase_network.get(source, ())):
                continue
            correlation = float(kinase_correlation_matrix.loc[source, target_name])
            edge_rows.append(
                {
                    "source_kinase": source,
                    "target_kinase": target_name,
                    "correlation": correlation,
                }
            )

    edge_table = pd.DataFrame.from_records(edge_rows)
    if edge_table.empty:
        edge_table = pd.DataFrame(
            columns=["source_kinase", "target_kinase", "correlation"]
        )
    edge_table = edge_table.astype(
        {
            "source_kinase": str,
            "target_kinase": str,
            "correlation": float,
        }
    )
    edge_table = edge_table.sort_values(
        ["source_kinase", "target_kinase"],
        ascending=[True, True],
        kind="stable",
    ).reset_index(drop=True)

    return SignalomeKinaseNetwork(
        correlation_matrix=kinase_correlation_matrix,
        node_table=node_table,
        edge_table=edge_table,
        neighbor_map={str(key): tuple(value) for key, value in kinase_network.items()},
    )


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
