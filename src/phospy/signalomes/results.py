from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from ..internal.pandas_copy import detached_frame_copy, detached_series_copy
from .constants import MODULE_ID_COLUMN
from .maps import SignalomeMapData
from .networks import SignalomeNetworkData
from .serialization import serialize_site_assignments_for_export

if TYPE_CHECKING:
    from .clustering import SignalomeModuleSelectionDiagnostics


@dataclass(slots=True, init=False)
class ExpandedSignalome:
    """One kinase-of-interest view over the global signalome state."""

    kinase: str
    linked_kinases: tuple[str, ...]
    regulated_module_ids: tuple[int, ...]
    _expression_matrix_source: pd.DataFrame
    _site_assignments_source: pd.DataFrame
    _row_positions: np.ndarray
    _expression_matrix_cache: pd.DataFrame | None
    _site_assignments_cache: pd.DataFrame | None

    def __init__(
        self,
        kinase: str,
        linked_kinases: tuple[str, ...],
        regulated_module_ids: tuple[int, ...],
        expression_matrix: pd.DataFrame,
        site_assignments: pd.DataFrame,
        row_positions: Sequence[int] | np.ndarray | None = None,
    ) -> None:
        self.kinase = str(kinase)
        self.linked_kinases = tuple(
            str(linked_kinase) for linked_kinase in linked_kinases
        )
        self.regulated_module_ids = tuple(
            int(module_id) for module_id in regulated_module_ids
        )

        if expression_matrix.shape[0] != site_assignments.shape[0]:
            msg = (
                "expression_matrix and site_assignments must contain the same number "
                "of rows for expanded signalome views"
            )
            raise ValueError(msg)

        self._expression_matrix_source = expression_matrix
        self._site_assignments_source = site_assignments

        if row_positions is None:
            self._row_positions = np.arange(site_assignments.shape[0], dtype=int)
            self._expression_matrix_cache = expression_matrix
            self._site_assignments_cache = site_assignments
            return

        row_positions_array = np.asarray(row_positions, dtype=int).reshape(-1)
        if row_positions_array.size > 0:
            if row_positions_array.min() < 0:
                msg = "row_positions must be non-negative"
                raise ValueError(msg)
            max_position = int(row_positions_array.max())
            if max_position >= site_assignments.shape[0]:
                msg = "row_positions are out of bounds for site_assignments"
                raise ValueError(msg)
            if max_position >= expression_matrix.shape[0]:
                msg = "row_positions are out of bounds for expression_matrix"
                raise ValueError(msg)
        self._row_positions = row_positions_array.copy()
        self._expression_matrix_cache = None
        self._site_assignments_cache = None

    @property
    def expression_matrix(self) -> pd.DataFrame:
        """Return the expanded expression matrix, materialising lazily."""

        if self._expression_matrix_cache is None:
            self._expression_matrix_cache = self._expression_matrix_source.take(
                self._row_positions
            )
        return self._expression_matrix_cache

    @property
    def site_assignments(self) -> pd.DataFrame:
        """Return expanded site assignments, materialising lazily."""

        if self._site_assignments_cache is None:
            self._site_assignments_cache = self._site_assignments_source.take(
                self._row_positions
            )
        return self._site_assignments_cache


@dataclass(slots=True, init=False)
class SignalomeModules:
    """Module-centric wide and long signalome views.

    Public accessors return detached copies by default. Use ``copy=False``
    or ``to_owned_frames()`` only when shared mutable access is required.
    """

    _module_table: pd.DataFrame
    _kinase_module_relationships: pd.DataFrame

    def __init__(
        self,
        module_table: pd.DataFrame,
        kinase_module_relationships: pd.DataFrame,
    ) -> None:
        self._module_table = module_table
        self._kinase_module_relationships = kinase_module_relationships

    @property
    def module_table(self) -> pd.DataFrame:
        """Return a detached module-by-kinase table."""

        return self.to_frame(copy=True)

    @property
    def module_table_live(self) -> pd.DataFrame:
        """Return the owned mutable module-by-kinase table."""

        return self._module_table

    @property
    def kinase_module_relationships(self) -> pd.DataFrame:
        """Return a detached kinase-to-module relationship table."""

        return self.to_relationship_table(copy=True)

    @property
    def kinase_module_relationships_live(self) -> pd.DataFrame:
        """Return the owned mutable kinase-to-module relationship table."""

        return self._kinase_module_relationships

    def to_frame(self, *, copy: bool = True) -> pd.DataFrame:
        """Return the wide signalome module table."""

        if copy:
            return detached_frame_copy(self._module_table)
        return self._module_table

    def to_relationship_table(self, *, copy: bool = True) -> pd.DataFrame:
        """Return the long kinase-to-module relationship table."""

        if copy:
            return detached_frame_copy(self._kinase_module_relationships)
        return self._kinase_module_relationships

    def to_owned_frames(self) -> tuple[pd.DataFrame, pd.DataFrame]:
        """Return owned mutable module and relationship tables."""

        return self._module_table, self._kinase_module_relationships


@dataclass(slots=True, init=False)
class SignalomeAssignments:
    """Site- and protein-level signalome assignments.

    Public accessors return detached copies by default. Use ``copy=False``
    or ``to_owned_frames()`` only when shared mutable access is required.
    """

    _site_assignments: pd.DataFrame
    _protein_assignments: pd.DataFrame

    def __init__(
        self,
        site_assignments: pd.DataFrame,
        protein_assignments: pd.DataFrame,
    ) -> None:
        self._site_assignments = site_assignments
        self._protein_assignments = protein_assignments

    @property
    def site_assignments(self) -> pd.DataFrame:
        """Return a detached site assignment table."""

        return self.sites(copy=True)

    @property
    def site_assignments_live(self) -> pd.DataFrame:
        """Return the owned mutable site assignment table."""

        return self._site_assignments

    @property
    def protein_assignments(self) -> pd.DataFrame:
        """Return a detached protein assignment table."""

        return self.proteins(copy=True)

    @property
    def protein_assignments_live(self) -> pd.DataFrame:
        """Return the owned mutable protein assignment table."""

        return self._protein_assignments

    @property
    def protein_modules(self) -> pd.Series:
        """Return a detached protein-to-module assignment series."""

        return detached_series_copy(self._protein_assignments.loc[:, MODULE_ID_COLUMN])

    @property
    def protein_modules_live(self) -> pd.Series:
        """Return the owned mutable protein-to-module assignment series."""

        return self._protein_assignments.loc[:, MODULE_ID_COLUMN]

    def sites(self, *, copy: bool = True) -> pd.DataFrame:
        """Return the site assignment table."""

        if copy:
            return detached_frame_copy(self._site_assignments)
        return self._site_assignments

    def proteins(self, *, copy: bool = True) -> pd.DataFrame:
        """Return the protein assignment table."""

        if copy:
            return detached_frame_copy(self._protein_assignments)
        return self._protein_assignments

    def to_owned_frames(self) -> tuple[pd.DataFrame, pd.DataFrame]:
        """Return owned mutable site and protein assignment tables."""

        return self._site_assignments, self._protein_assignments


@dataclass(slots=True, init=False)
class SignalomeKinaseNetwork:
    """Network-centric signalome outputs.

    Public accessors return detached copies by default. Use ``copy=False``
    or ``to_owned_frames()`` only when shared mutable access is required.
    """

    _adjacency_matrix: pd.DataFrame
    _correlation_matrix: pd.DataFrame
    _node_table: pd.DataFrame
    _edge_table: pd.DataFrame
    _neighbor_map: dict[str, tuple[str, ...]]

    def __init__(
        self,
        adjacency_matrix: pd.DataFrame,
        correlation_matrix: pd.DataFrame,
        node_table: pd.DataFrame,
        edge_table: pd.DataFrame,
        neighbor_map: dict[str, tuple[str, ...]],
    ) -> None:
        self._adjacency_matrix = adjacency_matrix
        self._correlation_matrix = correlation_matrix
        self._node_table = node_table
        self._edge_table = edge_table
        self._neighbor_map = neighbor_map

    @property
    def adjacency_matrix(self) -> pd.DataFrame:
        """Return a detached thresholded kinase adjacency matrix."""

        return self.adjacency(copy=True)

    @property
    def adjacency_matrix_live(self) -> pd.DataFrame:
        """Return the owned mutable thresholded kinase adjacency matrix."""

        return self._adjacency_matrix

    @property
    def correlation_matrix(self) -> pd.DataFrame:
        """Return a detached kinase correlation matrix."""

        return self.correlations(copy=True)

    @property
    def correlation_matrix_live(self) -> pd.DataFrame:
        """Return the owned mutable kinase correlation matrix."""

        return self.correlations(copy=False)

    @property
    def node_table(self) -> pd.DataFrame:
        """Return a detached kinase network node table."""

        return self.nodes(copy=True)

    @property
    def node_table_live(self) -> pd.DataFrame:
        """Return the owned mutable kinase network node table."""

        return self._node_table

    @property
    def edge_table(self) -> pd.DataFrame:
        """Return a detached kinase network edge table."""

        return self.edges(copy=True)

    @property
    def edge_table_live(self) -> pd.DataFrame:
        """Return the owned mutable kinase network edge table."""

        return self._edge_table

    @property
    def neighbor_map(self) -> dict[str, tuple[str, ...]]:
        """Return a detached kinase neighbor mapping."""

        return dict(self._neighbor_map)

    @property
    def neighbor_map_live(self) -> dict[str, tuple[str, ...]]:
        """Return the owned mutable kinase neighbor mapping."""

        return self._neighbor_map

    def adjacency(self, *, copy: bool = True) -> pd.DataFrame:
        """Return the thresholded kinase adjacency matrix."""

        if copy:
            return detached_frame_copy(self._adjacency_matrix)
        return self._adjacency_matrix

    def correlations(self, *, copy: bool = True) -> pd.DataFrame:
        """Return the raw kinase correlation matrix used to derive edge weights."""

        if copy:
            return detached_frame_copy(self._correlation_matrix)
        return self._correlation_matrix

    def nodes(self, *, copy: bool = True) -> pd.DataFrame:
        """Return the kinase network node table."""

        if copy:
            return detached_frame_copy(self._node_table)
        return self._node_table

    def edges(self, *, copy: bool = True) -> pd.DataFrame:
        """Return the kinase network edge table."""

        if copy:
            return detached_frame_copy(self._edge_table)
        return self._edge_table

    def to_owned_frames(self) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """Return owned mutable adjacency, node, and edge tables."""

        return self._adjacency_matrix, self._node_table, self._edge_table


@dataclass(slots=True, init=False)
class SignalomeResult:
    """Structured signalome outputs with stable access and export contracts.

    Public table properties return detached copies by default. Use
    ``to_owned_frames()`` or explicit ``*_live`` accessors when you intentionally
    need shared mutable state.
    """

    _scoring_matrix: pd.DataFrame
    _pred_mat: pd.DataFrame
    _expression_matrix: pd.DataFrame
    _modules: SignalomeModules
    _assignments: SignalomeAssignments
    _network: SignalomeKinaseNetwork
    _kinase_substrate_map: dict[str, tuple[str, ...]]
    _expanded_signalomes: dict[str, ExpandedSignalome]
    _module_selection_diagnostics: SignalomeModuleSelectionDiagnostics

    def __init__(
        self,
        scoring_matrix: pd.DataFrame,
        pred_mat: pd.DataFrame,
        expression_matrix: pd.DataFrame,
        modules: SignalomeModules,
        assignments: SignalomeAssignments,
        network: SignalomeKinaseNetwork,
        kinase_substrate_map: dict[str, tuple[str, ...]],
        expanded_signalomes: dict[str, ExpandedSignalome],
        module_selection_diagnostics: SignalomeModuleSelectionDiagnostics,
    ) -> None:
        self._scoring_matrix = scoring_matrix
        self._pred_mat = pred_mat
        self._expression_matrix = expression_matrix
        self._modules = modules
        self._assignments = assignments
        self._network = network
        self._kinase_substrate_map = kinase_substrate_map
        self._expanded_signalomes = expanded_signalomes
        self._module_selection_diagnostics = module_selection_diagnostics

    @property
    def scoring_matrix(self) -> pd.DataFrame:
        """Return a detached kinase scoring matrix."""

        return detached_frame_copy(self._scoring_matrix)

    @property
    def scoring_matrix_live(self) -> pd.DataFrame:
        """Return the owned mutable kinase scoring matrix."""

        return self._scoring_matrix

    @property
    def pred_mat(self) -> pd.DataFrame:
        """Return a detached prediction matrix."""

        return detached_frame_copy(self._pred_mat)

    @property
    def pred_mat_live(self) -> pd.DataFrame:
        """Return the owned mutable prediction matrix."""

        return self._pred_mat

    @property
    def expression_matrix(self) -> pd.DataFrame:
        """Return a detached phosphosite expression matrix."""

        return detached_frame_copy(self._expression_matrix)

    @property
    def expression_matrix_live(self) -> pd.DataFrame:
        """Return the owned mutable phosphosite expression matrix."""

        return self._expression_matrix

    @property
    def modules(self) -> SignalomeModules:
        """Return the module-centric signalome views."""

        return self._modules

    @property
    def assignments(self) -> SignalomeAssignments:
        """Return the site- and protein-level assignment views."""

        return self._assignments

    @property
    def network(self) -> SignalomeKinaseNetwork:
        """Return the network-centric signalome views."""

        return self._network

    @property
    def kinase_substrate_map(self) -> dict[str, tuple[str, ...]]:
        """Return a detached kinase-to-site mapping."""

        return dict(self._kinase_substrate_map)

    @property
    def kinase_substrate_map_live(self) -> dict[str, tuple[str, ...]]:
        """Return the owned mutable kinase-to-site mapping."""

        return self._kinase_substrate_map

    @property
    def expanded_signalomes(self) -> dict[str, ExpandedSignalome]:
        """Return detached kinase-of-interest signalome views."""

        return {
            kinase: ExpandedSignalome(
                kinase=expanded.kinase,
                linked_kinases=expanded.linked_kinases,
                regulated_module_ids=expanded.regulated_module_ids,
                expression_matrix=expanded._expression_matrix_source.take(
                    expanded._row_positions
                ),
                site_assignments=expanded._site_assignments_source.take(
                    expanded._row_positions
                ),
            )
            for kinase, expanded in self._expanded_signalomes.items()
        }

    @property
    def expanded_signalomes_live(self) -> dict[str, ExpandedSignalome]:
        """Return the owned mutable kinase-of-interest signalome views."""

        return self._expanded_signalomes

    @property
    def module_selection_diagnostics(self) -> SignalomeModuleSelectionDiagnostics:
        """Return module-selection diagnostics captured during clustering."""

        return self._module_selection_diagnostics

    @property
    def kinases_of_interest(self) -> tuple[str, ...]:
        """Return the kinases of interest included in ``expanded_signalomes``."""

        return tuple(self._expanded_signalomes)

    @property
    def signalome_modules(self) -> pd.DataFrame:
        """Return a detached module-by-kinase matrix."""

        return self.modules.to_frame(copy=True)

    @property
    def signalome_modules_live(self) -> pd.DataFrame:
        """Return the owned mutable module-by-kinase matrix."""

        return self.modules.to_frame(copy=False)

    @property
    def kinase_module_relationships(self) -> pd.DataFrame:
        """Return a detached kinase-to-module relationship table."""

        return self.modules.to_relationship_table(copy=True)

    @property
    def kinase_module_relationships_live(self) -> pd.DataFrame:
        """Return the owned mutable kinase-to-module relationship table."""

        return self.modules.to_relationship_table(copy=False)

    @property
    def site_assignments(self) -> pd.DataFrame:
        """Return a detached site assignment table."""

        return self.assignments.sites(copy=True)

    @property
    def site_assignments_live(self) -> pd.DataFrame:
        """Return the owned mutable site assignment table."""

        return self.assignments.sites(copy=False)

    @property
    def protein_assignments(self) -> pd.DataFrame:
        """Return a detached protein assignment table."""

        return self.assignments.proteins(copy=True)

    @property
    def protein_assignments_live(self) -> pd.DataFrame:
        """Return the owned mutable protein assignment table."""

        return self.assignments.proteins(copy=False)

    @property
    def protein_modules(self) -> pd.Series:
        """Return a detached protein-to-module assignment series."""

        return self.assignments.protein_modules

    @property
    def protein_modules_live(self) -> pd.Series:
        """Return the owned mutable protein-to-module assignment series."""

        return self.assignments.protein_modules_live

    @property
    def kinase_substrates(self) -> dict[str, tuple[str, ...]]:
        """Return kinase-to-site mappings derived from the relationship table."""

        return self.kinase_substrate_map

    @property
    def kinase_network(self) -> dict[str, tuple[str, ...]]:
        """Return the kinase neighbor mapping for compatibility."""

        return self.network.neighbor_map

    @property
    def kinase_network_live(self) -> dict[str, tuple[str, ...]]:
        """Return the owned mutable kinase neighbor mapping."""

        return self.network.neighbor_map_live

    @property
    def kinase_adjacency_matrix(self) -> pd.DataFrame:
        """Return a detached thresholded kinase adjacency matrix."""

        return self.network.adjacency(copy=True)

    @property
    def kinase_adjacency_matrix_live(self) -> pd.DataFrame:
        """Return the owned mutable thresholded kinase adjacency matrix."""

        return self.network.adjacency(copy=False)

    @property
    def kinase_correlation_matrix(self) -> pd.DataFrame:
        """Return a detached raw kinase correlation matrix."""

        return self.network.correlations(copy=True)

    @property
    def kinase_correlation_matrix_live(self) -> pd.DataFrame:
        """Return the owned mutable raw kinase correlation matrix."""

        return self.network.correlations(copy=False)

    @property
    def kinase_network_nodes(self) -> pd.DataFrame:
        """Return a detached kinase network node table."""

        return self.network.nodes(copy=True)

    @property
    def kinase_network_nodes_live(self) -> pd.DataFrame:
        """Return the owned mutable kinase network node table."""

        return self.network.nodes(copy=False)

    @property
    def kinase_network_edges(self) -> pd.DataFrame:
        """Return a detached kinase network edge table."""

        return self.network.edges(copy=True)

    @property
    def kinase_network_edges_live(self) -> pd.DataFrame:
        """Return the owned mutable kinase network edge table."""

        return self.network.edges(copy=False)

    def to_frames(
        self,
        *,
        include_inputs: bool = False,
        copy: bool = True,
    ) -> dict[str, pd.DataFrame]:
        """Return the canonical named signalome tables."""

        frames: dict[str, pd.DataFrame] = {
            "signalome_modules": self.modules.to_frame(copy=copy),
            "kinase_module_relationships": self.modules.to_relationship_table(
                copy=copy
            ),
            "site_assignments": self.assignments.sites(copy=copy),
            "protein_assignments": self.assignments.proteins(copy=copy),
            "kinase_network_nodes": self.network.nodes(copy=copy),
            "kinase_network_edges": self.network.edges(copy=copy),
            "kinase_adjacency_matrix": self.network.adjacency(copy=copy),
            "kinase_correlation_matrix": self.network.correlations(copy=copy),
        }
        if include_inputs:
            frames.update(
                {
                    "scoring_matrix": (
                        detached_frame_copy(self._scoring_matrix)
                        if copy
                        else self._scoring_matrix
                    ),
                    "pred_mat": detached_frame_copy(self._pred_mat)
                    if copy
                    else self._pred_mat,
                    "expression_matrix": (
                        detached_frame_copy(self._expression_matrix)
                        if copy
                        else self._expression_matrix
                    ),
                }
            )
        return frames

    def to_owned_frames(
        self, *, include_inputs: bool = False
    ) -> dict[str, pd.DataFrame]:
        """Return the owned mutable signalome tables for advanced workflows."""

        return self.to_frames(include_inputs=include_inputs, copy=False)

    def to_csv(self, directory: str | Path) -> dict[str, Path]:
        """Write the canonical signalome tables to CSV files."""

        output_dir = Path(directory)
        output_dir.mkdir(parents=True, exist_ok=True)
        written: dict[str, Path] = {}
        for name, frame in self.to_frames(copy=False).items():
            path = output_dir / f"{name}.csv"
            if name == "site_assignments":
                serialize_site_assignments_for_export(frame).to_csv(path)
            else:
                frame.to_csv(path)
            written[name] = path
        return written

    def to_map_data(self) -> SignalomeMapData:
        """Build serialisable map-ready plotting data from this result."""

        from .maps import build_signalome_map_data

        return build_signalome_map_data(self)

    def to_network_data(self) -> SignalomeNetworkData:
        """Build graph-ready plotting data from this result."""

        from .networks import build_signalome_network_data

        return build_signalome_network_data(self)
