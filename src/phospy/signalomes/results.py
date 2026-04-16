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
    from ..datasets.models import SiteToProteinResolutionDiagnostics
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

    def to_expression_matrix(self) -> pd.DataFrame:
        """Return a detached expanded expression matrix."""

        return detached_frame_copy(self.expression_matrix)

    def to_owned_expression_matrix(self) -> pd.DataFrame:
        """Return cheap shared owned expanded expression matrix state."""

        return self.expression_matrix

    def to_mutable_expression_matrix_unsafe(self) -> pd.DataFrame:
        """Return explicit mutable shared expanded expression matrix state."""

        return self.expression_matrix

    def to_site_assignments(self) -> pd.DataFrame:
        """Return detached expanded site assignments."""

        return detached_frame_copy(self.site_assignments)

    def to_owned_site_assignments(self) -> pd.DataFrame:
        """Return cheap shared owned expanded site assignment state."""

        return self.site_assignments

    def to_mutable_site_assignments_unsafe(self) -> pd.DataFrame:
        """Return explicit mutable shared expanded site assignment state."""

        return self.site_assignments


@dataclass(slots=True, init=False)
class SignalomeModules:
    """Module-centric wide and long signalome views.

    Public accessors return detached copies by default. Expert callers who
    intentionally need shared mutable access must opt in via
    ``to_mutable_tables_unsafe()``.
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

        return self.to_frame()

    @property
    def kinase_module_relationships(self) -> pd.DataFrame:
        """Return a detached kinase-to-module relationship table."""

        return self.to_relationship_table()

    def to_frame(self) -> pd.DataFrame:
        """Return the wide signalome module table."""

        return detached_frame_copy(self._module_table)

    def to_owned_frame(self) -> pd.DataFrame:
        """Return cheap shared owned module table state (no copy)."""

        return self._module_table

    def to_mutable_frame_unsafe(self) -> pd.DataFrame:
        """Return explicit mutable shared module table state."""

        return self._module_table

    def to_relationship_table(self) -> pd.DataFrame:
        """Return the long kinase-to-module relationship table."""

        return detached_frame_copy(self._kinase_module_relationships)

    def to_owned_relationship_table(self) -> pd.DataFrame:
        """Return cheap shared owned kinase-to-module relationship state."""

        return self._kinase_module_relationships

    def to_mutable_relationship_table_unsafe(self) -> pd.DataFrame:
        """Return explicit mutable shared kinase-to-module relationship state."""

        return self._kinase_module_relationships

    def to_mutable_tables_unsafe(self) -> tuple[pd.DataFrame, pd.DataFrame]:
        """Return owned mutable module and relationship tables.

        Warning: mutating these tables mutates the owning signalome result.
        """

        return (
            self.to_mutable_frame_unsafe(),
            self.to_mutable_relationship_table_unsafe(),
        )


@dataclass(slots=True, init=False)
class SignalomeAssignments:
    """Site- and protein-level signalome assignments.

    Public accessors return detached copies by default. Expert callers who
    intentionally need shared mutable access must opt in via
    ``to_mutable_tables_unsafe()``.
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

        return self.to_site_assignments()

    @property
    def protein_assignments(self) -> pd.DataFrame:
        """Return a detached protein assignment table."""

        return self.to_protein_assignments()

    @property
    def protein_modules(self) -> pd.Series:
        """Return a detached protein-to-module assignment series."""

        return self.to_protein_modules()

    def to_site_assignments(self) -> pd.DataFrame:
        """Return detached site assignments."""

        return detached_frame_copy(self._site_assignments)

    def to_owned_site_assignments(self) -> pd.DataFrame:
        """Return cheap shared owned site assignment state."""

        return self._site_assignments

    def to_mutable_site_assignments_unsafe(self) -> pd.DataFrame:
        """Return explicit mutable shared site assignment state."""

        return self._site_assignments

    def to_protein_assignments(self) -> pd.DataFrame:
        """Return detached protein assignments."""

        return detached_frame_copy(self._protein_assignments)

    def to_owned_protein_assignments(self) -> pd.DataFrame:
        """Return cheap shared owned protein assignment state."""

        return self._protein_assignments

    def to_mutable_protein_assignments_unsafe(self) -> pd.DataFrame:
        """Return explicit mutable shared protein assignment state."""

        return self._protein_assignments

    def to_protein_modules(self) -> pd.Series:
        """Return detached protein-to-module assignments."""

        return detached_series_copy(self._protein_assignments.loc[:, MODULE_ID_COLUMN])

    def to_owned_protein_modules(self) -> pd.Series:
        """Return cheap shared owned protein-to-module assignment state."""

        return self._protein_assignments.loc[:, MODULE_ID_COLUMN]

    def to_mutable_protein_modules_unsafe(self) -> pd.Series:
        """Return explicit mutable shared protein-to-module assignment state."""

        return self._protein_assignments.loc[:, MODULE_ID_COLUMN]

    def sites(self) -> pd.DataFrame:
        """Compatibility alias for ``to_site_assignments()``."""

        return self.to_site_assignments()

    def proteins(self) -> pd.DataFrame:
        """Compatibility alias for ``to_protein_assignments()``."""

        return self.to_protein_assignments()

    def to_mutable_tables_unsafe(self) -> tuple[pd.DataFrame, pd.DataFrame]:
        """Return owned mutable site and protein assignment tables.

        Warning: mutating these tables mutates the owning signalome result.
        """

        return (
            self.to_mutable_site_assignments_unsafe(),
            self.to_mutable_protein_assignments_unsafe(),
        )


@dataclass(slots=True, init=False)
class SignalomeKinaseNetwork:
    """Network-centric signalome outputs.

    Public accessors return detached copies by default. Expert callers who
    intentionally need shared mutable access must opt in via
    ``to_mutable_state_unsafe()``.
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

        return self.to_adjacency_matrix()

    @property
    def correlation_matrix(self) -> pd.DataFrame:
        """Return a detached kinase correlation matrix."""

        return self.to_correlation_matrix()

    @property
    def node_table(self) -> pd.DataFrame:
        """Return a detached kinase network node table."""

        return self.to_node_table()

    @property
    def edge_table(self) -> pd.DataFrame:
        """Return a detached kinase network edge table."""

        return self.to_edge_table()

    @property
    def neighbor_map(self) -> dict[str, tuple[str, ...]]:
        """Return a detached kinase neighbor mapping."""

        return self.to_neighbor_map()

    def to_adjacency_matrix(self) -> pd.DataFrame:
        """Return detached thresholded kinase adjacency matrix."""

        return detached_frame_copy(self._adjacency_matrix)

    def to_owned_adjacency_matrix(self) -> pd.DataFrame:
        """Return cheap shared owned adjacency matrix state."""

        return self._adjacency_matrix

    def to_mutable_adjacency_matrix_unsafe(self) -> pd.DataFrame:
        """Return explicit mutable shared adjacency matrix state."""

        return self._adjacency_matrix

    def to_correlation_matrix(self) -> pd.DataFrame:
        """Return detached raw kinase correlation matrix."""

        return detached_frame_copy(self._correlation_matrix)

    def to_owned_correlation_matrix(self) -> pd.DataFrame:
        """Return cheap shared owned correlation matrix state."""

        return self._correlation_matrix

    def to_mutable_correlation_matrix_unsafe(self) -> pd.DataFrame:
        """Return explicit mutable shared correlation matrix state."""

        return self._correlation_matrix

    def to_node_table(self) -> pd.DataFrame:
        """Return detached kinase network node table."""

        return detached_frame_copy(self._node_table)

    def to_owned_node_table(self) -> pd.DataFrame:
        """Return cheap shared owned kinase network node state."""

        return self._node_table

    def to_mutable_node_table_unsafe(self) -> pd.DataFrame:
        """Return explicit mutable shared kinase network node state."""

        return self._node_table

    def to_edge_table(self) -> pd.DataFrame:
        """Return detached kinase network edge table."""

        return detached_frame_copy(self._edge_table)

    def to_owned_edge_table(self) -> pd.DataFrame:
        """Return cheap shared owned kinase network edge state."""

        return self._edge_table

    def to_mutable_edge_table_unsafe(self) -> pd.DataFrame:
        """Return explicit mutable shared kinase network edge state."""

        return self._edge_table

    def to_neighbor_map(self) -> dict[str, tuple[str, ...]]:
        """Return detached kinase neighbor mapping."""

        return dict(self._neighbor_map)

    def to_owned_neighbor_map(self) -> dict[str, tuple[str, ...]]:
        """Return cheap shared owned kinase neighbor mapping state."""

        return self._neighbor_map

    def to_mutable_neighbor_map_unsafe(self) -> dict[str, tuple[str, ...]]:
        """Return explicit mutable shared kinase neighbor mapping state."""

        return self._neighbor_map

    def adjacency(self) -> pd.DataFrame:
        """Compatibility alias for ``to_adjacency_matrix()``."""

        return self.to_adjacency_matrix()

    def correlations(self) -> pd.DataFrame:
        """Compatibility alias for ``to_correlation_matrix()``."""

        return self.to_correlation_matrix()

    def nodes(self) -> pd.DataFrame:
        """Compatibility alias for ``to_node_table()``."""

        return self.to_node_table()

    def edges(self) -> pd.DataFrame:
        """Compatibility alias for ``to_edge_table()``."""

        return self.to_edge_table()

    def to_mutable_state_unsafe(
        self,
    ) -> tuple[
        pd.DataFrame,
        pd.DataFrame,
        pd.DataFrame,
        pd.DataFrame,
        dict[str, tuple[str, ...]],
    ]:
        """Return owned mutable network tables and neighbor map.

        Warning: mutating these objects mutates the owning signalome result.
        """

        return (
            self.to_mutable_adjacency_matrix_unsafe(),
            self.to_mutable_correlation_matrix_unsafe(),
            self.to_mutable_node_table_unsafe(),
            self.to_mutable_edge_table_unsafe(),
            self.to_mutable_neighbor_map_unsafe(),
        )


@dataclass(slots=True, init=False)
class SignalomeResult:
    """Structured signalome outputs with stable access and export contracts.

    Ownership convention:
    - ``to_*`` methods return detached safe copies
    - ``to_owned_*`` methods return cheap shared owned state
    - ``to_mutable_*_unsafe`` methods return explicit mutable shared state
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
    _site_to_protein_resolution_diagnostics: SiteToProteinResolutionDiagnostics | None

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
        site_to_protein_resolution_diagnostics: SiteToProteinResolutionDiagnostics
        | None = None,
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
        self._site_to_protein_resolution_diagnostics = (
            site_to_protein_resolution_diagnostics
        )

    @property
    def scoring_matrix(self) -> pd.DataFrame:
        """Return a detached kinase scoring matrix."""

        return self.to_scoring_matrix()

    @property
    def pred_mat(self) -> pd.DataFrame:
        """Return a detached prediction matrix."""

        return self.to_pred_mat()

    @property
    def expression_matrix(self) -> pd.DataFrame:
        """Return a detached phosphosite expression matrix."""

        return self.to_expression_matrix()

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

        return self.to_kinase_substrate_map()

    @property
    def expanded_signalomes(self) -> dict[str, ExpandedSignalome]:
        """Return detached kinase-of-interest signalome views."""

        return self.to_expanded_signalomes()

    def to_scoring_matrix(self) -> pd.DataFrame:
        """Return detached kinase scoring matrix."""

        return detached_frame_copy(self._scoring_matrix)

    def to_owned_scoring_matrix(self) -> pd.DataFrame:
        """Return cheap shared owned kinase scoring matrix state."""

        return self._scoring_matrix

    def to_mutable_scoring_matrix_unsafe(self) -> pd.DataFrame:
        """Return explicit mutable shared kinase scoring matrix state."""

        return self._scoring_matrix

    def to_pred_mat(self) -> pd.DataFrame:
        """Return detached signalome prediction matrix."""

        return detached_frame_copy(self._pred_mat)

    def to_owned_pred_mat(self) -> pd.DataFrame:
        """Return cheap shared owned signalome prediction matrix state."""

        return self._pred_mat

    def to_mutable_pred_mat_unsafe(self) -> pd.DataFrame:
        """Return explicit mutable shared signalome prediction matrix state."""

        return self._pred_mat

    def to_expression_matrix(self) -> pd.DataFrame:
        """Return detached phosphosite expression matrix."""

        return detached_frame_copy(self._expression_matrix)

    def to_owned_expression_matrix(self) -> pd.DataFrame:
        """Return cheap shared owned phosphosite expression matrix state."""

        return self._expression_matrix

    def to_mutable_expression_matrix_unsafe(self) -> pd.DataFrame:
        """Return explicit mutable shared phosphosite expression matrix state."""

        return self._expression_matrix

    def to_kinase_substrate_map(self) -> dict[str, tuple[str, ...]]:
        """Return detached kinase-to-site mappings."""

        return dict(self._kinase_substrate_map)

    def to_owned_kinase_substrate_map(self) -> dict[str, tuple[str, ...]]:
        """Return cheap shared owned kinase-to-site mapping state."""

        return self._kinase_substrate_map

    def to_mutable_kinase_substrate_map_unsafe(self) -> dict[str, tuple[str, ...]]:
        """Return explicit mutable shared kinase-to-site mapping state."""

        return self._kinase_substrate_map

    def to_expanded_signalomes(self) -> dict[str, ExpandedSignalome]:
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

    def to_owned_expanded_signalomes(self) -> dict[str, ExpandedSignalome]:
        """Return cheap shared owned kinase-of-interest signalome view state."""

        return self._expanded_signalomes

    def to_mutable_expanded_signalomes_unsafe(self) -> dict[str, ExpandedSignalome]:
        """Return explicit mutable shared kinase-of-interest signalome view state."""

        return self._expanded_signalomes

    @property
    def module_selection_diagnostics(self) -> SignalomeModuleSelectionDiagnostics:
        """Return module-selection diagnostics captured during clustering."""

        return self._module_selection_diagnostics

    @property
    def site_to_protein_resolution_diagnostics(
        self,
    ) -> SiteToProteinResolutionDiagnostics | None:
        """Return site-to-protein resolution diagnostics when available."""

        return self._site_to_protein_resolution_diagnostics

    def attach_site_to_protein_resolution_diagnostics(
        self,
        diagnostics: SiteToProteinResolutionDiagnostics | None,
    ) -> None:
        """Attach site-to-protein resolution diagnostics to this result."""

        self._site_to_protein_resolution_diagnostics = diagnostics

    @property
    def kinases_of_interest(self) -> tuple[str, ...]:
        """Return the kinases of interest included in ``expanded_signalomes``."""

        return tuple(self._expanded_signalomes)

    @property
    def signalome_modules(self) -> pd.DataFrame:
        """Return a detached module-by-kinase matrix."""

        return self.modules.to_frame()

    @property
    def kinase_module_relationships(self) -> pd.DataFrame:
        """Return a detached kinase-to-module relationship table."""

        return self.modules.to_relationship_table()

    @property
    def site_assignments(self) -> pd.DataFrame:
        """Return a detached site assignment table."""

        return self.assignments.to_site_assignments()

    @property
    def protein_assignments(self) -> pd.DataFrame:
        """Return a detached protein assignment table."""

        return self.assignments.to_protein_assignments()

    @property
    def protein_modules(self) -> pd.Series:
        """Return a detached protein-to-module assignment series."""

        return self.assignments.protein_modules

    @property
    def kinase_substrates(self) -> dict[str, tuple[str, ...]]:
        """Return kinase-to-site mappings derived from the relationship table."""

        return self.kinase_substrate_map

    @property
    def kinase_network(self) -> dict[str, tuple[str, ...]]:
        """Return the kinase neighbor mapping for compatibility."""

        return self.network.to_neighbor_map()

    @property
    def kinase_adjacency_matrix(self) -> pd.DataFrame:
        """Return a detached thresholded kinase adjacency matrix."""

        return self.network.to_adjacency_matrix()

    @property
    def kinase_correlation_matrix(self) -> pd.DataFrame:
        """Return a detached raw kinase correlation matrix."""

        return self.network.to_correlation_matrix()

    @property
    def kinase_network_nodes(self) -> pd.DataFrame:
        """Return a detached kinase network node table."""

        return self.network.to_node_table()

    @property
    def kinase_network_edges(self) -> pd.DataFrame:
        """Return a detached kinase network edge table."""

        return self.network.to_edge_table()

    def to_frames(
        self,
        *,
        include_inputs: bool = False,
    ) -> dict[str, pd.DataFrame]:
        """Return detached canonical named signalome tables."""

        frames: dict[str, pd.DataFrame] = {
            "signalome_modules": self.modules.to_frame(),
            "kinase_module_relationships": self.modules.to_relationship_table(),
            "site_assignments": self.assignments.to_site_assignments(),
            "protein_assignments": self.assignments.to_protein_assignments(),
            "kinase_network_nodes": self.network.to_node_table(),
            "kinase_network_edges": self.network.to_edge_table(),
            "kinase_adjacency_matrix": self.network.to_adjacency_matrix(),
            "kinase_correlation_matrix": self.network.to_correlation_matrix(),
        }
        if include_inputs:
            frames.update(
                {
                    "scoring_matrix": self.to_scoring_matrix(),
                    "pred_mat": self.to_pred_mat(),
                    "expression_matrix": self.to_expression_matrix(),
                }
            )
        return frames

    def to_owned_frames(
        self,
        *,
        include_inputs: bool = False,
    ) -> dict[str, pd.DataFrame]:
        """Return cheap shared owned canonical signalome tables."""

        frames: dict[str, pd.DataFrame] = {
            "signalome_modules": self.modules.to_owned_frame(),
            "kinase_module_relationships": self.modules.to_owned_relationship_table(),
            "site_assignments": self.assignments.to_owned_site_assignments(),
            "protein_assignments": self.assignments.to_owned_protein_assignments(),
            "kinase_network_nodes": self.network.to_owned_node_table(),
            "kinase_network_edges": self.network.to_owned_edge_table(),
            "kinase_adjacency_matrix": self.network.to_owned_adjacency_matrix(),
            "kinase_correlation_matrix": self.network.to_owned_correlation_matrix(),
        }
        if include_inputs:
            frames.update(
                {
                    "scoring_matrix": self.to_owned_scoring_matrix(),
                    "pred_mat": self.to_owned_pred_mat(),
                    "expression_matrix": self.to_owned_expression_matrix(),
                }
            )
        return frames

    def to_mutable_frames_unsafe(
        self, *, include_inputs: bool = False
    ) -> dict[str, pd.DataFrame]:
        """Return owned mutable signalome tables for advanced workflows.

        Warning: mutating these frames mutates this result's internal state and
        can invalidate assumptions in downstream code.
        """
        frames: dict[str, pd.DataFrame] = {
            "signalome_modules": self.modules.to_mutable_frame_unsafe(),
            "kinase_module_relationships": self.modules.to_mutable_relationship_table_unsafe(),
            "site_assignments": self.assignments.to_mutable_site_assignments_unsafe(),
            "protein_assignments": self.assignments.to_mutable_protein_assignments_unsafe(),
            "kinase_network_nodes": self.network.to_mutable_node_table_unsafe(),
            "kinase_network_edges": self.network.to_mutable_edge_table_unsafe(),
            "kinase_adjacency_matrix": self.network.to_mutable_adjacency_matrix_unsafe(),
            "kinase_correlation_matrix": self.network.to_mutable_correlation_matrix_unsafe(),
        }
        if include_inputs:
            frames.update(
                {
                    "scoring_matrix": self.to_mutable_scoring_matrix_unsafe(),
                    "pred_mat": self.to_mutable_pred_mat_unsafe(),
                    "expression_matrix": self.to_mutable_expression_matrix_unsafe(),
                }
            )
        return frames

    def expanded_signalomes_mutable_unsafe(self) -> dict[str, ExpandedSignalome]:
        """Return owned mutable kinase-of-interest signalome views.

        Warning: mutating returned objects mutates this result's internal state.
        """

        return self.to_mutable_expanded_signalomes_unsafe()

    def to_csv(self, directory: str | Path) -> dict[str, Path]:
        """Write the canonical signalome tables to CSV files."""

        output_dir = Path(directory)
        output_dir.mkdir(parents=True, exist_ok=True)
        written: dict[str, Path] = {}
        for name, frame in self.to_frames().items():
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
