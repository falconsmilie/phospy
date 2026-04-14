from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd

from .constants import MODULE_ID_COLUMN
from .maps import SignalomeMapData
from .networks import SignalomeNetworkData
from .serialization import serialize_site_assignments_for_export

if TYPE_CHECKING:
    from .clustering import SignalomeModuleSelectionDiagnostics


@dataclass(slots=True)
class ExpandedSignalome:
    """One kinase-of-interest view over the global signalome state."""

    kinase: str
    linked_kinases: tuple[str, ...]
    regulated_module_ids: tuple[int, ...]
    expression_matrix: pd.DataFrame
    site_assignments: pd.DataFrame


@dataclass(slots=True)
class SignalomeModules:
    """Module-centric wide and long signalome views.

    ``module_table`` is the traditional module-by-kinase percentage table.
    ``kinase_module_relationships`` is the graph-friendly long table derived from
    the non-zero cells of ``module_table``.

    The wrapped frames are mutable owned state. Accessors return the owned
    frames by default; pass ``copy=True`` for detached safe copies.
    """

    module_table: pd.DataFrame
    kinase_module_relationships: pd.DataFrame

    def to_frame(self, *, copy: bool = False) -> pd.DataFrame:
        """Return the wide signalome module table."""

        if copy:
            return self.module_table.copy(deep=True)
        return self.module_table

    def to_relationship_table(self, *, copy: bool = False) -> pd.DataFrame:
        """Return the long kinase-to-module relationship table."""

        if copy:
            return self.kinase_module_relationships.copy(deep=True)
        return self.kinase_module_relationships


@dataclass(slots=True)
class SignalomeAssignments:
    """Site- and protein-level signalome assignments.

    The wrapped frames are mutable owned state. Accessors return the owned
    frames by default; pass ``copy=True`` for detached safe copies.
    """

    site_assignments: pd.DataFrame
    protein_assignments: pd.DataFrame

    @property
    def protein_modules(self) -> pd.Series:
        """Return the protein-to-module assignment series for compatibility."""

        return self.protein_assignments.loc[:, MODULE_ID_COLUMN]

    def sites(self, *, copy: bool = False) -> pd.DataFrame:
        """Return the site assignment table."""

        if copy:
            return self.site_assignments.copy(deep=True)
        return self.site_assignments

    def proteins(self, *, copy: bool = False) -> pd.DataFrame:
        """Return the protein assignment table."""

        if copy:
            return self.protein_assignments.copy(deep=True)
        return self.protein_assignments


@dataclass(slots=True)
class SignalomeKinaseNetwork:
    """Network-centric signalome outputs.

    The wrapped frames are mutable owned state. Accessors return the owned
    frames by default; pass ``copy=True`` for detached safe copies.
    """

    correlation_matrix: pd.DataFrame
    node_table: pd.DataFrame
    edge_table: pd.DataFrame
    neighbor_map: dict[str, tuple[str, ...]]

    def adjacency(self, *, copy: bool = False) -> pd.DataFrame:
        """Return the kinase correlation matrix used to derive network edges."""

        if copy:
            return self.correlation_matrix.copy(deep=True)
        return self.correlation_matrix

    def nodes(self, *, copy: bool = False) -> pd.DataFrame:
        """Return the kinase network node table."""

        if copy:
            return self.node_table.copy(deep=True)
        return self.node_table

    def edges(self, *, copy: bool = False) -> pd.DataFrame:
        """Return the kinase network edge table."""

        if copy:
            return self.edge_table.copy(deep=True)
        return self.edge_table


@dataclass(slots=True)
class SignalomeResult:
    """Structured signalome outputs with stable access and export contracts.

    Access paths:

    ``modules``
        Wide module matrix plus long kinase-to-module relationships.
    ``assignments``
        Site-level and protein-level module assignments.
    ``network``
        Correlation matrix plus graph-friendly kinase node and edge tables.
    ``expanded_signalomes``
        Kinase-of-interest views derived from the module assignments.

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
    module_selection_diagnostics: SignalomeModuleSelectionDiagnostics

    @property
    def kinases_of_interest(self) -> tuple[str, ...]:
        """Return the kinases of interest included in ``expanded_signalomes``."""

        return tuple(self.expanded_signalomes)

    @property
    def signalome_modules(self) -> pd.DataFrame:
        """Return the module-by-kinase matrix."""

        return self.modules.module_table

    @property
    def kinase_module_relationships(self) -> pd.DataFrame:
        """Return the long kinase-to-module relationship table."""

        return self.modules.kinase_module_relationships

    @property
    def site_assignments(self) -> pd.DataFrame:
        """Return the site assignment table."""

        return self.assignments.site_assignments

    @property
    def protein_assignments(self) -> pd.DataFrame:
        """Return the protein assignment table."""

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
        """Return the kinase neighbor mapping for compatibility."""

        return dict(self.network.neighbor_map)

    @property
    def kinase_correlation_matrix(self) -> pd.DataFrame:
        """Return the kinase correlation matrix."""

        return self.network.correlation_matrix

    @property
    def kinase_network_nodes(self) -> pd.DataFrame:
        """Return the kinase network node table."""

        return self.network.node_table

    @property
    def kinase_network_edges(self) -> pd.DataFrame:
        """Return the kinase network edge table."""

        return self.network.edge_table

    def to_frames(
        self,
        *,
        include_inputs: bool = False,
        copy: bool = False,
    ) -> dict[str, pd.DataFrame]:
        """Return the canonical named signalome tables.

        Access is zero-copy by default for internal pipelines and plotting
        adapters. Pass ``copy=True`` when callers need detached frames.
        """

        frames: dict[str, pd.DataFrame] = {
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
                    "scoring_matrix": (
                        self.scoring_matrix.copy(deep=True)
                        if copy
                        else self.scoring_matrix
                    ),
                    "pred_mat": self.pred_mat.copy(deep=True)
                    if copy
                    else self.pred_mat,
                    "expression_matrix": (
                        self.expression_matrix.copy(deep=True)
                        if copy
                        else self.expression_matrix
                    ),
                }
            )
        return frames

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
