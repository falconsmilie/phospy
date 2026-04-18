"""Signalome domain models."""

from __future__ import annotations

from dataclasses import InitVar, dataclass

import pandas as pd

from phospy._frame_ownership import own_dataframe, own_optional_dataframe
from phospy.errors.validation import PhosPyValidationError


@dataclass(frozen=True, slots=True)
class SignalomeAssignments:
    """Signalome module assignment table."""

    table: pd.DataFrame
    _assume_owned: InitVar[bool] = False

    def __post_init__(self, _assume_owned: bool) -> None:
        object.__setattr__(
            self,
            "table",
            own_dataframe(
                self.table,
                field_name="signalome_result.module_assignments.table",
                error_type=PhosPyValidationError,
                assume_owned=_assume_owned,
            ),
        )

    @classmethod
    def _from_owned(cls, *, table: pd.DataFrame) -> SignalomeAssignments:
        return cls(table=table, _assume_owned=True)


@dataclass(frozen=True, slots=True)
class SignalomeModules:
    """Signalome module table."""

    table: pd.DataFrame
    _assume_owned: InitVar[bool] = False

    def __post_init__(self, _assume_owned: bool) -> None:
        object.__setattr__(
            self,
            "table",
            own_dataframe(
                self.table,
                field_name="signalome_result.signalome_modules.table",
                error_type=PhosPyValidationError,
                assume_owned=_assume_owned,
            ),
        )

    @classmethod
    def _from_owned(cls, *, table: pd.DataFrame) -> SignalomeModules:
        return cls(table=table, _assume_owned=True)


@dataclass(frozen=True, slots=True)
class KinaseNetwork:
    """Kinase network tables derived from signalome analysis."""

    edges: pd.DataFrame
    nodes: pd.DataFrame | None = None
    _assume_owned: InitVar[bool] = False

    def __post_init__(self, _assume_owned: bool) -> None:
        edges = own_dataframe(
            self.edges,
            field_name="signalome_result.kinase_network.edges",
            error_type=PhosPyValidationError,
            assume_owned=_assume_owned,
        )
        nodes = own_optional_dataframe(
            self.nodes,
            field_name="signalome_result.kinase_network.nodes",
            error_type=PhosPyValidationError,
            assume_owned=_assume_owned,
        )
        object.__setattr__(self, "edges", edges)
        object.__setattr__(self, "nodes", nodes)

    @classmethod
    def _from_owned(
        cls,
        *,
        edges: pd.DataFrame,
        nodes: pd.DataFrame | None = None,
    ) -> KinaseNetwork:
        return cls(edges=edges, nodes=nodes, _assume_owned=True)
