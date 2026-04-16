from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from .constants import (
    DEGREE_COLUMN,
    IS_KINASE_OF_INTEREST_COLUMN,
    KINASE_COLUMN,
    MODULE_ID_COLUMN,
    N_SUBSTRATES_COLUMN,
    PROTEIN_ID_COLUMN,
    SHARE_PERCENT_COLUMN,
    SITE_ID_COLUMN,
    TOP_KINASE_CANDIDATES_COLUMN,
    TOP_KINASE_IS_AMBIGUOUS_COLUMN,
    TOP_KINASE_TIE_COUNT_COLUMN,
    TOP_KINASE_WEIGHTS_COLUMN,
    TOP_SCORE_COLUMN,
    TOTAL_SHARE_PERCENT_COLUMN,
)
from .serialization import (
    serialize_top_kinase_candidates,
    serialize_top_kinase_weights,
)

if TYPE_CHECKING:
    from .results import SignalomeCoreResult

__all__ = [
    "SignalomeMapData",
    "build_signalome_map_data",
]

_SITE_POSITION_OFFSET_BASE = -0.3
_SITE_POSITION_OFFSET_SPAN = 0.6
_SITE_BASELINE_Y = -0.25
_SITE_ROW_Y_STEP = 0.15
_KINASE_COLLISION_AVOIDANCE_SPAN = 0.4
_KINASE_BASELINE_Y = 1.0


@dataclass(slots=True)
class SignalomeMapData:
    """Serialisable map-ready plotting data derived from a signalome result.

    This model does not render charts. It exposes deterministic coordinate and
    relationship tables that plotting or export layers can consume. Accessors
    return owned frames by default; pass ``copy=True`` for detached copies.
    """

    module_positions: pd.DataFrame
    site_positions: pd.DataFrame
    kinase_positions: pd.DataFrame
    kinase_module_links: pd.DataFrame

    def modules(self, *, copy: bool = False) -> pd.DataFrame:
        """Return the canonical module-position table."""

        if copy:
            return self.module_positions.copy(deep=True)
        return self.module_positions

    def sites(self, *, copy: bool = False) -> pd.DataFrame:
        """Return the canonical site-position table."""

        if copy:
            return self.site_positions.copy(deep=True)
        return self.site_positions

    def kinases(self, *, copy: bool = False) -> pd.DataFrame:
        """Return the canonical kinase-position table."""

        if copy:
            return self.kinase_positions.copy(deep=True)
        return self.kinase_positions

    def links(self, *, copy: bool = False) -> pd.DataFrame:
        """Return the canonical kinase-to-module link table."""

        if copy:
            return self.kinase_module_links.copy(deep=True)
        return self.kinase_module_links

    def to_frames(self, *, copy: bool = False) -> dict[str, pd.DataFrame]:
        """Return the named plotting tables that make up the map data."""

        return {
            "signalome_map_modules": self.modules(copy=copy),
            "signalome_map_sites": self.sites(copy=copy),
            "signalome_map_kinases": self.kinases(copy=copy),
            "signalome_map_links": self.links(copy=copy),
        }

    def to_csv(self, directory: str | Path) -> dict[str, Path]:
        """Write the canonical plotting tables to CSV files."""

        target_dir = Path(directory)
        target_dir.mkdir(parents=True, exist_ok=True)

        written_paths: dict[str, Path] = {}
        for name, frame in self.to_frames(copy=False).items():
            path = target_dir / f"{name}.csv"
            write_index = not isinstance(frame.index, pd.RangeIndex)
            frame.to_csv(
                path,
                encoding="utf-8",
                float_format="%.17g",
                lineterminator="\n",
                index=write_index,
            )
            written_paths[name] = path
        return written_paths


def build_signalome_map_data(
    signalome_result: SignalomeCoreResult,
) -> SignalomeMapData:
    """Build deterministic plotting data from a canonical signalome result."""

    module_positions = _build_module_positions(signalome_result)
    site_positions = _build_site_positions(signalome_result, module_positions)
    kinase_positions = _build_kinase_positions(signalome_result, module_positions)
    kinase_module_links = _build_kinase_module_links(
        signalome_result,
        module_positions,
        kinase_positions,
    )
    return SignalomeMapData(
        module_positions=module_positions,
        site_positions=site_positions,
        kinase_positions=kinase_positions,
        kinase_module_links=kinase_module_links,
    )


def _build_module_positions(signalome_result: SignalomeCoreResult) -> pd.DataFrame:
    module_table = signalome_result.modules.to_frame().sort_index()
    relationships = signalome_result.modules.to_relationship_table()
    site_assignments = signalome_result.assignments.sites()
    protein_assignments = signalome_result.assignments.proteins()

    module_positions = pd.DataFrame(index=module_table.index.copy())
    module_positions.index.name = MODULE_ID_COLUMN
    module_positions["x"] = np.arange(len(module_positions), dtype=float)
    module_positions["y"] = 0.0
    module_positions["n_sites"] = (
        site_assignments.groupby(MODULE_ID_COLUMN)
        .size()
        .reindex(module_positions.index, fill_value=0)
        .astype(int)
    )
    module_positions["n_proteins"] = (
        protein_assignments.groupby(MODULE_ID_COLUMN)
        .size()
        .reindex(module_positions.index, fill_value=0)
        .astype(int)
    )

    dominant_relationships = (
        relationships.sort_values(
            [MODULE_ID_COLUMN, SHARE_PERCENT_COLUMN, KINASE_COLUMN],
            ascending=[True, False, True],
            kind="stable",
        )
        .drop_duplicates(subset=[MODULE_ID_COLUMN], keep="first")
        .set_index(MODULE_ID_COLUMN)
    )
    dominant_relationships = dominant_relationships.reindex(module_positions.index)
    module_positions["dominant_kinase"] = dominant_relationships[KINASE_COLUMN].fillna(
        ""
    )
    module_positions["dominant_share_percent"] = (
        dominant_relationships[SHARE_PERCENT_COLUMN].fillna(0.0).astype(float)
    )
    return module_positions


def _build_site_positions(
    signalome_result: SignalomeCoreResult,
    module_positions: pd.DataFrame,
) -> pd.DataFrame:
    site_assignments = signalome_result.assignments.sites().reset_index()
    site_assignments = site_assignments.sort_values(
        [MODULE_ID_COLUMN, PROTEIN_ID_COLUMN, SITE_ID_COLUMN],
        ascending=[True, True, True],
        kind="stable",
    )
    expression_matrix = signalome_result.expression_matrix

    if site_assignments.empty:
        site_positions = pd.DataFrame(
            columns=[
                PROTEIN_ID_COLUMN,
                MODULE_ID_COLUMN,
                TOP_KINASE_CANDIDATES_COLUMN,
                TOP_KINASE_WEIGHTS_COLUMN,
                TOP_KINASE_TIE_COUNT_COLUMN,
                TOP_KINASE_IS_AMBIGUOUS_COLUMN,
                TOP_SCORE_COLUMN,
                "x",
                "y",
                "module_x",
                "module_y",
                "position_in_module",
                "expression_mean",
                "expression_std",
            ]
        )
        site_positions.index = pd.Index([], name=SITE_ID_COLUMN, dtype=object)
        return site_positions

    module_sizes = (
        site_assignments.groupby(MODULE_ID_COLUMN, sort=True)[SITE_ID_COLUMN]
        .transform("size")
        .to_numpy(dtype=int, copy=False)
    )
    position_zero_based = (
        site_assignments.groupby(MODULE_ID_COLUMN, sort=True)
        .cumcount()
        .to_numpy(dtype=int, copy=False)
    )
    offsets = np.zeros(len(site_assignments), dtype=float)
    multi_member_mask = module_sizes > 1
    offsets[multi_member_mask] = _SITE_POSITION_OFFSET_BASE + position_zero_based[
        multi_member_mask
    ] * (_SITE_POSITION_OFFSET_SPAN / (module_sizes[multi_member_mask] - 1))

    module_x = (
        site_assignments[MODULE_ID_COLUMN]
        .map(module_positions["x"])
        .to_numpy(dtype=float, copy=False)
    )
    module_y = (
        site_assignments[MODULE_ID_COLUMN]
        .map(module_positions["y"])
        .to_numpy(dtype=float, copy=False)
    )
    site_ids = site_assignments[SITE_ID_COLUMN].astype(str)
    expression_values = expression_matrix.loc[site_ids].to_numpy(
        dtype=float, copy=False
    )

    site_positions = pd.DataFrame(
        {
            PROTEIN_ID_COLUMN: site_assignments[PROTEIN_ID_COLUMN]
            .astype(str)
            .to_numpy(dtype=object, copy=False),
            MODULE_ID_COLUMN: site_assignments[MODULE_ID_COLUMN].to_numpy(
                dtype=int, copy=False
            ),
            TOP_KINASE_CANDIDATES_COLUMN: site_assignments[
                TOP_KINASE_CANDIDATES_COLUMN
            ].map(serialize_top_kinase_candidates),
            TOP_KINASE_WEIGHTS_COLUMN: site_assignments[TOP_KINASE_WEIGHTS_COLUMN].map(
                serialize_top_kinase_weights
            ),
            TOP_KINASE_TIE_COUNT_COLUMN: site_assignments[
                TOP_KINASE_TIE_COUNT_COLUMN
            ].to_numpy(dtype=int, copy=False),
            TOP_KINASE_IS_AMBIGUOUS_COLUMN: site_assignments[
                TOP_KINASE_IS_AMBIGUOUS_COLUMN
            ].to_numpy(dtype=bool, copy=False),
            TOP_SCORE_COLUMN: site_assignments[TOP_SCORE_COLUMN].to_numpy(
                dtype=float, copy=False
            ),
            "x": module_x + offsets,
            "y": _SITE_BASELINE_Y - position_zero_based * _SITE_ROW_Y_STEP,
            "module_x": module_x,
            "module_y": module_y,
            "position_in_module": position_zero_based + 1,
            "expression_mean": np.mean(expression_values, axis=1),
            "expression_std": np.std(expression_values, axis=1, ddof=0),
        },
        index=pd.Index(
            site_ids.to_numpy(dtype=object, copy=False), name=SITE_ID_COLUMN
        ),
    )
    site_positions.index.name = SITE_ID_COLUMN
    return site_positions


def _build_kinase_positions(
    signalome_result: SignalomeCoreResult,
    module_positions: pd.DataFrame,
) -> pd.DataFrame:
    relationships = signalome_result.modules.to_relationship_table()
    network_nodes = signalome_result.network.nodes()
    kinases_of_interest = set(signalome_result.kinases_of_interest)
    module_table = signalome_result.modules.to_frame()
    kinase_order = [str(kinase) for kinase in module_table.columns]

    records: list[dict[str, object]] = []
    for default_position, kinase in enumerate(kinase_order):
        kinase_links = relationships.loc[relationships.loc[:, KINASE_COLUMN] == kinase]
        total_share_percent = float(kinase_links.loc[:, SHARE_PERCENT_COLUMN].sum())
        if total_share_percent > 0.0:
            base_x = float(
                np.average(
                    module_positions.loc[
                        kinase_links.loc[:, MODULE_ID_COLUMN].to_numpy(dtype=int),
                        "x",
                    ].to_numpy(dtype=float),
                    weights=kinase_links.loc[:, SHARE_PERCENT_COLUMN].to_numpy(
                        dtype=float
                    ),
                )
            )
        else:
            base_x = float(default_position)

        records.append(
            {
                KINASE_COLUMN: kinase,
                "base_x": base_x,
                "module_count": int(kinase_links.shape[0]),
                TOTAL_SHARE_PERCENT_COLUMN: total_share_percent,
                DEGREE_COLUMN: int(network_nodes.loc[kinase, DEGREE_COLUMN]),
                N_SUBSTRATES_COLUMN: int(
                    network_nodes.loc[kinase, N_SUBSTRATES_COLUMN]
                ),
                IS_KINASE_OF_INTEREST_COLUMN: kinase in kinases_of_interest,
            }
        )

    kinase_positions = pd.DataFrame.from_records(records)
    kinase_positions = kinase_positions.sort_values(
        ["base_x", KINASE_COLUMN],
        ascending=[True, True],
        kind="stable",
    )

    adjusted_x: dict[str, float] = {}
    for base_x, group in kinase_positions.groupby("base_x", sort=True):
        offsets = _centered_offsets(len(group), span=_KINASE_COLLISION_AVOIDANCE_SPAN)
        for (_, row), offset in zip(group.iterrows(), offsets, strict=True):
            adjusted_x[str(row[KINASE_COLUMN])] = float(base_x) + float(offset)

    kinase_positions["x"] = kinase_positions.loc[:, KINASE_COLUMN].map(adjusted_x)
    kinase_positions["y"] = _KINASE_BASELINE_Y
    kinase_positions = kinase_positions.set_index(KINASE_COLUMN)
    kinase_positions.index.name = KINASE_COLUMN
    return kinase_positions.loc[
        kinase_order,
        [
            "x",
            "y",
            "base_x",
            "module_count",
            TOTAL_SHARE_PERCENT_COLUMN,
            DEGREE_COLUMN,
            N_SUBSTRATES_COLUMN,
            IS_KINASE_OF_INTEREST_COLUMN,
        ],
    ]


def _build_kinase_module_links(
    signalome_result: SignalomeCoreResult,
    module_positions: pd.DataFrame,
    kinase_positions: pd.DataFrame,
) -> pd.DataFrame:
    relationships = signalome_result.modules.to_relationship_table()
    if relationships.empty:
        return pd.DataFrame(
            columns=[
                KINASE_COLUMN,
                MODULE_ID_COLUMN,
                SHARE_PERCENT_COLUMN,
                "kinase_x",
                "kinase_y",
                "module_x",
                "module_y",
                IS_KINASE_OF_INTEREST_COLUMN,
            ]
        )

    links = relationships.sort_values(
        [MODULE_ID_COLUMN, SHARE_PERCENT_COLUMN, KINASE_COLUMN],
        ascending=[True, False, True],
        kind="stable",
    ).reset_index(drop=True)
    links["kinase_x"] = links.loc[:, KINASE_COLUMN].map(kinase_positions.loc[:, "x"])
    links["kinase_y"] = links.loc[:, KINASE_COLUMN].map(kinase_positions.loc[:, "y"])
    links["module_x"] = links.loc[:, MODULE_ID_COLUMN].map(module_positions.loc[:, "x"])
    links["module_y"] = links.loc[:, MODULE_ID_COLUMN].map(module_positions.loc[:, "y"])
    links[IS_KINASE_OF_INTEREST_COLUMN] = links.loc[:, KINASE_COLUMN].map(
        kinase_positions.loc[:, IS_KINASE_OF_INTEREST_COLUMN]
    )
    return links.loc[
        :,
        [
            KINASE_COLUMN,
            MODULE_ID_COLUMN,
            SHARE_PERCENT_COLUMN,
            "kinase_x",
            "kinase_y",
            "module_x",
            "module_y",
            IS_KINASE_OF_INTEREST_COLUMN,
        ],
    ]


def _centered_offsets(count: int, *, span: float) -> np.ndarray:
    if count <= 0:
        return np.array([], dtype=float)
    if count == 1:
        return np.array([0.0], dtype=float)
    return np.linspace(-span / 2.0, span / 2.0, num=count, dtype=float)
