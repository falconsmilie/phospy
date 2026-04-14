from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from .results import SignalomeResult

__all__ = [
    "SignalomeMapData",
    "build_signalome_map_data",
]


@dataclass(frozen=True, slots=True)
class SignalomeMapData:
    """Serialisable map-ready plotting data derived from a signalome result.

    This model does not render charts. It exposes deterministic coordinate and
    relationship tables that plotting or export layers can consume.
    """

    module_positions: pd.DataFrame
    site_positions: pd.DataFrame
    kinase_positions: pd.DataFrame
    kinase_module_links: pd.DataFrame

    def modules(self, *, copy: bool = True) -> pd.DataFrame:
        """Return the canonical module-position table."""

        if copy:
            return self.module_positions.copy(deep=True)
        return self.module_positions

    def sites(self, *, copy: bool = True) -> pd.DataFrame:
        """Return the canonical site-position table."""

        if copy:
            return self.site_positions.copy(deep=True)
        return self.site_positions

    def kinases(self, *, copy: bool = True) -> pd.DataFrame:
        """Return the canonical kinase-position table."""

        if copy:
            return self.kinase_positions.copy(deep=True)
        return self.kinase_positions

    def links(self, *, copy: bool = True) -> pd.DataFrame:
        """Return the canonical kinase-to-module link table."""

        if copy:
            return self.kinase_module_links.copy(deep=True)
        return self.kinase_module_links

    def to_frames(self, *, copy: bool = True) -> dict[str, pd.DataFrame]:
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
        for name, frame in self.to_frames(copy=True).items():
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


def build_signalome_map_data(signalome_result: SignalomeResult) -> SignalomeMapData:
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


def _build_module_positions(signalome_result: SignalomeResult) -> pd.DataFrame:
    module_table = signalome_result.modules.to_frame(copy=True).sort_index()
    relationships = signalome_result.modules.to_relationship_table(copy=True)
    site_assignments = signalome_result.assignments.sites(copy=True)
    protein_assignments = signalome_result.assignments.proteins(copy=True)

    module_positions = pd.DataFrame(index=module_table.index.copy())
    module_positions.index.name = "module_id"
    module_positions["x"] = np.arange(len(module_positions), dtype=float)
    module_positions["y"] = 0.0
    module_positions["n_sites"] = (
        site_assignments.groupby("module_id")
        .size()
        .reindex(module_positions.index, fill_value=0)
        .astype(int)
    )
    module_positions["n_proteins"] = (
        protein_assignments.groupby("module_id")
        .size()
        .reindex(module_positions.index, fill_value=0)
        .astype(int)
    )

    dominant_relationships = (
        relationships.sort_values(
            ["module_id", "share_percent", "kinase"],
            ascending=[True, False, True],
            kind="stable",
        )
        .drop_duplicates(subset=["module_id"], keep="first")
        .set_index("module_id")
    )
    dominant_relationships = dominant_relationships.reindex(module_positions.index)
    module_positions["dominant_kinase"] = dominant_relationships["kinase"].fillna("")
    module_positions["dominant_share_percent"] = (
        dominant_relationships["share_percent"].fillna(0.0).astype(float)
    )
    return module_positions


def _build_site_positions(
    signalome_result: SignalomeResult,
    module_positions: pd.DataFrame,
) -> pd.DataFrame:
    site_assignments = signalome_result.assignments.sites(copy=True).reset_index()
    site_assignments = site_assignments.sort_values(
        ["module_id", "protein_id", "site_id"],
        ascending=[True, True, True],
        kind="stable",
    )
    expression_matrix = signalome_result.expression_matrix

    records: list[dict[str, object]] = []
    for module_id, group in site_assignments.groupby("module_id", sort=True):
        module_x = float(module_positions.loc[module_id, "x"])
        module_y = float(module_positions.loc[module_id, "y"])
        offsets = _centered_offsets(len(group), span=0.6)

        for position_in_module, ((_, row), offset) in enumerate(
            zip(group.iterrows(), offsets, strict=True),
            start=1,
        ):
            site_id = str(row["site_id"])
            expression_values = expression_matrix.loc[site_id].to_numpy(dtype=float)
            records.append(
                {
                    "site_id": site_id,
                    "protein_id": str(row["protein_id"]),
                    "module_id": int(row["module_id"]),
                    "top_kinase_candidates": str(row["top_kinase_candidates"]),
                    "top_kinase_weights": str(row["top_kinase_weights"]),
                    "top_kinase_tie_count": int(row["top_kinase_tie_count"]),
                    "top_kinase_is_ambiguous": bool(row["top_kinase_is_ambiguous"]),
                    "top_score": float(row["top_score"]),
                    "x": module_x + float(offset),
                    "y": -0.25 - (position_in_module - 1) * 0.15,
                    "module_x": module_x,
                    "module_y": module_y,
                    "position_in_module": position_in_module,
                    "expression_mean": float(np.mean(expression_values)),
                    "expression_std": float(np.std(expression_values, ddof=0)),
                }
            )

    site_positions = pd.DataFrame.from_records(records).set_index("site_id")
    site_positions.index.name = "site_id"
    return site_positions


def _build_kinase_positions(
    signalome_result: SignalomeResult,
    module_positions: pd.DataFrame,
) -> pd.DataFrame:
    relationships = signalome_result.modules.to_relationship_table(copy=True)
    network_nodes = signalome_result.network.nodes(copy=True)
    kinases_of_interest = set(signalome_result.kinases_of_interest)
    kinase_order = [
        str(kinase) for kinase in signalome_result.signalome_modules.columns
    ]

    records: list[dict[str, object]] = []
    for default_position, kinase in enumerate(kinase_order):
        kinase_links = relationships.loc[relationships.loc[:, "kinase"] == kinase]
        total_share_percent = float(kinase_links.loc[:, "share_percent"].sum())
        if total_share_percent > 0.0:
            base_x = float(
                np.average(
                    module_positions.loc[
                        kinase_links.loc[:, "module_id"].to_numpy(dtype=int),
                        "x",
                    ].to_numpy(dtype=float),
                    weights=kinase_links.loc[:, "share_percent"].to_numpy(dtype=float),
                )
            )
        else:
            base_x = float(default_position)

        records.append(
            {
                "kinase": kinase,
                "base_x": base_x,
                "module_count": int(kinase_links.shape[0]),
                "total_share_percent": total_share_percent,
                "degree": int(network_nodes.loc[kinase, "degree"]),
                "n_substrates": int(network_nodes.loc[kinase, "n_substrates"]),
                "is_kinase_of_interest": kinase in kinases_of_interest,
            }
        )

    kinase_positions = pd.DataFrame.from_records(records)
    kinase_positions = kinase_positions.sort_values(
        ["base_x", "kinase"],
        ascending=[True, True],
        kind="stable",
    )

    adjusted_x: dict[str, float] = {}
    for base_x, group in kinase_positions.groupby("base_x", sort=True):
        offsets = _centered_offsets(len(group), span=0.4)
        for (_, row), offset in zip(group.iterrows(), offsets, strict=True):
            adjusted_x[str(row["kinase"])] = float(base_x) + float(offset)

    kinase_positions["x"] = kinase_positions.loc[:, "kinase"].map(adjusted_x)
    kinase_positions["y"] = 1.0
    kinase_positions = kinase_positions.set_index("kinase")
    kinase_positions.index.name = "kinase"
    return kinase_positions.loc[
        kinase_order,
        [
            "x",
            "y",
            "base_x",
            "module_count",
            "total_share_percent",
            "degree",
            "n_substrates",
            "is_kinase_of_interest",
        ],
    ]


def _build_kinase_module_links(
    signalome_result: SignalomeResult,
    module_positions: pd.DataFrame,
    kinase_positions: pd.DataFrame,
) -> pd.DataFrame:
    relationships = signalome_result.modules.to_relationship_table(copy=True)
    if relationships.empty:
        return pd.DataFrame(
            columns=[
                "kinase",
                "module_id",
                "share_percent",
                "kinase_x",
                "kinase_y",
                "module_x",
                "module_y",
                "is_kinase_of_interest",
            ]
        )

    links = relationships.sort_values(
        ["module_id", "share_percent", "kinase"],
        ascending=[True, False, True],
        kind="stable",
    ).reset_index(drop=True)
    links["kinase_x"] = links.loc[:, "kinase"].map(kinase_positions.loc[:, "x"])
    links["kinase_y"] = links.loc[:, "kinase"].map(kinase_positions.loc[:, "y"])
    links["module_x"] = links.loc[:, "module_id"].map(module_positions.loc[:, "x"])
    links["module_y"] = links.loc[:, "module_id"].map(module_positions.loc[:, "y"])
    links["is_kinase_of_interest"] = links.loc[:, "kinase"].map(
        kinase_positions.loc[:, "is_kinase_of_interest"]
    )
    return links.loc[
        :,
        [
            "kinase",
            "module_id",
            "share_percent",
            "kinase_x",
            "kinase_y",
            "module_x",
            "module_y",
            "is_kinase_of_interest",
        ],
    ]


def _centered_offsets(count: int, *, span: float) -> np.ndarray:
    if count <= 0:
        return np.array([], dtype=float)
    if count == 1:
        return np.array([0.0], dtype=float)
    return np.linspace(-span / 2.0, span / 2.0, num=count, dtype=float)
