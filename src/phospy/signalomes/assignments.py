from __future__ import annotations

from collections.abc import Mapping, Sequence

import numpy as np
import pandas as pd

from ..errors import InputCompatibilityError
from .results import ExpandedSignalome
from .site_ids import (
    parse_supported_signalome_site_ids,
    protein_id_from_supported_signalome_site_id,
    resolve_signalome_site_to_protein,
)

__all__ = [
    "build_expanded_signalomes",
    "build_kinase_module_relationship_table",
    "build_protein_assignment_table",
    "build_signalome_module_table",
    "build_site_assignments",
    "derive_protein_modules",
    "parse_supported_site_ids",
    "protein_id_from_supported_site_id",
    "resolve_site_to_protein",
    "select_kinase_substrates",
]


def _as_unique_string_index(
    index: pd.Index,
    *,
    context: str,
    label: str,
) -> pd.Index:
    resolved_index = pd.Index(index.astype(str), name=label)
    if resolved_index.has_duplicates:
        duplicates = sorted(
            {str(site_id) for site_id in resolved_index[resolved_index.duplicated()]}
        )
        preview = ", ".join(duplicates[:3])
        suffix = "..." if len(duplicates) > 3 else ""
        msg = f"{context} contains duplicate {label} values: {preview}{suffix}"
        raise InputCompatibilityError(msg)
    return resolved_index


def _as_unique_string_series_index(
    series: pd.Series,
    *,
    context: str,
    label: str,
) -> pd.Series:
    resolved = series.copy()
    resolved.index = _as_unique_string_index(
        resolved.index,
        context=context,
        label=label,
    )
    return resolved


def resolve_site_to_protein(
    *,
    site_ids: Sequence[str],
    site_to_protein: Mapping[str, str] | pd.Series | None,
) -> pd.Series:
    """Resolve aligned phosphosite IDs to validated protein IDs."""

    return resolve_signalome_site_to_protein(
        site_ids=site_ids,
        site_to_protein=site_to_protein,
        missing_mapping_context=(
            "site_to_protein must define a protein ID for every signalome site. "
            "Missing mappings for"
        ),
        invalid_mapping_context=(
            "site_to_protein must map every signalome site to a non-empty protein "
            "ID. Invalid mappings for"
        ),
        invalid_site_id_context=(
            "Signalome construction requires either an explicit site_to_protein "
            "mapping or phosphosite identifiers in the supported 'PROTEIN;SITE;...' "
            "format. Invalid site IDs"
        ),
    )


def parse_supported_site_ids(site_ids: Sequence[str]) -> pd.Series:
    """Parse supported phosphosite identifiers into protein IDs."""

    return parse_supported_signalome_site_ids(
        site_ids,
        invalid_site_id_context=(
            "Signalome construction requires either an explicit site_to_protein "
            "mapping or phosphosite identifiers in the supported 'PROTEIN;SITE;...' "
            "format. Invalid site IDs"
        ),
    )


def protein_id_from_supported_site_id(site_id: object) -> str | None:
    """Extract a protein ID from a supported ``PROTEIN;SITE;...`` site ID."""

    return protein_id_from_supported_signalome_site_id(site_id)


def derive_protein_modules(
    *,
    site_clusters: pd.Series,
    site_to_protein: pd.Series,
) -> pd.Series:
    """Collapse site clusters into protein-level module assignments."""

    aligned_site_to_protein = site_to_protein.loc[site_clusters.index.astype(str)]
    proteins = pd.Index(aligned_site_to_protein.tolist(), dtype=object)
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


def build_site_assignments(
    *,
    pred_mat: pd.DataFrame,
    protein_modules: pd.Series,
    site_to_protein: pd.Series,
) -> pd.DataFrame:
    """Build the site-level assignment table from prediction outputs.

    Tied top-kinase assignments are preserved as weighted multi-assignments.
    Kinases sharing the top score for a site receive equal fractional weight in
    ``top_kinase_weights`` as ordered ``(kinase, weight)`` tuples (e.g. two-way
    tie -> ``(("A", 0.5), ("B", 0.5))``), while the full tie set remains
    available in ``top_kinase_candidates``.
    """

    if pred_mat.shape[1] == 0:
        msg = "pred_mat must contain at least one kinase column"
        raise InputCompatibilityError(msg)

    site_index = _as_unique_string_index(
        pred_mat.index,
        context="pred_mat",
        label="site_id",
    )
    resolved_site_to_protein = _as_unique_string_series_index(
        site_to_protein,
        context="site_to_protein",
        label="site_id",
    )
    resolved_site_to_protein = resolved_site_to_protein.astype(str)
    missing_site_ids = sorted(
        {
            site_id
            for site_id in site_index
            if site_id not in resolved_site_to_protein.index
        }
    )
    if missing_site_ids:
        preview = ", ".join(missing_site_ids[:3])
        suffix = "..." if len(missing_site_ids) > 3 else ""
        msg = (
            "site_to_protein must define a protein ID for every pred_mat site. "
            f"Missing mappings for: {preview}{suffix}"
        )
        raise InputCompatibilityError(msg)

    resolved_protein_modules = _as_unique_string_series_index(
        protein_modules,
        context="protein_modules",
        label="protein_id",
    ).astype(int)

    protein_ids = resolved_site_to_protein.loc[site_index]
    missing_protein_ids = sorted(
        {
            str(protein_id)
            for protein_id in protein_ids.tolist()
            if str(protein_id) not in resolved_protein_modules.index
        }
    )
    if missing_protein_ids:
        preview = ", ".join(missing_protein_ids[:3])
        suffix = "..." if len(missing_protein_ids) > 3 else ""
        msg = (
            "protein_modules must define module IDs for every protein mapped "
            f"from pred_mat sites. Missing proteins: {preview}{suffix}"
        )
        raise InputCompatibilityError(msg)

    sorted_kinase_columns = sorted(str(kinase) for kinase in pred_mat.columns)
    sorted_pred_mat = pred_mat.loc[:, sorted_kinase_columns]

    top_scores = sorted_pred_mat.max(axis=1).astype(float)
    top_score_mask = sorted_pred_mat.eq(top_scores, axis=0)
    top_score_mask_values = top_score_mask.to_numpy(dtype=bool, copy=False)
    top_kinase_names = top_score_mask.columns.to_numpy(dtype=object, copy=False)

    top_kinase_tie_count = top_score_mask_values.sum(axis=1).astype(int)
    top_kinase_candidates: list[tuple[str, ...]] = []
    top_kinase_weights: list[tuple[tuple[str, float], ...]] = []
    for mask_row, tie_count in zip(
        top_score_mask_values, top_kinase_tie_count, strict=True
    ):
        tied_kinases = tuple(str(kinase) for kinase in top_kinase_names[mask_row])
        weight = 1.0 / float(tie_count)
        top_kinase_candidates.append(tied_kinases)
        top_kinase_weights.append(tuple((kinase, weight) for kinase in tied_kinases))

    module_ids = protein_ids.map(resolved_protein_modules).astype(int)

    site_assignments = pd.DataFrame(
        {
            "protein_id": protein_ids.to_numpy(dtype=object, copy=False),
            "module_id": module_ids.to_numpy(dtype=int, copy=False),
            "top_kinase_candidates": top_kinase_candidates,
            "top_kinase_weights": top_kinase_weights,
            "top_kinase_tie_count": top_kinase_tie_count,
            "top_kinase_is_ambiguous": top_kinase_tie_count > 1,
            "top_score": top_scores.to_numpy(dtype=float, copy=False),
        },
        index=site_index,
    )
    site_assignments = site_assignments.astype(
        {
            "protein_id": str,
            "module_id": int,
            "top_kinase_tie_count": int,
            "top_kinase_is_ambiguous": bool,
            "top_score": float,
        }
    )
    return site_assignments


def build_protein_assignment_table(*, site_assignments: pd.DataFrame) -> pd.DataFrame:
    """Build the protein-level assignment table from site assignments."""

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


def select_kinase_substrates(
    *,
    pred_mat: pd.DataFrame,
    cutoff: float,
) -> dict[str, tuple[str, ...]]:
    """Select site substrates per kinase from the prediction matrix."""

    kinase_names = pred_mat.columns.astype(str).to_numpy(dtype=object, copy=False)
    site_ids = pred_mat.index.astype(str).to_numpy(dtype=object, copy=False)
    substrate_mask = pred_mat.to_numpy(dtype=float, copy=False) > cutoff

    return {
        str(kinase): tuple(site_ids[substrate_mask[:, position]].tolist())
        for position, kinase in enumerate(kinase_names)
    }


def build_signalome_module_table(
    *,
    site_assignments: pd.DataFrame,
    kinase_substrates: Mapping[str, Sequence[str]],
) -> pd.DataFrame:
    """Build the wide module-by-kinase signalome table."""

    module_ids = sorted(
        {int(value) for value in site_assignments["module_id"].tolist()}
    )
    kinase_names = list(kinase_substrates)
    module_index = pd.Index(module_ids, name="module_id")
    kinase_index = pd.Index(kinase_names, name="kinase")

    protein_to_module = (
        site_assignments.loc[:, ["protein_id", "module_id"]]
        .drop_duplicates(subset=["protein_id"])
        .set_index("protein_id")
        .loc[:, "module_id"]
        .astype(int)
    )

    substrate_rows = [
        {"kinase": str(kinase), "site_id": str(site_id)}
        for kinase, substrates in kinase_substrates.items()
        for site_id in substrates
    ]
    if not substrate_rows:
        return pd.DataFrame(0.0, index=module_index, columns=kinase_index).round(3)

    substrate_table = pd.DataFrame.from_records(substrate_rows)
    substrate_table = substrate_table.astype({"kinase": str, "site_id": str})
    substrate_table = substrate_table.drop_duplicates(subset=["kinase", "site_id"])

    site_to_protein = site_assignments.loc[:, ["protein_id"]].copy()
    site_to_protein.index = pd.Index(site_to_protein.index.astype(str), name="site_id")
    site_to_protein["protein_id"] = site_to_protein["protein_id"].astype(str)

    protein_hits = substrate_table.join(site_to_protein, on="site_id", how="left")
    protein_hits = protein_hits.dropna(subset=["protein_id"])
    protein_hits = protein_hits.drop_duplicates(subset=["kinase", "protein_id"])
    protein_hits = protein_hits.join(
        protein_to_module.rename("module_id"),
        on="protein_id",
        how="inner",
    )

    if protein_hits.empty:
        return pd.DataFrame(0.0, index=module_index, columns=kinase_index).round(3)

    counts = (
        protein_hits.groupby(["module_id", "kinase"], sort=True)
        .size()
        .astype(float)
        .unstack(fill_value=0.0)
    )
    module_table = counts.reindex(
        index=module_index, columns=kinase_index, fill_value=0.0
    )
    module_table = module_table.astype(float)

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


def build_kinase_module_relationship_table(
    *,
    module_table: pd.DataFrame,
) -> pd.DataFrame:
    """Build the long kinase-to-module relationship table."""

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


def build_expanded_signalomes(
    *,
    kinases_of_interest: Sequence[str],
    kinase_network: Mapping[str, Sequence[str]],
    kinase_substrates: Mapping[str, Sequence[str]],
    signalome_modules: pd.DataFrame,
    site_assignments: pd.DataFrame,
    expression_matrix: pd.DataFrame,
    min_kinase_module_share_percent: float,
) -> dict[str, ExpandedSignalome]:
    """Build expanded kinase-specific signalome views."""

    expanded: dict[str, ExpandedSignalome] = {}
    available_sites = site_assignments.index.astype(str)
    site_positions = {
        str(site_id): position for position, site_id in enumerate(available_sites)
    }
    site_module_ids = site_assignments.loc[:, "module_id"].to_numpy(
        dtype=int, copy=False
    )
    signalome_module_values = signalome_modules.to_numpy(dtype=float, copy=False)
    signalome_module_ids = signalome_modules.index.to_numpy(dtype=int, copy=False)
    signalome_kinase_positions = {
        str(kinase): position
        for position, kinase in enumerate(signalome_modules.columns.astype(str))
    }

    for kinase in kinases_of_interest:
        linked_kinases = tuple(dict.fromkeys((kinase, *kinase_network.get(kinase, ()))))
        kinase_position = signalome_kinase_positions[str(kinase)]
        regulated_module_ids_array = signalome_module_ids[
            signalome_module_values[:, kinase_position]
            > min_kinase_module_share_percent
        ]
        regulated_module_ids = tuple(
            int(module_id) for module_id in regulated_module_ids_array
        )

        substrate_positions = np.fromiter(
            (
                site_positions[str(site_id)]
                for linked_kinase in linked_kinases
                for site_id in kinase_substrates.get(linked_kinase, ())
                if str(site_id) in site_positions
            ),
            dtype=int,
        )
        if substrate_positions.size > 0:
            unique_substrate_positions = np.unique(substrate_positions)
            selected_positions = unique_substrate_positions[
                np.isin(
                    site_module_ids[unique_substrate_positions],
                    regulated_module_ids_array,
                )
            ]
        else:
            selected_positions = np.asarray([], dtype=int)

        expanded[kinase] = ExpandedSignalome(
            kinase=kinase,
            linked_kinases=linked_kinases,
            regulated_module_ids=regulated_module_ids,
            expression_matrix=expression_matrix.take(selected_positions),
            site_assignments=site_assignments.take(selected_positions),
        )

    return expanded
