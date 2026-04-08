from __future__ import annotations

import json
from collections.abc import Mapping, Sequence

import pandas as pd

from .signalome_models import ExpandedSignalome
from .validation.errors import InputCompatibilityError

__all__ = [
    "build_expanded_signalomes",
    "build_kinase_module_relationship_table",
    "build_protein_assignment_table",
    "build_signalome_module_table",
    "build_site_assignments",
    "derive_protein_modules",
    "resolve_site_to_protein",
    "select_kinase_substrates",
]


def resolve_site_to_protein(
    *,
    site_ids: Sequence[str],
    site_to_protein: Mapping[str, str] | pd.Series | None,
) -> pd.Series:
    """Resolve aligned phosphosite IDs to validated protein IDs."""

    if site_to_protein is None:
        return parse_supported_site_ids(site_ids)

    if isinstance(site_to_protein, pd.Series):
        mapping: Mapping[str, str] = {
            str(site_id): str(protein_id)
            for site_id, protein_id in site_to_protein.items()
        }
    else:
        mapping = site_to_protein

    missing_site_ids = [site_id for site_id in site_ids if site_id not in mapping]
    if missing_site_ids:
        preview = ", ".join(missing_site_ids[:3])
        msg = (
            "site_to_protein must define a protein ID for every signalome site. "
            f"Missing mappings for: {preview}"
        )
        if len(missing_site_ids) > 3:
            msg += ", ..."
        raise InputCompatibilityError(msg)

    protein_ids = [str(mapping[site_id]).strip() for site_id in site_ids]
    invalid_site_ids = [
        site_id
        for site_id, protein_id in zip(site_ids, protein_ids, strict=True)
        if not protein_id
    ]
    if invalid_site_ids:
        preview = ", ".join(invalid_site_ids[:3])
        msg = (
            "site_to_protein must map every signalome site to a non-empty protein "
            f"ID. Invalid mappings for: {preview}"
        )
        if len(invalid_site_ids) > 3:
            msg += ", ..."
        raise InputCompatibilityError(msg)

    series = pd.Series(
        protein_ids, index=pd.Index(site_ids, dtype=object), dtype=object
    )
    series.index.name = "site_id"
    series.name = "protein_id"
    return series


def parse_supported_site_ids(site_ids: Sequence[str]) -> pd.Series:
    """Parse supported phosphosite identifiers into protein IDs."""

    protein_ids: list[str] = []
    invalid_site_ids: list[str] = []
    for site_id in site_ids:
        protein_id = protein_id_from_supported_site_id(site_id)
        if protein_id is None:
            invalid_site_ids.append(site_id)
            continue
        protein_ids.append(protein_id)

    if invalid_site_ids:
        preview = ", ".join(invalid_site_ids[:3])
        msg = (
            "Signalome construction requires either an explicit site_to_protein "
            "mapping or phosphosite identifiers in the supported 'PROTEIN;SITE;...' "
            f"format. Invalid site IDs: {preview}"
        )
        if len(invalid_site_ids) > 3:
            msg += ", ..."
        raise InputCompatibilityError(msg)

    series = pd.Series(
        protein_ids, index=pd.Index(site_ids, dtype=object), dtype=object
    )
    series.index.name = "site_id"
    series.name = "protein_id"
    return series


def protein_id_from_supported_site_id(site_id: str) -> str | None:
    """Extract a protein ID from a supported ``PROTEIN;SITE;...`` site ID."""

    parts = [part.strip() for part in str(site_id).split(";")]
    if len(parts) < 3:
        return None
    protein_id, residue = parts[0], parts[1]
    if not protein_id or not residue:
        return None
    return protein_id


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

    Tied top-kinase assignments are resolved deterministically by sorting the
    tied kinase names alphabetically. The full tie set is preserved in
    ``top_kinase_candidates`` and summarised by the tie-count diagnostics.
    """

    rows: list[dict[str, object]] = []
    for site_id in pred_mat.index.astype(str):
        protein_id = str(site_to_protein.loc[site_id])
        scores = pred_mat.loc[site_id]
        top_score = float(scores.max())
        top_kinases = sorted(
            str(kinase) for kinase in scores.index[scores.eq(top_score)]
        )
        rows.append(
            {
                "site_id": site_id,
                "protein_id": protein_id,
                "module_id": int(protein_modules.loc[protein_id]),
                "top_kinase": top_kinases[0],
                "top_kinase_candidates": json.dumps(top_kinases),
                "top_kinase_tie_count": len(top_kinases),
                "top_kinase_is_ambiguous": len(top_kinases) > 1,
                "top_score": top_score,
            }
        )

    site_assignments = pd.DataFrame.from_records(rows).set_index("site_id")
    site_assignments.index.name = "site_id"
    site_assignments = site_assignments.astype(
        {
            "protein_id": str,
            "module_id": int,
            "top_kinase": str,
            "top_kinase_candidates": str,
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

    selected: dict[str, tuple[str, ...]] = {}
    for kinase in pred_mat.columns.astype(str):
        mask = pred_mat.loc[:, kinase] > cutoff
        selected[kinase] = tuple(pred_mat.index[mask].astype(str).tolist())
    return selected


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
