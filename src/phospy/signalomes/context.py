"""Signalome context sidecar table builders."""

from __future__ import annotations

import json
from collections.abc import Mapping, Sequence

import numpy as np
import pandas as pd

from phospy.api.configs import (
    SIGNALOME_ASSIGNMENT_POLICY_CUTOFF_BINARY,
    SIGNALOME_ASSIGNMENT_POLICY_WEIGHTED_TOP,
    SignalomeAssignmentPolicy,
)
from phospy.errors.workflows import WorkflowStageError
from phospy.signalomes.assignments import _normalize_top_kinase_weights
from phospy.signalomes.constants import (
    MODULE_ID_COLUMN,
    PROTEIN_COLUMN,
    SITE_CLUSTER_COLUMN,
    SITE_ID_COLUMN,
    TOP_KINASE_COLUMN,
    TOP_KINASE_WEIGHTS_COLUMN,
    TOP_SCORE_COLUMN,
)

SITE_MEMBERSHIP_PROTEIN_MODULE_COLUMN = "protein_module_id"
SITE_MEMBERSHIP_INCLUDED_COLUMN = "included_in_module_table"
SITE_MEMBERSHIP_EXCLUDED_REASON_COLUMN = "excluded_reason"
SITE_MEMBERSHIP_GENE_SYMBOL_COLUMN = "gene_symbol"
SITE_MEMBERSHIP_TOP_KINASE_SCORE_COLUMN = "top_kinase_score"
SITE_MEMBERSHIP_TOP_KINASE_WEIGHT_COLUMN = "top_kinase_weight"
SITE_MEMBERSHIP_N_SUPPORTED_KINASES_COLUMN = "n_supported_kinases"

PROTEIN_SITE_CONTEXT_N_SITES_COLUMN = "n_sites"
PROTEIN_SITE_CONTEXT_SITE_IDS_COLUMN = "site_ids"
PROTEIN_SITE_CONTEXT_SITE_CLUSTERS_COLUMN = "site_clusters"
PROTEIN_SITE_CONTEXT_N_DISTINCT_SITE_CLUSTERS_COLUMN = "n_distinct_site_clusters"
PROTEIN_SITE_CONTEXT_PROTEIN_MODULE_COLUMN = "protein_module_id"
PROTEIN_SITE_CONTEXT_MULTI_SITE_COLUMN = "multi_site_protein"
PROTEIN_SITE_CONTEXT_AMBIGUOUS_CONTEXT_COLUMN = "ambiguous_module_context"
PROTEIN_SITE_CONTEXT_GENE_SYMBOL_COLUMN = "gene_symbol"
PROTEIN_SITE_CONTEXT_TOP_KINASES_BY_SITE_COLUMN = "top_kinases_by_site"
PROTEIN_SITE_CONTEXT_MODULE_IDS_BY_SITE_COLUMN = "module_ids_by_site"

EXCLUDED_REASON_DROPPED_ALL_MISSING_DOWNSTREAM_SCORES = (
    "dropped_all_missing_downstream_scores"
)
EXCLUDED_REASON_PROTEIN_NOT_ASSIGNED_TO_MODULE = "protein_not_assigned_to_module"
EXCLUDED_REASON_BELOW_SUBSTRATE_SUPPORT_CUTOFF = "below_substrate_support_cutoff"
EXCLUDED_REASON_NO_SUPPORTED_TOP_KINASE_WEIGHT = "no_supported_top_kinase_weight"
EXCLUDED_REASON_NOT_INCLUDED = "not_included_in_module_table"


def build_site_membership_table(
    *,
    module_assignments: pd.DataFrame,
    site_clusters: pd.Series,
    site_metadata: pd.DataFrame,
    prediction_matrix: pd.DataFrame,
    kinase_substrates: Mapping[str, Sequence[str]],
    substrate_support_cutoff: float,
    assignment_policy: SignalomeAssignmentPolicy,
) -> pd.DataFrame:
    """Build site-level signalome membership context sidecar table."""

    required_columns = {PROTEIN_COLUMN, MODULE_ID_COLUMN}
    if module_assignments.empty:
        return empty_site_membership_table()
    if not required_columns.issubset(module_assignments.columns):
        missing = sorted(required_columns.difference(module_assignments.columns))
        raise WorkflowStageError(
            "site membership context build requires module assignments with "
            f"columns {sorted(required_columns)}; missing columns: {missing}"
        )

    site_index = pd.Index(module_assignments.index.astype(str), name=SITE_ID_COLUMN)
    assignments = module_assignments.copy(deep=False)
    assignments.index = site_index

    cluster_series = site_clusters.copy(deep=False)
    cluster_series.index = pd.Index(
        cluster_series.index.astype(str), name=SITE_ID_COLUMN
    )
    aligned_clusters = cluster_series.reindex(site_index)
    aligned_clusters = aligned_clusters.astype("Int64")

    metadata = site_metadata.copy(deep=False)
    metadata.index = pd.Index(metadata.index.astype(str), name=SITE_ID_COLUMN)
    gene_symbols = (
        metadata.reindex(site_index)
        .loc[:, SITE_MEMBERSHIP_GENE_SYMBOL_COLUMN]
        .fillna("")
        if SITE_MEMBERSHIP_GENE_SYMBOL_COLUMN in metadata.columns
        else pd.Series("", index=site_index, dtype=object)
    )
    gene_symbols = gene_symbols.astype(str).str.strip()

    module_ids = assignments.loc[:, MODULE_ID_COLUMN].astype("int64")
    top_kinases = (
        assignments.loc[:, TOP_KINASE_COLUMN].astype(str)
        if TOP_KINASE_COLUMN in assignments.columns
        else pd.Series("", index=site_index, dtype=object)
    )
    top_scores = (
        assignments.loc[:, TOP_SCORE_COLUMN].astype(float)
        if TOP_SCORE_COLUMN in assignments.columns
        else pd.Series(np.nan, index=site_index, dtype=float)
    )

    aligned_prediction = prediction_matrix.reindex(site_index).astype(float)
    support_counts = (
        aligned_prediction.gt(float(substrate_support_cutoff))
        .sum(axis=1)
        .astype("int64")
    )
    supported_sites = _supported_site_set(kinase_substrates)

    has_top_kinase_support, top_kinase_weights = _resolve_top_kinase_context(
        assignments=assignments,
        top_kinases=top_kinases,
        site_index=site_index,
    )
    included_in_module = _resolve_included_in_module_table(
        assignment_policy=assignment_policy,
        module_ids=module_ids,
        site_index=site_index,
        supported_sites=supported_sites,
        has_top_kinase_support=has_top_kinase_support,
    )

    excluded_reasons = _resolve_excluded_reasons(
        assignment_policy=assignment_policy,
        module_ids=module_ids,
        site_index=site_index,
        included_in_module=included_in_module,
        aligned_clusters=aligned_clusters,
        supported_sites=supported_sites,
        has_top_kinase_support=has_top_kinase_support,
    )

    site_membership = pd.DataFrame(
        {
            SITE_ID_COLUMN: site_index.to_numpy(dtype=object, copy=False),
            PROTEIN_COLUMN: assignments.loc[:, PROTEIN_COLUMN]
            .astype(str)
            .to_numpy(
                dtype=object,
                copy=False,
            ),
            SITE_CLUSTER_COLUMN: aligned_clusters.to_numpy(copy=False),
            SITE_MEMBERSHIP_PROTEIN_MODULE_COLUMN: module_ids.to_numpy(
                dtype=np.int64,
                copy=False,
            ),
            SITE_MEMBERSHIP_INCLUDED_COLUMN: included_in_module.to_numpy(
                dtype=bool,
                copy=False,
            ),
            SITE_MEMBERSHIP_EXCLUDED_REASON_COLUMN: excluded_reasons.to_numpy(
                dtype=object,
                copy=False,
            ),
            SITE_MEMBERSHIP_GENE_SYMBOL_COLUMN: gene_symbols.to_numpy(
                dtype=object,
                copy=False,
            ),
            TOP_KINASE_COLUMN: top_kinases.to_numpy(dtype=object, copy=False),
            SITE_MEMBERSHIP_TOP_KINASE_SCORE_COLUMN: top_scores.to_numpy(
                dtype=float,
                copy=False,
            ),
            SITE_MEMBERSHIP_TOP_KINASE_WEIGHT_COLUMN: top_kinase_weights.to_numpy(
                dtype=float,
                copy=False,
            ),
            SITE_MEMBERSHIP_N_SUPPORTED_KINASES_COLUMN: support_counts.to_numpy(
                dtype=np.int64,
                copy=False,
            ),
        }
    )
    return site_membership.astype(
        {
            SITE_ID_COLUMN: str,
            PROTEIN_COLUMN: str,
            SITE_CLUSTER_COLUMN: "Int64",
            SITE_MEMBERSHIP_PROTEIN_MODULE_COLUMN: "int64",
            SITE_MEMBERSHIP_INCLUDED_COLUMN: bool,
            SITE_MEMBERSHIP_EXCLUDED_REASON_COLUMN: str,
            SITE_MEMBERSHIP_GENE_SYMBOL_COLUMN: str,
            TOP_KINASE_COLUMN: str,
            SITE_MEMBERSHIP_TOP_KINASE_SCORE_COLUMN: float,
            SITE_MEMBERSHIP_TOP_KINASE_WEIGHT_COLUMN: float,
            SITE_MEMBERSHIP_N_SUPPORTED_KINASES_COLUMN: "int64",
        }
    )


def build_protein_site_context_table(
    *,
    site_membership: pd.DataFrame,
) -> pd.DataFrame:
    """Build protein-level phosphosite context table."""

    required_columns = {
        SITE_ID_COLUMN,
        PROTEIN_COLUMN,
        SITE_CLUSTER_COLUMN,
        SITE_MEMBERSHIP_PROTEIN_MODULE_COLUMN,
        SITE_MEMBERSHIP_INCLUDED_COLUMN,
        SITE_MEMBERSHIP_GENE_SYMBOL_COLUMN,
        TOP_KINASE_COLUMN,
    }
    if site_membership.empty:
        return empty_protein_site_context_table()
    if not required_columns.issubset(site_membership.columns):
        missing = sorted(required_columns.difference(site_membership.columns))
        raise WorkflowStageError(
            "protein site-context build requires site-membership table columns "
            f"{sorted(required_columns)}; missing columns: {missing}"
        )

    membership = site_membership.copy(deep=False)
    protein_rows: list[dict[str, object]] = []
    for protein_id, group in membership.groupby(PROTEIN_COLUMN, sort=False):
        site_ids = [str(site_id) for site_id in group.loc[:, SITE_ID_COLUMN].tolist()]
        cluster_values = [
            None if pd.isna(value) else int(value)
            for value in group.loc[:, SITE_CLUSTER_COLUMN].tolist()
        ]
        module_values = [
            int(value)
            for value in group.loc[:, SITE_MEMBERSHIP_PROTEIN_MODULE_COLUMN]
            .astype("int64")
            .tolist()
        ]
        n_sites = len(site_ids)
        distinct_cluster_values = sorted(
            {
                int(value)
                for value in group.loc[:, SITE_CLUSTER_COLUMN].dropna().astype("int64")
            }
        )
        n_distinct_site_clusters = len(distinct_cluster_values)
        module_id = int(module_values[0]) if module_values else 0
        multi_site = n_sites > 1
        unique_module_ids = sorted(set(module_values))
        inclusion_states = set(
            bool(value)
            for value in group.loc[:, SITE_MEMBERSHIP_INCLUDED_COLUMN]
            .astype(bool)
            .tolist()
        )
        ambiguous_context = bool(
            multi_site
            and (
                n_distinct_site_clusters > 1
                or len(unique_module_ids) > 1
                or len(inclusion_states) > 1
            )
        )
        top_kinases_by_site = {
            site_id: str(top_kinase)
            for site_id, top_kinase in zip(
                site_ids,
                group.loc[:, TOP_KINASE_COLUMN].astype(str).tolist(),
                strict=True,
            )
        }
        module_ids_by_site = {
            site_id: int(module_id_value)
            for site_id, module_id_value in zip(site_ids, module_values, strict=True)
        }
        gene_symbols = [
            value
            for value in group.loc[:, SITE_MEMBERSHIP_GENE_SYMBOL_COLUMN]
            .astype(str)
            .str.strip()
            .tolist()
            if value != ""
        ]
        gene_symbol = gene_symbols[0] if gene_symbols else ""

        protein_rows.append(
            {
                PROTEIN_COLUMN: str(protein_id),
                PROTEIN_SITE_CONTEXT_N_SITES_COLUMN: int(n_sites),
                PROTEIN_SITE_CONTEXT_SITE_IDS_COLUMN: _serialize_json(site_ids),
                PROTEIN_SITE_CONTEXT_SITE_CLUSTERS_COLUMN: _serialize_json(
                    cluster_values
                ),
                PROTEIN_SITE_CONTEXT_N_DISTINCT_SITE_CLUSTERS_COLUMN: int(
                    n_distinct_site_clusters
                ),
                PROTEIN_SITE_CONTEXT_PROTEIN_MODULE_COLUMN: int(module_id),
                PROTEIN_SITE_CONTEXT_MULTI_SITE_COLUMN: bool(multi_site),
                PROTEIN_SITE_CONTEXT_AMBIGUOUS_CONTEXT_COLUMN: bool(ambiguous_context),
                PROTEIN_SITE_CONTEXT_GENE_SYMBOL_COLUMN: str(gene_symbol),
                PROTEIN_SITE_CONTEXT_TOP_KINASES_BY_SITE_COLUMN: _serialize_json(
                    top_kinases_by_site
                ),
                PROTEIN_SITE_CONTEXT_MODULE_IDS_BY_SITE_COLUMN: _serialize_json(
                    module_ids_by_site
                ),
            }
        )

    if not protein_rows:
        return empty_protein_site_context_table()
    protein_site_context = pd.DataFrame.from_records(protein_rows)
    return protein_site_context.astype(
        {
            PROTEIN_COLUMN: str,
            PROTEIN_SITE_CONTEXT_N_SITES_COLUMN: "int64",
            PROTEIN_SITE_CONTEXT_SITE_IDS_COLUMN: str,
            PROTEIN_SITE_CONTEXT_SITE_CLUSTERS_COLUMN: str,
            PROTEIN_SITE_CONTEXT_N_DISTINCT_SITE_CLUSTERS_COLUMN: "int64",
            PROTEIN_SITE_CONTEXT_PROTEIN_MODULE_COLUMN: "int64",
            PROTEIN_SITE_CONTEXT_MULTI_SITE_COLUMN: bool,
            PROTEIN_SITE_CONTEXT_AMBIGUOUS_CONTEXT_COLUMN: bool,
            PROTEIN_SITE_CONTEXT_GENE_SYMBOL_COLUMN: str,
            PROTEIN_SITE_CONTEXT_TOP_KINASES_BY_SITE_COLUMN: str,
            PROTEIN_SITE_CONTEXT_MODULE_IDS_BY_SITE_COLUMN: str,
        }
    )


def empty_site_membership_table() -> pd.DataFrame:
    """Return an empty site-membership table with stable schema."""

    return pd.DataFrame(
        {
            SITE_ID_COLUMN: pd.Series(dtype=str),
            PROTEIN_COLUMN: pd.Series(dtype=str),
            SITE_CLUSTER_COLUMN: pd.Series(dtype="Int64"),
            SITE_MEMBERSHIP_PROTEIN_MODULE_COLUMN: pd.Series(dtype="int64"),
            SITE_MEMBERSHIP_INCLUDED_COLUMN: pd.Series(dtype=bool),
            SITE_MEMBERSHIP_EXCLUDED_REASON_COLUMN: pd.Series(dtype=str),
            SITE_MEMBERSHIP_GENE_SYMBOL_COLUMN: pd.Series(dtype=str),
            TOP_KINASE_COLUMN: pd.Series(dtype=str),
            SITE_MEMBERSHIP_TOP_KINASE_SCORE_COLUMN: pd.Series(dtype=float),
            SITE_MEMBERSHIP_TOP_KINASE_WEIGHT_COLUMN: pd.Series(dtype=float),
            SITE_MEMBERSHIP_N_SUPPORTED_KINASES_COLUMN: pd.Series(dtype="int64"),
        }
    )


def empty_protein_site_context_table() -> pd.DataFrame:
    """Return an empty protein-site-context table with stable schema."""

    return pd.DataFrame(
        {
            PROTEIN_COLUMN: pd.Series(dtype=str),
            PROTEIN_SITE_CONTEXT_N_SITES_COLUMN: pd.Series(dtype="int64"),
            PROTEIN_SITE_CONTEXT_SITE_IDS_COLUMN: pd.Series(dtype=str),
            PROTEIN_SITE_CONTEXT_SITE_CLUSTERS_COLUMN: pd.Series(dtype=str),
            PROTEIN_SITE_CONTEXT_N_DISTINCT_SITE_CLUSTERS_COLUMN: pd.Series(
                dtype="int64"
            ),
            PROTEIN_SITE_CONTEXT_PROTEIN_MODULE_COLUMN: pd.Series(dtype="int64"),
            PROTEIN_SITE_CONTEXT_MULTI_SITE_COLUMN: pd.Series(dtype=bool),
            PROTEIN_SITE_CONTEXT_AMBIGUOUS_CONTEXT_COLUMN: pd.Series(dtype=bool),
            PROTEIN_SITE_CONTEXT_GENE_SYMBOL_COLUMN: pd.Series(dtype=str),
            PROTEIN_SITE_CONTEXT_TOP_KINASES_BY_SITE_COLUMN: pd.Series(dtype=str),
            PROTEIN_SITE_CONTEXT_MODULE_IDS_BY_SITE_COLUMN: pd.Series(dtype=str),
        }
    )


def _supported_site_set(
    kinase_substrates: Mapping[str, Sequence[str]],
) -> set[str]:
    supported: set[str] = set()
    for substrate_sites in kinase_substrates.values():
        for site_id in substrate_sites:
            supported.add(str(site_id))
    return supported


def _resolve_top_kinase_context(
    *,
    assignments: pd.DataFrame,
    top_kinases: pd.Series,
    site_index: pd.Index,
) -> tuple[pd.Series, pd.Series]:
    if TOP_KINASE_WEIGHTS_COLUMN not in assignments.columns:
        has_support = pd.Series(False, index=site_index, dtype=bool)
        weights = pd.Series(np.nan, index=site_index, dtype=float)
        return has_support, weights

    has_support_values: list[bool] = []
    top_kinase_weight_values: list[float] = []
    for site_id, top_kinase in zip(
        site_index.tolist(),
        top_kinases.astype(str).tolist(),
        strict=True,
    ):
        value = assignments.at[str(site_id), TOP_KINASE_WEIGHTS_COLUMN]
        normalized_weights = _normalize_top_kinase_weights(value, site_id=str(site_id))
        has_support_values.append(bool(normalized_weights))
        weight_lookup = {
            str(kinase): float(weight) for kinase, weight in normalized_weights
        }
        top_kinase_weight_values.append(
            float(weight_lookup.get(str(top_kinase), np.nan))
        )
    return (
        pd.Series(has_support_values, index=site_index, dtype=bool),
        pd.Series(top_kinase_weight_values, index=site_index, dtype=float),
    )


def _resolve_included_in_module_table(
    *,
    assignment_policy: SignalomeAssignmentPolicy,
    module_ids: pd.Series,
    site_index: pd.Index,
    supported_sites: set[str],
    has_top_kinase_support: pd.Series,
) -> pd.Series:
    module_assigned = module_ids.astype("int64") > 0
    if assignment_policy == SIGNALOME_ASSIGNMENT_POLICY_CUTOFF_BINARY:
        has_substrate_support = pd.Series(
            [str(site_id) in supported_sites for site_id in site_index.tolist()],
            index=site_index,
            dtype=bool,
        )
        return module_assigned & has_substrate_support
    if assignment_policy == SIGNALOME_ASSIGNMENT_POLICY_WEIGHTED_TOP:
        return module_assigned & has_top_kinase_support.astype(bool)
    return pd.Series(False, index=site_index, dtype=bool)


def _resolve_excluded_reasons(
    *,
    assignment_policy: SignalomeAssignmentPolicy,
    module_ids: pd.Series,
    site_index: pd.Index,
    included_in_module: pd.Series,
    aligned_clusters: pd.Series,
    supported_sites: set[str],
    has_top_kinase_support: pd.Series,
) -> pd.Series:
    reasons: list[str] = []
    for site_id in site_index.tolist():
        site_key = str(site_id)
        if bool(included_in_module.loc[site_key]):
            reasons.append("")
            continue
        if pd.isna(aligned_clusters.loc[site_key]):
            reasons.append(EXCLUDED_REASON_DROPPED_ALL_MISSING_DOWNSTREAM_SCORES)
            continue
        if int(module_ids.loc[site_key]) <= 0:
            reasons.append(EXCLUDED_REASON_PROTEIN_NOT_ASSIGNED_TO_MODULE)
            continue
        if assignment_policy == SIGNALOME_ASSIGNMENT_POLICY_CUTOFF_BINARY:
            if site_key not in supported_sites:
                reasons.append(EXCLUDED_REASON_BELOW_SUBSTRATE_SUPPORT_CUTOFF)
                continue
        elif assignment_policy == SIGNALOME_ASSIGNMENT_POLICY_WEIGHTED_TOP:
            if not bool(has_top_kinase_support.loc[site_key]):
                reasons.append(EXCLUDED_REASON_NO_SUPPORTED_TOP_KINASE_WEIGHT)
                continue
        reasons.append(EXCLUDED_REASON_NOT_INCLUDED)
    return pd.Series(reasons, index=site_index, dtype=object)


def _serialize_json(value: object) -> str:
    return json.dumps(value, separators=(",", ":"), ensure_ascii=True)


__all__ = [
    "PROTEIN_SITE_CONTEXT_AMBIGUOUS_CONTEXT_COLUMN",
    "PROTEIN_SITE_CONTEXT_GENE_SYMBOL_COLUMN",
    "PROTEIN_SITE_CONTEXT_MODULE_IDS_BY_SITE_COLUMN",
    "PROTEIN_SITE_CONTEXT_MULTI_SITE_COLUMN",
    "PROTEIN_SITE_CONTEXT_N_DISTINCT_SITE_CLUSTERS_COLUMN",
    "PROTEIN_SITE_CONTEXT_N_SITES_COLUMN",
    "PROTEIN_SITE_CONTEXT_PROTEIN_MODULE_COLUMN",
    "PROTEIN_SITE_CONTEXT_SITE_CLUSTERS_COLUMN",
    "PROTEIN_SITE_CONTEXT_SITE_IDS_COLUMN",
    "PROTEIN_SITE_CONTEXT_TOP_KINASES_BY_SITE_COLUMN",
    "SITE_MEMBERSHIP_EXCLUDED_REASON_COLUMN",
    "SITE_MEMBERSHIP_GENE_SYMBOL_COLUMN",
    "SITE_MEMBERSHIP_INCLUDED_COLUMN",
    "SITE_MEMBERSHIP_N_SUPPORTED_KINASES_COLUMN",
    "SITE_MEMBERSHIP_PROTEIN_MODULE_COLUMN",
    "SITE_MEMBERSHIP_TOP_KINASE_SCORE_COLUMN",
    "SITE_MEMBERSHIP_TOP_KINASE_WEIGHT_COLUMN",
    "build_protein_site_context_table",
    "build_site_membership_table",
    "empty_protein_site_context_table",
    "empty_site_membership_table",
]
