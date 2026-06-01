"""Signalome assignment domain services."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import pandas as pd

from phospy.errors.workflows import WorkflowStageError
from phospy.science.signalomes.constants import (
    DISPLAY_ID_COLUMN,
    GENE_SYMBOL_COLUMN,
    ISOFORM_ID_COLUMN,
    LEXICOGRAPHIC_TIE_BREAK_POLICY,
    MODULE_ID_COLUMN,
    MODULE_TOP_KINASE_CANDIDATES_COLUMN,
    MODULE_TOP_KINASE_COLUMN,
    MODULE_TOP_KINASE_IS_AMBIGUOUS_COLUMN,
    MODULE_TOP_KINASE_SELECTION_POLICY_COLUMN,
    MODULE_TOP_KINASE_TIE_COUNT_COLUMN,
    NO_SUPPORT_SELECTION_POLICY,
    PROTEIN_ACCESSION_COLUMN,
    PROTEIN_COLUMN,
    SITE_COLUMN,
    SITE_KEY_COLUMN,
    TOP_KINASE_CANDIDATES_COLUMN,
    TOP_KINASE_COLUMN,
    TOP_KINASE_IS_AMBIGUOUS_COLUMN,
    TOP_KINASE_SELECTION_POLICY_COLUMN,
    TOP_KINASE_TIE_COUNT_COLUMN,
    TOP_KINASE_WEIGHTS_COLUMN,
    TOP_SCORE_COLUMN,
    UNSUPPORTED_KINASE,
)


def build_module_assignments(
    *,
    prediction_matrix: pd.DataFrame,
    site_to_protein: pd.Series,
    site_metadata: pd.DataFrame | None = None,
    protein_modules: pd.Series | None = None,
) -> pd.DataFrame:
    """Build site-level module assignments with explicit tie-handling metadata."""

    if prediction_matrix.shape[1] == 0:
        raise WorkflowStageError("prediction matrix must contain at least one kinase")
    site_index = _as_unique_string_index(prediction_matrix.index, context="pred_mat")
    resolved_site_to_protein = _resolve_site_to_protein(
        site_index=site_index,
        site_to_protein=site_to_protein,
    )
    site_identity = _resolve_site_identity_columns(
        site_index=site_index,
        resolved_site_to_protein=resolved_site_to_protein,
        site_metadata=site_metadata,
    )

    sorted_kinase_columns = sorted(str(kinase) for kinase in prediction_matrix.columns)
    sorted_predictions = prediction_matrix.loc[:, sorted_kinase_columns].astype(float)
    top_scores = sorted_predictions.max(axis=1).astype(float)
    top_score_mask = sorted_predictions.eq(top_scores, axis=0).to_numpy(
        dtype=bool, copy=False
    )
    kinase_names = sorted_predictions.columns.to_numpy(dtype=object, copy=False)
    tie_counts = top_score_mask.sum(axis=1).astype("int64")

    top_kinase_candidates: list[tuple[str, ...]] = []
    top_kinase_weights: list[tuple[tuple[str, float], ...]] = []
    top_kinases: list[str] = []
    top_kinase_policies: list[str] = []
    for mask_row, tie_count in zip(top_score_mask, tie_counts, strict=True):
        candidates = tuple(str(kinase) for kinase in kinase_names[mask_row])
        top_kinase_candidates.append(candidates)
        if tie_count == 0:
            top_kinase_weights.append(())
            top_kinases.append(UNSUPPORTED_KINASE)
            top_kinase_policies.append(NO_SUPPORT_SELECTION_POLICY)
            continue
        weight = 1.0 / float(tie_count)
        top_kinase_weights.append(tuple((kinase, weight) for kinase in candidates))
        top_kinases.append(candidates[0])
        top_kinase_policies.append(LEXICOGRAPHIC_TIE_BREAK_POLICY)

    top_kinase_series = pd.Series(
        top_kinases,
        index=site_index.copy(),
        dtype=object,
        name=TOP_KINASE_COLUMN,
    )
    site_module_resolution = _resolve_site_module_resolution(
        top_kinases=top_kinase_series,
        top_kinase_weights=top_kinase_weights,
        site_to_protein=resolved_site_to_protein,
        protein_modules=protein_modules,
    )

    assignments = pd.DataFrame(
        {
            SITE_KEY_COLUMN: site_identity.loc[:, SITE_KEY_COLUMN].to_numpy(
                dtype=object,
                copy=False,
            ),
            DISPLAY_ID_COLUMN: site_identity.loc[:, DISPLAY_ID_COLUMN].to_numpy(
                dtype=object,
                copy=False,
            ),
            GENE_SYMBOL_COLUMN: site_identity.loc[:, GENE_SYMBOL_COLUMN].to_numpy(
                dtype=object,
                copy=False,
            ),
            SITE_COLUMN: site_identity.loc[:, SITE_COLUMN].to_numpy(
                dtype=object,
                copy=False,
            ),
            PROTEIN_COLUMN: resolved_site_to_protein.to_numpy(dtype=object, copy=False),
            PROTEIN_ACCESSION_COLUMN: site_identity.loc[
                :, PROTEIN_ACCESSION_COLUMN
            ].to_numpy(
                dtype=object,
                copy=False,
            ),
            ISOFORM_ID_COLUMN: site_identity.loc[:, ISOFORM_ID_COLUMN].to_numpy(
                dtype=object,
                copy=False,
            ),
            MODULE_ID_COLUMN: site_module_resolution.loc[:, MODULE_ID_COLUMN].to_numpy(
                dtype=np.int64, copy=False
            ),
            TOP_KINASE_COLUMN: top_kinase_series.to_numpy(dtype=object, copy=False),
            TOP_SCORE_COLUMN: top_scores.to_numpy(dtype=float, copy=False),
            TOP_KINASE_CANDIDATES_COLUMN: top_kinase_candidates,
            TOP_KINASE_WEIGHTS_COLUMN: top_kinase_weights,
            TOP_KINASE_TIE_COUNT_COLUMN: tie_counts,
            TOP_KINASE_IS_AMBIGUOUS_COLUMN: tie_counts > 1,
            TOP_KINASE_SELECTION_POLICY_COLUMN: top_kinase_policies,
            MODULE_TOP_KINASE_COLUMN: site_module_resolution.loc[
                :, MODULE_TOP_KINASE_COLUMN
            ].to_numpy(dtype=object, copy=False),
            MODULE_TOP_KINASE_CANDIDATES_COLUMN: site_module_resolution.loc[
                :, MODULE_TOP_KINASE_CANDIDATES_COLUMN
            ].to_numpy(dtype=object, copy=False),
            MODULE_TOP_KINASE_TIE_COUNT_COLUMN: site_module_resolution.loc[
                :, MODULE_TOP_KINASE_TIE_COUNT_COLUMN
            ].to_numpy(dtype=np.int64, copy=False),
            MODULE_TOP_KINASE_IS_AMBIGUOUS_COLUMN: site_module_resolution.loc[
                :, MODULE_TOP_KINASE_IS_AMBIGUOUS_COLUMN
            ].to_numpy(dtype=bool, copy=False),
            MODULE_TOP_KINASE_SELECTION_POLICY_COLUMN: site_module_resolution.loc[
                :, MODULE_TOP_KINASE_SELECTION_POLICY_COLUMN
            ].to_numpy(dtype=object, copy=False),
        },
        index=site_index.copy(),
    )
    return assignments.astype(
        {
            SITE_KEY_COLUMN: str,
            DISPLAY_ID_COLUMN: str,
            GENE_SYMBOL_COLUMN: str,
            SITE_COLUMN: str,
            PROTEIN_COLUMN: str,
            PROTEIN_ACCESSION_COLUMN: str,
            ISOFORM_ID_COLUMN: str,
            MODULE_ID_COLUMN: "int64",
            TOP_KINASE_COLUMN: str,
            TOP_SCORE_COLUMN: float,
            TOP_KINASE_TIE_COUNT_COLUMN: "int64",
            TOP_KINASE_IS_AMBIGUOUS_COLUMN: bool,
            TOP_KINASE_SELECTION_POLICY_COLUMN: str,
            MODULE_TOP_KINASE_COLUMN: str,
            MODULE_TOP_KINASE_TIE_COUNT_COLUMN: "int64",
            MODULE_TOP_KINASE_IS_AMBIGUOUS_COLUMN: bool,
            MODULE_TOP_KINASE_SELECTION_POLICY_COLUMN: str,
        }
    )


def select_kinase_substrates(
    *,
    prediction_matrix: pd.DataFrame,
    cutoff: float,
) -> dict[str, tuple[str, ...]]:
    """Select phosphosites supported per kinase above `cutoff`."""

    site_ids = prediction_matrix.index.astype(str).to_numpy(dtype=object, copy=False)
    kinase_names = prediction_matrix.columns.astype(str).to_numpy(
        dtype=object, copy=False
    )
    support_mask = prediction_matrix.to_numpy(dtype=float, copy=False) > float(cutoff)
    return {
        str(kinase): tuple(site_ids[support_mask[:, index]].tolist())
        for index, kinase in enumerate(kinase_names)
    }


def _normalize_top_kinase_weights(
    value: object,
    *,
    site_id: str,
) -> tuple[tuple[str, float], ...]:
    if isinstance(value, dict):
        pairs = tuple((str(kinase), float(weight)) for kinase, weight in value.items())
    elif isinstance(value, (tuple, list)):
        pairs = _normalize_top_kinase_weight_pairs(value, site_id=site_id)
    elif value is None:
        pairs = ()
    else:
        raise WorkflowStageError(
            "top_kinase_weights entries must be dicts or (kinase, weight) sequences; "
            f"received {type(value).__name__} at site_id='{site_id}'"
        )
    if not pairs:
        return ()
    positive_pairs = tuple((kinase, weight) for kinase, weight in pairs if weight > 0.0)
    if not positive_pairs:
        return ()
    return positive_pairs


def _normalize_top_kinase_weight_pairs(
    values: tuple[object, ...] | list[object],
    *,
    site_id: str,
) -> tuple[tuple[str, float], ...]:
    normalized: list[tuple[str, float]] = []
    for value in values:
        if not isinstance(value, (tuple, list)) or len(value) != 2:
            raise WorkflowStageError(
                "top_kinase_weights entries must be (kinase, weight) pairs; "
                f"received invalid entry at site_id='{site_id}'"
            )
        kinase, weight = value
        try:
            normalized.append((str(kinase), float(weight)))
        except (TypeError, ValueError) as exc:
            raise WorkflowStageError(
                "top_kinase_weights entries must contain float-compatible weights at "
                f"site_id='{site_id}'"
            ) from exc
    return tuple(normalized)


def _as_unique_string_index(index: pd.Index, *, context: str) -> pd.Index:
    resolved = pd.Index(index.astype(str), name=SITE_KEY_COLUMN)
    if not resolved.has_duplicates:
        return resolved
    duplicates = sorted(
        {str(site_id) for site_id in resolved[resolved.duplicated(keep=False)]}
    )
    preview = ", ".join(duplicates[:3])
    suffix = "..." if len(duplicates) > 3 else ""
    raise WorkflowStageError(
        f"{context} contains duplicate site identifiers: {preview}{suffix}"
    )


def _resolve_site_to_protein(
    *,
    site_index: pd.Index,
    site_to_protein: pd.Series,
) -> pd.Series:
    resolved = site_to_protein.copy()
    resolved.index = pd.Index(resolved.index.astype(str), name=SITE_KEY_COLUMN)
    missing = [site_id for site_id in site_index if site_id not in resolved.index]
    if missing:
        preview = ", ".join(missing[:3])
        suffix = "..." if len(missing) > 3 else ""
        raise WorkflowStageError(
            f"site-to-protein mapping is missing prediction sites: {preview}{suffix}"
        )
    aligned = resolved.loc[site_index].astype(str).str.strip()
    if (aligned == "").any():
        raise WorkflowStageError(
            "site-to-protein mapping contains empty protein identifiers"
        )
    aligned.index = site_index.copy()
    aligned.name = PROTEIN_COLUMN
    return aligned


def _resolve_site_identity_columns(
    *,
    site_index: pd.Index,
    resolved_site_to_protein: pd.Series,
    site_metadata: pd.DataFrame | None,
) -> pd.DataFrame:
    fallback_display = pd.Series(
        site_index.astype(str).tolist(),
        index=site_index.copy(),
        dtype=object,
        name=DISPLAY_ID_COLUMN,
    )
    empty = pd.Series(
        [""] * int(site_index.size),
        index=site_index.copy(),
        dtype=object,
    )
    if site_metadata is None:
        parsed_identity = [
            _parse_display_site_identity(display_id)
            for display_id in fallback_display.astype(str).tolist()
        ]
        parsed_gene_symbols = [gene_symbol for gene_symbol, _ in parsed_identity]
        parsed_sites = [site for _, site in parsed_identity]
        protein_values = resolved_site_to_protein.astype(str).tolist()
        identity = pd.DataFrame(
            {
                SITE_KEY_COLUMN: site_index.astype(str).tolist(),
                DISPLAY_ID_COLUMN: fallback_display.astype(str).tolist(),
                GENE_SYMBOL_COLUMN: [
                    gene_symbol if gene_symbol != "" else protein_id
                    for gene_symbol, protein_id in zip(
                        parsed_gene_symbols,
                        protein_values,
                        strict=True,
                    )
                ],
                SITE_COLUMN: parsed_sites,
                PROTEIN_ACCESSION_COLUMN: empty.astype(str).tolist(),
                ISOFORM_ID_COLUMN: empty.astype(str).tolist(),
            },
            index=site_index.copy(),
        )
    else:
        metadata = site_metadata.copy(deep=False)
        metadata.index = pd.Index(metadata.index.astype(str), name=SITE_KEY_COLUMN)
        if (
            not site_index.isin(metadata.index).all()
            and SITE_KEY_COLUMN in metadata.columns
        ):
            site_key_index = (
                metadata.loc[:, SITE_KEY_COLUMN].fillna("").astype(str).str.strip()
            )
            if not (site_key_index == "").any():
                remapped_metadata = metadata.copy(deep=False)
                remapped_metadata.index = pd.Index(
                    site_key_index.tolist(),
                    name=SITE_KEY_COLUMN,
                )
                if site_index.isin(remapped_metadata.index).all():
                    metadata = remapped_metadata
        aligned_metadata = metadata.reindex(site_index)

        def _resolve_column(
            column_name: str,
            *,
            fallback: pd.Series | None = None,
        ) -> pd.Series:
            if column_name in aligned_metadata.columns:
                return (
                    aligned_metadata.loc[:, column_name]
                    .fillna("")
                    .astype(str)
                    .str.strip()
                    .rename(column_name)
                )
            if fallback is not None:
                return fallback.astype(str).rename(column_name)
            return empty.astype(str).rename(column_name)

        identity = pd.DataFrame(
            {
                SITE_KEY_COLUMN: site_index.astype(str).tolist(),
                DISPLAY_ID_COLUMN: _resolve_column(
                    DISPLAY_ID_COLUMN, fallback=fallback_display
                ).tolist(),
                GENE_SYMBOL_COLUMN: _resolve_column(GENE_SYMBOL_COLUMN).tolist(),
                SITE_COLUMN: _resolve_column(SITE_COLUMN).tolist(),
                PROTEIN_ACCESSION_COLUMN: _resolve_column(
                    PROTEIN_ACCESSION_COLUMN
                ).tolist(),
                ISOFORM_ID_COLUMN: _resolve_column(ISOFORM_ID_COLUMN).tolist(),
            },
            index=site_index.copy(),
        )
    identity.loc[:, PROTEIN_COLUMN] = resolved_site_to_protein.astype(str).tolist()
    blank_display_mask = (
        identity.loc[:, DISPLAY_ID_COLUMN].astype(str).str.strip() == ""
    )
    if bool(blank_display_mask.any()):
        identity.loc[blank_display_mask, DISPLAY_ID_COLUMN] = (
            fallback_display.loc[blank_display_mask].astype(str).tolist()
        )
    identity = identity.astype(
        {
            SITE_KEY_COLUMN: str,
            DISPLAY_ID_COLUMN: str,
            GENE_SYMBOL_COLUMN: str,
            SITE_COLUMN: str,
            PROTEIN_COLUMN: str,
            PROTEIN_ACCESSION_COLUMN: str,
            ISOFORM_ID_COLUMN: str,
        }
    )
    return identity


def _parse_display_site_identity(display_id: str) -> tuple[str, str]:
    text = str(display_id).strip()
    if text.count(";") < 2:
        return "", ""
    tokens = [token.strip() for token in text.split(";")]
    if len(tokens) < 2:
        return "", ""
    gene_symbol = tokens[0]
    site = tokens[1]
    if gene_symbol == "" or site == "":
        return "", ""
    return gene_symbol, site


def _resolve_site_module_resolution(
    *,
    top_kinases: pd.Series,
    top_kinase_weights: Sequence[tuple[tuple[str, float], ...]],
    site_to_protein: pd.Series,
    protein_modules: pd.Series | None,
) -> pd.DataFrame:
    if protein_modules is None:
        protein_resolution = _derive_protein_modules_by_top_kinase(
            top_kinases=top_kinases,
            top_kinase_weights=top_kinase_weights,
            site_to_protein=site_to_protein,
        )
    else:
        protein_resolution = _build_protein_resolution_from_modules(
            top_kinases=top_kinases,
            top_kinase_weights=top_kinase_weights,
            site_to_protein=site_to_protein,
            protein_modules=protein_modules,
        )
    site_proteins = site_to_protein.to_numpy(dtype=object, copy=False)
    site_module_resolution = protein_resolution.loc[
        site_proteins,
        [
            MODULE_ID_COLUMN,
            MODULE_TOP_KINASE_COLUMN,
            MODULE_TOP_KINASE_CANDIDATES_COLUMN,
            MODULE_TOP_KINASE_TIE_COUNT_COLUMN,
            MODULE_TOP_KINASE_IS_AMBIGUOUS_COLUMN,
            MODULE_TOP_KINASE_SELECTION_POLICY_COLUMN,
        ],
    ]
    site_module_resolution.index = site_to_protein.index.copy()
    return site_module_resolution


def _build_protein_resolution_from_modules(
    *,
    top_kinases: pd.Series,
    top_kinase_weights: Sequence[tuple[tuple[str, float], ...]],
    site_to_protein: pd.Series,
    protein_modules: pd.Series,
) -> pd.DataFrame:
    if len(top_kinases) != len(top_kinase_weights):
        raise WorkflowStageError(
            "top kinase weights must align one-to-one with prediction sites"
        )
    normalized_modules = _normalize_protein_modules(protein_modules)
    protein_index = pd.Index(
        sorted(set(site_to_protein.astype(str).tolist())),
        name=PROTEIN_COLUMN,
        dtype=object,
    )
    module_ids = pd.Series(
        np.zeros(len(protein_index), dtype=np.int64),
        index=protein_index.copy(),
        dtype="int64",
    )
    shared_proteins = protein_index.intersection(normalized_modules.index)
    if not shared_proteins.empty:
        module_ids.loc[shared_proteins] = normalized_modules.loc[shared_proteins]

    top_table = pd.DataFrame(
        {
            PROTEIN_COLUMN: site_to_protein.to_numpy(dtype=object, copy=False),
            TOP_KINASE_COLUMN: top_kinases.to_numpy(dtype=object, copy=False),
            TOP_KINASE_WEIGHTS_COLUMN: list(top_kinase_weights),
        },
        index=site_to_protein.index.copy(),
    )
    top_table.loc[:, MODULE_ID_COLUMN] = (
        top_table.loc[:, PROTEIN_COLUMN]
        .astype(str)
        .map(module_ids)
        .fillna(0)
        .astype("int64")
    )
    module_resolution = _derive_module_top_kinase_resolution(top_table)
    protein_resolution = pd.DataFrame(
        {
            MODULE_ID_COLUMN: module_ids.to_numpy(dtype=np.int64, copy=False),
        },
        index=protein_index.copy(),
    )
    return protein_resolution.join(module_resolution, on=MODULE_ID_COLUMN).astype(
        {
            MODULE_ID_COLUMN: "int64",
            MODULE_TOP_KINASE_COLUMN: str,
            MODULE_TOP_KINASE_TIE_COUNT_COLUMN: "int64",
            MODULE_TOP_KINASE_IS_AMBIGUOUS_COLUMN: bool,
            MODULE_TOP_KINASE_SELECTION_POLICY_COLUMN: str,
        }
    )


def _normalize_protein_modules(protein_modules: pd.Series) -> pd.Series:
    resolved = protein_modules.copy()
    resolved.index = pd.Index(resolved.index.astype(str), name=PROTEIN_COLUMN)
    if resolved.index.has_duplicates:
        duplicated = sorted(
            {str(value) for value in resolved.index[resolved.index.duplicated()]}
        )
        preview = ", ".join(duplicated[:3])
        suffix = "..." if len(duplicated) > 3 else ""
        raise WorkflowStageError(
            f"protein_modules contains duplicate protein identifiers: {preview}{suffix}"
        )
    module_values = pd.to_numeric(resolved, errors="coerce")
    if module_values.isna().any():
        raise WorkflowStageError("protein_modules must contain integer module IDs")
    rounded = np.floor(module_values.to_numpy(dtype=float, copy=False))
    if not np.allclose(module_values.to_numpy(dtype=float, copy=False), rounded):
        raise WorkflowStageError("protein_modules must contain integer module IDs")
    integer_values = module_values.astype("int64")
    integer_values.loc[integer_values < 0] = 0
    return integer_values.astype("int64")


def _derive_protein_modules_by_top_kinase(
    *,
    top_kinases: pd.Series,
    top_kinase_weights: Sequence[tuple[tuple[str, float], ...]],
    site_to_protein: pd.Series,
) -> pd.DataFrame:
    if len(top_kinases) != len(top_kinase_weights):
        raise WorkflowStageError(
            "top kinase weights must align one-to-one with prediction sites"
        )
    top_table = pd.DataFrame(
        {
            PROTEIN_COLUMN: site_to_protein.to_numpy(dtype=object, copy=False),
            TOP_KINASE_COLUMN: top_kinases.to_numpy(dtype=object, copy=False),
            TOP_KINASE_WEIGHTS_COLUMN: list(top_kinase_weights),
        },
        index=site_to_protein.index.copy(),
    )

    protein_resolution_rows: list[dict[str, object]] = []
    for protein_id, group in top_table.groupby(PROTEIN_COLUMN, sort=True):
        supported_group = group.loc[
            group.loc[:, TOP_KINASE_COLUMN].astype(str) != UNSUPPORTED_KINASE
        ]
        counts = supported_group.loc[:, TOP_KINASE_COLUMN].astype(str).value_counts()
        if not counts.empty:
            max_count = int(counts.iloc[0])
            tied_kinases = tuple(
                sorted(kinase for kinase in counts[counts == max_count].index.to_list())
            )
            dominant_kinase = tied_kinases[0]
            selection_policy = LEXICOGRAPHIC_TIE_BREAK_POLICY
        else:
            tied_kinases = ()
            dominant_kinase = UNSUPPORTED_KINASE
            selection_policy = NO_SUPPORT_SELECTION_POLICY
        protein_resolution_rows.append(
            {
                PROTEIN_COLUMN: str(protein_id),
                MODULE_TOP_KINASE_COLUMN: dominant_kinase,
                MODULE_TOP_KINASE_CANDIDATES_COLUMN: tied_kinases,
                MODULE_TOP_KINASE_TIE_COUNT_COLUMN: len(tied_kinases),
                MODULE_TOP_KINASE_IS_AMBIGUOUS_COLUMN: len(tied_kinases) > 1,
                MODULE_TOP_KINASE_SELECTION_POLICY_COLUMN: selection_policy,
            }
        )

    protein_resolution = pd.DataFrame(protein_resolution_rows).set_index(PROTEIN_COLUMN)
    protein_resolution.index = pd.Index(
        protein_resolution.index.astype(str), name=PROTEIN_COLUMN
    )
    dominant_kinases = protein_resolution.loc[:, MODULE_TOP_KINASE_COLUMN].astype(str)
    supported_mask = dominant_kinases != UNSUPPORTED_KINASE
    module_by_kinase = {
        kinase: module_id
        for module_id, kinase in enumerate(
            sorted(set(dominant_kinases.loc[supported_mask].tolist())), start=1
        )
    }
    module_ids = pd.Series(
        np.zeros(len(protein_resolution), dtype=np.int64),
        index=protein_resolution.index.copy(),
        dtype="int64",
    )
    if supported_mask.any():
        module_ids.loc[supported_mask] = (
            dominant_kinases.loc[supported_mask].map(module_by_kinase).astype("int64")
        )
    protein_resolution.loc[:, MODULE_ID_COLUMN] = module_ids.to_numpy(
        dtype=np.int64, copy=False
    )
    return protein_resolution.astype(
        {
            MODULE_TOP_KINASE_COLUMN: str,
            MODULE_TOP_KINASE_TIE_COUNT_COLUMN: "int64",
            MODULE_TOP_KINASE_IS_AMBIGUOUS_COLUMN: bool,
            MODULE_TOP_KINASE_SELECTION_POLICY_COLUMN: str,
            MODULE_ID_COLUMN: "int64",
        }
    )


def _derive_module_top_kinase_resolution(top_table: pd.DataFrame) -> pd.DataFrame:
    module_resolution_rows: list[dict[str, object]] = []
    module_ids = sorted(
        {int(module_id) for module_id in top_table.loc[:, MODULE_ID_COLUMN].tolist()}
    )
    if 0 not in module_ids:
        module_ids = [0, *module_ids]
    for module_id in module_ids:
        module_group = top_table.loc[
            top_table.loc[:, MODULE_ID_COLUMN].astype("int64") == int(module_id)
        ]
        supported_group = module_group.loc[
            module_group.loc[:, TOP_KINASE_COLUMN].astype(str) != UNSUPPORTED_KINASE
        ]
        counts = supported_group.loc[:, TOP_KINASE_COLUMN].astype(str).value_counts()
        if not counts.empty and int(module_id) > 0:
            max_count = int(counts.iloc[0])
            tied_kinases = tuple(
                sorted(kinase for kinase in counts[counts == max_count].index.to_list())
            )
            dominant_kinase = tied_kinases[0]
            selection_policy = LEXICOGRAPHIC_TIE_BREAK_POLICY
        else:
            tied_kinases = ()
            dominant_kinase = UNSUPPORTED_KINASE
            selection_policy = NO_SUPPORT_SELECTION_POLICY
        module_resolution_rows.append(
            {
                MODULE_ID_COLUMN: int(module_id),
                MODULE_TOP_KINASE_COLUMN: dominant_kinase,
                MODULE_TOP_KINASE_CANDIDATES_COLUMN: tied_kinases,
                MODULE_TOP_KINASE_TIE_COUNT_COLUMN: len(tied_kinases),
                MODULE_TOP_KINASE_IS_AMBIGUOUS_COLUMN: len(tied_kinases) > 1,
                MODULE_TOP_KINASE_SELECTION_POLICY_COLUMN: selection_policy,
            }
        )
    module_resolution = pd.DataFrame(module_resolution_rows).set_index(MODULE_ID_COLUMN)
    module_resolution.index = pd.Index(
        module_resolution.index.astype("int64"),
        name=MODULE_ID_COLUMN,
    )
    return module_resolution


__all__ = [
    "build_module_assignments",
    "select_kinase_substrates",
]
