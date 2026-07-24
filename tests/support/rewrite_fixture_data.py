from __future__ import annotations

import ast
import csv
import json
from functools import lru_cache
from pathlib import Path

import pandas as pd

from phospy import (
    AnalysisReadyDatasetBuilder,
    AnalysisReadyPhosphoDataset,
)
from phospy.api import (
    DatasetBuildRequest,
    Organism,
)

ROOT = Path(__file__).resolve().parents[2]
REWRITE_PARITY_REFERENCE = (
    ROOT / "tests" / "fixtures" / "rewrite_parity" / "r_reference_l6"
)
REWRITE_PARITY_FRAGILE_SUPPORT_REFERENCE = (
    ROOT / "tests" / "fixtures" / "rewrite_parity" / "fragile_support_reference"
)
REWRITE_PARITY_FRAGILE_SUPPORT_MOTIF_MATRICES = (
    REWRITE_PARITY_FRAGILE_SUPPORT_REFERENCE / "motif_frequency_matrices"
)
REWRITE_PARITY_L6_PREDICTION_REFERENCE = (
    ROOT / "tests" / "fixtures" / "rewrite_parity" / "r_reference_l6_prediction"
)
REWRITE_PARITY_ADAPTIVE_SAMPLING_EDGE = (
    ROOT / "tests" / "fixtures" / "rewrite_parity" / "adaptive_sampling_edge"
)
REWRITE_PARITY_ADAPTIVE_SAMPLING_REPLAY = (
    ROOT / "tests" / "fixtures" / "rewrite_parity" / "adaptive_sampling_replay"
)
ADAPTIVE_SAMPLING_EDGE_TRACE_PREDICTIONS = (
    REWRITE_PARITY_ADAPTIVE_SAMPLING_EDGE / "trace_final_ensemble_predictions.csv"
)
ADAPTIVE_SAMPLING_EDGE_TRACE_TOP = (
    REWRITE_PARITY_ADAPTIVE_SAMPLING_EDGE / "trace_final_ensemble_top.csv"
)
ADAPTIVE_SAMPLING_EDGE_TRACE_CANDIDATES = (
    REWRITE_PARITY_ADAPTIVE_SAMPLING_EDGE / "trace_candidates.csv"
)
RAT_L6_PHOSPHO = REWRITE_PARITY_REFERENCE / "l6_phospho_matrix.csv"
RAT_L6_EXPECTED_PROFILE = REWRITE_PARITY_REFERENCE / "native_profile_scores.csv"
ACTIVITY_REFERENCE_PROVENANCE = REWRITE_PARITY_REFERENCE / "PROVENANCE.md"
ACTIVITY_REFERENCE_PREDMAT = REWRITE_PARITY_REFERENCE / "predMat.csv"
ACTIVITY_REFERENCE_ACTIVITY_MATRIX = (
    REWRITE_PARITY_REFERENCE / "kinase_activity_matrix.csv"
)
ACTIVITY_REFERENCE_THRESHOLDED_SUBSTRATE_MEAN_ACTIVITY = (
    REWRITE_PARITY_REFERENCE / "thresholded_substrate_mean_activity.csv"
)
ACTIVITY_REFERENCE_THRESHOLDED_SUBSTRATE_COUNTS = (
    REWRITE_PARITY_REFERENCE / "thresholded_substrate_counts.csv"
)
ACTIVITY_REFERENCE_TARGET_COUNTS = REWRITE_PARITY_REFERENCE / "kinase_target_counts.csv"
ACTIVITY_REFERENCE_TARGET_TABLE = REWRITE_PARITY_REFERENCE / "kinase_target_table.csv"
ACTIVITY_PARITY_FIXTURE_FILES: tuple[str, ...] = (
    "l6_phospho_matrix.csv",
    "native_profile_scores.csv",
    "predMat.csv",
    "kinase_activity_matrix.csv",
    "thresholded_substrate_mean_activity.csv",
    "thresholded_substrate_counts.csv",
    "kinase_target_counts.csv",
    "kinase_target_table.csv",
)
L6_PREDICTION_REFERENCE_PHOSPHO = (
    REWRITE_PARITY_L6_PREDICTION_REFERENCE / "l6_phospho_matrix.csv"
)
L6_PREDICTION_REFERENCE_PROFILE = (
    REWRITE_PARITY_L6_PREDICTION_REFERENCE / "native_profile_scores.csv"
)
L6_PREDICTION_REFERENCE_RANK_WEIGHTED_FUSION_SCORES = (
    REWRITE_PARITY_L6_PREDICTION_REFERENCE / "native_rank_weighted_fusion_scores.csv"
)
L6_PREDICTION_REFERENCE_SCORE_FUSION_WEIGHTS = (
    REWRITE_PARITY_L6_PREDICTION_REFERENCE / "native_score_fusion_weights.csv"
)
L6_PREDICTION_REFERENCE_CANDIDATES = (
    REWRITE_PARITY_L6_PREDICTION_REFERENCE / "native_candidate_substrates.csv"
)
L6_PREDICTION_REFERENCE_TOP30 = (
    REWRITE_PARITY_L6_PREDICTION_REFERENCE / "native_prediction_top30.csv"
)
L6_PREDICTION_REFERENCE_PREDMAT = REWRITE_PARITY_L6_PREDICTION_REFERENCE / "predMat.csv"
ADAPTIVE_SAMPLING_REPLAY_RANK_WEIGHTED_FUSION_SCORES = (
    REWRITE_PARITY_ADAPTIVE_SAMPLING_REPLAY / "rank_weighted_fusion_scores.csv"
)
REWRITE_PUBLIC_WORKFLOW_REFERENCE = (
    ROOT / "tests" / "fixtures" / "public_workflow_reference"
)
PUBLIC_PREDMAT_INPUT_PHOSPHO = (
    REWRITE_PUBLIC_WORKFLOW_REFERENCE / "predmat_input_phospho_matrix.csv"
)
PUBLIC_PREDMAT_INPUT_SUBSTRATE_MAP = (
    REWRITE_PUBLIC_WORKFLOW_REFERENCE / "predmat_input_substrate_map.json"
)
PUBLIC_PREDMAT_INPUT_SITE_SEQUENCES = (
    REWRITE_PUBLIC_WORKFLOW_REFERENCE / "predmat_input_site_sequences.json"
)
PUBLIC_PREDMAT_REWRITE_STABLE = (
    REWRITE_PUBLIC_WORKFLOW_REFERENCE / "predmat_rewrite_stable.csv"
)
PUBLIC_PREDMAT_REWRITE_R_PARITY = (
    REWRITE_PUBLIC_WORKFLOW_REFERENCE / "predmat_rewrite_r_parity.csv"
)
PUBLIC_PREDMAT_REWRITE_CONTRACT = (
    REWRITE_PUBLIC_WORKFLOW_REFERENCE / "predmat_rewrite_contract.json"
)
PUBLIC_PREDMAT_HISTORICAL_BASELINE_DEFAULT = (
    REWRITE_PUBLIC_WORKFLOW_REFERENCE / "predmat_historical_baseline_default.csv"
)
PUBLIC_PREDMAT_HISTORICAL_BASELINE_R_PARITY = (
    REWRITE_PUBLIC_WORKFLOW_REFERENCE / "predmat_historical_baseline_r_parity.csv"
)
RAT_L6_SITE_SEQUENCES = (
    ROOT
    / "src"
    / "phospy"
    / "data"
    / "reference_bundles"
    / "rat"
    / "l6_native"
    / "site_sequences.csv"
)
SIGNALOME_REWRITE_L6_ASSIGNMENTS = (
    REWRITE_PUBLIC_WORKFLOW_REFERENCE / "signalome_rewrite_l6_module_assignments.csv"
)
SIGNALOME_REWRITE_L6_MODULES = (
    REWRITE_PUBLIC_WORKFLOW_REFERENCE / "signalome_rewrite_l6_modules.csv"
)
SIGNALOME_REWRITE_L6_NETWORK_NODES = (
    REWRITE_PUBLIC_WORKFLOW_REFERENCE / "signalome_rewrite_l6_network_nodes.csv"
)
SIGNALOME_REWRITE_L6_NETWORK_EDGES = (
    REWRITE_PUBLIC_WORKFLOW_REFERENCE / "signalome_rewrite_l6_network_edges.csv"
)
SIGNALOME_REWRITE_L6_EXPANDED_SIGNALOME = (
    REWRITE_PUBLIC_WORKFLOW_REFERENCE / "signalome_rewrite_l6_expanded_signalome.csv"
)
SIGNALOME_REWRITE_L6_CONTRACT = (
    REWRITE_PUBLIC_WORKFLOW_REFERENCE / "signalome_rewrite_l6_contract.json"
)
KINASE_PUBLIC_PREDMAT_PROVENANCE_GOLDEN = (
    REWRITE_PUBLIC_WORKFLOW_REFERENCE / "kinase_public_predmat_provenance_golden.json"
)
SIGNALOME_L6_PROVENANCE_GOLDEN = (
    REWRITE_PUBLIC_WORKFLOW_REFERENCE / "signalome_l6_provenance_golden.json"
)


@lru_cache(maxsize=1)
def load_rat_l6_phospho() -> pd.DataFrame:
    return _read_float_matrix_csv(RAT_L6_PHOSPHO)


@lru_cache(maxsize=1)
def load_l6_prediction_reference_phospho() -> pd.DataFrame:
    return pd.read_csv(L6_PREDICTION_REFERENCE_PHOSPHO, index_col=0)


@lru_cache(maxsize=1)
def load_rat_l6_sequence_table() -> pd.Series:
    sequence_frame = pd.read_csv(RAT_L6_SITE_SEQUENCES)
    return sequence_frame.set_index("site_id").loc[:, "centralized_sequence"]


@lru_cache(maxsize=1)
def load_expected_profile_scores() -> pd.DataFrame:
    return pd.read_csv(RAT_L6_EXPECTED_PROFILE, index_col=0)


@lru_cache(maxsize=1)
def load_activity_reference_predmat() -> pd.DataFrame:
    return pd.read_csv(ACTIVITY_REFERENCE_PREDMAT, index_col=0)


@lru_cache(maxsize=1)
def load_activity_reference_activity_matrix() -> pd.DataFrame:
    frame = pd.read_csv(ACTIVITY_REFERENCE_ACTIVITY_MATRIX, index_col=0)
    frame.index = pd.Index(frame.index.astype(str), name="kinase")
    return frame.astype(float)


@lru_cache(maxsize=1)
def load_activity_reference_thresholded_substrate_mean_activity() -> pd.DataFrame:
    frame = pd.read_csv(
        ACTIVITY_REFERENCE_THRESHOLDED_SUBSTRATE_MEAN_ACTIVITY,
        index_col=0,
    )
    frame.index = pd.Index(frame.index.astype(str), name="kinase")
    return frame.astype(float)


@lru_cache(maxsize=1)
def load_activity_reference_thresholded_substrate_counts() -> pd.Series:
    frame = pd.read_csv(
        ACTIVITY_REFERENCE_THRESHOLDED_SUBSTRATE_COUNTS,
        index_col=0,
    )
    counts = frame.iloc[:, 0].astype(int)
    counts.index = pd.Index(counts.index.astype(str), name="kinase")
    counts.name = "n_substrates"
    return counts


@lru_cache(maxsize=1)
def load_activity_reference_target_counts() -> pd.Series:
    frame = pd.read_csv(ACTIVITY_REFERENCE_TARGET_COUNTS, index_col=0)
    counts = frame.iloc[:, 0].astype(int)
    counts.index = pd.Index(counts.index.astype(str), name="kinase")
    counts.name = "n_targets"
    return counts


@lru_cache(maxsize=1)
def load_activity_reference_target_table() -> pd.DataFrame:
    return pd.read_csv(ACTIVITY_REFERENCE_TARGET_TABLE).astype(
        {"site_id": str, "kinase": str, "score": float}
    )


@lru_cache(maxsize=1)
def load_activity_reference_provenance_text() -> str:
    return ACTIVITY_REFERENCE_PROVENANCE.read_text(encoding="utf-8")


def activity_parity_fixture_paths() -> tuple[Path, ...]:
    return tuple(
        REWRITE_PARITY_REFERENCE / name for name in ACTIVITY_PARITY_FIXTURE_FILES
    )


@lru_cache(maxsize=1)
def load_l6_prediction_reference_profile_scores() -> pd.DataFrame:
    return pd.read_csv(L6_PREDICTION_REFERENCE_PROFILE, index_col=0)


@lru_cache(maxsize=1)
def load_l6_prediction_reference_rank_weighted_fusion_scores() -> pd.DataFrame:
    return pd.read_csv(
        L6_PREDICTION_REFERENCE_RANK_WEIGHTED_FUSION_SCORES,
        index_col=0,
    )


@lru_cache(maxsize=1)
def load_l6_prediction_reference_score_fusion_weights() -> pd.DataFrame:
    frame = pd.read_csv(L6_PREDICTION_REFERENCE_SCORE_FUSION_WEIGHTS).set_index(
        "kinase"
    )
    frame.index = frame.index.astype(str)
    frame.index.name = "kinase"
    return frame


@lru_cache(maxsize=1)
def load_l6_prediction_reference_candidate_substrates() -> pd.DataFrame:
    return pd.read_csv(L6_PREDICTION_REFERENCE_CANDIDATES).astype(
        {"kinase": str, "site_id": str}
    )


@lru_cache(maxsize=1)
def load_l6_prediction_reference_top30() -> pd.DataFrame:
    return pd.read_csv(L6_PREDICTION_REFERENCE_TOP30).astype(
        {"kinase": str, "rank": int, "site_id": str}
    )


@lru_cache(maxsize=1)
def load_l6_prediction_reference_predmat() -> pd.DataFrame:
    return pd.read_csv(L6_PREDICTION_REFERENCE_PREDMAT, index_col=0)


@lru_cache(maxsize=1)
def load_fragile_support_profile_scores() -> pd.DataFrame:
    return pd.read_csv(
        REWRITE_PARITY_FRAGILE_SUPPORT_REFERENCE / "profile_scores.csv",
        index_col=0,
    )


@lru_cache(maxsize=1)
def load_fragile_support_motif_scores() -> pd.DataFrame:
    return pd.read_csv(
        REWRITE_PARITY_FRAGILE_SUPPORT_REFERENCE / "motif_scores.csv",
        index_col=0,
    )


@lru_cache(maxsize=1)
def load_fragile_support_motif_scores_full() -> pd.DataFrame:
    return pd.read_csv(
        REWRITE_PARITY_FRAGILE_SUPPORT_REFERENCE / "motif_scores_full.csv",
        index_col=0,
    )


@lru_cache(maxsize=1)
def load_fragile_support_motif_site_sequences_full() -> pd.Series:
    frame = pd.read_csv(
        REWRITE_PARITY_FRAGILE_SUPPORT_REFERENCE / "motif_site_sequences_full.csv"
    ).astype({"site_id": str, "centralized_sequence": str})
    sequence_series = frame.set_index("site_id").loc[:, "centralized_sequence"]
    sequence_series.index = pd.Index(sequence_series.index.astype(str), name="site_id")
    sequence_series.name = "centralized_sequence"
    return sequence_series


@lru_cache(maxsize=1)
def load_fragile_support_motif_frequency_matrices() -> dict[str, pd.DataFrame]:
    matrices: dict[str, pd.DataFrame] = {}
    for path in sorted(REWRITE_PARITY_FRAGILE_SUPPORT_MOTIF_MATRICES.glob("*.csv")):
        matrices[path.stem] = pd.read_csv(path, index_col=0).astype(float)
    return matrices


@lru_cache(maxsize=1)
def load_fragile_support_profile_sizes() -> pd.Series:
    frame = pd.read_csv(
        REWRITE_PARITY_FRAGILE_SUPPORT_REFERENCE / "profile_sizes.csv",
        index_col=0,
    )
    return frame.iloc[:, 0].astype(float)


@lru_cache(maxsize=1)
def load_fragile_support_motif_sizes() -> pd.Series:
    frame = pd.read_csv(
        REWRITE_PARITY_FRAGILE_SUPPORT_REFERENCE / "motif_sizes.csv",
        index_col=0,
    )
    return frame.iloc[:, 0].astype(float)


@lru_cache(maxsize=1)
def load_fragile_support_rank_weighted_fusion_scores() -> pd.DataFrame:
    return pd.read_csv(
        REWRITE_PARITY_FRAGILE_SUPPORT_REFERENCE / "rank_weighted_fusion_scores.csv",
        index_col=0,
    )


@lru_cache(maxsize=1)
def load_fragile_support_score_fusion_weights() -> pd.DataFrame:
    frame = pd.read_csv(
        REWRITE_PARITY_FRAGILE_SUPPORT_REFERENCE / "score_fusion_weights.csv",
        index_col=0,
    )
    frame.index = frame.index.astype(str)
    frame.index.name = "kinase"
    return frame


@lru_cache(maxsize=1)
def load_fragile_support_candidate_substrates() -> pd.DataFrame:
    return pd.read_csv(
        REWRITE_PARITY_FRAGILE_SUPPORT_REFERENCE / "candidate_substrates.csv"
    )


@lru_cache(maxsize=1)
def load_adaptive_sampling_edge_rank_weighted_fusion_scores() -> pd.DataFrame:
    return pd.read_csv(
        REWRITE_PARITY_ADAPTIVE_SAMPLING_EDGE / "rank_weighted_fusion_scores.csv",
        index_col=0,
    )


@lru_cache(maxsize=1)
def load_adaptive_sampling_replay_rank_weighted_fusion_scores() -> pd.DataFrame:
    return pd.read_csv(
        ADAPTIVE_SAMPLING_REPLAY_RANK_WEIGHTED_FUSION_SCORES,
        index_col=0,
    )


@lru_cache(maxsize=1)
def load_adaptive_sampling_replay_trace_candidates() -> pd.DataFrame:
    return pd.read_csv(REWRITE_PARITY_ADAPTIVE_SAMPLING_REPLAY / "trace_candidates.csv")


@lru_cache(maxsize=1)
def load_adaptive_sampling_replay_trace_initial_negatives() -> pd.DataFrame:
    return pd.read_csv(
        REWRITE_PARITY_ADAPTIVE_SAMPLING_REPLAY / "trace_initial_negatives.csv"
    )


@lru_cache(maxsize=1)
def load_adaptive_sampling_replay_trace_iteration_samples() -> pd.DataFrame:
    return pd.read_csv(
        REWRITE_PARITY_ADAPTIVE_SAMPLING_REPLAY / "trace_iteration_samples.csv"
    )


@lru_cache(maxsize=1)
def load_adaptive_sampling_replay_trace_final_predictions() -> pd.DataFrame:
    return pd.read_csv(
        REWRITE_PARITY_ADAPTIVE_SAMPLING_REPLAY / "trace_final_ensemble_predictions.csv"
    )


@lru_cache(maxsize=1)
def load_adaptive_sampling_replay_trace_final_top() -> pd.DataFrame:
    return pd.read_csv(
        REWRITE_PARITY_ADAPTIVE_SAMPLING_REPLAY / "trace_final_ensemble_top.csv"
    )


@lru_cache(maxsize=1)
def load_adaptive_sampling_edge_trace_predictions() -> pd.DataFrame:
    frame = pd.read_csv(ADAPTIVE_SAMPLING_EDGE_TRACE_PREDICTIONS).loc[
        :, ["kinase", "site", "prob_class_1"]
    ]
    return frame.astype({"kinase": str, "site": str, "prob_class_1": float})


@lru_cache(maxsize=1)
def load_adaptive_sampling_edge_trace_top() -> pd.DataFrame:
    frame = pd.read_csv(ADAPTIVE_SAMPLING_EDGE_TRACE_TOP).loc[
        :, ["kinase", "rank", "site", "prob_class_1"]
    ]
    return frame.astype(
        {
            "kinase": str,
            "rank": int,
            "site": str,
            "prob_class_1": float,
        }
    )


@lru_cache(maxsize=1)
def load_adaptive_sampling_edge_trace_candidates() -> pd.DataFrame:
    frame = pd.read_csv(ADAPTIVE_SAMPLING_EDGE_TRACE_CANDIDATES)
    expected_columns = {"kinase", "site", "selected_candidate"}
    if expected_columns.issubset(frame.columns):
        selected = frame.loc[:, "selected_candidate"].map(
            lambda value: str(value).strip().lower() == "true"
        )
        return pd.DataFrame(
            {
                "kinase": frame.loc[:, "kinase"].astype(str).values,
                "site": frame.loc[:, "site"].astype(str).values,
                "selected_candidate": selected.values,
            }
        )
    return frame


@lru_cache(maxsize=1)
def load_signalome_rewrite_l6_contract() -> dict[str, object]:
    return json.loads(SIGNALOME_REWRITE_L6_CONTRACT.read_text(encoding="utf-8"))


@lru_cache(maxsize=1)
def load_signalome_l6_provenance_golden() -> dict[str, object]:
    return json.loads(SIGNALOME_L6_PROVENANCE_GOLDEN.read_text(encoding="utf-8"))


@lru_cache(maxsize=1)
def load_public_predmat_input_phospho() -> pd.DataFrame:
    frame = _read_float_matrix_csv(PUBLIC_PREDMAT_INPUT_PHOSPHO)
    frame.index = pd.Index(frame.index.astype(str), name="phosphosite")
    return frame


@lru_cache(maxsize=1)
def load_public_predmat_input_substrate_map() -> dict[str, list[str]]:
    raw = json.loads(PUBLIC_PREDMAT_INPUT_SUBSTRATE_MAP.read_text(encoding="utf-8"))
    return {
        str(kinase): [str(site_id) for site_id in site_ids]
        for kinase, site_ids in raw.items()
    }


@lru_cache(maxsize=1)
def load_public_predmat_input_site_sequences() -> dict[str, str]:
    raw = json.loads(PUBLIC_PREDMAT_INPUT_SITE_SEQUENCES.read_text(encoding="utf-8"))
    return {str(site_id): str(sequence) for site_id, sequence in raw.items()}


@lru_cache(maxsize=1)
def load_public_predmat_rewrite_contract() -> dict[str, object]:
    return json.loads(PUBLIC_PREDMAT_REWRITE_CONTRACT.read_text(encoding="utf-8"))


@lru_cache(maxsize=1)
def load_kinase_public_predmat_provenance_golden() -> dict[str, object]:
    return json.loads(
        KINASE_PUBLIC_PREDMAT_PROVENANCE_GOLDEN.read_text(encoding="utf-8")
    )


@lru_cache(maxsize=1)
def load_public_predmat_rewrite_stable() -> pd.DataFrame:
    return pd.read_csv(PUBLIC_PREDMAT_REWRITE_STABLE, index_col=0)


@lru_cache(maxsize=1)
def load_public_predmat_rewrite_r_parity() -> pd.DataFrame:
    return pd.read_csv(PUBLIC_PREDMAT_REWRITE_R_PARITY, index_col=0)


@lru_cache(maxsize=1)
def load_public_predmat_historical_baseline_default() -> pd.DataFrame:
    return pd.read_csv(PUBLIC_PREDMAT_HISTORICAL_BASELINE_DEFAULT, index_col=0)


@lru_cache(maxsize=1)
def load_public_predmat_historical_baseline_r_parity() -> pd.DataFrame:
    return pd.read_csv(PUBLIC_PREDMAT_HISTORICAL_BASELINE_R_PARITY, index_col=0)


@lru_cache(maxsize=1)
def load_signalome_rewrite_l6_module_assignments() -> pd.DataFrame:
    frame = pd.read_csv(SIGNALOME_REWRITE_L6_ASSIGNMENTS, index_col=0)
    frame.index = pd.Index(frame.index.astype(str), name="site_id")
    return frame


@lru_cache(maxsize=1)
def load_signalome_rewrite_l6_modules() -> pd.DataFrame:
    frame = pd.read_csv(SIGNALOME_REWRITE_L6_MODULES, index_col=0)
    frame.index = pd.Index(frame.index.astype("int64"), name="module_id")
    frame.columns = pd.Index(frame.columns.astype(str), name="kinase")
    return frame.astype(float)


@lru_cache(maxsize=1)
def load_signalome_rewrite_l6_network_nodes() -> pd.DataFrame:
    frame = pd.read_csv(SIGNALOME_REWRITE_L6_NETWORK_NODES, index_col=0)
    frame.index = pd.Index(frame.index.astype(str), name="kinase")
    return frame.astype({"degree": "int64", "n_substrates": "int64"})


@lru_cache(maxsize=1)
def load_signalome_rewrite_l6_network_edges() -> pd.DataFrame:
    return pd.read_csv(SIGNALOME_REWRITE_L6_NETWORK_EDGES).astype(
        {
            "source_kinase": str,
            "target_kinase": str,
            "correlation": float,
            "valid_observations": "int64",
        }
    )


@lru_cache(maxsize=1)
def load_signalome_rewrite_l6_expanded_signalome() -> pd.DataFrame:
    return pd.read_csv(SIGNALOME_REWRITE_L6_EXPANDED_SIGNALOME).astype(
        {
            "kinase": str,
            "row_kind": str,
            "assignment_policy": str,
            "linked_kinases": str,
            "regulated_module_ids": str,
            "site_id": str,
            "site_order": "int64",
            "protein_id": str,
            "module_id": "int64",
            "support_kinases": str,
            "support_weight": float,
            "top_kinase": str,
            "top_score": float,
        }
    )


def _normalize_signalome_collection_value(value: object) -> str:
    parsed = value
    if value is None:
        return "[]"
    if isinstance(value, str):
        stripped = value.strip()
        if stripped == "":
            return "[]"
        try:
            parsed = ast.literal_eval(stripped)
        except (ValueError, SyntaxError):
            return stripped
    elif not isinstance(value, (tuple, list)) and pd.isna(value):
        return "[]"
    if isinstance(parsed, (tuple, list)):
        return json.dumps(
            [str(item) if item is not None else "" for item in parsed],
            ensure_ascii=True,
            separators=(",", ":"),
        )
    return str(parsed)


def _canonical_kinase_label(value: object) -> str:
    return str(value).strip().upper()


def _normalize_signalome_kinase_collection_value(value: object) -> str:
    parsed = value
    if value is None:
        return "[]"
    if isinstance(value, str):
        stripped = value.strip()
        if stripped == "":
            return "[]"
        try:
            parsed = ast.literal_eval(stripped)
        except (ValueError, SyntaxError):
            return _canonical_kinase_label(stripped)
    elif not isinstance(value, (tuple, list)) and pd.isna(value):
        return "[]"
    if isinstance(parsed, (tuple, list)):
        return json.dumps(
            [_canonical_kinase_label(item) for item in parsed if item is not None],
            ensure_ascii=True,
            separators=(",", ":"),
        )
    return _canonical_kinase_label(parsed)


def _normalize_signalome_kinase_weight_collection_value(value: object) -> str:
    parsed = value
    if value is None:
        return "[]"
    if isinstance(value, str):
        stripped = value.strip()
        if stripped == "":
            return "[]"
        try:
            parsed = ast.literal_eval(stripped)
        except (ValueError, SyntaxError):
            return stripped
    elif not isinstance(value, (tuple, list)) and pd.isna(value):
        return "[]"
    if isinstance(parsed, (tuple, list)):
        normalised_items: list[object] = []
        for item in parsed:
            if isinstance(item, tuple) and len(item) == 2:
                kinase, weight = item
                normalised_items.append((_canonical_kinase_label(kinase), weight))
            else:
                normalised_items.append(item)
        return json.dumps(
            [str(item) for item in normalised_items],
            ensure_ascii=True,
            separators=(",", ":"),
        )
    return str(parsed)


def normalize_signalome_modules_for_parity(frame: pd.DataFrame) -> pd.DataFrame:
    normalized = frame.copy(deep=True)
    normalized.index = pd.Index(normalized.index.astype("int64"), name="module_id")
    normalized.columns = pd.Index(
        [_canonical_kinase_label(value) for value in normalized.columns],
        name="kinase",
    )
    normalized = normalized.T.groupby(level=0, sort=False).mean().T
    return normalized.astype(float).sort_index().sort_index(axis=1)


def normalize_signalome_module_assignments_for_parity(
    frame: pd.DataFrame,
) -> pd.DataFrame:
    normalized = frame.copy(deep=True).astype("object")
    for column_name in ("site_key", "display_id", "site_id"):
        duplicate_name = f"{column_name}.1"
        if (
            column_name not in normalized.columns
            and duplicate_name in normalized.columns
        ):
            normalized = normalized.rename(columns={duplicate_name: column_name})
    normalized.index = pd.Index(normalized.index.astype(str), name="site_id")
    for column_name in (
        "site_key",
        "display_id",
        "gene_symbol",
        "site",
        "protein_id",
        "protein_accession",
        "isoform_id",
        "top_kinase",
        "top_kinase_selection_policy",
        "module_top_kinase",
        "module_top_kinase_selection_policy",
    ):
        if column_name in normalized.columns:
            normalized.loc[:, column_name] = (
                normalized.loc[:, column_name].fillna("").astype(str)
            )
    collection_columns = (
        "top_kinase_candidates",
        "top_kinase_weights",
        "module_top_kinase_candidates",
    )
    for column in collection_columns:
        if column in normalized.columns:
            normalized.loc[:, column] = normalized.loc[:, column].map(
                _normalize_signalome_kinase_weight_collection_value
                if column == "top_kinase_weights"
                else _normalize_signalome_kinase_collection_value
            )
    bool_columns = (
        "top_kinase_is_ambiguous",
        "module_top_kinase_is_ambiguous",
    )
    for column in bool_columns:
        if column in normalized.columns:
            normalized.loc[:, column] = normalized.loc[:, column].map(
                lambda value: (
                    value
                    if isinstance(value, bool)
                    else str(value).strip().lower() == "true"
                )
            )
    int_columns = (
        "module_id",
        "top_kinase_tie_count",
        "module_top_kinase_tie_count",
    )
    for column in int_columns:
        if column in normalized.columns:
            normalized.loc[:, column] = normalized.loc[:, column].astype("int64")
    if "top_score" in normalized.columns:
        normalized.loc[:, "top_score"] = normalized.loc[:, "top_score"].astype(float)
    kinase_scalar_columns = ("top_kinase", "module_top_kinase")
    for column in kinase_scalar_columns:
        if column in normalized.columns:
            normalized.loc[:, column] = normalized.loc[:, column].map(
                _canonical_kinase_label
            )
    return normalized.sort_index()


def normalize_signalome_network_nodes_for_parity(frame: pd.DataFrame) -> pd.DataFrame:
    normalized = frame.copy(deep=True)
    normalized.index = pd.Index(
        [_canonical_kinase_label(value) for value in normalized.index],
        name="kinase",
    )
    normalized = normalized.groupby(level=0, sort=False).sum()
    normalized = normalized.astype({"degree": "int64", "n_substrates": "int64"})
    return normalized.sort_index()


def normalize_signalome_network_edges_for_parity(frame: pd.DataFrame) -> pd.DataFrame:
    dtype_map: dict[str, object] = {
        "source_kinase": str,
        "target_kinase": str,
        "correlation": float,
    }
    if "valid_observations" in frame.columns:
        dtype_map["valid_observations"] = "int64"
    normalized = frame.copy(deep=True).astype(dtype_map)
    normalized.loc[:, "source_kinase"] = normalized.loc[:, "source_kinase"].map(
        _canonical_kinase_label
    )
    normalized.loc[:, "target_kinase"] = normalized.loc[:, "target_kinase"].map(
        _canonical_kinase_label
    )
    return normalized.sort_values(
        ["source_kinase", "target_kinase"],
        kind="mergesort",
    ).reset_index(drop=True)


def normalize_signalome_expanded_signalome_for_parity(
    frame: pd.DataFrame,
) -> pd.DataFrame:
    normalized = frame.copy(deep=True).astype("object")
    for column_name in (
        "kinase",
        "row_kind",
        "assignment_policy",
        "linked_kinases",
        "regulated_module_ids",
        "site_key",
        "display_id",
        "site_id",
        "gene_symbol",
        "site",
        "protein_id",
        "protein_accession",
        "isoform_id",
        "support_kinases",
        "top_kinase",
    ):
        if column_name in normalized.columns:
            normalized.loc[:, column_name] = (
                normalized.loc[:, column_name].fillna("").astype(str)
            )
    normalized = normalized.astype(
        {
            "site_order": "int64",
            "module_id": "int64",
            "support_weight": float,
            "top_score": float,
        }
    )
    normalized.loc[:, "kinase"] = normalized.loc[:, "kinase"].map(
        _canonical_kinase_label
    )
    normalized.loc[:, "top_kinase"] = normalized.loc[:, "top_kinase"].map(
        _canonical_kinase_label
    )
    kinase_collection_columns = ("linked_kinases", "support_kinases")
    for column in kinase_collection_columns:
        normalized.loc[:, column] = normalized.loc[:, column].map(
            _normalize_signalome_kinase_collection_value
        )
    return normalized.sort_values(
        [
            "kinase",
            "row_kind",
            "site_id",
            "site_order",
            "module_id",
            "protein_id",
            "top_kinase",
        ],
        kind="mergesort",
    ).reset_index(drop=True)


def site_metadata_for(phospho: pd.DataFrame) -> pd.DataFrame:
    split = phospho.index.to_series().astype(str).str.split(";", expand=True)
    site_sequences = load_rat_l6_sequence_table().reindex(phospho.index)
    if site_sequences.isna().any():
        missing = int(site_sequences.isna().sum())
        raise AssertionError(
            f"fixture missing site sequences for {missing} phosphosites"
        )
    return pd.DataFrame(
        {
            "gene_symbol": split.loc[:, 0].values,
            "site": split.loc[:, 1].values,
            "protein_id": split.loc[:, 0].values,
            "site_sequence": site_sequences.values,
            "localisation_confidence": [0.95] * int(len(phospho.index)),
        },
        index=phospho.index.copy(),
    )


def build_rat_l6_dataset(
    *,
    n_sites: int | None = 220,
    include_protein_id: bool = True,
) -> AnalysisReadyPhosphoDataset:
    phospho = load_rat_l6_phospho().copy(deep=True)
    if n_sites is not None:
        phospho = phospho.head(n_sites)
    site_metadata = site_metadata_for(phospho).copy(deep=True)
    if include_protein_id and "protein_id" not in site_metadata.columns:
        site_metadata.loc[:, "protein_id"] = (
            site_metadata.loc[:, "gene_symbol"].astype(str).tolist()
        )
    request = DatasetBuildRequest(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        input_intensity_scale="linear",
    )
    return AnalysisReadyDatasetBuilder().run(request)


def _read_float_matrix_csv(path: Path) -> pd.DataFrame:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle)
        try:
            header = next(reader)
        except StopIteration:
            return pd.DataFrame()
        if len(header) < 2:
            raise ValueError(f"expected indexed matrix CSV with >=2 columns: {path}")
        index_name = str(header[0]).strip() or None
        column_names = [str(name) for name in header[1:]]
        index_values: list[str] = []
        matrix_values: list[list[float]] = []
        for row in reader:
            if not row:
                continue
            if len(row) < len(column_names) + 1:
                row = row + [""] * ((len(column_names) + 1) - len(row))
            index_values.append(str(row[0]))
            matrix_values.append([_parse_float_cell(value) for value in row[1:]])
    frame = pd.DataFrame(
        matrix_values,
        index=pd.Index(index_values, name=index_name),
        columns=pd.Index(column_names),
        dtype="float64",
    )
    return frame


def _parse_float_cell(raw_value: object) -> float:
    text = str(raw_value).strip()
    if text == "":
        return float("nan")
    lowered = text.lower()
    if lowered in {"na", "nan", "null", "none"}:
        return float("nan")
    return float(text)
