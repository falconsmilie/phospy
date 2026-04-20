from __future__ import annotations

import ast
import json
from functools import lru_cache
from pathlib import Path

import pandas as pd

from phospy import (
    AnalysisReadyDatasetBuilder,
    AnalysisReadyPhosphoDataset,
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
ACTIVITY_REFERENCE_WEIGHTED = REWRITE_PARITY_REFERENCE / "kinase_activity_matrix.csv"
ACTIVITY_REFERENCE_KSEA_SCORES = REWRITE_PARITY_REFERENCE / "ksea_scores.csv"
ACTIVITY_REFERENCE_KSEA_COUNTS = REWRITE_PARITY_REFERENCE / "ksea_counts.csv"
ACTIVITY_REFERENCE_TARGET_COUNTS = REWRITE_PARITY_REFERENCE / "kinase_target_counts.csv"
ACTIVITY_REFERENCE_TARGET_TABLE = REWRITE_PARITY_REFERENCE / "kinase_target_table.csv"
ACTIVITY_PARITY_FIXTURE_FILES: tuple[str, ...] = (
    "l6_phospho_matrix.csv",
    "native_profile_scores.csv",
    "predMat.csv",
    "kinase_activity_matrix.csv",
    "ksea_scores.csv",
    "ksea_counts.csv",
    "kinase_target_counts.csv",
    "kinase_target_table.csv",
)
L6_PREDICTION_REFERENCE_PHOSPHO = (
    REWRITE_PARITY_L6_PREDICTION_REFERENCE / "l6_phospho_matrix.csv"
)
L6_PREDICTION_REFERENCE_PROFILE = (
    REWRITE_PARITY_L6_PREDICTION_REFERENCE / "native_profile_scores.csv"
)
L6_PREDICTION_REFERENCE_COMBINED = (
    REWRITE_PARITY_L6_PREDICTION_REFERENCE / "native_combined_scores.csv"
)
L6_PREDICTION_REFERENCE_WEIGHTS = (
    REWRITE_PARITY_L6_PREDICTION_REFERENCE / "native_combined_weights.csv"
)
L6_PREDICTION_REFERENCE_CANDIDATES = (
    REWRITE_PARITY_L6_PREDICTION_REFERENCE / "native_candidate_substrates.csv"
)
L6_PREDICTION_REFERENCE_TOP30 = (
    REWRITE_PARITY_L6_PREDICTION_REFERENCE / "native_prediction_top30.csv"
)
L6_PREDICTION_REFERENCE_PREDMAT = REWRITE_PARITY_L6_PREDICTION_REFERENCE / "predMat.csv"
ADAPTIVE_SAMPLING_REPLAY_COMBINED_SCORES = (
    REWRITE_PARITY_ADAPTIVE_SAMPLING_REPLAY / "combined_scores.csv"
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
PUBLIC_PREDMAT_LEGACY_DONOR_DEFAULT = (
    REWRITE_PUBLIC_WORKFLOW_REFERENCE / "predmat_legacy_default_donor.csv"
)
PUBLIC_PREDMAT_LEGACY_DONOR_R_PARITY = (
    REWRITE_PUBLIC_WORKFLOW_REFERENCE / "predmat_legacy_r_parity_donor.csv"
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
SIGNALOME_REWRITE_L6_ASSIGNMENTS_SELECTED = (
    REWRITE_PUBLIC_WORKFLOW_REFERENCE
    / "signalome_rewrite_l6_module_assignments_selected.csv"
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
SIGNALOME_REWRITE_L6_NETWORK_EDGES_SELECTED = (
    REWRITE_PUBLIC_WORKFLOW_REFERENCE
    / "signalome_rewrite_l6_network_edges_selected.csv"
)
SIGNALOME_REWRITE_L6_EXPANDED_SIGNALOME = (
    REWRITE_PUBLIC_WORKFLOW_REFERENCE / "signalome_rewrite_l6_expanded_signalome.csv"
)
SIGNALOME_REWRITE_L6_EXPANDED_AKT1_SELECTED = (
    REWRITE_PUBLIC_WORKFLOW_REFERENCE
    / "signalome_rewrite_l6_expanded_akt1_selected.csv"
)
SIGNALOME_REWRITE_L6_CONTRACT = (
    REWRITE_PUBLIC_WORKFLOW_REFERENCE / "signalome_rewrite_l6_contract.json"
)


@lru_cache(maxsize=1)
def load_rat_l6_phospho() -> pd.DataFrame:
    return pd.read_csv(RAT_L6_PHOSPHO, index_col=0)


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
def load_activity_reference_weighted_activity() -> pd.DataFrame:
    frame = pd.read_csv(ACTIVITY_REFERENCE_WEIGHTED, index_col=0)
    frame.index = pd.Index(frame.index.astype(str), name="kinase")
    return frame.astype(float)


@lru_cache(maxsize=1)
def load_activity_reference_ksea_scores() -> pd.DataFrame:
    frame = pd.read_csv(ACTIVITY_REFERENCE_KSEA_SCORES, index_col=0)
    frame.index = pd.Index(frame.index.astype(str), name="kinase")
    return frame.astype(float)


@lru_cache(maxsize=1)
def load_activity_reference_ksea_counts() -> pd.Series:
    frame = pd.read_csv(ACTIVITY_REFERENCE_KSEA_COUNTS, index_col=0)
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
def load_l6_prediction_reference_combined_scores() -> pd.DataFrame:
    return pd.read_csv(L6_PREDICTION_REFERENCE_COMBINED, index_col=0)


@lru_cache(maxsize=1)
def load_l6_prediction_reference_weights() -> pd.DataFrame:
    frame = pd.read_csv(L6_PREDICTION_REFERENCE_WEIGHTS).set_index("kinase")
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
def load_fragile_support_combined_scores() -> pd.DataFrame:
    return pd.read_csv(
        REWRITE_PARITY_FRAGILE_SUPPORT_REFERENCE / "combined_scores.csv",
        index_col=0,
    )


@lru_cache(maxsize=1)
def load_fragile_support_combined_weights() -> pd.DataFrame:
    frame = pd.read_csv(
        REWRITE_PARITY_FRAGILE_SUPPORT_REFERENCE / "combined_weights.csv",
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
def load_adaptive_sampling_edge_combined_scores() -> pd.DataFrame:
    return pd.read_csv(
        REWRITE_PARITY_ADAPTIVE_SAMPLING_EDGE / "combined_scores.csv",
        index_col=0,
    )


@lru_cache(maxsize=1)
def load_adaptive_sampling_replay_combined_scores() -> pd.DataFrame:
    return pd.read_csv(ADAPTIVE_SAMPLING_REPLAY_COMBINED_SCORES, index_col=0)


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
def load_adaptive_sampling_replay_trace_iteration_probabilities() -> pd.DataFrame:
    return pd.read_csv(
        REWRITE_PARITY_ADAPTIVE_SAMPLING_REPLAY / "trace_iteration_probabilities.csv"
    )


@lru_cache(maxsize=1)
def load_adaptive_sampling_replay_trace_iteration_decision_values() -> pd.DataFrame:
    return pd.read_csv(
        REWRITE_PARITY_ADAPTIVE_SAMPLING_REPLAY / "trace_iteration_decision_values.csv"
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
def load_public_predmat_input_phospho() -> pd.DataFrame:
    frame = pd.read_csv(PUBLIC_PREDMAT_INPUT_PHOSPHO).set_index("phosphosite")
    frame.index = pd.Index(frame.index.astype(str), name="phosphosite")
    return frame.astype(float)


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
def load_public_predmat_rewrite_stable() -> pd.DataFrame:
    return pd.read_csv(PUBLIC_PREDMAT_REWRITE_STABLE, index_col=0)


@lru_cache(maxsize=1)
def load_public_predmat_rewrite_r_parity() -> pd.DataFrame:
    return pd.read_csv(PUBLIC_PREDMAT_REWRITE_R_PARITY, index_col=0)


@lru_cache(maxsize=1)
def load_public_predmat_legacy_default_donor() -> pd.DataFrame:
    return pd.read_csv(PUBLIC_PREDMAT_LEGACY_DONOR_DEFAULT, index_col=0)


@lru_cache(maxsize=1)
def load_public_predmat_legacy_r_parity_donor() -> pd.DataFrame:
    return pd.read_csv(PUBLIC_PREDMAT_LEGACY_DONOR_R_PARITY, index_col=0)


@lru_cache(maxsize=1)
def load_signalome_rewrite_l6_module_assignments_selected() -> pd.DataFrame:
    frame = pd.read_csv(SIGNALOME_REWRITE_L6_ASSIGNMENTS_SELECTED, index_col=0)
    frame.index = pd.Index(frame.index.astype(str), name="site_id")
    return frame.astype(
        {
            "protein_id": str,
            "module_id": "int64",
            "top_kinase": str,
            "top_score": float,
            "top_kinase_tie_count": "int64",
            "top_kinase_is_ambiguous": bool,
        }
    )


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
def load_signalome_rewrite_l6_network_edges_selected() -> pd.DataFrame:
    return pd.read_csv(SIGNALOME_REWRITE_L6_NETWORK_EDGES_SELECTED).astype(
        {
            "source_kinase": str,
            "target_kinase": str,
            "correlation": float,
        }
    )


@lru_cache(maxsize=1)
def load_signalome_rewrite_l6_network_edges() -> pd.DataFrame:
    return pd.read_csv(SIGNALOME_REWRITE_L6_NETWORK_EDGES).astype(
        {
            "source_kinase": str,
            "target_kinase": str,
            "correlation": float,
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


@lru_cache(maxsize=1)
def load_signalome_rewrite_l6_expanded_akt1_selected() -> pd.DataFrame:
    return pd.read_csv(SIGNALOME_REWRITE_L6_EXPANDED_AKT1_SELECTED).astype(
        {
            "kinase": str,
            "assignment_policy": str,
            "linked_kinases": str,
            "regulated_module_ids": str,
            "site_id": str,
            "site_order": "int64",
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


def normalize_signalome_modules_for_parity(frame: pd.DataFrame) -> pd.DataFrame:
    normalized = frame.copy(deep=True)
    normalized.index = pd.Index(normalized.index.astype("int64"), name="module_id")
    normalized.columns = pd.Index(normalized.columns.astype(str), name="kinase")
    return normalized.astype(float).sort_index().sort_index(axis=1)


def normalize_signalome_module_assignments_for_parity(
    frame: pd.DataFrame,
) -> pd.DataFrame:
    normalized = frame.copy(deep=True)
    normalized.index = pd.Index(normalized.index.astype(str), name="site_id")
    collection_columns = (
        "top_kinase_candidates",
        "top_kinase_weights",
        "module_top_kinase_candidates",
    )
    for column in collection_columns:
        if column in normalized.columns:
            normalized.loc[:, column] = normalized.loc[:, column].map(
                _normalize_signalome_collection_value
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
    return normalized.sort_index()


def normalize_signalome_network_nodes_for_parity(frame: pd.DataFrame) -> pd.DataFrame:
    normalized = frame.copy(deep=True)
    normalized.index = pd.Index(normalized.index.astype(str), name="kinase")
    normalized = normalized.astype({"degree": "int64", "n_substrates": "int64"})
    return normalized.sort_index()


def normalize_signalome_network_edges_for_parity(frame: pd.DataFrame) -> pd.DataFrame:
    normalized = frame.copy(deep=True).astype(
        {
            "source_kinase": str,
            "target_kinase": str,
            "correlation": float,
        }
    )
    return normalized.sort_values(
        ["source_kinase", "target_kinase"],
        kind="mergesort",
    ).reset_index(drop=True)


def normalize_signalome_expanded_signalome_for_parity(
    frame: pd.DataFrame,
) -> pd.DataFrame:
    normalized = frame.copy(deep=True).astype(
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
            "site_sequence": site_sequences.values,
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
    )
    return AnalysisReadyDatasetBuilder().run(request)
