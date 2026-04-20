from __future__ import annotations

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
SIGNALOME_REWRITE_L6_NETWORK_EDGES_SELECTED = (
    REWRITE_PUBLIC_WORKFLOW_REFERENCE
    / "signalome_rewrite_l6_network_edges_selected.csv"
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


def build_rat_l6_dataset(*, n_sites: int | None = 220) -> AnalysisReadyPhosphoDataset:
    phospho = load_rat_l6_phospho().copy(deep=True)
    if n_sites is not None:
        phospho = phospho.head(n_sites)
    request = DatasetBuildRequest(
        phospho=phospho,
        site_metadata=site_metadata_for(phospho),
        organism=Organism.RAT,
    )
    return AnalysisReadyDatasetBuilder().run(request)
