from __future__ import annotations

import time
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from phospy import AnalysisReadyDatasetBuilder, KinaseWorkflow, SignalomeWorkflow
from phospy.api import (
    DatasetBuildRequest,
    KinasePredictionConfig,
    KinaseScoringConfig,
    KinaseWorkflowRequest,
    Organism,
    ReferenceBundle,
    ReferenceContextCompatibilityPolicy,
    SignalomeWorkflowRequest,
)
from phospy.api.results import KinaseScoringResult, KinaseWorkflowResult
from phospy.errors.workflows import SignalomeScaleError
from phospy.io.bundles._signalome.snapshots import SignalomeWorkflowConfigSnapshot
from phospy.io.bundles.signalome import save_signalome_workflow_bundle
from phospy.provenance.hashing import hash_table_tolerance
from phospy.science.datasets.preprocessing.models import (
    PreprocessingPlan,
    PreprocessingState,
)
from phospy.science.datasets.preprocessing.stages.normalisation import (
    NormalisationStage,
)
from phospy.science.prediction.execution import run_adaptive_ensemble_prediction
from phospy.science.prediction.motif_scoring import (
    DEFAULT_MOTIF_FLANK_SIZE,
    build_motif_library,
    score_phosphosite_motifs,
)
from phospy.science.signalomes.clustering import (
    MAX_FULL_CORRELATION_SITE_COUNT,
    SIGNALOME_CANDIDATE_SCORING_MODE_NOT_EVALUATED,
    SIGNALOME_CANDIDATE_SCORING_POLICY_FULL,
    SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED,
    SIGNALOME_CANDIDATE_SCORING_SKIP_REASON_EXPLICIT_MODULE_COUNT,
    SIGNALOME_CLUSTERING_ENGINE_EXACT_PYTHON,
    SIGNALOME_CLUSTERING_ENGINE_SCIPY_HIERARCHICAL,
    run_signalome_clustering_engine,
    select_module_count_with_diagnostics,
)
from phospy.science.signalomes.models import (
    SIGNALOME_MODULE_SELECTION_STRATEGY_EXPLICIT_MODULE_COUNT,
)
from tests.support.performance_contracts import (
    ADAPTIVE_PREDICTION_CONTRACT_CANDIDATE_KINASES,
    ADAPTIVE_PREDICTION_CONTRACT_N_KINASES,
    ADAPTIVE_PREDICTION_CONTRACT_N_SITES,
    ADAPTIVE_PREDICTION_CONTRACT_TOP_K,
    ADAPTIVE_PREDICTION_PEAK_MIB_MAX,
    ADAPTIVE_PREDICTION_RUNTIME_SECONDS_MAX,
    BUNDLE_PUBLISH_PEAK_MIB_MAX,
    BUNDLE_PUBLISH_RUNTIME_SECONDS_MAX,
    DIAGNOSTIC_RUNTIME_ABSOLUTE_SECONDS,
    DIAGNOSTIC_RUNTIME_RATIO_MULTIPLIER,
    KINASE_FILTERED_REFERENCE_PEAK_MIB_MAX,
    KINASE_FILTERED_REFERENCE_RUNTIME_SECONDS_MAX,
    MOTIF_PEAK_MIB_MAX,
    MOTIF_RUNTIME_SECONDS_MAX,
    PREPROCESSING_CONTRACT_N_SAMPLES,
    PREPROCESSING_CONTRACT_N_SITES,
    PROVENANCE_HASHING_PEAK_MIB_MAX,
    PROVENANCE_HASHING_RUNTIME_SECONDS_MAX,
    QUANTILE_PEAK_MIB_MAX,
    QUANTILE_RUNTIME_SECONDS_MAX,
    SIGNALOME_ABOVE_THRESHOLD_RUNTIME_SECONDS_MAX,
    SIGNALOME_BACKEND_MEDIUM_PEAK_MIB_MAX,
    SIGNALOME_BACKEND_MEDIUM_RUNTIME_SECONDS_MAX,
    SIGNALOME_BACKEND_SMALL_PEAK_MIB_MAX,
    SIGNALOME_BACKEND_SMALL_RUNTIME_SECONDS_MAX,
    SIGNALOME_BELOW_THRESHOLD_RUNTIME_SECONDS_MAX,
    SIGNALOME_CONTRACT_N_KINASES,
    SIGNALOME_CONTRACT_N_SITES,
    SIGNALOME_FULL_GUARD_RUNTIME_SECONDS_MAX,
    SIGNALOME_NEAR_THRESHOLD_RUNTIME_SECONDS_MAX,
    SIGNALOME_WORKFLOW_PEAK_MIB_MAX,
    SIGNALOME_WORKFLOW_PRECONDITIONED_PEAK_MIB_MAX,
    SIGNALOME_WORKFLOW_PRECONDITIONED_RUNTIME_SECONDS_MAX,
    SIGNALOME_WORKFLOW_RUNTIME_SECONDS_MAX,
    deterministic_kinase_substrate_map,
    deterministic_matrix,
    deterministic_site_ids,
    deterministic_site_metadata,
    deterministic_site_sequence_frame,
    deterministic_site_sequence_series,
    measure_runtime_and_peak_mib,
    median_runtime_seconds,
)
from tests.support.signalome_config import build_signalome_config

pytestmark = [pytest.mark.performance, pytest.mark.release_gate]

APPROXIMATION_REASON_TOKEN = "Used sampled within-cluster correlation estimates"


def _patch_cluster_tree_build_for_contract_scoring(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from phospy.science.signalomes.clustering import exact_python as exact_clustering
    from phospy.science.signalomes.clustering.backends import (
        exact_python as exact_tree_backend,
    )

    def _stub_build_cluster_tree(scoring_values: np.ndarray) -> object:
        n_sites = int(np.asarray(scoring_values, dtype=float).shape[0])
        return exact_clustering._WardClusterTree(n_sites=n_sites, merges=())

    def _stub_build_cluster_labels_from_tree(
        *,
        cluster_tree: object,
        cluster_counts: object,
    ) -> dict[int, np.ndarray]:
        n_sites = int(cluster_tree.n_sites)
        labels_by_count: dict[int, np.ndarray] = {}
        for requested_count in [int(value) for value in cluster_counts]:
            resolved_count = max(1, min(int(requested_count), n_sites))
            if resolved_count == 1:
                labels = np.zeros(n_sites, dtype=int)
            else:
                labels = np.arange(n_sites, dtype=int) % resolved_count
            labels_by_count[requested_count] = labels.astype(int, copy=False)
        return labels_by_count

    monkeypatch.setattr(
        exact_tree_backend,
        "build_cluster_tree",
        _stub_build_cluster_tree,
    )
    monkeypatch.setattr(
        exact_tree_backend,
        "build_cluster_labels_from_tree",
        _stub_build_cluster_labels_from_tree,
    )


def _build_signalome_scoring_matrix(
    *, n_sites: int, n_kinases: int, seed: int
) -> pd.DataFrame:
    matrix = deterministic_matrix(n_sites=n_sites, n_samples=n_kinases, seed=seed)
    matrix.columns = pd.Index(
        [f"KINASE_{index + 1:03d}" for index in range(matrix.shape[1])],
        name="kinase",
    )
    return matrix


def _build_signalome_realistic_scoring_matrix(
    *,
    n_sites: int,
    n_kinases: int,
    seed: int,
) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    n_modules = max(4, min(20, n_sites // 12))
    module_profiles = rng.normal(
        loc=10.0,
        scale=1.7,
        size=(n_modules, int(n_kinases)),
    )
    module_assignments = (np.arange(int(n_sites), dtype=int) * 5) % int(n_modules)
    noise = rng.normal(loc=0.0, scale=0.28, size=(int(n_sites), int(n_kinases)))
    values = np.round(module_profiles[module_assignments] + noise, decimals=6)
    site_ids = deterministic_site_ids(int(n_sites), start=30_000, gene_prefix="SIGSITE")
    return pd.DataFrame(
        values,
        index=site_ids,
        columns=pd.Index(
            [f"KINASE_{index + 1:03d}" for index in range(int(n_kinases))],
            name="kinase",
        ),
        dtype=float,
    )


def _build_signalome_site_to_protein(
    site_ids: pd.Index,
    *,
    sites_per_protein: int,
) -> pd.Series:
    if sites_per_protein < 1:
        raise ValueError("sites_per_protein must be >= 1")
    proteins = [
        f"PROT_{(index // int(sites_per_protein)) + 1:05d}"
        for index in range(int(site_ids.size))
    ]
    return pd.Series(
        proteins,
        index=pd.Index(site_ids.astype(str), name="site_id"),
        name="protein_id",
        dtype=str,
    )


def _build_adaptive_prediction_hot_path_inputs(
    *,
    n_sites: int,
    n_kinases: int,
    candidate_kinases: int,
    seed: int,
) -> tuple[pd.DataFrame, dict[str, list[str]]]:
    prediction_score_matrix = deterministic_matrix(
        n_sites=n_sites,
        n_samples=n_kinases,
        seed=seed,
        site_ids=deterministic_site_ids(n_sites, start=140_000, gene_prefix="PREDSITE"),
        sample_columns=pd.Index(
            [f"KINASE_{index + 1:03d}" for index in range(n_kinases)],
            name="kinase",
        ),
    )
    site_values = prediction_score_matrix.index.astype(str).tolist()
    candidate_substrates: dict[str, list[str]] = {}
    sites_per_kinase = min(96, len(site_values))
    for kinase_index in range(min(candidate_kinases, n_kinases)):
        kinase = f"KINASE_{kinase_index + 1:03d}"
        start = (kinase_index * 17) % len(site_values)
        candidate_substrates[kinase] = [
            site_values[(start + offset) % len(site_values)]
            for offset in range(sites_per_kinase)
        ]
    return prediction_score_matrix, candidate_substrates


def _collect_backend_contract_snapshot(
    *,
    backend_result: object,
    scoring_matrix: pd.DataFrame,
    runtime_seconds: float,
    peak_mib: float,
) -> dict[str, object]:
    diagnostics = backend_result.module_selection_diagnostics
    return {
        "backend_name": str(backend_result.backend_name),
        "site_count": int(scoring_matrix.shape[0]),
        "kinase_count": int(scoring_matrix.shape[1]),
        "selected_module_count": int(diagnostics.selected_module_count),
        "exact_tree_backend": str(backend_result.tree_implementation),
        "candidate_scoring_mode": str(backend_result.candidate_scoring_mode),
        "sampled_candidate_scoring_activated": bool(
            backend_result.candidate_scoring_evaluated
            and backend_result.candidate_scoring_mode
            == SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED
        ),
        "candidate_scoring_skipped": bool(
            not backend_result.candidate_scoring_evaluated
        ),
        "exact_tree_construction_occurred": bool(
            backend_result.exact_cluster_tree_built
        ),
        "runtime_seconds": float(runtime_seconds),
        "peak_mib": float(peak_mib),
    }


def _cluster_partitions_match(left: pd.Series, right: pd.Series) -> bool:
    if not left.index.equals(right.index):
        return False
    left_values = left.to_numpy(dtype=int, copy=False)
    right_values = right.to_numpy(dtype=int, copy=False)
    left_equal = left_values[:, None] == left_values[None, :]
    right_equal = right_values[:, None] == right_values[None, :]
    return bool(np.array_equal(left_equal, right_equal))


def _run_signalome_backend_contract(
    *,
    scoring_matrix: pd.DataFrame,
    site_to_protein: pd.Series,
    clustering_engine: str,
    max_clusters: int = 8,
    candidate_scoring_policy: str | None = None,
    max_exact_tree_sites: int | None = None,
    max_full_candidate_scoring_sites: int = MAX_FULL_CORRELATION_SITE_COUNT,
) -> tuple[object, float, float]:
    return measure_runtime_and_peak_mib(
        lambda: run_signalome_clustering_engine(
            scoring_matrix=scoring_matrix,
            site_to_protein=site_to_protein,
            requested_module_count=None,
            primary_threshold=0.45,
            fallback_threshold=0.15,
            max_clusters=max_clusters,
            candidate_scoring_policy=candidate_scoring_policy,
            max_exact_tree_sites=max_exact_tree_sites,
            max_full_candidate_scoring_sites=max_full_candidate_scoring_sites,
            clustering_engine=clustering_engine,
        ),
        warmup=True,
    )


def _build_kinase_workflow_inputs(
    *,
    n_sites: int,
    n_samples: int,
    eligible_kinases: int,
    substrates_per_kinase: int,
    offlane_kinases: int,
    offlane_sites_per_kinase: int,
) -> tuple[pd.Index, object, ReferenceBundle, set[str]]:
    site_ids = deterministic_site_ids(n_sites)
    phospho = deterministic_matrix(
        n_sites=n_sites,
        n_samples=n_samples,
        seed=7227,
        site_ids=site_ids,
    )
    site_metadata = deterministic_site_metadata(site_ids, include_protein_id=True)
    dataset = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
            input_intensity_scale="linear",
        )
    )

    offlane_site_ids = deterministic_site_ids(
        offlane_kinases * offlane_sites_per_kinase + 32,
        start=750_000,
        gene_prefix="OFFSITE",
    )
    kinase_substrate_map = deterministic_kinase_substrate_map(
        dataset_site_ids=site_ids,
        eligible_kinase_count=eligible_kinases,
        substrates_per_kinase=substrates_per_kinase,
        offlane_kinase_count=offlane_kinases,
        offlane_sites_per_kinase=offlane_sites_per_kinase,
        offlane_site_ids=offlane_site_ids,
    )
    reference_site_ids = site_ids.append(offlane_site_ids)
    site_sequences = deterministic_site_sequence_frame(
        reference_site_ids,
        sequence_width=(2 * DEFAULT_MOTIF_FLANK_SIZE) + 1,
    )
    references = ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=kinase_substrate_map,
        site_sequences=site_sequences,
    )
    eligible_kinase_names = {
        f"KINASE_{index + 1:03d}" for index in range(int(eligible_kinases))
    }
    return site_ids, dataset, references, eligible_kinase_names


def _build_kinase_result_for_signalome_performance(
    *,
    n_sites: int,
    n_samples: int,
    eligible_kinases: int,
    substrates_per_kinase: int,
    offlane_kinases: int,
    offlane_sites_per_kinase: int,
) -> tuple[KinaseWorkflowResult, set[str]]:
    _site_ids, dataset, references, eligible_kinase_names = (
        _build_kinase_workflow_inputs(
            n_sites=n_sites,
            n_samples=n_samples,
            eligible_kinases=eligible_kinases,
            substrates_per_kinase=substrates_per_kinase,
            offlane_kinases=offlane_kinases,
            offlane_sites_per_kinase=offlane_sites_per_kinase,
        )
    )
    request = KinaseWorkflowRequest(
        dataset=dataset,
        references=references,
        scoring_config=KinaseScoringConfig(
            min_substrates=2,
            include_diagnostic_scoring_tables=False,
            reference_context_compatibility_policy=(
                ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT
            ),
        ),
        prediction_config=KinasePredictionConfig(
            top_k=6,
            deterministic_max_selected_kinases=10,
            adaptive_ensemble_runs=10,
        ),
        activity_config=None,
    )
    return KinaseWorkflow().run(request), eligible_kinase_names


def test_signalome_full_vs_approximate_correlation_performance_contract(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _patch_cluster_tree_build_for_contract_scoring(monkeypatch)

    below_threshold_matrix = _build_signalome_scoring_matrix(
        n_sites=600,
        n_kinases=40,
        seed=1103,
    )
    below_diagnostics, below_runtime_seconds, _below_peak_mib = (
        measure_runtime_and_peak_mib(
            lambda: select_module_count_with_diagnostics(
                scoring_values=below_threshold_matrix,
                max_clusters=3,
            ),
            warmup=True,
        )
    )
    assert APPROXIMATION_REASON_TOKEN not in below_diagnostics.reason
    assert below_runtime_seconds < SIGNALOME_BELOW_THRESHOLD_RUNTIME_SECONDS_MAX

    above_threshold_sites = MAX_FULL_CORRELATION_SITE_COUNT + 50
    above_threshold_matrix = _build_signalome_scoring_matrix(
        n_sites=above_threshold_sites,
        n_kinases=12,
        seed=1107,
    )
    above_diagnostics, above_runtime_seconds, above_peak_mib = (
        measure_runtime_and_peak_mib(
            lambda: select_module_count_with_diagnostics(
                scoring_values=above_threshold_matrix,
                max_clusters=3,
                max_exact_tree_sites=above_threshold_sites,
            ),
            warmup=True,
        )
    )
    assert APPROXIMATION_REASON_TOKEN in above_diagnostics.reason
    assert above_runtime_seconds < SIGNALOME_ABOVE_THRESHOLD_RUNTIME_SECONDS_MAX

    full_matrix_mib = (above_threshold_sites * above_threshold_sites * 8) / (
        1024 * 1024
    )
    assert above_peak_mib < full_matrix_mib * 0.75


def test_signalome_backend_contracts_compare_exact_and_scipy_equivalent_small_fixture() -> (
    None
):
    scoring_matrix = _build_signalome_realistic_scoring_matrix(
        n_sites=48,
        n_kinases=14,
        seed=2301,
    )
    site_to_protein = _build_signalome_site_to_protein(
        scoring_matrix.index,
        sites_per_protein=2,
    )

    exact_result, exact_runtime, exact_peak_mib = _run_signalome_backend_contract(
        scoring_matrix=scoring_matrix,
        site_to_protein=site_to_protein,
        clustering_engine=SIGNALOME_CLUSTERING_ENGINE_EXACT_PYTHON,
        max_clusters=8,
        max_exact_tree_sites=96,
        max_full_candidate_scoring_sites=96,
    )
    scipy_result, scipy_runtime, scipy_peak_mib = _run_signalome_backend_contract(
        scoring_matrix=scoring_matrix,
        site_to_protein=site_to_protein,
        clustering_engine=SIGNALOME_CLUSTERING_ENGINE_SCIPY_HIERARCHICAL,
        max_clusters=8,
        max_exact_tree_sites=96,
        max_full_candidate_scoring_sites=96,
    )

    assert _cluster_partitions_match(
        exact_result.site_clusters,
        scipy_result.site_clusters,
    )
    assert _cluster_partitions_match(
        exact_result.protein_modules,
        scipy_result.protein_modules,
    )
    assert (
        exact_result.module_selection_diagnostics.selected_module_count
        == scipy_result.module_selection_diagnostics.selected_module_count
    )
    assert exact_result.candidate_scoring_mode == scipy_result.candidate_scoring_mode

    exact_snapshot = _collect_backend_contract_snapshot(
        backend_result=exact_result,
        scoring_matrix=scoring_matrix,
        runtime_seconds=exact_runtime,
        peak_mib=exact_peak_mib,
    )
    scipy_snapshot = _collect_backend_contract_snapshot(
        backend_result=scipy_result,
        scoring_matrix=scoring_matrix,
        runtime_seconds=scipy_runtime,
        peak_mib=scipy_peak_mib,
    )

    assert exact_snapshot["backend_name"] == SIGNALOME_CLUSTERING_ENGINE_EXACT_PYTHON
    assert (
        scipy_snapshot["backend_name"] == SIGNALOME_CLUSTERING_ENGINE_SCIPY_HIERARCHICAL
    )
    assert exact_snapshot["site_count"] == int(scoring_matrix.shape[0])
    assert scipy_snapshot["site_count"] == int(scoring_matrix.shape[0])
    assert exact_snapshot["kinase_count"] == int(scoring_matrix.shape[1])
    assert scipy_snapshot["kinase_count"] == int(scoring_matrix.shape[1])
    assert exact_snapshot["selected_module_count"] >= 1
    assert scipy_snapshot["selected_module_count"] >= 1
    assert exact_snapshot["exact_tree_backend"] == "exact_python_tree"
    assert scipy_snapshot["exact_tree_backend"] == "scipy_hierarchical_tree"
    assert (
        exact_snapshot["candidate_scoring_mode"]
        == SIGNALOME_CANDIDATE_SCORING_POLICY_FULL
    )
    assert (
        scipy_snapshot["candidate_scoring_mode"]
        == SIGNALOME_CANDIDATE_SCORING_POLICY_FULL
    )
    assert exact_snapshot["sampled_candidate_scoring_activated"] is False
    assert scipy_snapshot["sampled_candidate_scoring_activated"] is False
    assert exact_snapshot["candidate_scoring_skipped"] is False
    assert scipy_snapshot["candidate_scoring_skipped"] is False
    assert exact_snapshot["exact_tree_construction_occurred"] is True
    assert scipy_snapshot["exact_tree_construction_occurred"] is True
    assert (
        exact_snapshot["runtime_seconds"] < SIGNALOME_BACKEND_SMALL_RUNTIME_SECONDS_MAX
    )
    assert (
        scipy_snapshot["runtime_seconds"] < SIGNALOME_BACKEND_SMALL_RUNTIME_SECONDS_MAX
    )
    assert exact_snapshot["peak_mib"] < SIGNALOME_BACKEND_SMALL_PEAK_MIB_MAX
    assert scipy_snapshot["peak_mib"] < SIGNALOME_BACKEND_SMALL_PEAK_MIB_MAX


def test_signalome_backend_contracts_medium_fixture_activates_sampled_candidate_scoring() -> (
    None
):
    scoring_matrix = _build_signalome_realistic_scoring_matrix(
        n_sites=420,
        n_kinases=24,
        seed=2302,
    )
    site_to_protein = _build_signalome_site_to_protein(
        scoring_matrix.index,
        sites_per_protein=3,
    )

    exact_result, exact_runtime, exact_peak_mib = _run_signalome_backend_contract(
        scoring_matrix=scoring_matrix,
        site_to_protein=site_to_protein,
        clustering_engine=SIGNALOME_CLUSTERING_ENGINE_EXACT_PYTHON,
        max_clusters=10,
        max_exact_tree_sites=500,
        max_full_candidate_scoring_sites=180,
    )
    scipy_result, scipy_runtime, scipy_peak_mib = _run_signalome_backend_contract(
        scoring_matrix=scoring_matrix,
        site_to_protein=site_to_protein,
        clustering_engine=SIGNALOME_CLUSTERING_ENGINE_SCIPY_HIERARCHICAL,
        max_clusters=10,
        max_exact_tree_sites=500,
        max_full_candidate_scoring_sites=180,
    )

    assert (
        exact_result.candidate_scoring_mode
        == SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED
    )
    assert (
        scipy_result.candidate_scoring_mode
        == SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED
    )
    assert exact_result.candidate_scoring_evaluated is True
    assert scipy_result.candidate_scoring_evaluated is True
    assert exact_result.candidate_scoring_skip_reason is None
    assert scipy_result.candidate_scoring_skip_reason is None
    assert exact_result.candidate_scoring_sampling is not None
    assert scipy_result.candidate_scoring_sampling is not None
    assert (
        APPROXIMATION_REASON_TOKEN in exact_result.module_selection_diagnostics.reason
    )
    assert (
        APPROXIMATION_REASON_TOKEN in scipy_result.module_selection_diagnostics.reason
    )
    assert (
        exact_result.module_selection_diagnostics.selected_module_count
        == scipy_result.module_selection_diagnostics.selected_module_count
    )
    assert exact_result.exact_cluster_tree_built is True
    assert scipy_result.exact_cluster_tree_built is True

    exact_snapshot = _collect_backend_contract_snapshot(
        backend_result=exact_result,
        scoring_matrix=scoring_matrix,
        runtime_seconds=exact_runtime,
        peak_mib=exact_peak_mib,
    )
    scipy_snapshot = _collect_backend_contract_snapshot(
        backend_result=scipy_result,
        scoring_matrix=scoring_matrix,
        runtime_seconds=scipy_runtime,
        peak_mib=scipy_peak_mib,
    )
    assert exact_snapshot["sampled_candidate_scoring_activated"] is True
    assert scipy_snapshot["sampled_candidate_scoring_activated"] is True
    assert (
        exact_snapshot["runtime_seconds"] < SIGNALOME_BACKEND_MEDIUM_RUNTIME_SECONDS_MAX
    )
    assert (
        scipy_snapshot["runtime_seconds"] < SIGNALOME_BACKEND_MEDIUM_RUNTIME_SECONDS_MAX
    )
    assert exact_snapshot["peak_mib"] < SIGNALOME_BACKEND_MEDIUM_PEAK_MIB_MAX
    assert scipy_snapshot["peak_mib"] < SIGNALOME_BACKEND_MEDIUM_PEAK_MIB_MAX


def test_signalome_candidate_scoring_contract_full_vs_sampled_policy() -> None:
    scoring_matrix = _build_signalome_realistic_scoring_matrix(
        n_sites=420,
        n_kinases=24,
        seed=2312,
    )
    site_to_protein = _build_signalome_site_to_protein(
        scoring_matrix.index,
        sites_per_protein=3,
    )

    full_result, full_runtime, full_peak_mib = _run_signalome_backend_contract(
        scoring_matrix=scoring_matrix,
        site_to_protein=site_to_protein,
        clustering_engine=SIGNALOME_CLUSTERING_ENGINE_SCIPY_HIERARCHICAL,
        max_clusters=10,
        candidate_scoring_policy=SIGNALOME_CANDIDATE_SCORING_POLICY_FULL,
        max_exact_tree_sites=500,
        max_full_candidate_scoring_sites=500,
    )
    sampled_result, sampled_runtime, sampled_peak_mib = _run_signalome_backend_contract(
        scoring_matrix=scoring_matrix,
        site_to_protein=site_to_protein,
        clustering_engine=SIGNALOME_CLUSTERING_ENGINE_SCIPY_HIERARCHICAL,
        max_clusters=10,
        candidate_scoring_policy=SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED,
        max_exact_tree_sites=500,
        max_full_candidate_scoring_sites=500,
    )

    assert full_result.candidate_scoring_mode == SIGNALOME_CANDIDATE_SCORING_POLICY_FULL
    assert (
        sampled_result.candidate_scoring_mode
        == SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED
    )
    assert (
        APPROXIMATION_REASON_TOKEN
        not in full_result.module_selection_diagnostics.reason
    )
    assert (
        APPROXIMATION_REASON_TOKEN in sampled_result.module_selection_diagnostics.reason
    )
    assert full_result.candidate_scoring_evaluated is True
    assert sampled_result.candidate_scoring_evaluated is True
    assert full_result.exact_cluster_tree_built is True
    assert sampled_result.exact_cluster_tree_built is True
    assert full_runtime < SIGNALOME_BACKEND_MEDIUM_RUNTIME_SECONDS_MAX
    assert sampled_runtime < SIGNALOME_BACKEND_MEDIUM_RUNTIME_SECONDS_MAX
    assert full_peak_mib < SIGNALOME_BACKEND_MEDIUM_PEAK_MIB_MAX
    assert sampled_peak_mib < SIGNALOME_BACKEND_MEDIUM_PEAK_MIB_MAX


def test_signalome_exact_tree_guard_contract_near_threshold_fixture() -> None:
    near_limit = 72
    passing_matrix = _build_signalome_realistic_scoring_matrix(
        n_sites=near_limit,
        n_kinases=12,
        seed=2303,
    )
    site_to_protein = _build_signalome_site_to_protein(
        passing_matrix.index,
        sites_per_protein=2,
    )
    passing_result, passing_runtime, _passing_peak = _run_signalome_backend_contract(
        scoring_matrix=passing_matrix,
        site_to_protein=site_to_protein,
        clustering_engine=SIGNALOME_CLUSTERING_ENGINE_EXACT_PYTHON,
        max_clusters=7,
        max_exact_tree_sites=near_limit,
        max_full_candidate_scoring_sites=near_limit,
    )
    assert passing_result.exact_cluster_tree_built is True
    assert passing_result.module_selection_diagnostics.selected_module_count >= 1
    assert passing_runtime < SIGNALOME_NEAR_THRESHOLD_RUNTIME_SECONDS_MAX

    failing_matrix = _build_signalome_realistic_scoring_matrix(
        n_sites=near_limit + 1,
        n_kinases=12,
        seed=2304,
    )
    failing_site_to_protein = _build_signalome_site_to_protein(
        failing_matrix.index,
        sites_per_protein=2,
    )
    started = time.perf_counter()
    with pytest.raises(SignalomeScaleError) as exc_info:
        run_signalome_clustering_engine(
            scoring_matrix=failing_matrix,
            site_to_protein=failing_site_to_protein,
            requested_module_count=None,
            max_clusters=7,
            max_exact_tree_sites=near_limit,
            max_full_candidate_scoring_sites=near_limit + 20,
            clustering_engine=SIGNALOME_CLUSTERING_ENGINE_EXACT_PYTHON,
        )
    guard_runtime_seconds = time.perf_counter() - started
    message = str(exc_info.value).lower()
    assert "exact cluster-tree construction received" in message
    assert f"max_exact_tree_sites={near_limit}" in message
    assert "tree_implementation='exact_cluster_tree'" in message
    assert guard_runtime_seconds < SIGNALOME_FULL_GUARD_RUNTIME_SECONDS_MAX


def test_signalome_full_correlation_guard_contract_fixture() -> None:
    scoring_matrix = _build_signalome_realistic_scoring_matrix(
        n_sites=90,
        n_kinases=10,
        seed=2305,
    )
    site_to_protein = _build_signalome_site_to_protein(
        scoring_matrix.index,
        sites_per_protein=2,
    )
    started = time.perf_counter()
    with pytest.raises(SignalomeScaleError) as exc_info:
        run_signalome_clustering_engine(
            scoring_matrix=scoring_matrix,
            site_to_protein=site_to_protein,
            requested_module_count=None,
            max_clusters=6,
            candidate_scoring_policy=SIGNALOME_CANDIDATE_SCORING_POLICY_FULL,
            max_exact_tree_sites=120,
            max_full_candidate_scoring_sites=80,
            clustering_engine=SIGNALOME_CLUSTERING_ENGINE_EXACT_PYTHON,
        )
    guard_runtime_seconds = time.perf_counter() - started
    message = str(exc_info.value).lower()
    assert "full candidate-correlation scoring would evaluate" in message
    assert "exact cluster-tree construction has not been attempted" in message
    assert "use candidate_scoring_policy='sampled'" in message
    assert guard_runtime_seconds < SIGNALOME_FULL_GUARD_RUNTIME_SECONDS_MAX


def test_signalome_candidate_scoring_skip_contract_for_explicit_module_count() -> None:
    scoring_matrix = _build_signalome_realistic_scoring_matrix(
        n_sites=44,
        n_kinases=12,
        seed=2306,
    )
    site_to_protein = _build_signalome_site_to_protein(
        scoring_matrix.index,
        sites_per_protein=2,
    )
    backend_result, runtime_seconds, peak_mib = measure_runtime_and_peak_mib(
        lambda: run_signalome_clustering_engine(
            scoring_matrix=scoring_matrix,
            site_to_protein=site_to_protein,
            requested_module_count=4,
            max_clusters=8,
            clustering_engine=SIGNALOME_CLUSTERING_ENGINE_EXACT_PYTHON,
        ),
        warmup=True,
    )
    diagnostics = backend_result.module_selection_diagnostics
    snapshot = _collect_backend_contract_snapshot(
        backend_result=backend_result,
        scoring_matrix=scoring_matrix,
        runtime_seconds=runtime_seconds,
        peak_mib=peak_mib,
    )

    assert (
        diagnostics.strategy
        == SIGNALOME_MODULE_SELECTION_STRATEGY_EXPLICIT_MODULE_COUNT
    )
    assert (
        backend_result.candidate_scoring_mode
        == SIGNALOME_CANDIDATE_SCORING_MODE_NOT_EVALUATED
    )
    assert backend_result.candidate_scoring_evaluated is False
    assert (
        backend_result.candidate_scoring_skip_reason
        == SIGNALOME_CANDIDATE_SCORING_SKIP_REASON_EXPLICIT_MODULE_COUNT
    )
    assert snapshot["candidate_scoring_skipped"] is True
    assert snapshot["runtime_seconds"] < SIGNALOME_BACKEND_SMALL_RUNTIME_SECONDS_MAX
    assert snapshot["peak_mib"] < SIGNALOME_BACKEND_SMALL_PEAK_MIB_MAX


def test_signalome_workflow_performance_contract_reports_scale_guard_diagnostics() -> (
    None
):
    kinase_result, eligible_kinases = _build_kinase_result_for_signalome_performance(
        n_sites=220,
        n_samples=8,
        eligible_kinases=42,
        substrates_per_kinase=6,
        offlane_kinases=240,
        offlane_sites_per_kinase=8,
    )
    result, runtime_seconds, peak_mib = measure_runtime_and_peak_mib(
        lambda: SignalomeWorkflow().run(
            SignalomeWorkflowRequest(
                kinase_result=kinase_result,
                config=build_signalome_config(
                    substrate_support_cutoff=0.5,
                    module_selection_max_clusters=8,
                    candidate_scoring_policy=SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED,
                    max_exact_tree_sites=360,
                    max_full_candidate_scoring_sites=140,
                    clustering_engine=SIGNALOME_CLUSTERING_ENGINE_SCIPY_HIERARCHICAL,
                ),
            )
        ),
        warmup=True,
    )
    modules = result.signalome_modules.table
    assert not modules.empty
    assert set(modules.columns.astype(str).tolist()).issubset(eligible_kinases)
    assert result.provenance is not None
    workflow_parameters = result.provenance.workflow_parameters
    scale_guard = workflow_parameters["scale_guard"]
    backend_diagnostics = scale_guard["backend_diagnostics"]
    assert isinstance(backend_diagnostics, dict)
    assert (
        scale_guard["clustering_engine"]
        == SIGNALOME_CLUSTERING_ENGINE_SCIPY_HIERARCHICAL
    )
    assert (
        backend_diagnostics["backend_name"]
        == SIGNALOME_CLUSTERING_ENGINE_SCIPY_HIERARCHICAL
    )
    assert scale_guard["site_count"] == int(result.module_assignments.table.shape[0])
    assert scale_guard["input_protein_count"] >= 1
    assert scale_guard["input_kinase_count"] >= 1
    assert scale_guard["selected_module_count"] == int(
        result.module_selection_diagnostics.selected_module_count
    )
    assert scale_guard["candidate_module_counts_evaluated"] >= 1
    assert (
        scale_guard["candidate_module_count_upper_bound"]
        >= scale_guard["candidate_module_counts_evaluated"]
    )
    assert (
        scale_guard["tree_implementation"]
        == scale_guard["backend_diagnostics"]["tree_implementation"]
    )
    assert scale_guard["tree_generation_mode"] == "full_exact_tree_construction"
    assert scale_guard["tree_generation_is_approximate"] is False
    assert (
        scale_guard["tree_generation_scope"]
        == "module_count_selection_and_final_assignment"
    )
    assert scale_guard["tree_generation_guard_triggered"] is False
    assert (
        scale_guard["candidate_scoring_mode"]
        == SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED
    )
    assert scale_guard["candidate_scoring_strategy"] == (
        SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED
    )
    assert scale_guard["candidate_scoring_is_approximate"] is True
    assert scale_guard["candidate_scoring_guard_triggered"] is False
    assert scale_guard["candidate_scoring_evaluated"] is True
    assert scale_guard["candidate_scoring_skip_reason"] is None
    assert scale_guard["exact_cluster_tree_built"] is True
    assert isinstance(scale_guard["candidate_scoring_sampling"], dict)
    assert int(scale_guard["candidate_scoring_sampled_site_total"]) >= 0
    assert int(scale_guard["candidate_scoring_sampled_pair_count"]) >= 0
    assert workflow_parameters["score_preconditioning_diagnostics"][
        "retained_row_count"
    ] == int(result.module_assignments.table.shape[0])
    assert runtime_seconds < SIGNALOME_WORKFLOW_RUNTIME_SECONDS_MAX
    assert peak_mib < SIGNALOME_WORKFLOW_PEAK_MIB_MAX


def test_signalome_workflow_performance_contract_covers_all_missing_row_preconditioning() -> (
    None
):
    kinase_result, _eligible_kinases = _build_kinase_result_for_signalome_performance(
        n_sites=220,
        n_samples=8,
        eligible_kinases=42,
        substrates_per_kinase=6,
        offlane_kinases=240,
        offlane_sites_per_kinase=8,
    )
    rank_weighted_fusion_scores = (
        kinase_result.scoring_result.rank_weighted_fusion_scores
    )
    assert rank_weighted_fusion_scores is not None
    dropped_rows = 7
    sparse_scores = rank_weighted_fusion_scores.copy(deep=True)
    sparse_scores.iloc[:dropped_rows, :] = float("nan")
    sparse_kinase_result = KinaseWorkflowResult(
        dataset=kinase_result.dataset,
        references=kinase_result.references,
        scoring_result=KinaseScoringResult(
            profile_scores=kinase_result.scoring_result.profile_scores,
            motif_scores=kinase_result.scoring_result.motif_scores,
            rank_weighted_fusion_scores=sparse_scores,
            score_fusion_weights=kinase_result.scoring_result.score_fusion_weights,
        ),
        prediction_result=kinase_result.prediction_result,
        activity_result=kinase_result.activity_result,
        provenance=kinase_result.provenance,
    )
    result, runtime_seconds, peak_mib = measure_runtime_and_peak_mib(
        lambda: SignalomeWorkflow().run(
            SignalomeWorkflowRequest(
                kinase_result=sparse_kinase_result,
                config=build_signalome_config(
                    substrate_support_cutoff=0.5,
                    module_selection_max_clusters=8,
                    candidate_scoring_policy=SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED,
                    max_exact_tree_sites=360,
                    max_full_candidate_scoring_sites=140,
                    clustering_engine=SIGNALOME_CLUSTERING_ENGINE_EXACT_PYTHON,
                    score_preconditioning_policy="allow_and_report",
                ),
            )
        ),
        warmup=True,
    )
    assert result.score_preconditioning_diagnostics.input_row_count == int(
        kinase_result.prediction_result.pred_mat.shape[0]
    )
    assert (
        result.score_preconditioning_diagnostics.dropped_all_missing_row_count
        == dropped_rows
    )
    assert result.score_preconditioning_diagnostics.retained_row_count == (
        int(kinase_result.prediction_result.pred_mat.shape[0]) - dropped_rows
    )
    assert result.module_assignments.table.shape[0] == (
        int(kinase_result.prediction_result.pred_mat.shape[0]) - dropped_rows
    )
    assert result.provenance is not None
    preconditioning = result.provenance.workflow_parameters[
        "score_preconditioning_diagnostics"
    ]
    assert preconditioning["dropped_all_missing_row_count"] == dropped_rows
    assert preconditioning["retained_row_count"] == int(
        result.module_assignments.table.shape[0]
    )
    assert preconditioning["policy"] == "allow_and_report"
    assert runtime_seconds < SIGNALOME_WORKFLOW_PRECONDITIONED_RUNTIME_SECONDS_MAX
    assert peak_mib < SIGNALOME_WORKFLOW_PRECONDITIONED_PEAK_MIB_MAX


def test_quantile_normalisation_performance_contract() -> None:
    phospho = deterministic_matrix(
        n_sites=PREPROCESSING_CONTRACT_N_SITES,
        n_samples=PREPROCESSING_CONTRACT_N_SAMPLES,
        seed=181,
    )
    state = PreprocessingState(
        phospho=phospho,
        site_metadata=deterministic_site_metadata(
            phospho.index, include_protein_id=False
        ),
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan(normalisation_policy="quantile"),
    )
    stage = NormalisationStage()
    normalized_result, runtime_seconds, peak_mib = measure_runtime_and_peak_mib(
        lambda: stage.run(state),
        warmup=True,
    )
    normalized_phospho = normalized_result.state.phospho

    assert normalized_phospho.shape == phospho.shape
    assert int(normalized_phospho.isna().sum().sum()) == 0
    assert runtime_seconds < QUANTILE_RUNTIME_SECONDS_MAX
    assert peak_mib < QUANTILE_PEAK_MIB_MAX


def test_motif_scoring_contract_scales_with_eligible_overlap() -> None:
    dataset_sites = deterministic_site_ids(SIGNALOME_CONTRACT_N_SITES)
    offlane_sites = deterministic_site_ids(
        360 * 10 + 64,
        start=850_000,
        gene_prefix="OFFSITE",
    )
    full_reference_map = deterministic_kinase_substrate_map(
        dataset_site_ids=dataset_sites,
        eligible_kinase_count=SIGNALOME_CONTRACT_N_KINASES,
        substrates_per_kinase=12,
        offlane_kinase_count=360,
        offlane_sites_per_kinase=10,
        offlane_site_ids=offlane_sites,
    )
    eligible_kinases = {
        f"KINASE_{index + 1:03d}" for index in range(SIGNALOME_CONTRACT_N_KINASES)
    }
    filtered_map = full_reference_map.loc[
        full_reference_map.loc[:, "kinase"].astype(str).isin(eligible_kinases)
    ]
    sequence_series = deterministic_site_sequence_series(
        dataset_sites.append(offlane_sites),
        window_width=(2 * DEFAULT_MOTIF_FLANK_SIZE) + 1,
    )

    def _run_motif_scoring():
        motif_frequency_matrices, motif_sizes = build_motif_library(
            kinase_substrate_map=filtered_map,
            site_sequences=sequence_series,
            flank_size=DEFAULT_MOTIF_FLANK_SIZE,
        )
        return score_phosphosite_motifs(
            site_sequences=sequence_series.loc[dataset_sites],
            motif_frequency_matrices=motif_frequency_matrices,
            motif_sizes=motif_sizes,
            site_index=dataset_sites,
            min_motif_size=2,
            flank_size=DEFAULT_MOTIF_FLANK_SIZE,
        )

    motif_result, runtime_seconds, peak_mib = measure_runtime_and_peak_mib(
        _run_motif_scoring,
        warmup=True,
    )
    scored_kinases = set(motif_result.motif_scores.columns.astype(str).tolist())

    assert scored_kinases == eligible_kinases
    assert not any(name.startswith("OFFLANE_") for name in scored_kinases)
    assert runtime_seconds < MOTIF_RUNTIME_SECONDS_MAX
    assert peak_mib < MOTIF_PEAK_MIB_MAX


def test_adaptive_prediction_hot_path_performance_contract_with_fixed_seed() -> None:
    prediction_score_matrix, candidate_substrates = (
        _build_adaptive_prediction_hot_path_inputs(
            n_sites=ADAPTIVE_PREDICTION_CONTRACT_N_SITES,
            n_kinases=ADAPTIVE_PREDICTION_CONTRACT_N_KINASES,
            candidate_kinases=ADAPTIVE_PREDICTION_CONTRACT_CANDIDATE_KINASES,
            seed=9417,
        )
    )
    prediction_config = KinasePredictionConfig(
        top_k=ADAPTIVE_PREDICTION_CONTRACT_TOP_K,
        deterministic_max_selected_kinases=12,
        adaptive_ensemble_runs=8,
        mode="adaptive_ensemble",
        n_iterations=12,
        random_state=93103,
    )

    adaptive_scores, runtime_seconds, peak_mib = measure_runtime_and_peak_mib(
        lambda: run_adaptive_ensemble_prediction(
            prediction_score_matrix=prediction_score_matrix,
            candidate_substrates=candidate_substrates,
            prediction_config=prediction_config,
            random_state=93103,
        ),
        warmup=True,
    )
    repeated_scores = run_adaptive_ensemble_prediction(
        prediction_score_matrix=prediction_score_matrix,
        candidate_substrates=candidate_substrates,
        prediction_config=prediction_config,
        random_state=93103,
    )

    pd.testing.assert_frame_equal(adaptive_scores, repeated_scores, check_dtype=False)
    assert adaptive_scores.shape[0] == ADAPTIVE_PREDICTION_CONTRACT_N_SITES
    assert adaptive_scores.shape[1] == ADAPTIVE_PREDICTION_CONTRACT_CANDIDATE_KINASES
    assert runtime_seconds < ADAPTIVE_PREDICTION_RUNTIME_SECONDS_MAX
    assert peak_mib < ADAPTIVE_PREDICTION_PEAK_MIB_MAX


def test_large_reference_map_contract_keeps_filtered_scoring_bounded() -> None:
    _site_ids, dataset, references, eligible_kinases = _build_kinase_workflow_inputs(
        n_sites=250,
        n_samples=8,
        eligible_kinases=42,
        substrates_per_kinase=6,
        offlane_kinases=320,
        offlane_sites_per_kinase=8,
    )
    request = KinaseWorkflowRequest(
        dataset=dataset,
        references=references,
        scoring_config=KinaseScoringConfig(
            min_substrates=2,
            include_diagnostic_scoring_tables=False,
            reference_context_compatibility_policy=(
                ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT
            ),
        ),
        prediction_config=KinasePredictionConfig(
            top_k=6,
            deterministic_max_selected_kinases=10,
            adaptive_ensemble_runs=10,
        ),
        activity_config=None,
    )

    workflow_result, runtime_seconds, peak_mib = measure_runtime_and_peak_mib(
        lambda: KinaseWorkflow().run(request),
        warmup=True,
    )
    rank_weighted_fusion_scores = (
        workflow_result.scoring_result.rank_weighted_fusion_scores
    )
    assert rank_weighted_fusion_scores is not None

    downstream_kinases = set(rank_weighted_fusion_scores.columns.astype(str).tolist())
    assert downstream_kinases.issubset(eligible_kinases)
    assert not any(name.startswith("OFFLANE_") for name in downstream_kinases)
    assert runtime_seconds < KINASE_FILTERED_REFERENCE_RUNTIME_SECONDS_MAX
    assert peak_mib < KINASE_FILTERED_REFERENCE_PEAK_MIB_MAX


def test_diagnostic_scoring_tables_contract_has_bounded_runtime_overhead() -> None:
    _site_ids, dataset, references, _eligible_kinases = _build_kinase_workflow_inputs(
        n_sites=250,
        n_samples=8,
        eligible_kinases=42,
        substrates_per_kinase=6,
        offlane_kinases=320,
        offlane_sites_per_kinase=8,
    )
    workflow = KinaseWorkflow()
    default_request = KinaseWorkflowRequest(
        dataset=dataset,
        references=references,
        scoring_config=KinaseScoringConfig(
            min_substrates=2,
            include_diagnostic_scoring_tables=False,
            reference_context_compatibility_policy=(
                ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT
            ),
        ),
        prediction_config=KinasePredictionConfig(
            top_k=6,
            deterministic_max_selected_kinases=10,
            adaptive_ensemble_runs=10,
        ),
        activity_config=None,
    )
    diagnostic_request = KinaseWorkflowRequest(
        dataset=dataset,
        references=references,
        scoring_config=KinaseScoringConfig(
            min_substrates=2,
            include_diagnostic_scoring_tables=True,
            reference_context_compatibility_policy=(
                ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT
            ),
        ),
        prediction_config=KinasePredictionConfig(
            top_k=6,
            deterministic_max_selected_kinases=10,
            adaptive_ensemble_runs=10,
        ),
        activity_config=None,
    )

    default_runtime = median_runtime_seconds(
        lambda: workflow.run(default_request),
        repeats=3,
        warmup=True,
    )
    diagnostic_runtime = median_runtime_seconds(
        lambda: workflow.run(diagnostic_request),
        repeats=3,
        warmup=True,
    )

    assert diagnostic_runtime <= (
        default_runtime * DIAGNOSTIC_RUNTIME_RATIO_MULTIPLIER
        + DIAGNOSTIC_RUNTIME_ABSOLUTE_SECONDS
    )

    default_result = workflow.run(default_request)
    diagnostic_result = workflow.run(diagnostic_request)
    assert default_result.scoring_result.motif_scores is None
    assert default_result.scoring_result.score_fusion_weights is None
    assert diagnostic_result.scoring_result.motif_scores is not None
    assert diagnostic_result.scoring_result.score_fusion_weights is not None


def test_bundle_publishing_performance_contract_with_representative_tables(
    tmp_path: Path,
) -> None:
    kinase_result, _eligible_kinases = _build_kinase_result_for_signalome_performance(
        n_sites=220,
        n_samples=8,
        eligible_kinases=42,
        substrates_per_kinase=6,
        offlane_kinases=240,
        offlane_sites_per_kinase=8,
    )
    signalome_request = SignalomeWorkflowRequest(
        kinase_result=kinase_result,
        config=build_signalome_config(
            substrate_support_cutoff=0.5,
            module_selection_max_clusters=8,
            candidate_scoring_policy=SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED,
            max_exact_tree_sites=360,
            max_full_candidate_scoring_sites=140,
            clustering_engine=SIGNALOME_CLUSTERING_ENGINE_SCIPY_HIERARCHICAL,
        ),
    )
    signalome_result = SignalomeWorkflow().run(signalome_request)
    config_snapshot = SignalomeWorkflowConfigSnapshot.from_request(signalome_request)
    call_index = [0]

    def _write_bundle():
        bundle_root = tmp_path / f"signalome_contract_bundle_{call_index[0]}"
        call_index[0] += 1
        return save_signalome_workflow_bundle(
            signalome_result,
            bundle_root,
            config_snapshot=config_snapshot,
            output_format="csv",
        )

    written_paths, runtime_seconds, peak_mib = measure_runtime_and_peak_mib(
        _write_bundle,
        warmup=True,
    )
    assert "manifest" in written_paths
    assert "signalome.module_assignments" in written_paths
    assert "signalome.signalome_modules" in written_paths
    assert "prediction.pred_mat" in written_paths
    for written_path in written_paths.values():
        assert Path(written_path).exists()
    assert runtime_seconds < BUNDLE_PUBLISH_RUNTIME_SECONDS_MAX
    assert peak_mib < BUNDLE_PUBLISH_PEAK_MIB_MAX


def test_provenance_hashing_performance_contract_large_dataframes() -> None:
    phospho = deterministic_matrix(
        n_sites=PREPROCESSING_CONTRACT_N_SITES,
        n_samples=PREPROCESSING_CONTRACT_N_SAMPLES,
        seed=9511,
    )
    site_metadata = deterministic_site_metadata(phospho.index, include_protein_id=True)
    site_metadata.loc[:, "category"] = (
        np.arange(int(site_metadata.shape[0]), dtype=int) % 9
    ).astype(str)

    (phospho_hash, metadata_hash), runtime_seconds, peak_mib = (
        measure_runtime_and_peak_mib(
            lambda: (
                hash_table_tolerance(
                    phospho,
                    name="performance_contracts.phospho",
                ),
                hash_table_tolerance(
                    site_metadata,
                    name="performance_contracts.site_metadata",
                ),
            ),
            warmup=True,
        )
    )
    repeated_phospho_hash = hash_table_tolerance(
        phospho,
        name="performance_contracts.phospho",
    )
    repeated_metadata_hash = hash_table_tolerance(
        site_metadata, name="performance_contracts.site_metadata"
    )

    assert phospho_hash == repeated_phospho_hash
    assert metadata_hash == repeated_metadata_hash
    assert len(phospho_hash) == 64
    assert len(metadata_hash) == 64
    assert runtime_seconds < PROVENANCE_HASHING_RUNTIME_SECONDS_MAX
    assert peak_mib < PROVENANCE_HASHING_PEAK_MIB_MAX
