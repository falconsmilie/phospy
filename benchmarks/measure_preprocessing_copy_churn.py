from __future__ import annotations

import sys
import time
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))


@dataclass(slots=True)
class CopyCounts:
    dataframe_deep_copies: int = 0
    dataframe_shallow_copies: int = 0
    series_deep_copies: int = 0
    series_shallow_copies: int = 0

    def add(self, other: CopyCounts) -> None:
        self.dataframe_deep_copies += other.dataframe_deep_copies
        self.dataframe_shallow_copies += other.dataframe_shallow_copies
        self.series_deep_copies += other.series_deep_copies
        self.series_shallow_copies += other.series_shallow_copies


@contextmanager
def count_copies():
    counts = CopyCounts()

    original_dataframe_copy = pd.DataFrame.copy
    original_series_copy = pd.Series.copy

    def dataframe_copy(self, *args, **kwargs):
        deep = bool(kwargs.get("deep", True))
        if deep:
            counts.dataframe_deep_copies += 1
        else:
            counts.dataframe_shallow_copies += 1
        return original_dataframe_copy(self, *args, **kwargs)

    def series_copy(self, *args, **kwargs):
        deep = bool(kwargs.get("deep", True))
        if deep:
            counts.series_deep_copies += 1
        else:
            counts.series_shallow_copies += 1
        return original_series_copy(self, *args, **kwargs)

    pd.DataFrame.copy = dataframe_copy
    pd.Series.copy = series_copy
    try:
        yield counts
    finally:
        pd.DataFrame.copy = original_dataframe_copy
        pd.Series.copy = original_series_copy


def _run_with_measurements(*, repeats: int, operation) -> tuple[float, CopyCounts]:
    durations: list[float] = []
    aggregate_counts = CopyCounts()
    for _ in range(repeats):
        with count_copies() as counts:
            start = time.perf_counter()
            operation()
            durations.append(time.perf_counter() - start)
        aggregate_counts.add(counts)
    mean_seconds = sum(durations) / max(len(durations), 1)
    return mean_seconds, aggregate_counts


def _format_counts(counts: CopyCounts, *, repeats: int) -> str:
    return (
        f"df_deep={counts.dataframe_deep_copies // repeats}, "
        f"df_shallow={counts.dataframe_shallow_copies // repeats}, "
        f"series_deep={counts.series_deep_copies // repeats}, "
        f"series_shallow={counts.series_shallow_copies // repeats}"
    )


def _build_signalome_result_for_access_bench():
    from phospy.internal.types import SIGNALOME_MODULE_SELECTION_STRATEGY_SINGLE_MODULE
    from phospy.signalomes.clustering import SignalomeModuleSelectionDiagnostics
    from phospy.signalomes.results import (
        ExpandedSignalome,
        SignalomeAssignments,
        SignalomeKinaseNetwork,
        SignalomeModules,
        SignalomeResult,
    )

    module_table = pd.DataFrame(
        {"KINASE_A": [100.0]},
        index=pd.Index([1], name="module_id"),
    )
    kinase_module_relationships = pd.DataFrame(
        {"kinase": ["KINASE_A"], "module_id": [1], "module_percent": [100.0]}
    )
    site_assignments = pd.DataFrame(
        {
            "protein_id": ["P1"],
            "module_id": [1],
            "top_kinase_candidates": ['["KINASE_A"]'],
            "top_kinase_weights": ['{"KINASE_A": 1.0}'],
            "top_kinase_tie_count": [1],
            "top_kinase_is_ambiguous": [False],
            "top_score": [0.95],
        },
        index=pd.Index(["SITE_1"], name="site_id"),
    )
    protein_assignments = pd.DataFrame(
        {"module_id": [1]},
        index=pd.Index(["P1"], name="protein_id"),
    )
    adjacency = pd.DataFrame(
        {"KINASE_A": [0.0]},
        index=pd.Index(["KINASE_A"], name="kinase"),
    )
    correlation = pd.DataFrame(
        {"KINASE_A": [1.0]},
        index=pd.Index(["KINASE_A"], name="kinase"),
    )
    node_table = pd.DataFrame(
        {"degree": [0], "n_substrates": [1]},
        index=pd.Index(["KINASE_A"], name="kinase"),
    )
    edge_table = pd.DataFrame(columns=["source_kinase", "target_kinase", "correlation"])
    expression_matrix = pd.DataFrame(
        {"sample_1": [1.0]},
        index=pd.Index(["SITE_1"], name="site_id"),
    )
    expanded = {
        "KINASE_A": ExpandedSignalome(
            kinase="KINASE_A",
            linked_kinases=("KINASE_A",),
            regulated_module_ids=(1,),
            expression_matrix=expression_matrix,
            site_assignments=site_assignments,
            row_positions=[0],
        )
    }

    return SignalomeResult(
        scoring_matrix=pd.DataFrame({"KINASE_A": [0.5]}, index=expression_matrix.index),
        pred_mat=pd.DataFrame({"KINASE_A": [0.8]}, index=expression_matrix.index),
        expression_matrix=expression_matrix,
        modules=SignalomeModules(
            module_table=module_table,
            kinase_module_relationships=kinase_module_relationships,
        ),
        assignments=SignalomeAssignments(
            site_assignments=site_assignments,
            protein_assignments=protein_assignments,
        ),
        network=SignalomeKinaseNetwork(
            adjacency_matrix=adjacency,
            correlation_matrix=correlation,
            node_table=node_table,
            edge_table=edge_table,
            neighbor_map={"KINASE_A": tuple()},
        ),
        kinase_substrate_map={"KINASE_A": ("SITE_1",)},
        expanded_signalomes=expanded,
        module_selection_diagnostics=SignalomeModuleSelectionDiagnostics(
            strategy=SIGNALOME_MODULE_SELECTION_STRATEGY_SINGLE_MODULE,
            selected_module_count=1,
            requested_module_count=1,
            threshold_used=None,
            max_clusters_evaluated=1,
            candidate_scores={},
            reason="benchmark",
        ),
    )


def main() -> None:
    from phospy.datasets import (
        AnalysisReadyPhosphoDataset,
        DatasetLoader,
        DatasetSchema,
        PhosphoDataset,
    )
    from phospy.preprocessing import CorePreprocessingConfig

    repeats = 5
    fixture_dir = ROOT / "examples" / "data" / "simple_workflow"
    total_path = fixture_dir / "total.tsv"
    phospho_path = fixture_dir / "phospho.tsv"
    total_df = pd.read_csv(total_path, sep="\t")
    phospho_df = pd.read_csv(phospho_path, sep="\t")
    loader = DatasetLoader(schema=DatasetSchema())

    load_from_files_seconds, load_from_files_counts = _run_with_measurements(
        repeats=repeats,
        operation=lambda: loader.resolve_inputs(total=total_path, phospho=phospho_path),
    )
    print(
        "loader.resolve_inputs(file, file): "
        f"mean_seconds={load_from_files_seconds:.6f}, "
        f"{_format_counts(load_from_files_counts, repeats=repeats)}"
    )

    load_from_memory_seconds, load_from_memory_counts = _run_with_measurements(
        repeats=repeats,
        operation=lambda: loader.resolve_inputs(total=total_df, phospho=phospho_df),
    )
    print(
        "loader.resolve_inputs(df, df): "
        f"mean_seconds={load_from_memory_seconds:.6f}, "
        f"{_format_counts(load_from_memory_counts, repeats=repeats)}"
    )

    dataset = PhosphoDataset.from_files(
        phospho_path=phospho_path,
        total_path=total_path,
    )
    preprocessing_seconds, preprocessing_counts = _run_with_measurements(
        repeats=repeats,
        operation=lambda: dataset.preprocessing.run(config=CorePreprocessingConfig()),
    )
    print(
        "dataset.preprocessing.run(config): "
        f"mean_seconds={preprocessing_seconds:.6f}, "
        f"{_format_counts(preprocessing_counts, repeats=repeats)}"
    )

    core_result = dataset.preprocessing.run(config=CorePreprocessingConfig())
    analysis_ready_seconds, analysis_ready_counts = _run_with_measurements(
        repeats=repeats,
        operation=lambda: AnalysisReadyPhosphoDataset.from_core_processing_result(
            core_result,
            schema=dataset.schema,
            comparisons=dataset.comparisons,
        ),
    )
    print(
        "AnalysisReadyPhosphoDataset.from_core_processing_result(...): "
        f"mean_seconds={analysis_ready_seconds:.6f}, "
        f"{_format_counts(analysis_ready_counts, repeats=repeats)}"
    )

    signalome_result = _build_signalome_result_for_access_bench()
    signalome_access_seconds, signalome_access_counts = _run_with_measurements(
        repeats=repeats,
        operation=lambda: (
            signalome_result.to_frames(include_inputs=True),
            signalome_result.expanded_signalomes,
        ),
    )
    print(
        "SignalomeResult.to_frames(include_inputs=True) + expanded_signalomes: "
        f"mean_seconds={signalome_access_seconds:.6f}, "
        f"{_format_counts(signalome_access_counts, repeats=repeats)}"
    )


if __name__ == "__main__":
    main()
