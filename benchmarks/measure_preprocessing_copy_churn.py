#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import platform
import sys
import time
from collections.abc import Callable
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path

import numpy as np
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

    @property
    def dataframe_total(self) -> int:
        return self.dataframe_deep_copies + self.dataframe_shallow_copies

    @property
    def series_total(self) -> int:
        return self.series_deep_copies + self.series_shallow_copies

    def per_run(self, repeats: int) -> dict[str, int]:
        safe_repeats = max(1, repeats)
        return {
            "dataframe_deep": self.dataframe_deep_copies // safe_repeats,
            "dataframe_shallow": self.dataframe_shallow_copies // safe_repeats,
            "dataframe_total": self.dataframe_total // safe_repeats,
            "series_deep": self.series_deep_copies // safe_repeats,
            "series_shallow": self.series_shallow_copies // safe_repeats,
            "series_total": self.series_total // safe_repeats,
        }


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


def _run_with_measurements(
    *,
    repeats: int,
    warmups: int,
    operation_factory: Callable[[], Callable[[], object]],
) -> tuple[float, CopyCounts]:
    durations: list[float] = []
    aggregate_counts = CopyCounts()

    for _ in range(max(0, warmups)):
        operation_factory()()

    for _ in range(repeats):
        operation = operation_factory()
        with count_copies() as counts:
            start = time.perf_counter()
            operation()
            durations.append(time.perf_counter() - start)
        aggregate_counts.add(counts)

    mean_seconds = sum(durations) / max(len(durations), 1)
    return mean_seconds, aggregate_counts


def _build_preprocessing_inputs(
    *,
    schema,
    config,
    n_genes: int,
    sites_per_gene: int,
    random_state: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(random_state)

    genes = np.asarray([f"GENE{i + 1}" for i in range(n_genes)], dtype=object)
    total_values = rng.normal(
        loc=25.0,
        scale=5.0,
        size=(n_genes, len(schema.total_cols)),
    )
    total_df = pd.DataFrame(total_values, columns=schema.total_cols)
    total_df.insert(0, "genes", genes)

    n_sites = n_genes * sites_per_gene
    site_gene = np.repeat(genes, sites_per_gene)
    site_ids = np.asarray(
        [
            f"{gene}_S{(index % sites_per_gene) + 1}"
            for index, gene in enumerate(site_gene)
        ],
        dtype=object,
    )
    phospho_values = rng.normal(
        loc=12.0,
        scale=3.0,
        size=(n_sites, len(schema.phospho_cols)),
    )
    phospho_values[rng.random(phospho_values.shape) < 0.03] = config.phospho_sentinel
    localization_prob = rng.uniform(0.5, 1.0, size=n_sites)
    low_mask = rng.random(n_sites) < 0.1
    localization_prob[low_mask] = rng.uniform(0.0, 0.7, size=low_mask.sum())

    phospho_df = pd.DataFrame(phospho_values, columns=schema.phospho_cols)
    phospho_df.insert(
        0, "centralized_sequence", ["RARTSSFAEPGGGGGGGGGPGGSASPARPAR"] * n_sites
    )
    phospho_df.insert(0, "localization_prob", localization_prob)
    phospho_df.insert(0, "gene_p_site", site_ids)
    phospho_df.insert(0, "gene_names", site_gene)
    phospho_df.insert(0, "uid", [f"UID_{index + 1}" for index in range(n_sites)])
    return total_df, phospho_df


def _build_markdown(report: dict[str, object]) -> str:
    benchmarks = report["benchmarks"]
    guards = report["guards"]
    lines = [
        "# Preprocessing owned-path benchmark",
        "",
        "This report tracks the preprocessing ownership model:",
        "",
        "- public safe path (`CoreProcessor.process`)",
        "- owned fast path (`CoreProcessor.process_owned`)",
        "- large-matrix owned preprocessing",
        "",
        "## Benchmark cases",
        "",
        "| Case | Mean seconds | Rows (total/phospho) | DataFrame copies/run | Series copies/run |",
        "| --- | ---: | ---: | ---: | ---: |",
    ]
    for benchmark_name in (
        "public_safe_path",
        "owned_fast_path",
        "large_matrix_owned_path",
    ):
        entry = benchmarks[benchmark_name]
        lines.append(
            "| "
            f"`{benchmark_name}` | "
            f"{entry['mean_seconds']:.6f} | "
            f"{entry['rows']['total']}/{entry['rows']['phospho']} | "
            f"{entry['copies_per_run']['dataframe_total']} "
            f"(deep={entry['copies_per_run']['dataframe_deep']}, "
            f"shallow={entry['copies_per_run']['dataframe_shallow']}) | "
            f"{entry['copies_per_run']['series_total']} "
            f"(deep={entry['copies_per_run']['series_deep']}, "
            f"shallow={entry['copies_per_run']['series_shallow']}) |"
        )

    lines.extend(
        [
            "",
            "## Guard checks",
            "",
            f"- all checks passed: `{guards['all_passed']}`",
            f"- public boundary adds copy overhead: `{guards['public_boundary_adds_expected_dataframe_copies']}`",
            f"- owned path dataframe copy budget: `{guards['owned_path_within_dataframe_copy_budget']}`",
            f"- large-matrix path dataframe copy budget: `{guards['large_path_within_dataframe_copy_budget']}`",
            f"- owned path series copy budget: `{guards['owned_path_within_series_copy_budget']}`",
            "",
            "## Run locally",
            "",
            "```bash",
            "python benchmarks/measure_preprocessing_copy_churn.py --check",
            "```",
            "",
        ]
    )
    return "\n".join(lines)


def _run_benchmarks(
    *,
    repeats: int,
    warmups: int,
    small_genes: int,
    small_sites_per_gene: int,
    large_genes: int,
    large_sites_per_gene: int,
    min_public_minus_owned_dataframe_copies: int,
    max_owned_dataframe_copies_per_run: int,
    max_large_dataframe_copies_per_run: int,
    max_owned_series_copies_per_run: int,
) -> dict[str, object]:
    from phospy.datasets import DatasetSchema
    from phospy.preprocessing import CorePreprocessingConfig, CoreProcessor

    schema = DatasetSchema()
    config = CorePreprocessingConfig()
    processor = CoreProcessor(schema=schema)

    total_small, phospho_small = _build_preprocessing_inputs(
        schema=schema,
        config=config,
        n_genes=small_genes,
        sites_per_gene=small_sites_per_gene,
        random_state=7,
    )
    total_large, phospho_large = _build_preprocessing_inputs(
        schema=schema,
        config=config,
        n_genes=large_genes,
        sites_per_gene=large_sites_per_gene,
        random_state=11,
    )

    def public_operation_factory() -> Callable[[], object]:
        total_df = total_small.copy(deep=True)
        phospho_df = phospho_small.copy(deep=True)

        def operation() -> object:
            return processor.process(total_df, phospho_df, config=config)

        return operation

    def owned_operation_factory() -> Callable[[], object]:
        total_df = total_small.copy(deep=True)
        phospho_df = phospho_small.copy(deep=True)

        def operation() -> object:
            return processor.process_owned(total_df, phospho_df, config=config)

        return operation

    def large_operation_factory() -> Callable[[], object]:
        total_df = total_large.copy(deep=True)
        phospho_df = phospho_large.copy(deep=True)

        def operation() -> object:
            return processor.process_owned(total_df, phospho_df, config=config)

        return operation

    public_seconds, public_counts = _run_with_measurements(
        repeats=repeats,
        warmups=warmups,
        operation_factory=public_operation_factory,
    )
    owned_seconds, owned_counts = _run_with_measurements(
        repeats=repeats,
        warmups=warmups,
        operation_factory=owned_operation_factory,
    )
    large_seconds, large_counts = _run_with_measurements(
        repeats=repeats,
        warmups=warmups,
        operation_factory=large_operation_factory,
    )

    public_per_run = public_counts.per_run(repeats)
    owned_per_run = owned_counts.per_run(repeats)
    large_per_run = large_counts.per_run(repeats)

    benchmark_results = {
        "public_safe_path": {
            "mean_seconds": public_seconds,
            "rows": {"total": len(total_small), "phospho": len(phospho_small)},
            "copies_per_run": public_per_run,
        },
        "owned_fast_path": {
            "mean_seconds": owned_seconds,
            "rows": {"total": len(total_small), "phospho": len(phospho_small)},
            "copies_per_run": owned_per_run,
        },
        "large_matrix_owned_path": {
            "mean_seconds": large_seconds,
            "rows": {"total": len(total_large), "phospho": len(phospho_large)},
            "copies_per_run": large_per_run,
        },
    }

    guards = {
        "public_boundary_adds_expected_dataframe_copies": (
            public_per_run["dataframe_total"] - owned_per_run["dataframe_total"]
            >= min_public_minus_owned_dataframe_copies
        ),
        "owned_path_within_dataframe_copy_budget": (
            owned_per_run["dataframe_total"] <= max_owned_dataframe_copies_per_run
        ),
        "large_path_within_dataframe_copy_budget": (
            large_per_run["dataframe_total"] <= max_large_dataframe_copies_per_run
        ),
        "owned_path_within_series_copy_budget": (
            owned_per_run["series_total"] <= max_owned_series_copies_per_run
        ),
    }
    guards["all_passed"] = all(guards.values())

    return {
        "repeats": repeats,
        "warmups": warmups,
        "environment": {
            "python_version": platform.python_version(),
            "platform": platform.platform(),
            "pandas_version": pd.__version__,
        },
        "inputs": {
            "small": {
                "genes": small_genes,
                "sites_per_gene": small_sites_per_gene,
            },
            "large": {
                "genes": large_genes,
                "sites_per_gene": large_sites_per_gene,
            },
        },
        "thresholds": {
            "min_public_minus_owned_dataframe_copies": min_public_minus_owned_dataframe_copies,
            "max_owned_dataframe_copies_per_run": max_owned_dataframe_copies_per_run,
            "max_large_dataframe_copies_per_run": max_large_dataframe_copies_per_run,
            "max_owned_series_copies_per_run": max_owned_series_copies_per_run,
        },
        "benchmarks": benchmark_results,
        "guards": guards,
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Benchmark preprocessing copy churn for public safe paths and owned fast paths."
        )
    )
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--warmups", type=int, default=1)
    parser.add_argument("--small-genes", type=int, default=1000)
    parser.add_argument("--small-sites-per-gene", type=int, default=2)
    parser.add_argument("--large-genes", type=int, default=6000)
    parser.add_argument("--large-sites-per-gene", type=int, default=3)
    parser.add_argument(
        "--min-public-minus-owned-dataframe-copies",
        type=int,
        default=2,
    )
    parser.add_argument("--max-owned-dataframe-copies-per-run", type=int, default=40)
    parser.add_argument("--max-large-dataframe-copies-per-run", type=int, default=40)
    parser.add_argument("--max-owned-series-copies-per-run", type=int, default=8)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("benchmarks/reports/latest"),
        help=(
            "Directory for preprocessing_copy_churn.json and "
            "preprocessing_copy_churn.md"
        ),
    )
    parser.add_argument(
        "--stdout-only",
        action="store_true",
        help="Print the JSON report to stdout instead of writing files.",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="Exit with status 1 when copy-churn guard checks fail.",
    )
    args = parser.parse_args()

    report = _run_benchmarks(
        repeats=max(1, args.repeats),
        warmups=max(0, args.warmups),
        small_genes=max(1, args.small_genes),
        small_sites_per_gene=max(1, args.small_sites_per_gene),
        large_genes=max(1, args.large_genes),
        large_sites_per_gene=max(1, args.large_sites_per_gene),
        min_public_minus_owned_dataframe_copies=max(
            0, args.min_public_minus_owned_dataframe_copies
        ),
        max_owned_dataframe_copies_per_run=max(
            0, args.max_owned_dataframe_copies_per_run
        ),
        max_large_dataframe_copies_per_run=max(
            0, args.max_large_dataframe_copies_per_run
        ),
        max_owned_series_copies_per_run=max(0, args.max_owned_series_copies_per_run),
    )

    if args.stdout_only:
        print(json.dumps(report, indent=2, sort_keys=True))
    else:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        json_path = args.output_dir / "preprocessing_copy_churn.json"
        markdown_path = args.output_dir / "preprocessing_copy_churn.md"
        json_path.write_text(
            json.dumps(report, indent=2, sort_keys=True), encoding="utf-8"
        )
        markdown_path.write_text(_build_markdown(report), encoding="utf-8")
        print(f"Wrote {json_path}")
        print(f"Wrote {markdown_path}")

    if args.check and not bool(report["guards"]["all_passed"]):
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
