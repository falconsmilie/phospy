from __future__ import annotations

import sys
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
    dataframe_copies: int = 0
    series_copies: int = 0


@contextmanager
def count_copies():
    counts = CopyCounts()

    original_dataframe_copy = pd.DataFrame.copy
    original_series_copy = pd.Series.copy

    def dataframe_copy(self, *args, **kwargs):
        counts.dataframe_copies += 1
        return original_dataframe_copy(self, *args, **kwargs)

    def series_copy(self, *args, **kwargs):
        counts.series_copies += 1
        return original_series_copy(self, *args, **kwargs)

    pd.DataFrame.copy = dataframe_copy
    pd.Series.copy = series_copy
    try:
        yield counts
    finally:
        pd.DataFrame.copy = original_dataframe_copy
        pd.Series.copy = original_series_copy


def main() -> None:
    from phospy import AnalysisReadyPhosphoDataset, PhosphoDataset

    dataset = PhosphoDataset.from_files(
        phospho_path=ROOT / "examples" / "data" / "phospho.tsv",
        total_path=ROOT / "examples" / "data" / "total.tsv",
        comparisons=["Basal", "Treat"],
    )

    with count_copies() as preprocessing_counts:
        core = dataset.preprocessing.run()

    print(
        "dataset.preprocessing.run(): "
        f"DataFrame.copy={preprocessing_counts.dataframe_copies}, "
        f"Series.copy={preprocessing_counts.series_copies}"
    )

    with count_copies() as analysis_ready_counts:
        AnalysisReadyPhosphoDataset.from_core_processing_result(core)

    print(
        "AnalysisReadyPhosphoDataset.from_core_processing_result(): "
        f"DataFrame.copy={analysis_ready_counts.dataframe_copies}, "
        f"Series.copy={analysis_ready_counts.series_copies}"
    )


if __name__ == "__main__":
    main()
