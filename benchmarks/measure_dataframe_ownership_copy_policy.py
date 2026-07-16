#!/usr/bin/env python3
"""Benchmark DataFrame ownership borrow-copy policy.

Targets:
- `phospy.frames.ownership._borrow_dataframe`
"""

from __future__ import annotations

import sys
import time
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))


@dataclass(slots=True)
class _CopyCounts:
    dataframe_deep: int = 0
    dataframe_shallow: int = 0


@contextmanager
def _count_dataframe_copies() -> Iterator[_CopyCounts]:
    counts = _CopyCounts()
    original_copy = pd.DataFrame.copy

    def wrapped_copy(self, *args, **kwargs):  # type: ignore[no-untyped-def]
        deep = kwargs.get("deep", args[0] if args else True)
        if bool(deep):
            counts.dataframe_deep += 1
        else:
            counts.dataframe_shallow += 1
        return original_copy(self, *args, **kwargs)

    pd.DataFrame.copy = wrapped_copy  # pyright: ignore[reportAttributeAccessIssue] - benchmark instrumentation monkeypatches pandas copy counts.
    try:
        yield counts
    finally:
        pd.DataFrame.copy = original_copy  # pyright: ignore[reportAttributeAccessIssue] - restore benchmark instrumentation.


@dataclass(frozen=True, slots=True)
class _BorrowMetrics:
    runtime_seconds: float
    deep_copies: int
    shallow_copies: int
    owner_mutation_leaks: int


def _numeric_frame() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample_a": [1.0, 2.0, 3.0, 4.0],
            "sample_b": [2.0, 3.0, 4.0, 5.0],
            "sample_c": [3.0, 4.0, 5.0, 6.0],
        },
        index=pd.Index(["site_1", "site_2", "site_3", "site_4"], name="site_id"),
    )


def _extension_frame() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "site_id": pd.Series(
                ["site_1", "site_2", "site_3", "site_4"], dtype="string"
            ),
            "gene_symbol": pd.Series(
                ["MAPK14", "GSK3B", "AKT1", "MTOR"], dtype="string"
            ),
        }
    )


def _exercise_borrow(
    frame: pd.DataFrame,
    *,
    repeats: int,
    replacement: object,
) -> _BorrowMetrics:
    from phospy.frames.ownership import _borrow_dataframe

    original_value = frame.iloc[0, 0]
    owner_mutation_leaks = 0
    with _count_dataframe_copies() as counts:
        started = time.perf_counter()
        for _ in range(repeats):
            borrowed = _borrow_dataframe(frame)
            try:
                borrowed.iloc[0, 0] = replacement  # pyright: ignore[reportCallIssue,reportArgumentType] - benchmark intentionally probes pandas mutation behavior.
            except ValueError:
                pass
            if frame.iloc[0, 0] != original_value:
                owner_mutation_leaks += 1
        runtime_seconds = time.perf_counter() - started

    return _BorrowMetrics(
        runtime_seconds=runtime_seconds,
        deep_copies=counts.dataframe_deep,
        shallow_copies=counts.dataframe_shallow,
        owner_mutation_leaks=owner_mutation_leaks,
    )


def main() -> None:
    repeats = 200
    numeric_metrics = _exercise_borrow(
        _numeric_frame(),
        repeats=repeats,
        replacement=999.0,
    )
    extension_metrics = _exercise_borrow(
        _extension_frame(),
        repeats=repeats,
        replacement="changed",
    )

    print("benchmark_suite=dataframe_ownership_copy_policy_v1")
    print(f"pandas_version={pd.__version__}")
    print(f"repeats={repeats}")
    print(f"numeric_borrow_runtime_seconds={numeric_metrics.runtime_seconds:.6f}")
    print(f"numeric_borrow_deep_copies={numeric_metrics.deep_copies}")
    print(f"numeric_borrow_shallow_copies={numeric_metrics.shallow_copies}")
    print(f"numeric_owner_mutation_leaks={numeric_metrics.owner_mutation_leaks}")
    print(f"extension_borrow_runtime_seconds={extension_metrics.runtime_seconds:.6f}")
    print(f"extension_borrow_deep_copies={extension_metrics.deep_copies}")
    print(f"extension_borrow_shallow_copies={extension_metrics.shallow_copies}")
    print(f"extension_owner_mutation_leaks={extension_metrics.owner_mutation_leaks}")


if __name__ == "__main__":
    main()
