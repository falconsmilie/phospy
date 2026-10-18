from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from tests.support import performance_contracts
from tests.support.performance_contracts import (
    measure_runtime_and_peak_mib,
    with_missing_fraction,
)


def test_with_missing_fraction_sets_requested_cells_for_fortran_backed_frame() -> None:
    values = np.asfortranarray(np.arange(24, dtype=float).reshape(4, 6))
    matrix = pd.DataFrame(values)

    result = with_missing_fraction(matrix, missing_fraction=0.25, seed=123)

    assert int(result.isna().sum().sum()) == 6
    assert int(matrix.isna().sum().sum()) == 0
    assert result.shape == matrix.shape
    assert result.index.equals(matrix.index)
    assert result.columns.equals(matrix.columns)


def test_runtime_and_memory_measure_same_callable_in_independent_executions(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls: list[tuple[int, bool]] = []
    perf_counter_values = iter((10.0, 12.5))
    tracing = False

    def workload() -> int:
        calls.append((len(calls) + 1, tracing))
        return len(calls)

    def start_tracing() -> None:
        nonlocal tracing
        tracing = True

    def stop_tracing() -> None:
        nonlocal tracing
        tracing = False

    monkeypatch.setattr(
        performance_contracts.time,
        "perf_counter",
        lambda: next(perf_counter_values),
    )
    monkeypatch.setattr(
        performance_contracts.tracemalloc,
        "start",
        start_tracing,
    )
    monkeypatch.setattr(
        performance_contracts.tracemalloc,
        "stop",
        stop_tracing,
    )
    monkeypatch.setattr(
        performance_contracts.tracemalloc,
        "get_traced_memory",
        lambda: (128, 3 * 1024 * 1024),
    )

    result, runtime_seconds, peak_mib = measure_runtime_and_peak_mib(
        workload,
        warmup=True,
    )

    assert calls == [(1, False), (2, False), (3, True)]
    assert result == 2
    assert runtime_seconds == 2.5
    assert peak_mib == 3.0
