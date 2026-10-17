from __future__ import annotations

import json
import tracemalloc

import numpy as np
import pandas as pd
import pandas.testing as pdt
import pytest

from phospy.api import KinasePredictionResult, KinaseScoringResult
from tests.support.site_keys import protein_site_key_index

pytestmark = [pytest.mark.performance, pytest.mark.release_gate]

_ROWS = 10_000
_COLUMNS = 100
_FLOAT64_BYTES = 8


def test_representative_kinase_exports_count_copies_and_separate_memory_measures(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    site_index = protein_site_key_index(
        protein_identifiers=[f"P{row:05d}" for row in range(_ROWS)],
        sites=["S1"] * _ROWS,
    )
    columns = pd.Index(
        [f"K{column:03d}" for column in range(_COLUMNS)],
        name="kinase",
    )
    values = np.arange(_ROWS * _COLUMNS, dtype=np.float64).reshape(
        _ROWS,
        _COLUMNS,
    )
    values /= float(values.size)
    prediction_input = pd.DataFrame(values, index=site_index, columns=columns)
    score_input = pd.DataFrame(1.0 - values, index=site_index, columns=columns)
    prediction_result = KinasePredictionResult(pred_mat=prediction_input)
    scoring_result = KinaseScoringResult(profile_scores=score_input)

    prediction_getter = KinasePredictionResult.pred_mat.fget
    score_getter = KinaseScoringResult.authoritative_scores.fget
    assert prediction_getter is not None
    assert score_getter is not None
    calls = {"prediction": 0, "score": 0}

    def _prediction_export(result: KinasePredictionResult) -> pd.DataFrame:
        calls["prediction"] += 1
        return prediction_getter(result)

    def _score_export(result: KinaseScoringResult) -> pd.DataFrame:
        calls["score"] += 1
        return score_getter(result)

    monkeypatch.setattr(
        KinasePredictionResult,
        "pred_mat",
        property(_prediction_export),
    )
    monkeypatch.setattr(
        KinaseScoringResult,
        "authoritative_scores",
        property(_score_export),
    )

    tracemalloc.start()
    prediction_export = prediction_result.pred_mat
    score_export = scoring_result.authoritative_scores
    current_bytes, peak_bytes = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    raw_payload_bytes_each = _FLOAT64_BYTES * _ROWS * _COLUMNS
    prediction_pandas_bytes = int(
        prediction_export.memory_usage(index=True, deep=True).sum()
    )
    score_pandas_bytes = int(score_export.memory_usage(index=True, deep=True).sum())
    report = {
        "shape": [_ROWS, _COLUMNS],
        "dtype": "float64",
        "snapshot_calls": calls,
        "raw_payload_bytes_each": raw_payload_bytes_each,
        "raw_payload_bytes_pair": 2 * raw_payload_bytes_each,
        "pandas_bytes": {
            "prediction": prediction_pandas_bytes,
            "score": score_pandas_bytes,
        },
        "pandas_overhead_bytes": {
            "prediction": prediction_pandas_bytes - raw_payload_bytes_each,
            "score": score_pandas_bytes - raw_payload_bytes_each,
        },
        "tracemalloc_current_bytes": int(current_bytes),
        "tracemalloc_peak_bytes": int(peak_bytes),
    }
    print(json.dumps(report, sort_keys=True))

    assert calls == {"prediction": 1, "score": 1}
    assert raw_payload_bytes_each == 8_000_000
    assert prediction_export.to_numpy(copy=False).nbytes == raw_payload_bytes_each
    assert score_export.to_numpy(copy=False).nbytes == raw_payload_bytes_each
    assert prediction_pandas_bytes >= raw_payload_bytes_each
    assert score_pandas_bytes >= raw_payload_bytes_each
    assert report["pandas_overhead_bytes"]["prediction"] >= 0
    assert report["pandas_overhead_bytes"]["score"] >= 0
    assert peak_bytes >= 2 * raw_payload_bytes_each
    pdt.assert_frame_equal(prediction_export, prediction_input)
    pdt.assert_frame_equal(score_export, score_input)
    assert prediction_export is not prediction_input
    assert score_export is not score_input
