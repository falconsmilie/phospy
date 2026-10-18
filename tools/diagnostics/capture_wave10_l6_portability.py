"""Capture bounded Wave 10 L6 portability evidence.

This module is diagnostic instrumentation.  It wraps production seams only for
the lifetime of this process and executes the normal public L6 workflow entry
point from the committed fixtures.
"""

from __future__ import annotations

import argparse
import contextlib
import dataclasses
import hashlib
import importlib.metadata
import io
import json
import os
import platform
import subprocess
import sys
import traceback
from collections.abc import Iterator, Mapping, Sequence
from pathlib import Path
from typing import Any, cast
from unittest.mock import patch

import numpy as np
import pandas as pd
from scipy import stats
from threadpoolctl import threadpool_info

from phospy.kinase import executor as kinase_executor
from phospy.kinase.prediction import execution as prediction_execution
from phospy.kinase.prediction import runner as prediction_runner
from phospy.kinase.prediction import sampling_core
from phospy.run_lineage.serialization import to_payload as provenance_to_payload
from tests.support import l6_prediction_parity_metrics as l6_metrics
from tests.support.l6_prediction_sensitivity import (
    derive_l6_prediction_sensitivity_budget,
)
from tests.support.rewrite_fixture_data import (
    load_l6_prediction_reference_candidate_substrates,
    load_l6_prediction_reference_predmat,
    load_l6_prediction_reference_profile_scores,
    load_l6_prediction_reference_rank_weighted_fusion_scores,
    load_l6_prediction_reference_score_fusion_weights,
    load_l6_prediction_reference_top30,
)

SCHEMA = "phospy-wave10-l6-portability-diagnostic-v1"
FLOAT_ENCODING = "IEEE-754 binary64 hexadecimal via float.hex"
TIE_WIDTHS = (0.0, 1e-15, 1e-12, 1e-9)
VALUE_THRESHOLDS = (0.5,)
REQUIRED_SKX_FEATURES = ("AVX512F", "AVX512CD", "AVX512_SKX")


def _json_bytes(value: object) -> bytes:
    return (
        json.dumps(
            value,
            ensure_ascii=False,
            allow_nan=False,
            sort_keys=True,
            separators=(",", ":"),
        )
        + "\n"
    ).encode("utf-8")


def _write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(_json_bytes(value))


def _sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def _git(*args: str) -> str:
    return subprocess.check_output(["git", *args], text=True, encoding="utf-8").strip()


def _source_identity() -> dict[str, object]:
    return {
        "commit": _git("rev-parse", "HEAD"),
        "tree": _git("show", "-s", "--format=%T", "HEAD"),
        "parents": _git("show", "-s", "--format=%P", "HEAD").split(),
        "branch": _git("branch", "--show-current"),
        "github": {
            key.lower(): os.environ.get(key)
            for key in (
                "GITHUB_REPOSITORY",
                "GITHUB_RUN_ID",
                "GITHUB_RUN_ATTEMPT",
                "GITHUB_JOB",
                "GITHUB_WORKFLOW",
                "GITHUB_REF",
                "RUNNER_NAME",
                "RUNNER_OS",
                "RUNNER_ARCH",
            )
        },
    }


def _public_source(source: Mapping[str, object]) -> dict[str, object]:
    github = cast(Mapping[str, object], source.get("github", {}))
    return {key: source.get(key) for key in ("commit", "tree", "parents", "branch")} | {
        "github": {
            key: github.get(key)
            for key in (
                "github_repository",
                "github_run_id",
                "github_run_attempt",
                "github_job",
                "github_workflow",
                "github_ref",
            )
        }
    }


def _public_environment(environment: Mapping[str, object]) -> dict[str, object]:
    operating_system = cast(Mapping[str, object], environment["operating_system"])
    python = cast(Mapping[str, object], environment["python"])
    numpy_runtime = cast(Mapping[str, object], environment["numpy_runtime"])
    cpu_features = cast(Mapping[str, object], numpy_runtime["cpu_features"])
    threadpools = cast(list[Mapping[str, object]], environment["threadpools"])
    return {
        "source": _public_source(cast(Mapping[str, object], environment["source"])),
        "operating_system": {
            key: operating_system.get(key)
            for key in ("system", "release", "version", "platform", "architecture")
        },
        "python": {key: python.get(key) for key in ("version", "implementation")},
        "dependencies": environment["dependencies"],
        "threadpools": [
            {
                key: pool.get(key)
                for key in (
                    "internal_api",
                    "user_api",
                    "prefix",
                    "version",
                    "architecture",
                    "threading_layer",
                    "num_threads",
                )
            }
            for pool in threadpools
        ],
        "numpy_cpu_features": {
            key: bool(cpu_features.get(key, False))
            for key in ("AVX", "AVX2", "AVX512F", "AVX512CD", "AVX512_SKX", "FMA3")
        },
    }


def _numpy_runtime() -> dict[str, object]:
    config: object
    try:
        config = np.__config__.show(mode="dicts")
    except TypeError:
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            np.__config__.show()
        config = {"text": buffer.getvalue()}
    try:
        cpu_features = dict(np._core._multiarray_umath.__cpu_features__)  # type: ignore[attr-defined]
    except AttributeError:
        cpu_features = {}
    runtime_buffer = io.StringIO()
    with contextlib.redirect_stdout(runtime_buffer):
        np.show_runtime()
    return {
        "configuration": config,
        "show_runtime": runtime_buffer.getvalue().splitlines(),
        "cpu_features": {
            str(key): bool(value) for key, value in sorted(cpu_features.items())
        },
    }


def runtime_identity() -> dict[str, object]:
    distributions = (
        "numpy",
        "pandas",
        "scipy",
        "scikit-learn",
        "threadpoolctl",
    )
    return {
        "schema": SCHEMA,
        "metadata": {"encoding": "UTF-8", "line_endings": "LF"},
        "source": _source_identity(),
        "operating_system": {
            "system": platform.system(),
            "release": platform.release(),
            "version": platform.version(),
            "platform": platform.platform(),
            "architecture": platform.machine(),
            "processor": platform.processor(),
        },
        "python": {
            "version": platform.python_version(),
            "full_version": sys.version,
            "implementation": platform.python_implementation(),
            "executable": sys.executable,
        },
        "dependencies": {
            name: importlib.metadata.version(name) for name in distributions
        },
        "threadpools": threadpool_info(),
        "numpy_runtime": _numpy_runtime(),
        "environment": {
            name: os.environ.get(name)
            for name in (
                "OPENBLAS_CORETYPE",
                "OMP_NUM_THREADS",
                "OPENBLAS_NUM_THREADS",
                "MKL_NUM_THREADS",
                "VECLIB_MAXIMUM_THREADS",
                "NUMEXPR_NUM_THREADS",
            )
        },
    }


def _normalise_scalar(value: object) -> object:
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        if not np.isfinite(value):
            return {"non_finite": str(value)}
        return value
    if isinstance(value, np.generic):
        return _normalise_scalar(value.item())
    if dataclasses.is_dataclass(value):
        return {
            field.name: _normalise_scalar(getattr(value, field.name))
            for field in dataclasses.fields(value)
        }
    if isinstance(value, Mapping):
        return {str(key): _normalise_scalar(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_normalise_scalar(item) for item in value]
    return str(value)


def array_payload(values: object) -> dict[str, object]:
    array = np.asarray(values)
    if np.issubdtype(array.dtype, np.number):
        floats = np.asarray(array, dtype="<f8", order="C")
        flat = floats.reshape(-1)
        hex_values = [float(item).hex() for item in flat]
        canonical = _json_bytes(
            {
                "shape": list(floats.shape),
                "values_hex": hex_values,
            }
        )
        return {
            "shape": list(floats.shape),
            "dtype": "<f8",
            "value_encoding": FLOAT_ENCODING,
            "values_hex": hex_values,
            "binary64_bytes_sha256": _sha256_bytes(floats.tobytes(order="C")),
            "canonical_sha256": _sha256_bytes(canonical),
        }
    strings = np.asarray(array, dtype=str, order="C")
    values_list = strings.reshape(-1).tolist()
    canonical = _json_bytes({"shape": list(strings.shape), "values": values_list})
    return {
        "shape": list(strings.shape),
        "dtype": "str",
        "values": values_list,
        "canonical_sha256": _sha256_bytes(canonical),
    }


def frame_payload(frame: pd.DataFrame) -> dict[str, object]:
    payload = {
        "row_ids": [str(value) for value in frame.index],
        "column_ids": [str(value) for value in frame.columns],
        "values": array_payload(frame.to_numpy(dtype=float, copy=False)),
    }
    return {
        **payload,
        "matrix_sha256": _sha256_bytes(_json_bytes(payload)),
    }


def frame_from_payload(payload: Mapping[str, object]) -> pd.DataFrame:
    values_payload = cast(Mapping[str, object], payload["values"])
    shape = tuple(cast(list[int], values_payload["shape"]))
    values = np.asarray(
        [
            float.fromhex(value)
            for value in cast(list[str], values_payload["values_hex"])
        ],
        dtype=float,
    ).reshape(shape)
    return pd.DataFrame(
        values,
        index=cast(list[str], payload["row_ids"]),
        columns=cast(list[str], payload["column_ids"]),
    )


def _row_identities(values: np.ndarray, universe: pd.DataFrame | None) -> list[str]:
    universe_lookup: dict[str, list[str]] = {}
    if universe is not None:
        for row_id, row in zip(
            universe.index,
            universe.to_numpy(dtype="<f8", copy=False),
            strict=True,
        ):
            digest = _sha256_bytes(np.asarray(row, dtype="<f8").tobytes())
            universe_lookup.setdefault(digest, []).append(str(row_id))
    identities: list[str] = []
    for row in np.asarray(values, dtype="<f8"):
        digest = _sha256_bytes(np.asarray(row, dtype="<f8").tobytes())
        matches = universe_lookup.get(digest, [])
        suffix = "|".join(matches) if matches else "unmapped"
        identities.append(f"sha256:{digest}:{suffix}")
    return identities


class TraceState:
    def __init__(self) -> None:
        self.checkpoints: list[dict[str, object]] = []
        self.arrays: dict[str, dict[str, object]] = {}
        self.identity_sets: dict[str, list[str]] = {}
        self.current_kinase = "unassigned"
        self.current_ensemble = -1
        self.ensemble_counts: dict[str, int] = {}
        self.current_model: dict[str, object] | None = None
        self.universe: pd.DataFrame | None = None
        self.candidates: dict[str, list[str]] = {}
        self.kinase_ranking: list[dict[str, object]] = []
        self.selected_kinases: list[str] = []
        self.execution_config: dict[str, object] = {}

    def array(self, values: object) -> dict[str, object]:
        payload = array_payload(values)
        digest = cast(str, payload["canonical_sha256"])
        self.arrays.setdefault(digest, payload)
        return {
            "array_ref": digest,
            "shape": payload["shape"],
            "dtype": payload["dtype"],
        }

    def identities(self, values: np.ndarray) -> dict[str, object]:
        identities = _row_identities(values, self.universe)
        digest = _sha256_bytes(_json_bytes(identities))
        self.identity_sets.setdefault(digest, identities)
        return {"identity_set_ref": digest, "count": len(identities)}

    def checkpoint(
        self,
        kind: str,
        payload: object,
        *,
        iteration: int | None = None,
    ) -> dict[str, object]:
        record = {
            "sequence": len(self.checkpoints),
            "kind": kind,
            "kinase": self.current_kinase,
            "ensemble": self.current_ensemble,
            "iteration": iteration,
            "payload": payload,
            "payload_sha256": _sha256_bytes(_json_bytes(payload)),
        }
        self.checkpoints.append(record)
        return record


class TracingModel:
    def __init__(self, model: Any, trace: TraceState) -> None:
        self._model = model
        self._trace = trace
        self._iteration = sum(
            1
            for item in trace.checkpoints
            if item["kind"] == "svc_fit_input"
            and item["kinase"] == trace.current_kinase
            and item["ensemble"] == trace.current_ensemble
        )

    @property
    def classes_(self) -> np.ndarray:
        return np.asarray(self._model.classes_)  # type: ignore[attr-defined]

    def fit(self, x: np.ndarray, y: np.ndarray) -> None:
        input_payload = {
            "values": self._trace.array(x),
            "labels": self._trace.array(y),
            "sample_identities": self._trace.identities(x),
        }
        record = self._trace.checkpoint(
            "svc_fit_input", input_payload, iteration=self._iteration
        )
        self._trace.current_model = record
        self._model.fit(x, y)  # type: ignore[attr-defined]
        svc = self._model.steps[-1][1]  # type: ignore[attr-defined]
        support_indices = np.asarray(svc.support_, dtype=int)
        identity_ref = cast(dict[str, object], input_payload["sample_identities"])
        row_ids = self._trace.identity_sets[cast(str, identity_ref["identity_set_ref"])]
        fitted_payload = {
            "classes": self._trace.array(svc.classes_),
            "support_indices": support_indices.tolist(),
            "support_identities": [row_ids[int(index)] for index in support_indices],
            "support_vectors": self._trace.array(svc.support_vectors_),
            "support_counts": self._trace.array(svc.n_support_),
            "dual_coef": self._trace.array(svc.dual_coef_),
            "intercept": self._trace.array(svc.intercept_),
            "prob_a": self._trace.array(svc.probA_),
            "prob_b": self._trace.array(svc.probB_),
        }
        self._trace.checkpoint("svc_fitted", fitted_payload, iteration=self._iteration)

    def predict_proba(self, x: np.ndarray) -> np.ndarray:
        probabilities = np.asarray(
            self._model.predict_proba(x),
            dtype=float,  # type: ignore[attr-defined]
        )
        decisions = np.asarray(
            self._model.decision_function(x),
            dtype=float,  # type: ignore[attr-defined]
        )
        self._trace.checkpoint(
            "svc_evaluation",
            {
                "input": self._trace.array(x),
                "sample_identities": self._trace.identities(x),
                "decision_function": self._trace.array(decisions),
                "predict_proba": self._trace.array(probabilities),
            },
            iteration=self._iteration,
        )
        return probabilities

    def decision_function(self, x: np.ndarray) -> np.ndarray:
        return np.asarray(
            self._model.decision_function(x),
            dtype=float,  # type: ignore[attr-defined]
        )


@contextlib.contextmanager
def _instrument_prediction(trace: TraceState) -> Iterator[None]:
    original_make_svm = sampling_core.make_svm
    original_resample = sampling_core._resample_training_rows
    original_generate = (
        prediction_execution.PredictionSamplingRandomSource.generators_for_kinase
    )

    def traced_make_svm(**kwargs: Any) -> TracingModel:
        return TracingModel(original_make_svm(**kwargs), trace)

    def traced_resample(**kwargs: Any) -> tuple[np.ndarray, np.ndarray]:
        resampled_x, resampled_y = original_resample(**kwargs)
        trace.checkpoint(
            "iterative_resampling",
            {
                "values": trace.array(resampled_x),
                "labels": trace.array(resampled_y),
                "sample_identities": trace.identities(resampled_x),
                "rng_state_after": _normalise_scalar(
                    kwargs["resampling_rng"].bit_generator.state  # type: ignore[union-attr]
                ),
            },
            iteration=(
                cast(int, trace.current_model["iteration"])
                if trace.current_model is not None
                else None
            ),
        )
        return resampled_x, resampled_y

    def traced_generate(
        source: Any, *, kinase: str
    ) -> tuple[np.random.Generator, np.random.Generator]:
        trace.current_kinase = str(kinase)
        trace.current_ensemble = -1
        generators = original_generate(source, kinase=kinase)
        trace.checkpoint(
            "random_streams",
            {
                "configured_random_state": source.random_state,  # type: ignore[attr-defined]
                "negative_sampling_state": _normalise_scalar(
                    generators[0].bit_generator.state
                ),
                "resampling_state": _normalise_scalar(
                    generators[1].bit_generator.state
                ),
            },
        )
        return generators

    def traced_sampling(**kwargs: Any) -> np.ndarray:
        kinase = trace.current_kinase
        ensemble = trace.ensemble_counts.get(kinase, 0)
        trace.ensemble_counts[kinase] = ensemble + 1
        trace.current_ensemble = ensemble
        trace.checkpoint(
            "ensemble_training_input",
            {
                "train_values": trace.array(kwargs["train_values"]),
                "train_labels": trace.array(kwargs["train_labels"]),
                "training_sample_identities": trace.identities(
                    np.asarray(kwargs["train_values"])
                ),
                "test_values": trace.array(kwargs["test_values"]),
                "test_sample_identities": trace.identities(
                    np.asarray(kwargs["test_values"])
                ),
                "kernel": kwargs["kernel"],
                "n_iterations": kwargs["n_iterations"],
                "sampling_policy": _normalise_scalar(kwargs["sampling_policy"]),
                "resampling_rng_state_before": _normalise_scalar(
                    kwargs["resampling_rng"].bit_generator.state  # type: ignore[union-attr]
                ),
            },
        )
        result = sampling_core.run_adaptive_sampling_ensemble(**kwargs)
        trace.checkpoint("ensemble_result", trace.array(result))
        return result

    def traced_candidates(**kwargs: Any) -> dict[str, list[str]]:
        candidates = prediction_runner.build_candidate_substrate_list(**kwargs)
        trace.candidates = {
            str(key): [str(value) for value in values]
            for key, values in candidates.items()
        }
        scores = cast(pd.DataFrame, kwargs["scores"])
        trace.universe = scores.copy(deep=True)
        trace.checkpoint(
            "candidate_selection",
            {
                "score_matrix": frame_payload(scores),
                "candidates": trace.candidates,
                "top": kwargs["top"],
                "score_threshold": kwargs["score_threshold"],
                "inclusion": kwargs["inclusion"],
            },
        )
        return candidates

    def traced_adaptive(**kwargs: Any) -> pd.DataFrame:
        trace.execution_config = {
            "prediction_config": _normalise_scalar(kwargs["prediction_config"]),
            "random_state": kwargs["random_state"],
            "kernel": kwargs.get("kernel", "rbf"),
        }
        output = prediction_execution.run_adaptive_ensemble_prediction(**kwargs)
        trace.checkpoint("aggregated_prediction_values", frame_payload(output))
        return output

    def traced_rank(**kwargs: Any) -> pd.Series:
        ranking = prediction_runner.rank_kinases_for_prediction(**kwargs)
        trace.kinase_ranking = [
            {"kinase": str(key), "score": float(cast(Any, value)).hex()}
            for key, value in ranking.items()
        ]
        trace.checkpoint("kinase_assignment", trace.kinase_ranking)
        return ranking

    def traced_outputs(**kwargs: Any) -> tuple[pd.DataFrame, pd.DataFrame]:
        selected = kwargs["selected_kinases"]
        trace.selected_kinases = [str(value) for value in selected]
        pred_mat, substrate_list = prediction_runner.build_prediction_outputs(**kwargs)
        trace.checkpoint(
            "published_prediction_outputs",
            {
                "selected_kinases": trace.selected_kinases,
                "prediction_matrix": frame_payload(pred_mat),
                "substrate_list": _table_records(substrate_list),
            },
        )
        return pred_mat, substrate_list

    with (
        patch.object(sampling_core, "make_svm", traced_make_svm),
        patch.object(sampling_core, "_resample_training_rows", traced_resample),
        patch.object(
            prediction_execution,
            "run_adaptive_sampling_ensemble",
            traced_sampling,
        ),
        patch.object(
            prediction_execution.PredictionSamplingRandomSource,
            "generators_for_kinase",
            traced_generate,
        ),
        patch.object(
            kinase_executor, "build_candidate_substrate_list", traced_candidates
        ),
        patch.object(
            kinase_executor, "run_adaptive_ensemble_prediction", traced_adaptive
        ),
        patch.object(kinase_executor, "rank_kinases_for_prediction", traced_rank),
        patch.object(kinase_executor, "build_prediction_outputs", traced_outputs),
    ):
        yield


def _table_records(frame: pd.DataFrame | None) -> list[dict[str, object]] | None:
    if frame is None:
        return None
    return [
        {str(key): _normalise_scalar(value) for key, value in row.items()}
        for row in frame.to_dict(orient="records")
    ]


def _align_prediction_matrices(
    observed: pd.DataFrame, expected: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    observed_canonical = l6_metrics._canonicalise_kinase_columns(observed)
    expected_canonical = l6_metrics._canonicalise_kinase_columns(expected)
    common_rows = expected_canonical.index.intersection(observed_canonical.index)
    common_columns = expected_canonical.columns.intersection(observed_canonical.columns)
    expected_aligned = expected_canonical.loc[common_rows, common_columns].astype(float)
    observed_aligned = observed_canonical.loc[common_rows, common_columns].astype(float)
    return expected_aligned, observed_aligned, observed_aligned - expected_aligned


def _safe_correlation(left: np.ndarray, right: np.ndarray, method: str) -> float:
    finite = np.isfinite(left) & np.isfinite(right)
    if int(finite.sum()) < 2:
        return 0.0
    x = left[finite]
    y = right[finite]
    if np.unique(x).size < 2 or np.unique(y).size < 2:
        return 0.0
    result = stats.pearsonr(x, y) if method == "pearson" else stats.spearmanr(x, y)
    return float(result.statistic)


def matrix_metrics(observed: pd.DataFrame, expected: pd.DataFrame) -> dict[str, object]:
    expected_aligned, observed_aligned, difference = _align_prediction_matrices(
        observed, expected
    )
    expected_values = expected_aligned.to_numpy(dtype=float).reshape(-1)
    observed_values = observed_aligned.to_numpy(dtype=float).reshape(-1)
    diff_values = difference.to_numpy(dtype=float)
    per_kinase: dict[str, object] = {}
    for kinase in expected_aligned.columns:
        left = expected_aligned.loc[:, kinase].to_numpy(dtype=float)
        right = observed_aligned.loc[:, kinase].to_numpy(dtype=float)
        delta = right - left
        per_kinase[str(kinase)] = {
            "mae": float(np.nanmean(np.abs(delta))),
            "max_abs_difference": float(np.nanmax(np.abs(delta))),
            "pearson": _safe_correlation(left, right, "pearson"),
            "spearman": _safe_correlation(left, right, "spearman"),
        }
    largest: list[dict[str, object]] = []
    positions = np.dstack(
        np.unravel_index(
            np.argsort(np.abs(diff_values), axis=None)[::-1], diff_values.shape
        )
    )[0]
    for row_position, column_position in positions[:100]:
        expected_value = float(
            cast(Any, expected_aligned.iat[row_position, column_position])
        )
        observed_value = float(
            cast(Any, observed_aligned.iat[row_position, column_position])
        )
        difference_value = observed_value - expected_value
        largest.append(
            {
                "row_id": str(expected_aligned.index[row_position]),
                "kinase": str(expected_aligned.columns[column_position]),
                "expected_hex": expected_value.hex(),
                "observed_hex": observed_value.hex(),
                "difference_hex": difference_value.hex(),
                "absolute_difference": abs(difference_value),
            }
        )
    absolute = np.abs(observed_values - expected_values)
    return {
        "shape": list(expected_aligned.shape),
        "mae": float(np.nanmean(absolute)),
        "max_abs_difference": float(np.nanmax(absolute)),
        "pearson": _safe_correlation(expected_values, observed_values, "pearson"),
        "spearman": _safe_correlation(expected_values, observed_values, "spearman"),
        "per_kinase": per_kinase,
        "largest_100_absolute_differences": largest,
    }


def _rankings(frame: pd.DataFrame) -> dict[str, object]:
    rankings: dict[str, object] = {}
    for kinase in frame.columns:
        series = frame.loc[:, kinase].sort_values(
            ascending=False, kind="mergesort", na_position="last"
        )
        ordered = [str(value) for value in series.index]
        tie_groups: dict[str, list[list[str]]] = {}
        for width in TIE_WIDTHS:
            groups: list[list[str]] = []
            active: list[str] = []
            anchor: float | None = None
            for row_id, score in series.items():
                numeric = float(cast(Any, score))
                if anchor is None or abs(numeric - anchor) <= width:
                    active.append(str(row_id))
                    if anchor is None:
                        anchor = numeric
                else:
                    if len(active) > 1:
                        groups.append(active)
                    active = [str(row_id)]
                    anchor = numeric
            if len(active) > 1:
                groups.append(active)
            tie_groups[float(width).hex()] = groups
        rankings[str(kinase)] = {
            "top_1": ordered[:1],
            "top_5": ordered[:5],
            "top_10": ordered[:10],
            "tie_groups": tie_groups,
        }
    return {
        "candidate_tie_widths": [float(width).hex() for width in TIE_WIDTHS],
        "tie_width_policy": "fixed before observing hosted-runner differences",
        "by_kinase": rankings,
    }


def _threshold_crossings(
    observed: pd.DataFrame, expected: pd.DataFrame
) -> dict[str, object]:
    expected_aligned, observed_aligned, _ = _align_prediction_matrices(
        observed, expected
    )
    crossings: dict[str, object] = {}
    for threshold in VALUE_THRESHOLDS:
        changed: list[dict[str, object]] = []
        for row_id in expected_aligned.index:
            for kinase in expected_aligned.columns:
                left = float(cast(Any, expected_aligned.at[row_id, kinase]))
                right = float(cast(Any, observed_aligned.at[row_id, kinase]))
                if (left >= threshold) != (right >= threshold):
                    changed.append(
                        {
                            "row_id": str(row_id),
                            "kinase": str(kinase),
                            "expected_hex": left.hex(),
                            "observed_hex": right.hex(),
                        }
                    )
        crossings[float(threshold).hex()] = changed
    return {
        "candidate_thresholds": [float(value).hex() for value in VALUE_THRESHOLDS],
        "crossings": crossings,
    }


def _ordered_memberships(frame: pd.DataFrame) -> dict[str, list[str]]:
    kinase_column = "kinase"
    site_column = "site_id" if "site_id" in frame.columns else "substrate_site"
    return {
        str(kinase): [str(value) for value in group.loc[:, site_column]]
        for kinase, group in frame.groupby(kinase_column, sort=False)
    }


def _membership_difference(
    observed: Mapping[str, Sequence[str]], expected: Mapping[str, Sequence[str]]
) -> dict[str, object]:
    differences: dict[str, object] = {}
    for kinase in sorted(set(observed) | set(expected)):
        observed_values = list(observed.get(kinase, ()))
        expected_values = list(expected.get(kinase, ()))
        differences[kinase] = {
            "observed_order": observed_values,
            "expected_order": expected_values,
            "included_only_observed": sorted(
                set(observed_values) - set(expected_values)
            ),
            "excluded_from_observed": sorted(
                set(expected_values) - set(observed_values)
            ),
            "ordering_equal": observed_values == expected_values,
        }
    return differences


def _result_diagnostics(result: Any, trace: TraceState) -> dict[str, object]:
    prediction = result.prediction_result
    substrate_list = prediction.substrate_list
    observed_memberships = (
        _ordered_memberships(
            l6_metrics._normalize_prediction_topk_frame(substrate_list)
        )
        if substrate_list is not None
        else {}
    )
    expected_top = l6_metrics._normalize_prediction_topk_frame(
        load_l6_prediction_reference_top30()
    )
    expected_memberships = _ordered_memberships(expected_top)
    expected_candidates_frame = l6_metrics._canonicalise_kinase_column(
        load_l6_prediction_reference_candidate_substrates(), column_name="kinase"
    )
    expected_candidates = _ordered_memberships(expected_candidates_frame)
    site_metadata = result.dataset.site_metadata
    display_lookup = (
        {
            str(site_key): str(display_id)
            for site_key, display_id in site_metadata.loc[:, ["display_id"]]
            .reset_index()
            .itertuples(index=False)
        }
        if "display_id" in site_metadata.columns
        else {}
    )
    observed_candidates = {
        l6_metrics._canonical_kinase_label(kinase): [
            display_lookup.get(site, site) for site in sites
        ]
        for kinase, sites in trace.candidates.items()
    }
    return {
        "candidate_inclusion_exclusion": _membership_difference(
            observed_candidates, expected_candidates
        ),
        "kinase_assignment": {
            "ranking": trace.kinase_ranking,
            "selected": trace.selected_kinases,
            "expected_published_kinases": list(expected_memberships),
            "selected_only_observed": sorted(
                set(trace.selected_kinases) - set(expected_memberships)
            ),
            "selected_only_expected": sorted(
                set(expected_memberships) - set(trace.selected_kinases)
            ),
        },
        "selected_site_differences": _membership_difference(
            observed_memberships, expected_memberships
        ),
        "mask_and_attrition": {
            "eligibility_report": _normalise_scalar(result.eligibility_report),
            "site_attrition_summary": _normalise_scalar(result.site_attrition_summary),
            "attrition_provenance": _normalise_scalar(result.attrition_provenance),
            "dataset_site_count": int(len(result.dataset.phospho)),
        },
        "published_workflow_outputs": {
            "prediction_matrix": frame_payload(prediction.pred_mat),
            "substrate_list": _table_records(substrate_list),
        },
        "provenance_and_diagnostics": {
            "provenance": (
                provenance_to_payload(result.provenance)
                if result.provenance is not None
                else None
            ),
            "caveats": _normalise_scalar(result.caveats),
        },
    }


def capture(output: Path, lane: str, process_number: int) -> None:
    output.mkdir(parents=True, exist_ok=True)
    environment = runtime_identity()
    _write_json(output / "environment.json", environment)
    trace = TraceState()
    try:
        l6_metrics._run_l6_workflow.cache_clear()
        with _instrument_prediction(trace):
            result = l6_metrics._run_l6_workflow("stable")

        observed = l6_metrics._canonicalise_kinase_columns(
            l6_metrics._with_dataset_display_index(
                result.prediction_result.pred_mat, result
            )
        )
        expected = l6_metrics._canonicalise_kinase_columns(
            load_l6_prediction_reference_predmat()
        )
        expected_aligned, observed_aligned, difference = _align_prediction_matrices(
            observed, expected
        )
        _write_json(
            output / "expected-prediction-matrix.json", frame_payload(expected_aligned)
        )
        _write_json(
            output / "observed-prediction-matrix.json", frame_payload(observed_aligned)
        )
        _write_json(
            output / "difference-prediction-matrix.json", frame_payload(difference)
        )
        _write_json(output / "checkpoints.json", {"checkpoints": trace.checkpoints})
        _write_json(output / "array-store.json", {"arrays": trace.arrays})
        _write_json(
            output / "identity-store.json", {"identity_sets": trace.identity_sets}
        )

        scoring = result.scoring_result
        rank_weighted_fusion_scores = scoring.rank_weighted_fusion_scores
        score_fusion_weights = scoring.score_fusion_weights
        if rank_weighted_fusion_scores is None or score_fusion_weights is None:
            raise AssertionError("L6 diagnostic scoring surfaces are incomplete")
        scoring_surfaces = {
            "observed_profile": frame_payload(
                l6_metrics._canonicalise_kinase_columns(
                    l6_metrics._with_dataset_display_index(
                        scoring.profile_scores, result
                    )
                )
            ),
            "expected_profile": frame_payload(
                l6_metrics._canonicalise_kinase_columns(
                    load_l6_prediction_reference_profile_scores()
                )
            ),
            "observed_rank_weighted_fusion": frame_payload(
                l6_metrics._canonicalise_kinase_columns(
                    l6_metrics._with_dataset_display_index(
                        rank_weighted_fusion_scores, result
                    )
                )
            ),
            "expected_rank_weighted_fusion": frame_payload(
                l6_metrics._canonicalise_kinase_columns(
                    load_l6_prediction_reference_rank_weighted_fusion_scores()
                )
            ),
            "observed_fusion_weights": frame_payload(
                l6_metrics._canonicalise_kinase_index(score_fusion_weights)
            ),
            "expected_fusion_weights": frame_payload(
                l6_metrics._canonicalise_kinase_index(
                    load_l6_prediction_reference_score_fusion_weights()
                )
            ),
        }
        sensitivity = derive_l6_prediction_sensitivity_budget()
        summary = {
            "schema": SCHEMA,
            "metadata": {
                "encoding": "UTF-8",
                "line_endings": "LF",
                "float_encoding": FLOAT_ENCODING,
            },
            "lane": lane,
            "fresh_process": process_number,
            "source": environment["source"],
            "execution_config": trace.execution_config,
            "workflow_inputs": {
                "analysis_ready_phospho": frame_payload(result.dataset.phospho)
            },
            "fixed_promoted_reference": {
                "source": "tests/fixtures/rewrite_parity/r_reference_l6_prediction/predMat.csv",
                "sensitivity_budget": _normalise_scalar(sensitivity),
            },
            "matrix_metrics": matrix_metrics(observed, expected),
            "rankings": {
                "observed": _rankings(observed_aligned),
                "expected": _rankings(expected_aligned),
            },
            "threshold_crossings": _threshold_crossings(observed, expected),
            "downstream": _result_diagnostics(result, trace),
            "scoring_surfaces": scoring_surfaces,
            "checkpoint_count": len(trace.checkpoints),
            "status": "complete",
        }
        _write_json(output / "summary.json", summary)
    except BaseException as exc:
        _write_json(
            output / "error.json",
            {
                "schema": SCHEMA,
                "lane": lane,
                "fresh_process": process_number,
                "source": environment["source"],
                "status": "failed",
                "exception_type": type(exc).__name__,
                "message": str(exc),
                "traceback": traceback.format_exc().splitlines(),
                "checkpoints_captured": len(trace.checkpoints),
            },
        )
        _write_json(output / "checkpoints.json", {"checkpoints": trace.checkpoints})
        _write_json(output / "array-store.json", {"arrays": trace.arrays})
        _write_json(
            output / "identity-store.json", {"identity_sets": trace.identity_sets}
        )
        write_manifest(output)
        raise
    write_manifest(output)


def _load_json(path: Path) -> dict[str, object]:
    return cast(dict[str, object], json.loads(path.read_text(encoding="utf-8")))


def _checkpoint_differences(left: Path, right: Path) -> dict[str, object]:
    left_file = left / "checkpoints.json"
    right_file = right / "checkpoints.json"
    if not left_file.exists() or not right_file.exists():
        return {"comparable": False, "reason": "checkpoint file unavailable"}
    left_items = cast(list[dict[str, object]], _load_json(left_file)["checkpoints"])
    right_items = cast(list[dict[str, object]], _load_json(right_file)["checkpoints"])
    left_arrays_path = left / "array-store.json"
    right_arrays_path = right / "array-store.json"
    left_arrays = (
        cast(dict[str, dict[str, object]], _load_json(left_arrays_path)["arrays"])
        if left_arrays_path.exists()
        else {}
    )
    right_arrays = (
        cast(dict[str, dict[str, object]], _load_json(right_arrays_path)["arrays"])
        if right_arrays_path.exists()
        else {}
    )
    differing: list[dict[str, object]] = []
    for position, (left_item, right_item) in enumerate(
        zip(left_items, right_items, strict=False)
    ):
        if left_item["payload_sha256"] != right_item["payload_sha256"]:
            difference: dict[str, object] = {
                "position": position,
                "left": {
                    key: left_item[key]
                    for key in (
                        "kind",
                        "kinase",
                        "ensemble",
                        "iteration",
                        "payload_sha256",
                    )
                },
                "right": {
                    key: right_item[key]
                    for key in (
                        "kind",
                        "kinase",
                        "ensemble",
                        "iteration",
                        "payload_sha256",
                    )
                },
            }
            array_difference = _first_array_difference(
                left_item.get("payload"),
                right_item.get("payload"),
                left_arrays=left_arrays,
                right_arrays=right_arrays,
            )
            if array_difference is not None:
                difference["first_array_difference"] = array_difference
            differing.append(difference)
    if len(left_items) != len(right_items):
        differing.append(
            {
                "position": min(len(left_items), len(right_items)),
                "reason": "checkpoint count differs",
                "left_count": len(left_items),
                "right_count": len(right_items),
            }
        )
    return {
        "comparable": True,
        "left_count": len(left_items),
        "right_count": len(right_items),
        "first_divergence": differing[0] if differing else None,
        "difference_count": len(differing),
        "all_checkpoint_payloads_equal": not differing,
    }


def _array_values(payload: Mapping[str, object]) -> np.ndarray:
    shape = tuple(cast(list[int], payload["shape"]))
    return np.asarray(
        [float.fromhex(value) for value in cast(list[str], payload["values_hex"])],
        dtype=float,
    ).reshape(shape)


def _first_array_difference(
    left: object,
    right: object,
    *,
    left_arrays: Mapping[str, Mapping[str, object]],
    right_arrays: Mapping[str, Mapping[str, object]],
    path: tuple[str, ...] = (),
) -> dict[str, object] | None:
    if isinstance(left, dict) and isinstance(right, dict):
        if "array_ref" in left and "array_ref" in right:
            left_ref = str(left["array_ref"])
            right_ref = str(right["array_ref"])
            if left_ref == right_ref:
                return None
            result: dict[str, object] = {
                "path": list(path),
                "left_array_ref": left_ref,
                "right_array_ref": right_ref,
            }
            left_payload = left_arrays.get(left_ref)
            right_payload = right_arrays.get(right_ref)
            if left_payload is None or right_payload is None:
                result["values_available"] = False
                return result
            left_values = _array_values(left_payload)
            right_values = _array_values(right_payload)
            result["values_available"] = True
            result["left_shape"] = list(left_values.shape)
            result["right_shape"] = list(right_values.shape)
            if left_values.shape != right_values.shape:
                return result
            delta = right_values - left_values
            absolute = np.abs(delta)
            positions = np.dstack(
                np.unravel_index(np.argsort(absolute, axis=None)[::-1], absolute.shape)
            )[0]
            result["mae"] = float(np.nanmean(absolute))
            result["max_abs_difference"] = float(np.nanmax(absolute))
            result["largest_100_absolute_differences"] = [
                {
                    "index": [int(value) for value in position],
                    "left_hex": float(left_values[tuple(position)]).hex(),
                    "right_hex": float(right_values[tuple(position)]).hex(),
                    "difference_hex": float(delta[tuple(position)]).hex(),
                    "absolute_difference": float(absolute[tuple(position)]),
                }
                for position in positions[:100]
            ]
            return result
        for key in sorted(set(left) | set(right)):
            nested_result = _first_array_difference(
                left.get(key),
                right.get(key),
                left_arrays=left_arrays,
                right_arrays=right_arrays,
                path=(*path, str(key)),
            )
            if nested_result is not None:
                return nested_result
    elif isinstance(left, list) and isinstance(right, list):
        for index, (left_item, right_item) in enumerate(zip(left, right, strict=False)):
            nested_result = _first_array_difference(
                left_item,
                right_item,
                left_arrays=left_arrays,
                right_arrays=right_arrays,
                path=(*path, str(index)),
            )
            if nested_result is not None:
                return nested_result
    return None


def _prefix_checkpoint_differences(
    left_summary: Mapping[str, object], right_summary: Mapping[str, object]
) -> dict[str, object]:
    left_scoring = cast(Mapping[str, object], left_summary["scoring_surfaces"])
    right_scoring = cast(Mapping[str, object], right_summary["scoring_surfaces"])
    left_inputs = cast(Mapping[str, object], left_summary["workflow_inputs"])
    right_inputs = cast(Mapping[str, object], right_summary["workflow_inputs"])
    surfaces = [
        (
            "analysis_ready_phospho",
            left_inputs["analysis_ready_phospho"],
            right_inputs["analysis_ready_phospho"],
        ),
        (
            "profile_scores",
            left_scoring["observed_profile"],
            right_scoring["observed_profile"],
        ),
        (
            "score_fusion_weights",
            left_scoring["observed_fusion_weights"],
            right_scoring["observed_fusion_weights"],
        ),
        (
            "rank_weighted_fusion_scores",
            left_scoring["observed_rank_weighted_fusion"],
            right_scoring["observed_rank_weighted_fusion"],
        ),
    ]
    differences = [
        {
            "kind": kind,
            "left_matrix_sha256": cast(Mapping[str, object], left)["matrix_sha256"],
            "right_matrix_sha256": cast(Mapping[str, object], right)["matrix_sha256"],
        }
        for kind, left, right in surfaces
        if cast(Mapping[str, object], left)["matrix_sha256"]
        != cast(Mapping[str, object], right)["matrix_sha256"]
    ]
    return {
        "all_equal": not differences,
        "difference_count": len(differences),
        "first_divergence": differences[0] if differences else None,
    }


def _ranking_differences(
    left: Mapping[str, object], right: Mapping[str, object]
) -> dict[str, object]:
    left_by_kinase = cast(Mapping[str, Mapping[str, object]], left["by_kinase"])
    right_by_kinase = cast(Mapping[str, Mapping[str, object]], right["by_kinase"])
    changes: list[dict[str, object]] = []
    tie_changes: list[str] = []
    for kinase in sorted(set(left_by_kinase) | set(right_by_kinase)):
        left_ranking = left_by_kinase.get(kinase, {})
        right_ranking = right_by_kinase.get(kinase, {})
        for width in ("top_1", "top_5", "top_10"):
            left_order = list(cast(Sequence[str], left_ranking.get(width, ())))
            right_order = list(cast(Sequence[str], right_ranking.get(width, ())))
            if left_order != right_order:
                changes.append(
                    {
                        "kinase": kinase,
                        "width": width,
                        "left_order": left_order,
                        "right_order": right_order,
                        "membership_equal": set(left_order) == set(right_order),
                    }
                )
        if left_ranking.get("tie_groups") != right_ranking.get("tie_groups"):
            tie_changes.append(kinase)
    return {
        "candidate_tie_widths": left.get("candidate_tie_widths"),
        "top_order_or_membership_change_count": len(changes),
        "top_order_or_membership_changes": changes,
        "tie_group_changed_kinases": tie_changes,
    }


def _downstream_comparison(
    left_summary: Mapping[str, object], right_summary: Mapping[str, object]
) -> dict[str, object]:
    left_downstream = cast(Mapping[str, object], left_summary["downstream"])
    right_downstream = cast(Mapping[str, object], right_summary["downstream"])
    section_names = (
        "candidate_inclusion_exclusion",
        "kinase_assignment",
        "selected_site_differences",
        "mask_and_attrition",
        "published_workflow_outputs",
        "provenance_and_diagnostics",
    )
    return {
        "sections": {
            name: {
                "equal": left_downstream[name] == right_downstream[name],
                "left_sha256": _sha256_bytes(_json_bytes(left_downstream[name])),
                "right_sha256": _sha256_bytes(_json_bytes(right_downstream[name])),
            }
            for name in section_names
        },
        "rankings": _ranking_differences(
            cast(Mapping[str, object], left_summary["rankings"])["observed"],  # type: ignore[arg-type]
            cast(Mapping[str, object], right_summary["rankings"])["observed"],  # type: ignore[arg-type]
        ),
        "threshold_crossings_equal": (
            left_summary["threshold_crossings"] == right_summary["threshold_crossings"]
        ),
    }


def compare(left: Path, right: Path, output: Path, label: str) -> None:
    left_summary = left / "summary.json"
    right_summary = right / "summary.json"
    matrix_files = (
        "expected-prediction-matrix.json",
        "observed-prediction-matrix.json",
        "difference-prediction-matrix.json",
    )
    left_summary_payload = _load_json(left_summary) if left_summary.exists() else None
    right_summary_payload = (
        _load_json(right_summary) if right_summary.exists() else None
    )
    numerical_prefix = (
        _prefix_checkpoint_differences(
            cast(dict[str, object], left_summary_payload),
            cast(dict[str, object], right_summary_payload),
        )
        if left_summary_payload is not None and right_summary_payload is not None
        else {"all_equal": False, "reason": "summary unavailable"}
    )
    detailed_checkpoints = _checkpoint_differences(left, right)
    first_divergence: object | None = numerical_prefix.get("first_divergence")
    if first_divergence is None:
        first_divergence = detailed_checkpoints.get("first_divergence")
    matrix_comparison: dict[str, object] | None = None
    left_observed = left / "observed-prediction-matrix.json"
    right_observed = right / "observed-prediction-matrix.json"
    if left_observed.exists() and right_observed.exists():
        left_frame = frame_from_payload(_load_json(left_observed))
        right_frame = frame_from_payload(_load_json(right_observed))
        matrix_comparison = {
            "metrics": matrix_metrics(right_frame, left_frame),
            "left_rankings": _rankings(left_frame),
            "right_rankings": _rankings(right_frame),
        }
    comparison = {
        "schema": SCHEMA,
        "metadata": {"encoding": "UTF-8", "line_endings": "LF"},
        "label": label,
        "left": str(left),
        "right": str(right),
        "source": _source_identity(),
        "complete": left_summary.exists() and right_summary.exists(),
        "matrix_file_sha256": {
            name: {
                "left": _sha256_bytes((left / name).read_bytes())
                if (left / name).exists()
                else None,
                "right": _sha256_bytes((right / name).read_bytes())
                if (right / name).exists()
                else None,
                "equal": (
                    (left / name).exists()
                    and (right / name).exists()
                    and (left / name).read_bytes() == (right / name).read_bytes()
                ),
            }
            for name in matrix_files
        },
        "numerical_prefix_checkpoints": numerical_prefix,
        "detailed_checkpoints": detailed_checkpoints,
        "first_numerical_divergence": first_divergence,
        "observed_matrix_comparison": matrix_comparison,
        "downstream_comparison": (
            _downstream_comparison(
                cast(Mapping[str, object], left_summary_payload),
                cast(Mapping[str, object], right_summary_payload),
            )
            if left_summary_payload is not None and right_summary_payload is not None
            else None
        ),
        "downstream_sha256": {
            "left": (
                _sha256_bytes(
                    _json_bytes(
                        cast(dict[str, object], left_summary_payload)["downstream"]
                    )
                )
                if left_summary_payload is not None
                else None
            ),
            "right": (
                _sha256_bytes(
                    _json_bytes(
                        cast(dict[str, object], right_summary_payload)["downstream"]
                    )
                )
                if right_summary_payload is not None
                else None
            ),
        },
    }
    _write_json(output, comparison)


def write_manifest(root: Path) -> None:
    root.mkdir(parents=True, exist_ok=True)
    manifest_path = root / "manifest.json"
    files: list[dict[str, object]] = []
    for path in sorted(item for item in root.rglob("*") if item.is_file()):
        if path.resolve() == manifest_path.resolve():
            continue
        content = path.read_bytes()
        files.append(
            {
                "path": path.relative_to(root).as_posix(),
                "size_bytes": len(content),
                "sha256": _sha256_bytes(content),
            }
        )
    _write_json(
        manifest_path,
        {
            "schema": SCHEMA,
            "metadata": {"encoding": "UTF-8", "line_endings": "LF"},
            "source": _source_identity(),
            "manifest_policy": "self-excluding; manifest.json is not listed or hashed",
            "files": files,
        },
    )


def write_redacted_manifest(root: Path) -> None:
    root.mkdir(parents=True, exist_ok=True)
    manifest_path = root / "manifest.json"
    files: list[dict[str, object]] = []
    for path in sorted(item for item in root.rglob("*") if item.is_file()):
        if path.resolve() == manifest_path.resolve():
            continue
        content = path.read_bytes()
        files.append(
            {
                "path": path.relative_to(root).as_posix(),
                "size_bytes": len(content),
                "sha256": _sha256_bytes(content),
            }
        )
    _write_json(
        manifest_path,
        {
            "schema": SCHEMA,
            "metadata": {"encoding": "UTF-8", "line_endings": "LF"},
            "source": _public_source(_source_identity()),
            "manifest_policy": "self-excluding; manifest.json is not listed or hashed",
            "files": files,
        },
    )


def verify_manifest(root: Path) -> list[dict[str, object]]:
    manifest_path = root / "manifest.json"
    manifest = _load_json(manifest_path)
    if manifest.get("manifest_policy") != (
        "self-excluding; manifest.json is not listed or hashed"
    ):
        raise AssertionError(f"invalid manifest policy: {manifest_path}")
    listed = {
        cast(str, item["path"]): item
        for item in cast(list[dict[str, object]], manifest["files"])
    }
    actual = {
        path.relative_to(root).as_posix(): path
        for path in root.rglob("*")
        if path.is_file() and path.resolve() != manifest_path.resolve()
    }
    if set(listed) != set(actual):
        raise AssertionError(
            f"manifest membership mismatch: missing={sorted(set(actual) - set(listed))}, "
            f"unexpected={sorted(set(listed) - set(actual))}"
        )
    for relative, path in actual.items():
        content = path.read_bytes()
        expected = listed[relative]
        if len(content) != expected["size_bytes"]:
            raise AssertionError(f"manifest size mismatch: {relative}")
        if _sha256_bytes(content) != expected["sha256"]:
            raise AssertionError(f"manifest SHA-256 mismatch: {relative}")
    return [listed[key] for key in sorted(listed)]


def _decision_difference_counts(summary: Mapping[str, object]) -> dict[str, int]:
    downstream = cast(Mapping[str, object], summary["downstream"])
    candidates = cast(
        Mapping[str, Mapping[str, object]],
        downstream["candidate_inclusion_exclusion"],
    )
    selected_sites = cast(
        Mapping[str, Mapping[str, object]], downstream["selected_site_differences"]
    )
    assignment = cast(Mapping[str, object], downstream["kinase_assignment"])
    thresholds = cast(Mapping[str, object], summary["threshold_crossings"])
    crossings = cast(Mapping[str, list[object]], thresholds["crossings"])
    return {
        "candidate_inclusions": sum(
            len(cast(list[object], item["included_only_observed"]))
            for item in candidates.values()
        ),
        "candidate_exclusions": sum(
            len(cast(list[object], item["excluded_from_observed"]))
            for item in candidates.values()
        ),
        "selected_site_inclusions": sum(
            len(cast(list[object], item["included_only_observed"]))
            for item in selected_sites.values()
        ),
        "selected_site_exclusions": sum(
            len(cast(list[object], item["excluded_from_observed"]))
            for item in selected_sites.values()
        ),
        "selected_kinases_only_observed": len(
            cast(list[object], assignment["selected_only_observed"])
        ),
        "selected_kinases_only_expected": len(
            cast(list[object], assignment["selected_only_expected"])
        ),
        "threshold_crossings": sum(len(items) for items in crossings.values()),
    }


def _changed_memberships(value: object) -> dict[str, object]:
    memberships = cast(Mapping[str, Mapping[str, object]], value)
    changed: dict[str, object] = {}
    for kinase, item in memberships.items():
        if (
            not item["included_only_observed"]
            and not item["excluded_from_observed"]
            and item["ordering_equal"]
        ):
            continue
        observed = cast(list[str], item["observed_order"])
        expected = cast(list[str], item["expected_order"])
        first_mismatch: dict[str, object] | None = None
        for position, (observed_id, expected_id) in enumerate(
            zip(observed, expected, strict=False)
        ):
            if observed_id != expected_id:
                first_mismatch = {
                    "position": position,
                    "observed": observed_id,
                    "expected": expected_id,
                }
                break
        if first_mismatch is None and len(observed) != len(expected):
            first_mismatch = {
                "position": min(len(observed), len(expected)),
                "observed": observed[min(len(observed), len(expected)) :][:1],
                "expected": expected[min(len(observed), len(expected)) :][:1],
            }
        changed[kinase] = {
            "included_only_observed": item["included_only_observed"],
            "excluded_from_observed": item["excluded_from_observed"],
            "ordering_equal": item["ordering_equal"],
            "first_order_difference": first_mismatch,
        }
    return changed


def _promoted_reference_effects(summary: Mapping[str, object]) -> dict[str, object]:
    downstream = cast(Mapping[str, object], summary["downstream"])
    rankings = cast(Mapping[str, object], summary["rankings"])
    assignment = cast(Mapping[str, object], downstream["kinase_assignment"])
    return {
        "rankings": _ranking_differences(
            cast(Mapping[str, object], rankings["expected"]),
            cast(Mapping[str, object], rankings["observed"]),
        ),
        "candidate_set_changes": _changed_memberships(
            downstream["candidate_inclusion_exclusion"]
        ),
        "selected_site_changes": _changed_memberships(
            downstream["selected_site_differences"]
        ),
        "kinase_assignment_changes": {
            "selected_only_observed": assignment["selected_only_observed"],
            "selected_only_expected": assignment["selected_only_expected"],
        },
        "threshold_crossings": summary["threshold_crossings"],
    }


def _scientific_effect_annotation_records(
    *, group: str, label: str, fresh_process: int, summary: Mapping[str, object]
) -> list[dict[str, object]]:
    effects = _promoted_reference_effects(summary)
    rankings = cast(Mapping[str, object], effects["rankings"])
    ranking_changes = cast(list[object], rankings["top_order_or_membership_changes"])
    candidates = cast(Mapping[str, object], effects["candidate_set_changes"])
    selected = cast(Mapping[str, object], effects["selected_site_changes"])
    base: dict[str, object] = {
        "annotation_group": group,
        "label": label,
        "fresh_process": fresh_process,
    }
    records: list[dict[str, object]] = [
        {
            **base,
            "kind": "promoted_reference_scientific_effects",
            "ranking_change_count": rankings["top_order_or_membership_change_count"],
            "tie_group_changed_kinases": rankings["tie_group_changed_kinases"],
            "candidate_set_changed_kinase_count": len(candidates),
            "selected_site_changed_kinase_count": len(selected),
            "kinase_assignment_changes": effects["kinase_assignment_changes"],
            "threshold_crossings": effects["threshold_crossings"],
        }
    ]
    for offset in range(0, len(ranking_changes), 10):
        records.append(
            {
                **base,
                "kind": "promoted_reference_ranking_changes",
                "changes": ranking_changes[offset : offset + 10],
            }
        )
    for kind, values in (
        ("promoted_reference_candidate_changes", candidates),
        ("promoted_reference_selected_site_changes", selected),
    ):
        names = sorted(values)
        for offset in range(0, len(names), 6):
            chunk = names[offset : offset + 6]
            records.append(
                {
                    **base,
                    "kind": kind,
                    "kinases": {name: values[name] for name in chunk},
                }
            )
    return records


def _environment_markdown(environment: Mapping[str, object]) -> str:
    public = _public_environment(environment)
    operating_system = cast(Mapping[str, object], public["operating_system"])
    python = cast(Mapping[str, object], public["python"])
    return (
        f"- OS: `{operating_system.get('platform')}`; architecture: "
        f"`{operating_system.get('architecture')}`\n"
        f"- Python: `{python.get('version')}` ({python.get('implementation')})\n"
        f"- Dependencies: `{json.dumps(public['dependencies'], sort_keys=True)}`\n"
        f"- Threadpool/BLAS: `{json.dumps(public['threadpools'], sort_keys=True)}`\n"
        f"- CPU features: `{json.dumps(public['numpy_cpu_features'], sort_keys=True)}`"
    )


def inspection_markdown(root: Path, label: str) -> str:
    files = verify_manifest(root)
    lines = [
        f"## Wave 10 diagnostic inspection: {label}",
        "",
        "Manifest verification: **passed** (self-excluding, exact membership, "
        "sizes, and SHA-256 hashes).",
        "",
    ]
    process_summaries = sorted(root.glob("process-*/summary.json"))
    if process_summaries:
        summaries = [_load_json(path) for path in process_summaries]
        environment = _load_json(process_summaries[0].parent / "environment.json")
        lines.extend(["### Environment", "", _environment_markdown(environment)])
        for index, summary in enumerate(summaries, start=1):
            metrics = cast(Mapping[str, object], summary["matrix_metrics"])
            counts = _decision_difference_counts(summary)
            lines.extend(
                [
                    "",
                    f"### Fresh process {index}",
                    "",
                    f"- Checkpoints: `{summary['checkpoint_count']}`",
                    f"- Promoted-reference MAE: `{metrics['mae']}`",
                    f"- Promoted-reference maximum absolute difference: "
                    f"`{metrics['max_abs_difference']}`",
                    f"- Pearson: `{metrics['pearson']}`; Spearman: "
                    f"`{metrics['spearman']}`",
                    f"- Decision-difference counts: "
                    f"`{json.dumps(counts, sort_keys=True)}`",
                ]
            )
        comparison_path = root / "fresh-process-comparison.json"
        comparison = _load_json(comparison_path)
        detailed = cast(Mapping[str, object], comparison["detailed_checkpoints"])
        observed_comparison = cast(
            Mapping[str, object], comparison["observed_matrix_comparison"]
        )
        repeat_metrics = cast(Mapping[str, object], observed_comparison["metrics"])
        lines.extend(
            [
                "",
                "### Fresh-process repeatability",
                "",
                f"- First numerical divergence: "
                f"`{json.dumps(comparison['first_numerical_divergence'], sort_keys=True)}`",
                f"- All detailed checkpoint payloads equal: "
                f"`{detailed['all_checkpoint_payloads_equal']}`",
                f"- Matrix MAE: `{repeat_metrics['mae']}`; maximum absolute "
                f"difference: `{repeat_metrics['max_abs_difference']}`",
            ]
        )
    preflight_path = root / "windows-hardware-preflight.json"
    if preflight_path.exists():
        preflight = _load_json(preflight_path)
        lines.extend(
            [
                "### Windows hardware preflight",
                "",
                _environment_markdown(preflight),
                "",
                f"- SkylakeX capabilities: "
                f"`{json.dumps(preflight['capabilities'], sort_keys=True)}`",
                f"- Provenance-test disposition: "
                f"`{preflight['provenance_test_disposition']}`",
            ]
        )
    cross_paths = sorted((root / "cross-lane-comparisons").glob("*.json"))
    if cross_paths:
        lines.extend(["", "### Cross-lane comparisons", ""])
        for path in cross_paths:
            comparison = _load_json(path)
            observed_comparison = cast(
                Mapping[str, object], comparison["observed_matrix_comparison"]
            )
            metrics = cast(Mapping[str, object], observed_comparison["metrics"])
            largest = cast(
                list[Mapping[str, object]],
                metrics.get(
                    "largest_100_absolute_differences",
                    metrics.get("largest_10_absolute_differences", []),
                ),
            )
            top_five = [
                {
                    key: item[key]
                    for key in (
                        "row_id",
                        "kinase",
                        "absolute_difference",
                        "expected_hex",
                        "observed_hex",
                    )
                }
                for item in largest[:5]
            ]
            lines.extend(
                [
                    f"#### {comparison['label']}",
                    "",
                    f"- First numerical divergence: "
                    f"`{json.dumps(comparison['first_numerical_divergence'], sort_keys=True)}`",
                    f"- Matrix MAE: `{metrics['mae']}`; maximum absolute "
                    f"difference: `{metrics['max_abs_difference']}`",
                    f"- Pearson: `{metrics['pearson']}`; Spearman: "
                    f"`{metrics['spearman']}`",
                    f"- Five largest differing cells: "
                    f"`{json.dumps(top_five, sort_keys=True)}`",
                    f"- Downstream SHA-256: "
                    f"`{json.dumps(comparison['downstream_sha256'], sort_keys=True)}`",
                    f"- Downstream scientific effects: "
                    f"`{json.dumps(comparison.get('downstream_comparison'), sort_keys=True)}`",
                    "",
                ]
            )
    lines.extend(
        [
            "### Artifact files",
            "",
            "| File | Size (bytes) | SHA-256 |",
            "|---|---:|---|",
            *[
                f"| `{item['path']}` | {item['size_bytes']} | `{item['sha256']}` |"
                for item in files
            ],
        ]
    )
    return "\n".join(lines) + "\n"


def _bounded_matrix_metrics(metrics: Mapping[str, object]) -> dict[str, object]:
    largest = metrics.get(
        "largest_100_absolute_differences",
        metrics.get("largest_10_absolute_differences", []),
    )
    return {
        key: value
        for key, value in metrics.items()
        if key
        not in (
            "largest_100_absolute_differences",
            "largest_10_absolute_differences",
        )
    } | {"largest_10_absolute_differences": cast(list[object], largest)[:10]}


def _bounded_first_divergence(value: object) -> object:
    if not isinstance(value, dict):
        return value
    bounded = dict(value)
    array_difference = bounded.get("first_array_difference")
    if isinstance(array_difference, dict):
        bounded_array = dict(array_difference)
        largest = bounded_array.pop("largest_100_absolute_differences", None)
        if isinstance(largest, list):
            bounded_array["largest_10_absolute_differences"] = largest[:10]
        bounded["first_array_difference"] = bounded_array
    return bounded


def _public_cross_comparison(comparison: Mapping[str, object]) -> dict[str, object]:
    observed = cast(Mapping[str, object], comparison["observed_matrix_comparison"])
    detailed = cast(Mapping[str, object], comparison["detailed_checkpoints"])
    return {
        "schema": comparison["schema"],
        "metadata": comparison["metadata"],
        "label": comparison["label"],
        "source": _public_source(cast(Mapping[str, object], comparison["source"])),
        "complete": comparison["complete"],
        "matrix_file_sha256": comparison["matrix_file_sha256"],
        "numerical_prefix_checkpoints": comparison["numerical_prefix_checkpoints"],
        "detailed_checkpoints": {
            key: detailed.get(key)
            for key in (
                "comparable",
                "left_count",
                "right_count",
                "difference_count",
                "all_checkpoint_payloads_equal",
            )
        },
        "first_numerical_divergence": _bounded_first_divergence(
            comparison["first_numerical_divergence"]
        ),
        "observed_matrix_comparison": {
            "metrics": _bounded_matrix_metrics(
                cast(Mapping[str, object], observed["metrics"])
            )
        },
        "downstream_comparison": comparison["downstream_comparison"],
        "downstream_sha256": comparison["downstream_sha256"],
    }


def bounded_compare(left: Path, right: Path, output: Path, label: str) -> None:
    compare(left, right, output, label)
    _write_json(output, _public_cross_comparison(_load_json(output)))


def _metric_annotation_records(
    *,
    group: str,
    kind_prefix: str,
    base: Mapping[str, object],
    metrics: Mapping[str, object],
) -> list[dict[str, object]]:
    per_kinase = cast(Mapping[str, object], metrics["per_kinase"])
    names = sorted(per_kinase)
    largest = cast(
        list[object],
        metrics.get(
            "largest_100_absolute_differences",
            metrics.get("largest_10_absolute_differences", []),
        ),
    )[:10]
    core_metrics = {
        key: value
        for key, value in metrics.items()
        if key
        not in (
            "per_kinase",
            "largest_100_absolute_differences",
            "largest_10_absolute_differences",
        )
    }
    records: list[dict[str, object]] = [
        {
            "annotation_group": group,
            "kind": f"{kind_prefix}_metrics",
            **base,
            "matrix_metrics": core_metrics,
        }
    ]
    for offset in range(0, len(names), 6):
        chunk_names = names[offset : offset + 6]
        records.append(
            {
                "annotation_group": group,
                "kind": f"{kind_prefix}_per_kinase_metrics",
                **base,
                "kinases": {name: per_kinase[name] for name in chunk_names},
            }
        )
    records.append(
        {
            "annotation_group": group,
            "kind": f"{kind_prefix}_largest_cells",
            **base,
            "largest_10_absolute_differences": largest,
        }
    )
    return records


def inspection_annotations(root: Path, label: str) -> list[dict[str, object]]:
    files = verify_manifest(root)
    annotations: list[dict[str, object]] = [
        {
            "annotation_group": "identity",
            "kind": "artifact_manifest",
            "label": label,
            "verification": "passed",
            "manifest_policy": "self-excluding",
            "files": files,
        }
    ]
    process_summaries = sorted(root.glob("process-*/summary.json"))
    if process_summaries:
        environment = _load_json(process_summaries[0].parent / "environment.json")
        annotations.append(
            {
                "annotation_group": "identity",
                "kind": "environment",
                "label": label,
                **_public_environment(environment),
            }
        )
        for index, path in enumerate(process_summaries, start=1):
            summary = _load_json(path)
            metrics = cast(Mapping[str, object], summary["matrix_metrics"])
            downstream = cast(Mapping[str, object], summary["downstream"])
            group = f"process-{index}"
            base = {
                "label": label,
                "fresh_process": index,
                "checkpoint_count": summary["checkpoint_count"],
                "matrix_sha256": {
                    name.removesuffix("-prediction-matrix.json"): _load_json(
                        path.parent / name
                    )["matrix_sha256"]
                    for name in (
                        "expected-prediction-matrix.json",
                        "observed-prediction-matrix.json",
                        "difference-prediction-matrix.json",
                    )
                },
                "decision_difference_counts": _decision_difference_counts(summary),
                "mask_and_attrition_sha256": _sha256_bytes(
                    _json_bytes(downstream["mask_and_attrition"])
                ),
                "provenance_and_diagnostics_sha256": _sha256_bytes(
                    _json_bytes(downstream["provenance_and_diagnostics"])
                ),
                "published_workflow_outputs_sha256": _sha256_bytes(
                    _json_bytes(downstream["published_workflow_outputs"])
                ),
            }
            annotations.extend(
                _metric_annotation_records(
                    group=group,
                    kind_prefix="promoted_reference",
                    base=base,
                    metrics=metrics,
                )
            )
            annotations.extend(
                _scientific_effect_annotation_records(
                    group=f"{group}-effects",
                    label=label,
                    fresh_process=index,
                    summary=summary,
                )
            )
        comparison = _load_json(root / "fresh-process-comparison.json")
        observed = cast(Mapping[str, object], comparison["observed_matrix_comparison"])
        detailed = cast(Mapping[str, object], comparison["detailed_checkpoints"])
        repeat_base = {
            "label": label,
            "numerical_prefix_checkpoints": comparison["numerical_prefix_checkpoints"],
            "detailed_checkpoints": {
                key: detailed.get(key)
                for key in (
                    "comparable",
                    "left_count",
                    "right_count",
                    "difference_count",
                    "all_checkpoint_payloads_equal",
                )
            },
            "downstream_sha256": comparison["downstream_sha256"],
        }
        annotations.extend(
            _metric_annotation_records(
                group="repeatability",
                kind_prefix="fresh_process_repeatability",
                base=repeat_base,
                metrics=cast(Mapping[str, object], observed["metrics"]),
            )
        )
        annotations.append(
            {
                "annotation_group": "repeatability",
                "kind": "fresh_process_repeatability_first_divergence",
                **repeat_base,
                "first_numerical_divergence": _bounded_first_divergence(
                    comparison["first_numerical_divergence"]
                ),
            }
        )
    preflight_path = root / "windows-hardware-preflight.json"
    if preflight_path.exists():
        preflight = _load_json(preflight_path)
        annotations.append(
            {
                "annotation_group": "windows",
                "kind": "windows_hardware_preflight",
                "label": label,
                **_public_environment(preflight),
                "capabilities": preflight["capabilities"],
                "provenance_test_disposition": preflight["provenance_test_disposition"],
            }
        )
    for path in sorted((root / "cross-lane-comparisons").glob("*.json")):
        comparison = _load_json(path)
        observed = cast(Mapping[str, object], comparison["observed_matrix_comparison"])
        cross_label = str(comparison["label"])
        group = f"cross:{cross_label}"
        effect_group = f"cross-effects:{cross_label}"
        downstream_comparison = cast(
            Mapping[str, object], comparison["downstream_comparison"]
        )
        sections = cast(Mapping[str, object], downstream_comparison["sections"])
        rankings = cast(Mapping[str, object], downstream_comparison["rankings"])
        metric_base = {"label": cross_label}
        effect_base = {
            "annotation_group": effect_group,
            "label": cross_label,
        }
        annotations.append(
            {
                **effect_base,
                "kind": "cross_lane_downstream_effects",
                "numerical_prefix_checkpoints": comparison[
                    "numerical_prefix_checkpoints"
                ],
                "detailed_checkpoint_difference_count": cast(
                    Mapping[str, object], comparison["detailed_checkpoints"]
                )["difference_count"],
                "downstream_sha256": comparison["downstream_sha256"],
                "threshold_crossings_equal": downstream_comparison[
                    "threshold_crossings_equal"
                ],
                "ranking_change_count": rankings[
                    "top_order_or_membership_change_count"
                ],
                "tie_group_changed_kinases": rankings["tie_group_changed_kinases"],
            }
        )
        section_names = sorted(sections)
        for offset in range(0, len(section_names), 3):
            chunk = section_names[offset : offset + 3]
            annotations.append(
                {
                    **effect_base,
                    "kind": "cross_lane_downstream_sections",
                    "sections": {name: sections[name] for name in chunk},
                }
            )
        ranking_changes = cast(
            list[object], rankings["top_order_or_membership_changes"]
        )
        for offset in range(0, len(ranking_changes), 10):
            annotations.append(
                {
                    **effect_base,
                    "kind": "cross_lane_ranking_changes",
                    "changes": ranking_changes[offset : offset + 10],
                }
            )
        checkpoint_base = {
            "label": cross_label,
            "numerical_prefix_checkpoints": comparison["numerical_prefix_checkpoints"],
            "detailed_checkpoint_difference_count": cast(
                Mapping[str, object], comparison["detailed_checkpoints"]
            )["difference_count"],
        }
        annotations.extend(
            _metric_annotation_records(
                group=group,
                kind_prefix="cross_lane",
                base=metric_base,
                metrics=cast(Mapping[str, object], observed["metrics"]),
            )
        )
        annotations.append(
            {
                "annotation_group": group,
                "kind": "cross_lane_first_divergence",
                **checkpoint_base,
                "first_numerical_divergence": _bounded_first_divergence(
                    comparison["first_numerical_divergence"]
                ),
            }
        )
    return annotations


def emit_inspection_annotations(root: Path, label: str, groups: Sequence[str]) -> None:
    for record in inspection_annotations(root, label):
        if groups and record["annotation_group"] not in groups:
            continue
        kind = str(record["kind"])
        message = json.dumps(
            record,
            ensure_ascii=False,
            allow_nan=False,
            sort_keys=True,
            separators=(",", ":"),
        )
        if len(message) > 3900:
            raise AssertionError(
                f"annotation exceeds 3900 characters ({len(message)}): {kind}"
            )
        escaped = message.replace("%", "%25").replace("\r", "%0D").replace("\n", "%0A")
        print(f"::notice title=Wave10 {label} {kind}::{escaped}")


def windows_preflight(output: Path) -> None:
    identity = runtime_identity()
    try:
        command = (
            "Get-CimInstance Win32_Processor | "
            "Select-Object Manufacturer,Name,ProcessorId | ConvertTo-Json -Compress"
        )
        completed = subprocess.run(
            ["powershell", "-NoProfile", "-Command", command],
            check=True,
            capture_output=True,
            text=True,
            encoding="utf-8",
        )
        processors = json.loads(completed.stdout)
    except (OSError, subprocess.SubprocessError, json.JSONDecodeError) as exc:
        processors = {"error": f"{type(exc).__name__}: {exc}"}
    cpu_features = cast(
        dict[str, bool],
        identity["numpy_runtime"]["cpu_features"],  # type: ignore[index]
    )
    feature_availability = {
        name: bool(cpu_features.get(name, False)) for name in REQUIRED_SKX_FEATURES
    }
    all_required = all(feature_availability.values())
    identity["windows_processor"] = processors
    identity["capabilities"] = {
        "required": list(REQUIRED_SKX_FEATURES),
        "availability": feature_availability,
        "all_required": all_required,
        "feature_source": "NumPy runtime __cpu_features__",
    }
    identity["provenance_test_disposition"] = (
        "eligible_for_controlled_exact_byte_test"
        if all_required
        else "controlled_hardware_unavailable"
    )
    _write_json(output, identity)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    capture_parser = subparsers.add_parser("capture")
    capture_parser.add_argument("--output", type=Path, required=True)
    capture_parser.add_argument("--lane", required=True)
    capture_parser.add_argument("--process", type=int, required=True)
    compare_parser = subparsers.add_parser("compare")
    compare_parser.add_argument("--left", type=Path, required=True)
    compare_parser.add_argument("--right", type=Path, required=True)
    compare_parser.add_argument("--output", type=Path, required=True)
    compare_parser.add_argument("--label", required=True)
    bounded_compare_parser = subparsers.add_parser("bounded-compare")
    bounded_compare_parser.add_argument("--left", type=Path, required=True)
    bounded_compare_parser.add_argument("--right", type=Path, required=True)
    bounded_compare_parser.add_argument("--output", type=Path, required=True)
    bounded_compare_parser.add_argument("--label", required=True)
    manifest_parser = subparsers.add_parser("manifest")
    manifest_parser.add_argument("--root", type=Path, required=True)
    redacted_manifest_parser = subparsers.add_parser("redacted-manifest")
    redacted_manifest_parser.add_argument("--root", type=Path, required=True)
    inspect_parser = subparsers.add_parser("inspect")
    inspect_parser.add_argument("--root", type=Path, required=True)
    inspect_parser.add_argument("--label", required=True)
    annotate_parser = subparsers.add_parser("annotate")
    annotate_parser.add_argument("--root", type=Path, required=True)
    annotate_parser.add_argument("--label", required=True)
    annotate_parser.add_argument("--group", action="append", default=[])
    runtime_parser = subparsers.add_parser("runtime")
    runtime_parser.add_argument("--output", type=Path, required=True)
    preflight_parser = subparsers.add_parser("windows-preflight")
    preflight_parser.add_argument("--output", type=Path, required=True)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    if args.command == "capture":
        capture(args.output, args.lane, args.process)
    elif args.command == "compare":
        compare(args.left, args.right, args.output, args.label)
    elif args.command == "bounded-compare":
        bounded_compare(args.left, args.right, args.output, args.label)
    elif args.command == "manifest":
        write_manifest(args.root)
    elif args.command == "redacted-manifest":
        write_redacted_manifest(args.root)
    elif args.command == "inspect":
        sys.stdout.write(inspection_markdown(args.root, args.label))
    elif args.command == "annotate":
        emit_inspection_annotations(args.root, args.label, args.group)
    elif args.command == "runtime":
        _write_json(args.output, runtime_identity())
    elif args.command == "windows-preflight":
        windows_preflight(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
