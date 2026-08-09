"""Fast schema validation for retained release-scale benchmark evidence.

The validator intentionally checks only report structure and invariant metadata.
It does not enforce runtime or memory thresholds because retained benchmark
evidence is a dated machine-specific observation.
"""

from __future__ import annotations

import json
import math
import re
from collections.abc import Mapping, Sequence
from datetime import date, datetime
from pathlib import Path
from typing import Any, NoReturn

EVIDENCE_CATEGORY = "benchmark_evidence"
BENCHMARK_SCHEMA_VERSION = "release_scale_builder_differential_report_v1"
BENCHMARK_NAME = "release_scale_builder_differential"
EXPECTED_SITES = 50_000
EXPECTED_SAMPLES = 48
RSS_UNAVAILABLE = "unavailable"
REQUIRED_DEPENDENCIES = ("numpy", "pandas", "scipy", "phospy")
REQUIRED_TIMING_KEYS = (
    "request_preparation_seconds",
    "builder_execution_seconds",
    "preprocessing_execution_seconds",
    "preprocessing_report_assembly_seconds",
    "provenance_fingerprinting_seconds",
    "differential_analysis_seconds",
    "result_table_assembly_seconds",
    "result_table_fingerprinting_seconds",
    "scientific_summary_assembly_seconds",
    "single_workload_runtime_seconds",
    "total_runtime_seconds",
)
REQUIRED_OUTPUT_FINGERPRINT_KEYS = (
    "input_table_fingerprints_digest",
    "output_table_fingerprints_digest",
    "preprocessing_trace_digest",
    "differential_result_table_digest",
    "scientific_summary_digest",
)

_SHA256_HEX_RE = re.compile(r"^[0-9a-f]{64}$")


class BenchmarkReportSchemaError(AssertionError):
    """Raised when retained release-scale benchmark evidence is malformed."""


def validate_report_path(path: str | Path) -> dict[str, Any]:
    """Load and validate one retained release-scale benchmark report."""

    report_path = Path(path)
    payload = json.loads(report_path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        _fail("report root must be a JSON object")
    validate_release_scale_benchmark_report(payload)
    return payload


def validate_release_scale_benchmark_report(payload: Mapping[str, Any]) -> None:
    """Validate retained benchmark-evidence report structure and invariants."""

    _require_equal(payload, "evidence_category", EVIDENCE_CATEGORY)
    _require_equal(payload, "benchmark_schema_version", BENCHMARK_SCHEMA_VERSION)
    _require_equal(payload, "benchmark", BENCHMARK_NAME)
    _require_equal(payload, "benchmark_name", BENCHMARK_NAME)
    _require_non_empty_string(payload, "description")
    generated_at = _require_non_empty_string(payload, "generated_at_utc")
    generation_date = _require_non_empty_string(payload, "generation_date_utc")
    _validate_utc_timestamp(generated_at, field="generated_at_utc")
    _validate_generation_date(
        generation_date,
        generated_at_utc=generated_at,
        field="generation_date_utc",
    )

    command = _require_mapping(payload, "command")
    _require_non_empty_string(command, "executable")
    _require_string_sequence(command, "argv", minimum_length=1)
    _require_string_sequence(command, "full_argv", minimum_length=2)
    _require_non_empty_string(command, "command_line")

    config = _require_mapping(payload, "config")
    workload = _require_mapping(payload, "workload")
    dataset_dimensions = _require_mapping(payload, "dataset_dimensions")
    _require_equal(dataset_dimensions, "sites", EXPECTED_SITES)
    _require_equal(dataset_dimensions, "samples", EXPECTED_SAMPLES)
    _require_equal(config, "n_sites", EXPECTED_SITES)
    _require_equal(config, "n_samples", EXPECTED_SAMPLES)
    _require_equal(workload, "sites", EXPECTED_SITES)
    _require_equal(workload, "samples", EXPECTED_SAMPLES)
    _require_int(config, "seed")
    _require_int(workload, "seed")
    _require_number(config, "missing_fraction")
    _require_number(workload, "missing_fraction")

    output_dimensions = _require_mapping(payload, "output_dimensions")
    matrix_dimensions = _require_mapping(output_dimensions, "analysis_ready_matrix")
    _require_equal(matrix_dimensions, "rows", EXPECTED_SITES)
    _require_equal(matrix_dimensions, "columns", EXPECTED_SAMPLES)
    result_table_dimensions = _require_mapping(
        output_dimensions,
        "differential_result_table",
    )
    _require_equal(result_table_dimensions, "rows", EXPECTED_SITES)
    _require_positive_int(result_table_dimensions, "columns")
    _require_equal(output_dimensions, "tested_feature_count", EXPECTED_SITES)
    _require_equal(payload, "tested_feature_count", EXPECTED_SITES)

    metrics = _require_mapping(payload, "metrics")
    _require_equal(metrics, "sites", EXPECTED_SITES)
    _require_equal(metrics, "samples", EXPECTED_SAMPLES)
    _require_equal(metrics, "output_rows", EXPECTED_SITES)
    _require_positive_int(metrics, "output_columns")
    _require_equal(metrics, "tested_feature_count", EXPECTED_SITES)

    timings = _require_mapping(payload, "timings")
    runtime = _require_mapping(payload, "runtime")
    runtime_timings = _require_mapping(runtime, "timings_seconds")
    for timing_key in REQUIRED_TIMING_KEYS:
        _require_finite_non_negative_number(timings, timing_key)
        _require_finite_non_negative_number(runtime_timings, timing_key)
    _require_finite_non_negative_number(runtime, "total_runtime_seconds")
    _require_finite_non_negative_number(metrics, "total_runtime_seconds")

    peak_memory = _require_mapping(payload, "peak_memory")
    _validate_peak_memory_value(metrics, "process_rss_peak_mib")
    _validate_peak_memory_value(peak_memory, "process_rss_peak_mib")

    python = _require_mapping(payload, "python")
    _require_non_empty_string(python, "version")
    _require_non_empty_string(python, "executable")
    environment = _require_mapping(payload, "environment")
    environment_python = _require_mapping(environment, "python")
    _require_non_empty_string(environment_python, "version")
    _require_non_empty_string(environment_python, "executable")
    machine = _require_mapping(payload, "machine")
    for field in ("platform", "system", "release", "version", "machine", "os_name"):
        _require_non_empty_string(machine, field)
    _require_string(machine, "processor")

    dependencies = _require_mapping(payload, "dependencies")
    for dependency_name in REQUIRED_DEPENDENCIES:
        value = _require_non_empty_string(dependencies, dependency_name)
        if value == "unavailable":
            _fail(f"dependencies.{dependency_name} must record an installed version")

    output_fingerprints = _require_mapping(payload, "output_fingerprints")
    for fingerprint_key in REQUIRED_OUTPUT_FINGERPRINT_KEYS:
        _require_digest(output_fingerprints, fingerprint_key)
    _require_digest(payload, "scientific_summary_digest")
    _require_digest(metrics, "scientific_summary_digest")

    scientific_summary = _require_mapping(payload, "scientific_summary")
    summary_matrix_dimensions = _require_mapping(
        scientific_summary,
        "matrix_dimensions",
    )
    _require_equal(summary_matrix_dimensions, "rows", EXPECTED_SITES)
    _require_equal(summary_matrix_dimensions, "columns", EXPECTED_SAMPLES)
    _require_equal(scientific_summary, "tested_feature_count", EXPECTED_SITES)
    for digest_key in (
        "input_table_fingerprints_digest",
        "output_table_fingerprints_digest",
        "preprocessing_trace_digest",
        "differential_result_table_digest",
    ):
        _require_digest(scientific_summary, digest_key)
    result_table_fingerprint = _require_mapping(
        scientific_summary,
        "differential_result_table_fingerprint",
    )
    _require_equal(result_table_fingerprint, "rows", EXPECTED_SITES)
    _require_positive_int(result_table_fingerprint, "columns")
    _require_non_empty_string(result_table_fingerprint, "hash_algorithm")
    _require_digest(result_table_fingerprint, "hash_value")
    _validate_table_fingerprint_list(
        scientific_summary,
        "input_table_fingerprints",
    )
    _validate_table_fingerprint_list(
        scientific_summary,
        "output_table_fingerprints",
    )


def _validate_table_fingerprint_list(
    payload: Mapping[str, Any],
    field: str,
) -> None:
    fingerprints = _require_sequence(payload, field, minimum_length=1)
    for index, fingerprint in enumerate(fingerprints):
        if not isinstance(fingerprint, Mapping):
            _fail(f"{field}[{index}] must be an object")
        _require_non_empty_string(fingerprint, "name")
        _require_positive_int(fingerprint, "rows")
        _require_positive_int(fingerprint, "columns")
        _require_non_empty_string(fingerprint, "tolerance_hash_algorithm")
        _require_digest(fingerprint, "tolerance_hash_value")


def _validate_peak_memory_value(payload: Mapping[str, Any], field: str) -> None:
    value = _require_present(payload, field)
    if value == RSS_UNAVAILABLE:
        return
    if isinstance(value, bool) or not isinstance(value, int | float):
        _fail(f"{field} must be a finite non-negative number or {RSS_UNAVAILABLE!r}")
    if not math.isfinite(float(value)) or float(value) < 0.0:
        _fail(f"{field} must be finite and non-negative")


def _validate_utc_timestamp(value: str, *, field: str) -> None:
    if not value.endswith("Z"):
        _fail(f"{field} must be a UTC ISO-8601 timestamp ending with 'Z'")
    try:
        datetime.fromisoformat(value.removesuffix("Z") + "+00:00")
    except ValueError:
        _fail(f"{field} must be a valid UTC ISO-8601 timestamp")


def _validate_generation_date(
    value: str,
    *,
    generated_at_utc: str,
    field: str,
) -> None:
    try:
        date.fromisoformat(value)
    except ValueError:
        _fail(f"{field} must be a valid ISO date")
    if value != generated_at_utc.split("T", maxsplit=1)[0]:
        _fail(f"{field} must match the date portion of generated_at_utc")


def _require_present(payload: Mapping[str, Any], field: str) -> Any:
    if field not in payload:
        _fail(f"missing required field: {field}")
    return payload[field]


def _require_mapping(payload: Mapping[str, Any], field: str) -> Mapping[str, Any]:
    value = _require_present(payload, field)
    if not isinstance(value, Mapping):
        _fail(f"{field} must be an object")
    return value


def _require_sequence(
    payload: Mapping[str, Any],
    field: str,
    *,
    minimum_length: int,
) -> Sequence[Any]:
    value = _require_present(payload, field)
    if isinstance(value, str) or not isinstance(value, Sequence):
        _fail(f"{field} must be an array")
    if len(value) < minimum_length:
        _fail(f"{field} must contain at least {minimum_length} item(s)")
    return value


def _require_string_sequence(
    payload: Mapping[str, Any],
    field: str,
    *,
    minimum_length: int,
) -> tuple[str, ...]:
    value = _require_sequence(payload, field, minimum_length=minimum_length)
    strings: list[str] = []
    for index, item in enumerate(value):
        if not isinstance(item, str) or item == "":
            _fail(f"{field}[{index}] must be a non-empty string")
        strings.append(item)
    return tuple(strings)


def _require_string(payload: Mapping[str, Any], field: str) -> str:
    value = _require_present(payload, field)
    if not isinstance(value, str):
        _fail(f"{field} must be a string")
    return value


def _require_non_empty_string(payload: Mapping[str, Any], field: str) -> str:
    value = _require_string(payload, field)
    if value == "":
        _fail(f"{field} must be a non-empty string")
    return value


def _require_equal(payload: Mapping[str, Any], field: str, expected: object) -> None:
    value = _require_present(payload, field)
    if value != expected:
        _fail(f"{field} must be {expected!r}; observed {value!r}")


def _require_int(payload: Mapping[str, Any], field: str) -> int:
    value = _require_present(payload, field)
    if isinstance(value, bool) or not isinstance(value, int):
        _fail(f"{field} must be an integer")
    return int(value)


def _require_positive_int(payload: Mapping[str, Any], field: str) -> int:
    value = _require_int(payload, field)
    if value <= 0:
        _fail(f"{field} must be positive")
    return value


def _require_number(payload: Mapping[str, Any], field: str) -> float:
    value = _require_present(payload, field)
    if isinstance(value, bool) or not isinstance(value, int | float):
        _fail(f"{field} must be a number")
    numeric = float(value)
    if not math.isfinite(numeric):
        _fail(f"{field} must be finite")
    return numeric


def _require_finite_non_negative_number(
    payload: Mapping[str, Any],
    field: str,
) -> float:
    numeric = _require_number(payload, field)
    if numeric < 0.0:
        _fail(f"{field} must be non-negative")
    return numeric


def _require_digest(payload: Mapping[str, Any], field: str) -> str:
    value = _require_non_empty_string(payload, field)
    if _SHA256_HEX_RE.fullmatch(value) is None:
        _fail(f"{field} must be a 64-character lowercase SHA-256 hex digest")
    return value


def _fail(message: str) -> NoReturn:
    raise BenchmarkReportSchemaError(message)
