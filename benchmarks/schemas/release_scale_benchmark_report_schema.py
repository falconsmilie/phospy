"""Fast schema validation for retained release-scale benchmark evidence.

The validator intentionally checks only report structure and invariant metadata.
It does not enforce runtime or memory thresholds because retained benchmark
evidence is a dated machine-specific observation.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
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
SOURCE_PROVENANCE_SCHEMA_VERSION = "release_scale_source_provenance_v1"
REQUIRED_SOURCE_DIGEST_KEYS = (
    ("source_tree", True),
    ("benchmark_script", False),
    ("pyproject_toml", False),
)

_SHA256_HEX_RE = re.compile(r"^[0-9a-f]{64}$")
_PYTHON_REQUIREMENT_SPECIFIER_RE = re.compile(
    r"^\s*(?P<operator>===|!=|==|<=|>=|<|>|~=)\s*"
    r"(?P<version>[0-9]+(?:\.[0-9]+){0,2}(?:\.\*)?)\s*$"
)


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

    source_provenance = _require_mapping(payload, "source_provenance")

    command = _require_mapping(payload, "command")
    _require_non_empty_string(command, "executable")
    _require_string_sequence(command, "argv", minimum_length=1)
    _require_string_sequence(command, "full_argv", minimum_length=2)
    _require_non_empty_string(command, "command_line")

    config = _require_mapping(payload, "config")
    workload = _require_mapping(payload, "workload")
    workload_configuration = _require_mapping(workload, "configuration")
    dataset_dimensions = _require_mapping(payload, "dataset_dimensions")
    dataset_sites = _require_int(dataset_dimensions, "sites")
    dataset_samples = _require_int(dataset_dimensions, "samples")
    config_sites = _require_int(config, "n_sites")
    config_samples = _require_int(config, "n_samples")
    workload_sites = _require_int(workload, "sites")
    workload_samples = _require_int(workload, "samples")
    workload_config_sites = _require_int(workload_configuration, "n_sites")
    workload_config_samples = _require_int(workload_configuration, "n_samples")
    _require_equal_value("dataset_dimensions.sites", dataset_sites, EXPECTED_SITES)
    _require_equal_value(
        "dataset_dimensions.samples", dataset_samples, EXPECTED_SAMPLES
    )
    for field, observed in (
        ("config.n_sites", config_sites),
        ("workload.sites", workload_sites),
        ("workload.configuration.n_sites", workload_config_sites),
    ):
        _require_equal_value(field, observed, dataset_sites)
    for field, observed in (
        ("config.n_samples", config_samples),
        ("workload.samples", workload_samples),
        ("workload.configuration.n_samples", workload_config_samples),
    ):
        _require_equal_value(field, observed, dataset_samples)
    config_seed = _require_int(config, "seed")
    workload_seed = _require_int(workload, "seed")
    workload_config_seed = _require_int(workload_configuration, "seed")
    _require_equal_value("workload.seed", workload_seed, config_seed)
    _require_equal_value(
        "workload.configuration.seed", workload_config_seed, config_seed
    )
    config_missing_fraction = _require_number(config, "missing_fraction")
    workload_missing_fraction = _require_number(workload, "missing_fraction")
    workload_config_missing_fraction = _require_number(
        workload_configuration,
        "missing_fraction",
    )
    _require_equal_value(
        "workload.missing_fraction",
        workload_missing_fraction,
        config_missing_fraction,
    )
    _require_equal_value(
        "workload.configuration.missing_fraction",
        workload_config_missing_fraction,
        config_missing_fraction,
    )

    output_dimensions = _require_mapping(payload, "output_dimensions")
    matrix_dimensions = _require_mapping(output_dimensions, "analysis_ready_matrix")
    analysis_ready_rows = _require_positive_int(matrix_dimensions, "rows")
    analysis_ready_columns = _require_positive_int(matrix_dimensions, "columns")
    _require_equal_value(
        "output_dimensions.analysis_ready_matrix.rows",
        analysis_ready_rows,
        dataset_sites,
    )
    _require_equal_value(
        "output_dimensions.analysis_ready_matrix.columns",
        analysis_ready_columns,
        dataset_samples,
    )
    result_table_dimensions = _require_mapping(
        output_dimensions,
        "differential_result_table",
    )
    result_table_rows = _require_positive_int(result_table_dimensions, "rows")
    result_table_columns = _require_positive_int(result_table_dimensions, "columns")
    output_tested_feature_count = _require_int(
        output_dimensions,
        "tested_feature_count",
    )
    payload_tested_feature_count = _require_int(payload, "tested_feature_count")

    metrics = _require_mapping(payload, "metrics")
    _require_equal_value("metrics.sites", _require_int(metrics, "sites"), dataset_sites)
    _require_equal_value(
        "metrics.samples",
        _require_int(metrics, "samples"),
        dataset_samples,
    )
    metrics_output_rows = _require_positive_int(metrics, "output_rows")
    metrics_output_columns = _require_positive_int(metrics, "output_columns")
    metrics_tested_feature_count = _require_int(metrics, "tested_feature_count")
    for field, observed in (
        ("output_dimensions.tested_feature_count", output_tested_feature_count),
        ("tested_feature_count", payload_tested_feature_count),
        ("metrics.tested_feature_count", metrics_tested_feature_count),
    ):
        _require_equal_value(field, observed, EXPECTED_SITES)

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

    python_payload = _require_mapping(payload, "python")
    _require_non_empty_string(python_payload, "version")
    _require_non_empty_string(python_payload, "executable")
    environment = _require_mapping(payload, "environment")
    environment_python = _require_mapping(environment, "python")
    _require_non_empty_string(environment_python, "version")
    _require_non_empty_string(environment_python, "executable")
    environment_python_version = _require_python_version_info(
        environment_python,
        "version_info",
    )
    machine = _require_mapping(payload, "machine")
    for field in ("platform", "system", "release", "version", "machine", "os_name"):
        _require_non_empty_string(machine, field)
    _require_string(machine, "processor")

    dependencies = _require_mapping(payload, "dependencies")
    environment_dependencies = _require_mapping(environment, "dependencies")
    for dependency_name in REQUIRED_DEPENDENCIES:
        value = _require_non_empty_string(dependencies, dependency_name)
        if value == "unavailable":
            _fail(f"dependencies.{dependency_name} must record an installed version")
        environment_value = _require_non_empty_string(
            environment_dependencies,
            dependency_name,
        )
        _require_equal_value(
            f"environment.dependencies.{dependency_name}",
            environment_value,
            value,
        )

    source_project_version, requires_python, source_python_version = (
        _validate_source_provenance(
            source_provenance,
            dependencies=dependencies,
        )
    )
    _require_equal_value(
        "environment.python.version_info",
        environment_python_version,
        source_python_version,
    )
    if not _python_version_satisfies_requirement(
        source_python_version,
        requires_python,
    ):
        _fail(
            "source_provenance.runtime.python_version_info must satisfy "
            "source_provenance.project.requires_python"
        )
    _require_equal_value(
        "dependencies.phospy",
        _require_non_empty_string(dependencies, "phospy"),
        source_project_version,
    )

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
    summary_matrix_rows = _require_positive_int(summary_matrix_dimensions, "rows")
    summary_matrix_columns = _require_positive_int(summary_matrix_dimensions, "columns")
    _require_equal_value(
        "scientific_summary.matrix_dimensions.rows",
        summary_matrix_rows,
        analysis_ready_rows,
    )
    _require_equal_value(
        "scientific_summary.matrix_dimensions.columns",
        summary_matrix_columns,
        analysis_ready_columns,
    )
    summary_tested_feature_count = _require_int(
        scientific_summary,
        "tested_feature_count",
    )
    _require_equal_value(
        "scientific_summary.tested_feature_count",
        summary_tested_feature_count,
        payload_tested_feature_count,
    )
    for digest_key in (
        "input_table_fingerprints_digest",
        "output_table_fingerprints_digest",
        "preprocessing_trace_digest",
        "differential_result_table_digest",
        "differential_policy_provenance_digest",
        "differential_workflow_provenance_digest",
        "dataset_workflow_provenance_digest",
    ):
        _require_digest(scientific_summary, digest_key)
    result_table_fingerprint = _require_mapping(
        scientific_summary,
        "differential_result_table_fingerprint",
    )
    fingerprint_rows = _require_positive_int(result_table_fingerprint, "rows")
    fingerprint_columns = _require_positive_int(result_table_fingerprint, "columns")
    for field, observed in (
        ("metrics.output_rows", metrics_output_rows),
        ("output_dimensions.differential_result_table.rows", result_table_rows),
        (
            "scientific_summary.differential_result_table_fingerprint.rows",
            fingerprint_rows,
        ),
    ):
        _require_equal_value(field, observed, EXPECTED_SITES)
    _require_equal_value(
        "metrics.output_columns",
        metrics_output_columns,
        result_table_columns,
    )
    _require_equal_value(
        "scientific_summary.differential_result_table_fingerprint.columns",
        fingerprint_columns,
        result_table_columns,
    )
    _require_non_empty_string(result_table_fingerprint, "hash_algorithm")
    result_table_hash = _require_digest(result_table_fingerprint, "hash_value")
    _require_equal_value(
        "scientific_summary.differential_result_table_digest",
        _require_digest(scientific_summary, "differential_result_table_digest"),
        result_table_hash,
    )
    _require_equal_value(
        "output_fingerprints.differential_result_table_digest",
        _require_digest(output_fingerprints, "differential_result_table_digest"),
        result_table_hash,
    )
    for digest_key in (
        "input_table_fingerprints_digest",
        "output_table_fingerprints_digest",
        "preprocessing_trace_digest",
    ):
        _require_equal_value(
            f"output_fingerprints.{digest_key}",
            _require_digest(output_fingerprints, digest_key),
            _require_digest(scientific_summary, digest_key),
        )
    scientific_summary_digest = _require_digest(payload, "scientific_summary_digest")
    _require_equal_value(
        "metrics.scientific_summary_digest",
        _require_digest(metrics, "scientific_summary_digest"),
        scientific_summary_digest,
    )
    _require_equal_value(
        "output_fingerprints.scientific_summary_digest",
        _require_digest(output_fingerprints, "scientific_summary_digest"),
        scientific_summary_digest,
    )
    _validate_table_fingerprint_list(
        scientific_summary,
        "input_table_fingerprints",
    )
    _validate_table_fingerprint_list(
        scientific_summary,
        "output_table_fingerprints",
    )


def _validate_source_provenance(
    payload: Mapping[str, Any],
    *,
    dependencies: Mapping[str, Any],
) -> tuple[str, str, tuple[int, int, int]]:
    _require_equal(payload, "schema", SOURCE_PROVENANCE_SCHEMA_VERSION)
    project = _require_mapping(payload, "project")
    imported = _require_mapping(payload, "imported_phospy")
    runtime = _require_mapping(payload, "runtime")
    digests = _require_mapping(payload, "digests")

    _require_equal(project, "name", "phospy")
    project_version = _require_non_empty_string(project, "version")
    requires_python = _require_non_empty_string(project, "requires_python")
    _require_non_empty_string(project, "pyproject_path")
    project_pyproject_digest = _require_digest(project, "pyproject_sha256")

    _require_non_empty_string(imported, "module_path")
    _require_non_empty_string(imported, "resolved_path")
    _require_non_empty_string(imported, "expected_package_root")

    runtime_phospy_version = _require_non_empty_string(runtime, "phospy_version")
    distribution_phospy_version = _require_non_empty_string(
        runtime,
        "distribution_phospy_version",
    )
    _require_equal_value(
        "source_provenance.runtime.phospy_version",
        runtime_phospy_version,
        project_version,
    )
    _require_equal_value(
        "source_provenance.runtime.distribution_phospy_version",
        distribution_phospy_version,
        project_version,
    )
    _require_equal_value(
        "dependencies.phospy",
        _require_non_empty_string(dependencies, "phospy"),
        project_version,
    )
    _require_non_empty_string(runtime, "python_version")
    _require_non_empty_string(runtime, "python_executable")
    source_python_version = _require_python_version_info(
        runtime,
        "python_version_info",
    )

    for section_name, require_file_count in REQUIRED_SOURCE_DIGEST_KEYS:
        section = _require_mapping(digests, section_name)
        _require_non_empty_string(section, "path")
        _require_non_empty_string(section, "algorithm")
        section_digest = _require_digest(section, "sha256")
        if section_name == "pyproject_toml":
            _require_equal_value(
                "source_provenance.digests.pyproject_toml.sha256",
                section_digest,
                project_pyproject_digest,
            )
        if require_file_count:
            _require_positive_int(section, "file_count")

    git = _require_mapping(payload, "git")
    git_available = _require_bool(git, "available")
    if git_available:
        _require_non_empty_string(git, "commit")
        _require_non_empty_string(git, "tree")
        _require_bool(git, "dirty")
        _require_string_sequence(git, "status_porcelain_v1", minimum_length=0)
        _require_digest(git, "status_porcelain_v1_sha256")

    return project_version, requires_python, source_python_version


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


def _python_version_satisfies_requirement(
    version_info: Sequence[int],
    requires_python: str,
) -> bool:
    specifiers = tuple(
        item.strip() for item in requires_python.split(",") if item.strip()
    )
    if not specifiers:
        _fail("source_provenance.project.requires_python must not be empty")
    version = (int(version_info[0]), int(version_info[1]), int(version_info[2]))
    for specifier in specifiers:
        match = _PYTHON_REQUIREMENT_SPECIFIER_RE.fullmatch(specifier)
        if match is None:
            _fail(
                "source_provenance.project.requires_python contains an unsupported "
                f"specifier: {specifier!r}"
            )
        operator = match.group("operator")
        expected_raw = match.group("version")
        if not _version_satisfies_specifier(version, operator, expected_raw):
            return False
    return True


def _version_satisfies_specifier(
    version: tuple[int, int, int],
    operator: str,
    expected_raw: str,
) -> bool:
    if expected_raw.endswith(".*"):
        expected_prefix = tuple(
            int(part) for part in expected_raw.removesuffix(".*").split(".")
        )
        prefix_matches = version[: len(expected_prefix)] == expected_prefix
        if operator == "==":
            return prefix_matches
        if operator == "!=":
            return not prefix_matches
        _fail("wildcard requires-python specifiers only support == and !=")
    expected = _version_tuple(expected_raw)
    if operator == "==":
        return _compare_versions(version, expected) == 0
    if operator == "!=":
        return _compare_versions(version, expected) != 0
    if operator == ">=":
        return _compare_versions(version, expected) >= 0
    if operator == ">":
        return _compare_versions(version, expected) > 0
    if operator == "<=":
        return _compare_versions(version, expected) <= 0
    if operator == "<":
        return _compare_versions(version, expected) < 0
    if operator == "~=":
        upper_bound = _compatible_release_upper_bound(expected_raw)
        return (
            _compare_versions(version, expected) >= 0
            and _compare_versions(version, upper_bound) < 0
        )
    if operator == "===":
        return ".".join(str(part) for part in version) == expected_raw
    _fail(f"unsupported requires-python operator: {operator!r}")


def _version_tuple(value: str) -> tuple[int, int, int]:
    parts = tuple(int(part) for part in value.split("."))
    if len(parts) > 3:
        _fail(f"unsupported Python version specifier: {value!r}")
    normalized = (*parts, *(0 for _ in range(3 - len(parts))))
    return (normalized[0], normalized[1], normalized[2])


def _compare_versions(left: Sequence[int], right: Sequence[int]) -> int:
    left_tuple = tuple(int(part) for part in left[:3])
    right_tuple = tuple(int(part) for part in right[:3])
    if left_tuple < right_tuple:
        return -1
    if left_tuple > right_tuple:
        return 1
    return 0


def _compatible_release_upper_bound(value: str) -> tuple[int, int, int]:
    parts = [int(part) for part in value.split(".")]
    if len(parts) == 1:
        return (parts[0] + 1, 0, 0)
    if len(parts) == 2:
        return (parts[0] + 1, 0, 0)
    return (parts[0], parts[1] + 1, 0)


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


def _require_bool(payload: Mapping[str, Any], field: str) -> bool:
    value = _require_present(payload, field)
    if not isinstance(value, bool):
        _fail(f"{field} must be a boolean")
    return value


def _require_equal(payload: Mapping[str, Any], field: str, expected: object) -> None:
    value = _require_present(payload, field)
    if value != expected:
        _fail(f"{field} must be {expected!r}; observed {value!r}")


def _require_equal_value(field: str, observed: object, expected: object) -> None:
    if observed != expected:
        _fail(f"{field} must equal {expected!r}; observed {observed!r}")


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


def _require_python_version_info(
    payload: Mapping[str, Any],
    field: str,
) -> tuple[int, int, int]:
    value = _require_sequence(payload, field, minimum_length=3)
    parts = value[:3]
    if not all(isinstance(part, int) and not isinstance(part, bool) for part in parts):
        _fail(f"{field} must start with integer major/minor/patch values")
    return int(parts[0]), int(parts[1]), int(parts[2])


def _fail(message: str) -> NoReturn:
    raise BenchmarkReportSchemaError(message)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Validate retained release-scale benchmark evidence JSON.",
    )
    parser.add_argument("report_path", nargs="+", type=Path)
    args = parser.parse_args(argv)
    for report_path in args.report_path:
        validate_report_path(report_path)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except BenchmarkReportSchemaError as error:
        print(f"release-scale benchmark report schema error: {error}", file=sys.stderr)
        raise SystemExit(1) from error
