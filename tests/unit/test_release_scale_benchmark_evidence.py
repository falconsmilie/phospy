from __future__ import annotations

import importlib.util
import math
import re
import sys
from pathlib import Path
from types import ModuleType

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

BENCHMARK_SCRIPT = (
    REPO_ROOT / "benchmarks" / ("measure_release_scale_builder_differential.py")
)
SCHEMA_MODULE_PATH = (
    REPO_ROOT / "benchmarks" / "schemas" / "release_scale_benchmark_report_schema.py"
)
EVIDENCE_DIRECTORY = REPO_ROOT / "benchmarks" / "evidence"
EVIDENCE_REPORT_RE = re.compile(
    r"^release-scale-builder-differential-\d{4}-\d{2}-\d{2}\.json$"
)
DIGEST_RE = re.compile(r"^[0-9a-f]{64}$")
_SCHEMA: ModuleType | None = None


def _load_schema() -> ModuleType:
    global _SCHEMA
    if _SCHEMA is not None:
        return _SCHEMA
    spec = importlib.util.spec_from_file_location(
        "release_scale_benchmark_report_schema_under_test",
        SCHEMA_MODULE_PATH,
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    try:
        spec.loader.exec_module(module)
    finally:
        sys.modules.pop(spec.name, None)
    _SCHEMA = module
    return module


BENCHMARK_SCHEMA_VERSION = _load_schema().BENCHMARK_SCHEMA_VERSION
EXPECTED_SAMPLES = _load_schema().EXPECTED_SAMPLES
EXPECTED_SITES = _load_schema().EXPECTED_SITES
REQUIRED_DEPENDENCIES = _load_schema().REQUIRED_DEPENDENCIES
REQUIRED_OUTPUT_FINGERPRINT_KEYS = _load_schema().REQUIRED_OUTPUT_FINGERPRINT_KEYS
REQUIRED_TIMING_KEYS = _load_schema().REQUIRED_TIMING_KEYS
RSS_UNAVAILABLE = _load_schema().RSS_UNAVAILABLE
validate_report_path = _load_schema().validate_report_path


def _retained_evidence_reports() -> list[Path]:
    if not EVIDENCE_DIRECTORY.exists():
        return []
    return sorted(EVIDENCE_DIRECTORY.glob("release-scale-builder-differential-*.json"))


def _retained_report_path() -> Path:
    reports = _retained_evidence_reports()
    assert reports, "expected at least one retained release-scale benchmark report"
    return reports[-1]


def _retained_report_payload() -> dict[str, object]:
    return validate_report_path(_retained_report_path())


def _load_benchmark_script() -> ModuleType:
    spec = importlib.util.spec_from_file_location(
        "release_scale_benchmark_script_under_test",
        BENCHMARK_SCRIPT,
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    try:
        spec.loader.exec_module(module)
    finally:
        sys.modules.pop(spec.name, None)
    return module


def _dummy_release_scale_result(module: ModuleType) -> object:
    digest = "a" * 64
    timings = {key: 0.01 for key in REQUIRED_TIMING_KEYS}
    metrics = {
        "sites": EXPECTED_SITES,
        "samples": EXPECTED_SAMPLES,
        "missing_fraction": 0.03,
        "original_missing_cell_count": 72_000,
        "final_missing_cell_count": 0,
        "output_rows": EXPECTED_SITES,
        "output_columns": 12,
        "tested_feature_count": EXPECTED_SITES,
        "process_rss_peak_mib": RSS_UNAVAILABLE,
        "scientific_summary_digest": digest,
        **timings,
    }
    scientific_summary = {
        "matrix_dimensions": {
            "rows": EXPECTED_SITES,
            "columns": EXPECTED_SAMPLES,
        },
        "original_missing_cell_count": 72_000,
        "expected_original_missing_cell_count": 72_000,
        "final_missing_cell_count": 0,
        "input_table_fingerprints": [
            {
                "name": "dataset.phospho",
                "rows": EXPECTED_SITES,
                "columns": EXPECTED_SAMPLES,
                "tolerance_hash_algorithm": "sha256-float-round-8dp-v1",
                "tolerance_hash_value": digest,
            }
        ],
        "output_table_fingerprints": [
            {
                "name": "dataset.phospho",
                "rows": EXPECTED_SITES,
                "columns": EXPECTED_SAMPLES,
                "tolerance_hash_algorithm": "sha256-float-round-8dp-v1",
                "tolerance_hash_value": digest,
            }
        ],
        "input_table_fingerprints_digest": digest,
        "output_table_fingerprints_digest": digest,
        "preprocessing_trace_digest": digest,
        "tested_feature_count": EXPECTED_SITES,
        "differential_result_table_fingerprint": {
            "name": "differential.C2_vs_C1.result_table",
            "rows": EXPECTED_SITES,
            "columns": 12,
            "hash_algorithm": "sha256-float-round-8dp-v1",
            "hash_value": digest,
        },
        "differential_result_table_digest": digest,
    }
    return module.ReleaseScaleBenchmarkResult(
        config=module.default_config(),
        timings=timings,
        metrics=metrics,
        scientific_summary=scientific_summary,
    )


def test_release_scale_benchmark_evidence_retained_report_exists() -> None:
    reports = _retained_evidence_reports()

    assert reports
    for report in reports:
        assert EVIDENCE_REPORT_RE.fullmatch(report.name)


def test_release_scale_benchmark_evidence_schema_version() -> None:
    payload = _retained_report_payload()

    assert payload["benchmark_schema_version"] == BENCHMARK_SCHEMA_VERSION


def test_release_scale_benchmark_evidence_default_dimensions() -> None:
    payload = _retained_report_payload()

    assert payload["dataset_dimensions"] == {
        "sites": EXPECTED_SITES,
        "samples": EXPECTED_SAMPLES,
    }
    assert payload["config"]["n_sites"] == EXPECTED_SITES
    assert payload["config"]["n_samples"] == EXPECTED_SAMPLES


def test_release_scale_benchmark_evidence_required_environment_fields() -> None:
    payload = _retained_report_payload()

    assert payload["python"]["version"]
    assert payload["python"]["executable"]
    assert payload["environment"]["python"]["version"]
    assert payload["environment"]["python"]["executable"]
    for field in ("platform", "system", "release", "version", "machine", "os_name"):
        assert payload["machine"][field]
    assert isinstance(payload["machine"]["processor"], str)


def test_release_scale_benchmark_evidence_timings_are_finite_non_negative() -> None:
    payload = _retained_report_payload()

    for timing_key in REQUIRED_TIMING_KEYS:
        value = payload["timings"][timing_key]
        assert isinstance(value, int | float)
        assert math.isfinite(float(value))
        assert float(value) >= 0.0


def test_release_scale_benchmark_evidence_peak_memory_contract() -> None:
    payload = _retained_report_payload()

    peak_rss = payload["peak_memory"]["process_rss_peak_mib"]
    assert peak_rss == RSS_UNAVAILABLE or (
        isinstance(peak_rss, int | float)
        and math.isfinite(float(peak_rss))
        and float(peak_rss) >= 0.0
    )


def test_release_scale_benchmark_evidence_dependency_versions() -> None:
    payload = _retained_report_payload()

    for dependency_name in REQUIRED_DEPENDENCIES:
        value = payload["dependencies"][dependency_name]
        assert isinstance(value, str)
        assert value
        assert value != "unavailable"


def test_release_scale_benchmark_evidence_fingerprint_formats() -> None:
    payload = _retained_report_payload()

    for fingerprint_key in REQUIRED_OUTPUT_FINGERPRINT_KEYS:
        assert DIGEST_RE.fullmatch(payload["output_fingerprints"][fingerprint_key])


def test_release_scale_benchmark_evidence_scientific_summary_digest_format() -> None:
    payload = _retained_report_payload()

    assert DIGEST_RE.fullmatch(payload["scientific_summary_digest"])
    assert DIGEST_RE.fullmatch(payload["metrics"]["scientific_summary_digest"])


def test_release_scale_benchmark_report_path_handling_writes_explicit_path(
    tmp_path: Path,
) -> None:
    module = _load_benchmark_script()
    report_path = tmp_path / "release-scale-builder-differential-test.json"
    result = _dummy_release_scale_result(module)

    written_path = module.write_report(
        result,
        report_path=report_path,
        command={
            "executable": sys.executable,
            "argv": [str(BENCHMARK_SCRIPT), "--report-path", str(report_path)],
            "full_argv": [
                sys.executable,
                str(BENCHMARK_SCRIPT),
                "--report-path",
                str(report_path),
            ],
            "command_line": (
                f"{sys.executable} {BENCHMARK_SCRIPT} --report-path {report_path}"
            ),
        },
    )

    assert written_path == report_path.resolve()
    payload = validate_report_path(written_path)
    assert payload["report_path"] == str(report_path.resolve())
    assert module.default_report_path() == REPO_ROOT / "benchmarks" / "reports" / (
        "release-scale-builder-differential.json"
    )


def test_release_scale_benchmark_evidence_ci_policy_does_not_invoke_workload() -> None:
    workflow_text = "\n".join(
        path.read_text(encoding="utf-8")
        for path in sorted((REPO_ROOT / ".github" / "workflows").glob("*.yml"))
    )

    assert "measure_release_scale_builder_differential.py" not in workflow_text
    assert "benchmark-release-scale" not in workflow_text
    assert "pytest tests/unit/test_benchmark_scripts_smoke.py" in workflow_text


def test_release_scale_benchmark_evidence_makefile_policy_is_explicit_only() -> None:
    makefile = (REPO_ROOT / "Makefile").read_text(encoding="utf-8")
    release_check_line = next(
        line for line in makefile.splitlines() if line.startswith("release-check:")
    )
    test_performance_recipe = makefile.split("test-performance:", maxsplit=1)[1].split(
        "test-release-gates:", maxsplit=1
    )[0]

    assert "benchmark-release-scale" not in release_check_line
    assert (
        "measure_release_scale_builder_differential.py" not in test_performance_recipe
    )
    assert "benchmark-release-scale:" in makefile
