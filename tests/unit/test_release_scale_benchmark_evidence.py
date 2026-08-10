from __future__ import annotations

import importlib.util
import json
import math
import re
import sys
import tomllib
from pathlib import Path
from types import ModuleType

import pandas as pd
import pytest

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


def _pyproject() -> dict[str, object]:
    with (REPO_ROOT / "pyproject.toml").open("rb") as handle:
        return tomllib.load(handle)


def _project_version() -> str:
    project = _pyproject()["project"]
    assert isinstance(project, dict)
    version = project["version"]
    assert isinstance(version, str)
    return version


def _project_requires_python() -> str:
    project = _pyproject()["project"]
    assert isinstance(project, dict)
    requires_python = project["requires-python"]
    assert isinstance(requires_python, str)
    return requires_python


def _current_supported_python_version_info() -> list[object]:
    major, minor, patch = sys.version_info[:3]
    if not ((major, minor) >= (3, 11) and (major, minor) < (3, 13)):
        major, minor, patch = 3, 12, 10
    return [major, minor, patch, "final", 0]


def _current_supported_python_version_text() -> str:
    version_info = _current_supported_python_version_info()
    return ".".join(str(part) for part in version_info[:3])


def _dummy_source_provenance() -> dict[str, object]:
    project_version = _project_version()
    python_version_info = _current_supported_python_version_info()
    return {
        "schema": "release_scale_source_provenance_v1",
        "project": {
            "name": "phospy",
            "version": project_version,
            "requires_python": _project_requires_python(),
            "pyproject_path": "pyproject.toml",
            "pyproject_sha256": "b" * 64,
        },
        "imported_phospy": {
            "module_path": "src/phospy/__init__.py",
            "resolved_path": str((SRC / "phospy" / "__init__.py").resolve()),
            "expected_package_root": "src/phospy",
        },
        "runtime": {
            "python_version": ".".join(str(part) for part in python_version_info[:3]),
            "python_version_info": python_version_info,
            "python_executable": sys.executable,
            "phospy_version": project_version,
            "distribution_phospy_version": project_version,
        },
        "digests": {
            "source_tree": {
                "path": "src/phospy",
                "sha256": "c" * 64,
                "algorithm": "sha256-sorted-relative-paths-and-file-bytes-v1",
                "file_count": 1,
            },
            "benchmark_script": {
                "path": "benchmarks/measure_release_scale_builder_differential.py",
                "sha256": "d" * 64,
                "algorithm": "sha256-file-bytes-v1",
            },
            "pyproject_toml": {
                "path": "pyproject.toml",
                "sha256": "b" * 64,
                "algorithm": "sha256-file-bytes-v1",
            },
        },
        "git": {
            "available": True,
            "commit": "1" * 40,
            "tree": "2" * 40,
            "dirty": False,
            "status_porcelain_v1": [],
            "status_porcelain_v1_sha256": "e" * 64,
        },
    }


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
        "differential_policy_provenance_digest": digest,
        "differential_workflow_provenance_digest": digest,
        "dataset_workflow_provenance_digest": digest,
    }
    return module.ReleaseScaleBenchmarkResult(
        config=module.default_config(),
        timings=timings,
        metrics=metrics,
        scientific_summary=scientific_summary,
        source_provenance=_dummy_source_provenance(),
    )


def _valid_schema_payload() -> dict[str, object]:
    module = _load_benchmark_script()
    payload = module.report_payload(
        _dummy_release_scale_result(module),
        generated_at_utc="2026-08-10T00:00:00Z",
        command={
            "executable": sys.executable,
            "argv": [str(BENCHMARK_SCRIPT)],
            "full_argv": [sys.executable, str(BENCHMARK_SCRIPT)],
            "command_line": f"{sys.executable} {BENCHMARK_SCRIPT}",
        },
        source_provenance=_dummy_source_provenance(),
    )
    version_text = _current_supported_python_version_text()
    version_info = _current_supported_python_version_info()
    project_version = _project_version()
    dependencies = {
        "numpy": "2.5.1",
        "pandas": "3.0.3",
        "phospy": project_version,
        "scipy": "1.18.0",
    }
    payload["python"] = {
        "version": f"{version_text} test",
        "executable": sys.executable,
    }
    environment = payload["environment"]
    assert isinstance(environment, dict)
    environment["python"] = {
        "version": f"{version_text} test",
        "version_info": version_info,
        "implementation": "CPython",
        "executable": sys.executable,
    }
    environment["dependencies"] = dict(dependencies)
    payload["dependencies"] = dict(dependencies)
    return payload


def _assert_invalid_report(payload: dict[str, object], expected: str) -> None:
    schema = _load_schema()
    with pytest.raises(schema.BenchmarkReportSchemaError, match=expected):
        schema.validate_release_scale_benchmark_report(payload)


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


def test_release_scale_benchmark_evidence_records_current_project_version() -> None:
    payload = _retained_report_payload()
    project_version = _project_version()

    assert payload["dependencies"]["phospy"] == project_version
    assert payload["source_provenance"]["project"]["version"] == project_version
    assert (
        payload["source_provenance"]["runtime"]["distribution_phospy_version"]
        == project_version
    )


def test_release_scale_benchmark_evidence_records_supported_python_version() -> None:
    payload = _retained_report_payload()
    schema = _load_schema()

    assert schema._python_version_satisfies_requirement(
        tuple(payload["source_provenance"]["runtime"]["python_version_info"][:3]),
        payload["source_provenance"]["project"]["requires_python"],
    )


def test_release_scale_benchmark_evidence_records_source_provenance() -> None:
    payload = _retained_report_payload()
    source_provenance = payload["source_provenance"]

    assert source_provenance["imported_phospy"]["module_path"]
    for section_name in ("source_tree", "benchmark_script", "pyproject_toml"):
        assert DIGEST_RE.fullmatch(source_provenance["digests"][section_name]["sha256"])


def test_release_scale_benchmark_evidence_fingerprint_formats() -> None:
    payload = _retained_report_payload()

    for fingerprint_key in REQUIRED_OUTPUT_FINGERPRINT_KEYS:
        assert DIGEST_RE.fullmatch(payload["output_fingerprints"][fingerprint_key])


def test_release_scale_benchmark_evidence_scientific_summary_digest_format() -> None:
    payload = _retained_report_payload()

    assert DIGEST_RE.fullmatch(payload["scientific_summary_digest"])
    assert DIGEST_RE.fullmatch(payload["metrics"]["scientific_summary_digest"])


def test_release_scale_schema_rejects_differential_columns_fingerprint_mismatch() -> (
    None
):
    payload = _valid_schema_payload()
    payload["metrics"]["output_columns"] = 48
    payload["output_dimensions"]["differential_result_table"]["columns"] = 48
    payload["scientific_summary"]["differential_result_table_fingerprint"][
        "columns"
    ] = 28

    _assert_invalid_report(payload, "fingerprint\\.columns")


def test_release_scale_schema_rejects_metrics_and_output_dimension_mismatch() -> None:
    payload = _valid_schema_payload()
    payload["metrics"]["output_columns"] = 13

    _assert_invalid_report(payload, "metrics\\.output_columns")


def test_release_scale_schema_rejects_project_runtime_version_mismatch() -> None:
    payload = _valid_schema_payload()
    payload["source_provenance"]["runtime"]["distribution_phospy_version"] = "1.5.2"

    _assert_invalid_report(payload, "distribution_phospy_version")


def test_release_scale_schema_rejects_unsupported_recorded_python() -> None:
    payload = _valid_schema_payload()
    payload["source_provenance"]["runtime"]["python_version"] = "3.13.0"
    payload["source_provenance"]["runtime"]["python_version_info"] = [
        3,
        13,
        0,
        "final",
        0,
    ]
    payload["environment"]["python"]["version"] = "3.13.0 test"
    payload["environment"]["python"]["version_info"] = [3, 13, 0, "final", 0]

    _assert_invalid_report(payload, "requires_python")


def test_release_scale_schema_rejects_bad_source_digest_format() -> None:
    payload = _valid_schema_payload()
    payload["source_provenance"]["digests"]["source_tree"]["sha256"] = "A" * 64

    _assert_invalid_report(payload, "sha256")


def test_release_scale_schema_rejects_tested_feature_count_mismatch() -> None:
    payload = _valid_schema_payload()
    payload["scientific_summary"]["tested_feature_count"] = EXPECTED_SITES - 1

    _assert_invalid_report(payload, "tested_feature_count")


def test_release_scale_schema_rejects_scientific_summary_digest_mismatch() -> None:
    payload = _valid_schema_payload()
    payload["metrics"]["scientific_summary_digest"] = "f" * 64

    _assert_invalid_report(payload, "scientific_summary_digest")


def test_release_scale_schema_rejects_repeated_table_digest_mismatch() -> None:
    payload = _valid_schema_payload()
    payload["output_fingerprints"]["differential_result_table_digest"] = "f" * 64

    _assert_invalid_report(payload, "differential_result_table_digest")


def test_release_scale_source_tree_digest_is_deterministic(tmp_path: Path) -> None:
    module = _load_benchmark_script()
    root = tmp_path / "src" / "phospy"
    root.mkdir(parents=True)
    (root / "b.py").write_text("B = 2\n", encoding="utf-8")
    (root / "a.py").write_text("A = 1\n", encoding="utf-8")
    cache = root / "__pycache__"
    cache.mkdir()
    (cache / "ignored.pyc").write_bytes(b"ignored")

    first = module._source_tree_digest(root, repo_root=tmp_path)
    second = module._source_tree_digest(root, repo_root=tmp_path)

    assert first["sha256"] == second["sha256"]
    assert first["file_count"] == 2
    assert DIGEST_RE.fullmatch(first["sha256"])


def test_release_scale_generation_uses_actual_differential_table_dimensions(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    module = _load_benchmark_script()
    config = module.ReleaseScaleBenchmarkConfig(n_sites=3, n_samples=4)
    table = pd.DataFrame(
        [[float(row + column) for column in range(7)] for row in range(3)]
    )
    digest = "a" * 64
    timings = {
        "request_preparation_seconds": 0.01,
        "builder_execution_seconds": 0.02,
        "differential_analysis_seconds": 0.03,
        "result_table_assembly_seconds": 0.04,
    }
    summary_payload = {
        "matrix_dimensions": {"rows": 3, "columns": 4},
        "original_missing_cell_count": 0,
        "final_missing_cell_count": 0,
        "tested_feature_count": 3,
        "differential_result_table_fingerprint": {
            "rows": 3,
            "columns": 7,
            "hash_value": digest,
        },
    }

    monkeypatch.setattr(
        module,
        "validate_benchmark_source_environment",
        _dummy_source_provenance,
    )
    monkeypatch.setattr(module, "_assert_release_scale_outputs", lambda **_: None)
    monkeypatch.setattr(
        module,
        "_run_release_scale_workflow",
        lambda _: module._ReleaseScaleWorkflowMeasurement(
            dataset=object(),
            result=object(),
            table=table,
            timings=timings,
            original_missing_cell_count=0,
        ),
    )
    monkeypatch.setattr(
        module,
        "_build_release_scale_scientific_summary",
        lambda **_: module._ScientificSummaryMeasurement(
            summary=module._ReleaseScaleScientificSummary.from_payload(summary_payload),
            timings={"result_table_fingerprinting_seconds": 0.05},
        ),
    )

    result = module.run_benchmark(config)

    assert result.metrics["output_rows"] == 3
    assert result.metrics["output_columns"] == 7


def test_release_scale_generation_fails_stale_metadata_before_workload(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    module = _load_benchmark_script()
    called = False
    project_version = _project_version()

    def fake_run_benchmark(**_: object) -> object:
        nonlocal called
        called = True
        raise AssertionError("workload should not run")

    def fake_distribution_version(name: str) -> str:
        if name == "phospy":
            return "1.5.2"
        return "0.0"

    monkeypatch.setattr(
        module, "_current_python_version_info", lambda: (3, 12, 10, "final", 0)
    )
    monkeypatch.setattr(module.phospy, "__version__", project_version)
    monkeypatch.setattr(module.importlib.metadata, "version", fake_distribution_version)
    monkeypatch.setattr(module, "run_benchmark", fake_run_benchmark)

    with pytest.raises(module.BenchmarkSourceProvenanceError, match="distribution"):
        module.main([])

    assert called is False


def test_release_scale_generation_fails_import_origin_before_workload(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    module = _load_benchmark_script()
    called = False
    project_version = _project_version()

    def fake_run_benchmark(**_: object) -> object:
        nonlocal called
        called = True
        raise AssertionError("workload should not run")

    def fake_distribution_version(name: str) -> str:
        if name == "phospy":
            return project_version
        return "0.0"

    wrong_origin = tmp_path / "site-packages" / "phospy" / "__init__.py"
    wrong_origin.parent.mkdir(parents=True)
    wrong_origin.write_text("__version__ = '1.6.0'\n", encoding="utf-8")

    monkeypatch.setattr(
        module, "_current_python_version_info", lambda: (3, 12, 10, "final", 0)
    )
    monkeypatch.setattr(module.phospy, "__version__", project_version)
    monkeypatch.setattr(module.phospy, "__file__", str(wrong_origin))
    monkeypatch.setattr(module.importlib.metadata, "version", fake_distribution_version)
    monkeypatch.setattr(module, "run_benchmark", fake_run_benchmark)

    with pytest.raises(
        module.BenchmarkSourceProvenanceError, match="intended checkout"
    ):
        module.main([])

    assert called is False


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
    payload = json.loads(written_path.read_text(encoding="utf-8"))
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
    benchmark_recipe = makefile.split("benchmark-release-scale:", maxsplit=1)[1].split(
        "dataset-builder-demo:",
        maxsplit=1,
    )[0]
    test_performance_recipe = makefile.split("test-performance:", maxsplit=1)[1].split(
        "test-release-gates:", maxsplit=1
    )[0]

    assert "benchmark-release-scale" not in release_check_line
    assert (
        "measure_release_scale_builder_differential.py" not in test_performance_recipe
    )
    assert "benchmark-release-scale:" in makefile
    assert "measure_release_scale_builder_differential.py" in benchmark_recipe


def test_release_scale_benchmark_docs_use_observation_only_language() -> None:
    performance_doc = (REPO_ROOT / "docs" / "performance.md").read_text(
        encoding="utf-8"
    )
    release_scale_section = performance_doc.split(
        "## Optional Release-Scale Local Benchmark",
        maxsplit=1,
    )[1]
    normalized = release_scale_section.lower()
    compact = re.sub(r"\s+", " ", normalized)

    assert "dated observation" in normalized
    assert "not independent proof" in compact
    assert "not scientific validation" in normalized
    assert "external parity" in normalized
    assert "portable guarantee" in normalized
    assert "not a release failure" in normalized
