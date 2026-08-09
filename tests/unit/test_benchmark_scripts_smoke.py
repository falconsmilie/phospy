from __future__ import annotations

import ast
import importlib.util
import sys
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

BENCHMARK_DIR = REPO_ROOT / "benchmarks"
BENCHMARK_SCRIPTS = sorted(BENCHMARK_DIR.glob("*.py"))
LOCAL_IMPORT_ROOTS = {"phospy", "tests"}
DISALLOWED_SEAM_TOKENS = (
    "phospy.science.prediction.aggregation",
    "phospy.science.prediction.sampling",
    "phospy.preprocessing",
    "phospy.science.signalomes.analysis",
    "tests/test_parity-with_metrics.py",
    "tests/test_end_to_end_parity.py",
)
WRITE_CALL_TOKENS = (
    ".to_csv(",
    ".to_json(",
    ".to_parquet(",
    ".to_pickle(",
    ".write_text(",
    ".write_bytes(",
    "open(",
)
REPORTS_DIRECTORY_TOKEN = "benchmarks/reports"
RELEASE_SCALE_BENCHMARK_SCRIPT = (
    BENCHMARK_DIR / "measure_release_scale_builder_differential.py"
)
REPEATED_WORKFLOW_BENCHMARK_SCRIPT = (
    BENCHMARK_DIR / "measure_repeated_workflow_dataset_snapshot_reuse.py"
)


def _read_source(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def _parse(path: Path) -> ast.Module:
    return ast.parse(_read_source(path), filename=str(path))


def _load_script_module(script_path: Path) -> ModuleType:
    module_name = f"benchmark_smoke_{script_path.stem}"
    spec = importlib.util.spec_from_file_location(module_name, script_path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    try:
        spec.loader.exec_module(module)
    finally:
        sys.modules.pop(module_name, None)
    return module


def _iter_local_import_modules(tree: ast.Module) -> list[str]:
    modules: list[str] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                root = alias.name.split(".", maxsplit=1)[0]
                if root in LOCAL_IMPORT_ROOTS:
                    modules.append(alias.name)
        elif isinstance(node, ast.ImportFrom):
            if node.level != 0 or node.module is None:
                continue
            root = node.module.split(".", maxsplit=1)[0]
            if root in LOCAL_IMPORT_ROOTS:
                modules.append(node.module)
    return sorted(set(modules))


def _find_main_function(
    tree: ast.Module,
) -> ast.FunctionDef | ast.AsyncFunctionDef | None:
    for node in tree.body:
        if (
            isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
            and node.name == "main"
        ):
            return node
    return None


def _script_has_key_value_print_in_main(script_path: Path) -> bool:
    source = _read_source(script_path)
    tree = _parse(script_path)
    main_fn = _find_main_function(tree)
    if main_fn is None:
        return False
    for node in ast.walk(main_fn):
        if not isinstance(node, ast.Call):
            continue
        if not isinstance(node.func, ast.Name) or node.func.id != "print":
            continue
        for argument in node.args:
            source_segment = ast.get_source_segment(source, argument)
            if source_segment is not None and "=" in source_segment:
                return True
    return False


def test_benchmark_script_inventory_is_non_empty() -> None:
    assert BENCHMARK_SCRIPTS, "expected at least one benchmark script"


@pytest.mark.parametrize("script_path", BENCHMARK_SCRIPTS, ids=lambda path: path.name)
def test_benchmark_scripts_have_module_header(script_path: Path) -> None:
    docstring = ast.get_docstring(_parse(script_path))
    assert docstring and docstring.strip(), (
        f"{script_path.name} must include a short module header describing scope"
    )


@pytest.mark.parametrize("script_path", BENCHMARK_SCRIPTS, ids=lambda path: path.name)
def test_benchmark_scripts_import_smoke(script_path: Path) -> None:
    _load_script_module(script_path)


@pytest.mark.parametrize("script_path", BENCHMARK_SCRIPTS, ids=lambda path: path.name)
def test_benchmark_scripts_local_imports_resolve(script_path: Path) -> None:
    for module_name in _iter_local_import_modules(_parse(script_path)):
        assert importlib.util.find_spec(module_name) is not None, (
            f"{script_path.name} references missing local module '{module_name}'"
        )


@pytest.mark.parametrize("script_path", BENCHMARK_SCRIPTS, ids=lambda path: path.name)
def test_benchmark_scripts_do_not_reference_retired_seams(script_path: Path) -> None:
    source = _read_source(script_path)
    for token in DISALLOWED_SEAM_TOKENS:
        assert token not in source, (
            f"{script_path.name} references retired benchmark seam '{token}'"
        )


@pytest.mark.parametrize("script_path", BENCHMARK_SCRIPTS, ids=lambda path: path.name)
def test_benchmark_scripts_emit_key_value_metric_from_main(script_path: Path) -> None:
    assert _script_has_key_value_print_in_main(script_path), (
        f"{script_path.name} must print at least one key=value metric from main()"
    )


@pytest.mark.parametrize("script_path", BENCHMARK_SCRIPTS, ids=lambda path: path.name)
def test_benchmark_scripts_write_only_to_reports_directory(script_path: Path) -> None:
    source = _read_source(script_path)
    if any(token in source for token in WRITE_CALL_TOKENS):
        assert REPORTS_DIRECTORY_TOKEN in source, (
            f"{script_path.name} may only write report artefacts under "
            f"'{REPORTS_DIRECTORY_TOKEN}/'"
        )


def test_release_scale_benchmark_defaults_are_declared_without_dataset_build() -> None:
    module = _load_script_module(RELEASE_SCALE_BENCHMARK_SCRIPT)
    config = module.default_config()

    assert config.n_sites == 50_000
    assert config.n_samples == 48
    assert config.missing_fraction == 0.03


def test_release_scale_benchmark_declares_linear_imputation_input_scale() -> None:
    module = _load_script_module(RELEASE_SCALE_BENCHMARK_SCRIPT)
    request = module._build_release_scale_dataset_request(
        module.ReleaseScaleBenchmarkConfig(n_sites=8, n_samples=4)
    )

    input_scale = request.preprocessing_config.missing_data.input_scale

    assert getattr(input_scale, "value", input_scale) == "linear"


def test_release_scale_benchmark_documents_main_key_value_metrics() -> None:
    module = _load_script_module(RELEASE_SCALE_BENCHMARK_SCRIPT)

    assert {
        "total_runtime_seconds",
        "request_preparation_seconds",
        "builder_execution_seconds",
        "preprocessing_execution_seconds",
        "provenance_fingerprinting_seconds",
        "differential_analysis_seconds",
        "result_table_assembly_seconds",
        "tested_feature_count",
        "process_rss_peak_mib",
    }.issubset(set(module.DOCUMENTED_MAIN_METRIC_KEYS))


def test_release_scale_benchmark_metric_formatter_emits_key_value_pairs() -> None:
    module = _load_script_module(RELEASE_SCALE_BENCHMARK_SCRIPT)
    metrics = {key: "value" for key in module.DOCUMENTED_MAIN_METRIC_KEYS}

    formatted = module.format_metric_values(metrics)

    for key in module.DOCUMENTED_MAIN_METRIC_KEYS:
        assert f"{key}=" in formatted


def test_release_scale_benchmark_report_path_is_under_benchmarks_reports() -> None:
    module = _load_script_module(RELEASE_SCALE_BENCHMARK_SCRIPT)
    report_path = module.default_report_path()

    assert report_path == BENCHMARK_DIR / "reports" / (
        "release-scale-builder-differential.json"
    )
    assert "tests" not in report_path.relative_to(REPO_ROOT).parts


def test_release_scale_benchmark_report_payload_records_environment_and_outputs() -> (
    None
):
    module = _load_script_module(RELEASE_SCALE_BENCHMARK_SCRIPT)
    result = module.ReleaseScaleBenchmarkResult(
        config=module.default_config(),
        timings={"total_runtime_seconds": 1.25, "builder_execution_seconds": 0.5},
        metrics={
            "output_rows": 50_000,
            "output_columns": 12,
            "tested_feature_count": 50_000,
            "total_runtime_seconds": 1.25,
            "process_rss_peak_mib": 256.0,
            "scientific_summary_digest": "summary-digest",
        },
        scientific_summary={
            "input_table_fingerprints_digest": "input-digest",
            "output_table_fingerprints_digest": "output-digest",
            "preprocessing_trace_digest": "trace-digest",
            "differential_result_table_digest": "result-digest",
        },
    )

    payload = module.report_payload(result)

    assert payload["python"]["version"]
    assert payload["environment"]["python"]["executable"]
    assert {"numpy", "pandas", "scipy", "phospy"} <= set(payload["dependencies"])
    assert payload["machine"]["platform"]
    assert payload["runtime"]["total_runtime_seconds"] == 1.25
    assert payload["peak_memory"]["process_rss_peak_mib"] == 256.0
    assert payload["output_fingerprints"] == {
        "input_table_fingerprints_digest": "input-digest",
        "output_table_fingerprints_digest": "output-digest",
        "preprocessing_trace_digest": "trace-digest",
        "differential_result_table_digest": "result-digest",
        "scientific_summary_digest": "summary-digest",
    }


def test_repeated_workflow_snapshot_reuse_benchmark_schema() -> None:
    module = _load_script_module(REPEATED_WORKFLOW_BENCHMARK_SCRIPT)
    measurement = module.RunMeasurement(
        runtime_seconds=0.25,
        peak_tracemalloc_mib=12.5,
        full_frame_deep_copies={"dataset.phospho": 1},
        projected_frame_deep_copies={"dataset.site_metadata projected": 1},
        snapshot_constructions={"dataset.phospho internal snapshot": 1},
    )
    result = module.BenchmarkResult(
        config=module.BenchmarkConfig(
            n_sites=40,
            n_samples=4,
            n_kinases=3,
            substrates_per_kinase=8,
        ),
        dataset_dimensions={
            "phospho_rows": 40,
            "phospho_columns": 4,
            "site_metadata_rows": 40,
            "site_metadata_columns": 10,
        },
        frame_dtypes={
            "phospho": {"A_1": "float64"},
            "site_metadata": {"site_key": "object"},
        },
        setup=measurement,
        differential_first=measurement,
        differential_repeated=measurement,
        kinase_first=measurement,
        kinase_repeated=measurement,
        total_full_frame_deep_copies={"dataset.phospho": 1},
        total_projected_frame_deep_copies={"dataset.site_metadata projected": 1},
        total_snapshot_constructions={"dataset.phospho internal snapshot": 1},
        environment={"python_version": "3.x", "platform": "test"},
        dependencies={"phospy": "test", "numpy": "test", "pandas": "test"},
    )

    payload = module.report_payload(result)

    assert {
        "dataset_phospho_rows",
        "dataset_phospho_columns",
        "frame_dtypes_json",
        "setup_run_seconds",
        "differential_first_run_seconds",
        "differential_repeated_run_seconds",
        "kinase_first_run_seconds",
        "kinase_repeated_run_seconds",
        "setup_full_frame_deep_copy_counts_json",
        "per_run_full_frame_deep_copy_counts_json",
        "full_frame_deep_copy_counts_json",
        "projected_frame_deep_copy_counts_json",
        "per_run_projected_frame_deep_copy_counts_json",
        "setup_snapshot_construction_counts_json",
        "per_run_snapshot_construction_counts_json",
        "snapshot_construction_counts_json",
        "environment_json",
        "dependency_versions_json",
    }.issubset(set(module.DOCUMENTED_MAIN_METRIC_KEYS))
    assert payload["benchmark"] == "repeated_workflow_dataset_snapshot_reuse"
    assert payload["dataset"]["dimensions"]["phospho_rows"] == 40
    assert payload["dataset"]["frame_dtypes"]["phospho"]["A_1"] == "float64"
    assert payload["setup"]["full_frame_deep_copies"] == {"dataset.phospho": 1}
    assert payload["runs"]["differential"]["first"]["runtime_seconds"] == 0.25
    assert payload["runs"]["kinase"]["repeated"]["peak_tracemalloc_mib"] == 12.5
    assert payload["totals"]["full_frame_deep_copies"] == {"dataset.phospho": 1}
    assert payload["totals"]["projected_frame_deep_copies"] == {
        "dataset.site_metadata projected": 1
    }
    assert payload["totals"]["snapshot_constructions"] == {
        "dataset.phospho internal snapshot": 1
    }
    assert payload["environment"]["python_version"] == "3.x"
    assert payload["dependencies"]["pandas"] == "test"
    assert "machine-specific" in payload["observation_scope"].lower()


def test_repeated_workflow_benchmark_copy_attribution_uses_frame_origin() -> None:
    module = _load_script_module(REPEATED_WORKFLOW_BENCHMARK_SCRIPT)
    source = pd.DataFrame({"value": [1.0, 2.0]}, index=pd.Index(["a", "b"]))
    unrelated_same_shape = pd.DataFrame(
        {"value": [3.0, 4.0]},
        index=pd.Index(["a", "b"]),
    )

    with module._instrument_copy_accounting() as counts:
        module._mark_frame_source(counts, source, label="dataset.phospho")
        source.copy(deep=True)
        unrelated_same_shape.copy(deep=True)

    assert dict(counts.full_frame_deep_copies) == {"dataset.phospho": 1}
