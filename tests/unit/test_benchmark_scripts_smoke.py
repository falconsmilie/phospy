from __future__ import annotations

import ast
import importlib.util
import sys
from pathlib import Path

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
    "phospy.prediction.aggregation",
    "phospy.prediction.sampling",
    "phospy.preprocessing",
    "phospy.signalomes.analysis",
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


def _read_source(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def _parse(path: Path) -> ast.Module:
    return ast.parse(_read_source(path), filename=str(path))


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
    module_name = f"benchmark_smoke_{script_path.stem}"
    spec = importlib.util.spec_from_file_location(module_name, script_path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    try:
        spec.loader.exec_module(module)
    finally:
        sys.modules.pop(module_name, None)


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
