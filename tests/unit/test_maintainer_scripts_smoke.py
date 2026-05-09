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

ACTIVE_SCRIPT_DIR = REPO_ROOT / "scripts" / "active"
SUPPORT_SCRIPT_DIR = REPO_ROOT / "scripts" / "support"

ACTIVE_PYTHON_SCRIPTS = sorted(ACTIVE_SCRIPT_DIR.glob("*.py"))
MAINTAINER_PYTHON_SCRIPTS = ACTIVE_PYTHON_SCRIPTS + [
    REPO_ROOT / "scripts" / "run_pyright.py",
    SUPPORT_SCRIPT_DIR / "public_workflow_reference.py",
]
ACTIVE_R_SCRIPT = ACTIVE_SCRIPT_DIR / "generate_r_l6_fixtures.R"
SCRIPTS_README = REPO_ROOT / "scripts" / "README.md"

LOCAL_IMPORT_ROOTS = {"phospy", "tests", "scripts"}
GENERATOR_OUTPUT_TOKENS: dict[str, tuple[str, ...]] = {
    "generate_l6_prediction_parity_fixtures.py": (
        "tests/fixtures/rewrite_parity/r_reference_l6_prediction",
        "native_profile_scores.csv",
        "predMat.csv",
    ),
    "generate_provenance_goldens.py": (
        "tests/fixtures/public_workflow_reference",
        "kinase_public_predmat_provenance_golden.json",
        "signalome_l6_provenance_golden.json",
    ),
    "generate_public_predmat_rewrite_reference.py": (
        "tests/fixtures/public_workflow_reference",
        "predmat_rewrite_stable.csv",
        "predmat_rewrite_r_parity.csv",
    ),
    "generate_signalome_public_workflow_reference.py": (
        "tests/fixtures/public_workflow_reference",
        "signalome_rewrite_l6_contract.json",
        "signalome_rewrite_l6_module_assignments.csv",
    ),
}
GENERATOR_SOURCE_MARKERS: dict[str, tuple[str, ...]] = {
    "generate_l6_prediction_parity_fixtures.py": (
        "DEFAULT_OUTPUT_DIR",
        "r_reference_l6_prediction",
        "predMat.csv",
    ),
    "generate_provenance_goldens.py": (
        "PUBLIC_WORKFLOW_REFERENCE",
        "KINASE_GOLDEN_PATH",
        "SIGNALOME_GOLDEN_PATH",
    ),
    "generate_public_predmat_rewrite_reference.py": (
        "DEFAULT_OUTPUT_DIR",
        "public_workflow_reference",
        "predmat_rewrite_stable.csv",
    ),
    "generate_signalome_public_workflow_reference.py": (
        "DEFAULT_OUTPUT_DIR",
        "contract_path",
        "signalome_rewrite_l6_contract.json",
    ),
}
GENERATOR_OUTPUT_PATHS: dict[str, tuple[Path, ...]] = {
    "generate_l6_prediction_parity_fixtures.py": (
        REPO_ROOT
        / "tests"
        / "fixtures"
        / "rewrite_parity"
        / "r_reference_l6_prediction"
        / "native_profile_scores.csv",
        REPO_ROOT
        / "tests"
        / "fixtures"
        / "rewrite_parity"
        / "r_reference_l6_prediction"
        / "predMat.csv",
    ),
    "generate_provenance_goldens.py": (
        REPO_ROOT
        / "tests"
        / "fixtures"
        / "public_workflow_reference"
        / "kinase_public_predmat_provenance_golden.json",
        REPO_ROOT
        / "tests"
        / "fixtures"
        / "public_workflow_reference"
        / "signalome_l6_provenance_golden.json",
    ),
    "generate_public_predmat_rewrite_reference.py": (
        REPO_ROOT
        / "tests"
        / "fixtures"
        / "public_workflow_reference"
        / "predmat_rewrite_stable.csv",
        REPO_ROOT
        / "tests"
        / "fixtures"
        / "public_workflow_reference"
        / "predmat_rewrite_r_parity.csv",
    ),
    "generate_signalome_public_workflow_reference.py": (
        REPO_ROOT
        / "tests"
        / "fixtures"
        / "public_workflow_reference"
        / "signalome_rewrite_l6_contract.json",
        REPO_ROOT
        / "tests"
        / "fixtures"
        / "public_workflow_reference"
        / "signalome_rewrite_l6_module_assignments.csv",
    ),
}


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


def _has_argparse_description(tree: ast.Module) -> bool:
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        func = node.func
        is_argparse_ctor = (
            isinstance(func, ast.Attribute)
            and isinstance(func.value, ast.Name)
            and func.value.id == "argparse"
            and func.attr == "ArgumentParser"
        ) or (isinstance(func, ast.Name) and func.id == "ArgumentParser")
        if not is_argparse_ctor:
            continue
        for keyword in node.keywords:
            if keyword.arg != "description":
                continue
            value = keyword.value
            if isinstance(value, ast.Constant) and isinstance(value.value, str):
                return bool(value.value.strip())
            return True
    return False


def _file_fingerprint(path: Path) -> tuple[bool, int | None, int | None]:
    if not path.exists():
        return (False, None, None)
    stat = path.stat()
    return (True, int(stat.st_size), int(stat.st_mtime_ns))


def _import_script(script_path: Path) -> None:
    module_name = f"maintainer_script_smoke_{script_path.stem}"
    spec = importlib.util.spec_from_file_location(module_name, script_path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    try:
        spec.loader.exec_module(module)
    finally:
        sys.modules.pop(module_name, None)


def test_active_script_inventory_is_non_empty() -> None:
    assert ACTIVE_PYTHON_SCRIPTS, "expected active maintainer Python scripts"


@pytest.mark.parametrize(
    "script_path",
    MAINTAINER_PYTHON_SCRIPTS,
    ids=lambda path: path.relative_to(REPO_ROOT).as_posix(),
)
def test_maintainer_python_scripts_parse(script_path: Path) -> None:
    _parse(script_path)


@pytest.mark.parametrize(
    "script_path",
    MAINTAINER_PYTHON_SCRIPTS,
    ids=lambda path: path.relative_to(REPO_ROOT).as_posix(),
)
def test_maintainer_python_scripts_import_smoke(script_path: Path) -> None:
    _import_script(script_path)


@pytest.mark.parametrize(
    "script_path",
    MAINTAINER_PYTHON_SCRIPTS,
    ids=lambda path: path.relative_to(REPO_ROOT).as_posix(),
)
def test_maintainer_python_local_imports_resolve(script_path: Path) -> None:
    for module_name in _iter_local_import_modules(_parse(script_path)):
        assert importlib.util.find_spec(module_name) is not None, (
            f"{script_path.name} references missing local module '{module_name}'"
        )


@pytest.mark.parametrize(
    "script_path",
    MAINTAINER_PYTHON_SCRIPTS,
    ids=lambda path: path.relative_to(REPO_ROOT).as_posix(),
)
def test_maintainer_python_scripts_have_docstring_or_argparse_description(
    script_path: Path,
) -> None:
    tree = _parse(script_path)
    docstring = ast.get_docstring(tree)
    has_docstring = bool(docstring and docstring.strip())
    assert has_docstring or _has_argparse_description(tree), (
        f"{script_path.name} must include a module docstring or argparse description"
    )


@pytest.mark.parametrize(
    "script_path",
    ACTIVE_PYTHON_SCRIPTS,
    ids=lambda path: path.relative_to(REPO_ROOT).as_posix(),
)
def test_active_generators_do_not_write_fixtures_on_import(script_path: Path) -> None:
    output_paths = GENERATOR_OUTPUT_PATHS.get(script_path.name)
    assert output_paths is not None, f"missing output path map for {script_path.name}"
    before = {path: _file_fingerprint(path) for path in output_paths}
    _import_script(script_path)
    after = {path: _file_fingerprint(path) for path in output_paths}
    assert before == after, (
        f"{script_path.name} changed fixture outputs during import; generation must "
        "stay inside explicit CLI entrypoints"
    )


def test_active_r_script_is_present_and_documented() -> None:
    assert ACTIVE_R_SCRIPT.exists(), "expected active R maintainer script to exist"
    source = _read_source(ACTIVE_R_SCRIPT)
    assert "required_pkgs" in source, "R script must declare required package contract"
    readme = _read_source(SCRIPTS_README)
    assert ACTIVE_R_SCRIPT.name in readme, (
        f"scripts/README.md must document {ACTIVE_R_SCRIPT.name}"
    )


def test_active_generator_output_locations_are_explicit_and_documented() -> None:
    readme = _read_source(SCRIPTS_README)
    for script_name, tokens in GENERATOR_OUTPUT_TOKENS.items():
        source = _read_source(ACTIVE_SCRIPT_DIR / script_name)
        assert script_name in readme, f"scripts/README.md must list {script_name}"
        for marker in GENERATOR_SOURCE_MARKERS[script_name]:
            assert marker in source, (
                f"{script_name} must keep explicit output-location markers ({marker})"
            )
        for token in tokens:
            assert token in readme, (
                f"scripts/README.md must document '{token}' for {script_name}"
            )
