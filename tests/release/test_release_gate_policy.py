from __future__ import annotations

import ast
import re
from pathlib import Path
from typing import Any

import pytest

try:
    import tomllib  # pyright: ignore[reportMissingImports]
except ModuleNotFoundError:  # pragma: no cover - Python 3.10 fallback
    import tomli as tomllib  # pyright: ignore[reportMissingImports]


ROOT = Path(__file__).resolve().parents[2]

pytestmark = pytest.mark.release_gate

AUTHORITATIVE_RELEASE_COMMAND = "make release-check"
PUBLIC_RELEASE_INSTRUCTION_DOCS = (
    Path("README.md"),
    Path("docs/maintenance.md"),
    Path("docs/testing/README.md"),
    Path("docs/testing/pytest_markers.md"),
    Path("docs/contributing.md"),
    Path(".github/CONTRIBUTING.md"),
)
RELEASE_CHECK_DEPENDENCIES = (
    "lint",
    "type-check",
    "test-unit",
    "test-parity",
    "test-performance",
    "validate-reference-bundles",
    "test-release-gates",
    "build",
    "verify-installed-distributions",
)
PROHIBITED_RELEASE_SCALE_ROWS = 50_000
PROHIBITED_RELEASE_SCALE_SAMPLES = 48
REQUIRED_TEST_SOURCE_ROOTS = (
    Path("tests/unit"),
    Path("tests/integration"),
    Path("tests/parity"),
    Path("tests/workflows"),
    Path("tests/validation"),
    Path("tests/science"),
    Path("tests/architecture"),
    Path("tests/performance"),
    Path("tests/release"),
    Path("tests/golden"),
    Path("tests/support"),
)
STATIC_RELEASE_SCALE_POLICY_EXCLUSIONS = {
    Path("tests/release/test_release_gate_policy.py"),
}
RELEASE_SCALE_BENCHMARK_PATH = Path(
    "benchmarks/measure_release_scale_builder_differential.py"
)
RELEASE_SCALE_REMOVED_PYTEST_PATH = Path(
    "tests/performance/test_end_to_end_release_scale_contract.py"
)
RELEASE_SCALE_BENCHMARK_TARGET = "benchmark-release-scale"
RELEASE_SCALE_BENCHMARK_TOKENS = (
    RELEASE_SCALE_BENCHMARK_TARGET,
    RELEASE_SCALE_BENCHMARK_PATH.as_posix(),
    str(RELEASE_SCALE_BENCHMARK_PATH).replace("/", "\\"),
    RELEASE_SCALE_REMOVED_PYTEST_PATH.as_posix(),
    str(RELEASE_SCALE_REMOVED_PYTEST_PATH).replace("/", "\\"),
)


def _read(relative_path: str | Path) -> str:
    return (ROOT / relative_path).read_text(encoding="utf-8").replace("\r\n", "\n")


def _load_pyproject() -> dict[str, Any]:
    with (ROOT / "pyproject.toml").open("rb") as handle:
        return tomllib.load(handle)


def _pytest_config() -> dict[str, Any]:
    pyproject = _load_pyproject()
    return pyproject["tool"]["pytest"]["ini_options"]


def _pytest_markers() -> set[str]:
    markers = _pytest_config()["markers"]
    return {marker.split(":", maxsplit=1)[0] for marker in markers}


def _make_target_line(target_name: str) -> str:
    lines = _read("Makefile").splitlines()
    return next(line for line in lines if line.startswith(f"{target_name}:"))


def _make_target_body(target_name: str) -> str:
    lines = _read("Makefile").splitlines()
    start_index = next(
        index for index, line in enumerate(lines) if line.startswith(f"{target_name}:")
    )
    body: list[str] = []
    for line in lines[start_index + 1 :]:
        if line.startswith("\t"):
            body.append(line[1:])
            continue
        if line.strip() or body:
            break
    return "\n".join(body)


def _workflow_job_block(workflow_text: str, job_name: str) -> str:
    match = re.search(
        rf"(?ms)^  {re.escape(job_name)}:\n(?P<body>.*?)(?=^  [A-Za-z0-9_-]+:|\Z)",
        workflow_text,
    )
    assert match is not None, f"Workflow job is missing: {job_name}"
    return match.group("body")


def _has_needs(job_block: str, dependency: str) -> bool:
    list_match = re.search(
        r"(?ms)^\s+needs:\s*\n(?P<items>(?:\s+- .+\n)+)",
        job_block,
    )
    if list_match is not None:
        return (
            re.search(
                rf"(?m)^\s+- {re.escape(dependency)}\s*$",
                list_match.group("items"),
            )
            is not None
        )
    return (
        re.search(rf"(?m)^\s+needs:\s*{re.escape(dependency)}\s*$", job_block)
        is not None
    )


def _make_target_dependencies(target_name: str) -> tuple[str, ...]:
    target_line = _make_target_line(target_name)
    return tuple(target_line.split(":", maxsplit=1)[1].split())


def _make_dependency_closure(target_name: str) -> set[str]:
    seen: set[str] = set()
    pending = list(_make_target_dependencies(target_name))
    while pending:
        dependency = pending.pop()
        if dependency in seen:
            continue
        seen.add(dependency)
        try:
            pending.extend(_make_target_dependencies(dependency))
        except StopIteration:
            continue
    return seen


def _assert_supported_python_matrix(job_block: str) -> None:
    assert "python-version: ['3.10', '3.11', '3.12']" in job_block


def _iter_required_test_python_sources() -> tuple[Path, ...]:
    paths: set[Path] = set()
    for relative_root in REQUIRED_TEST_SOURCE_ROOTS:
        root = ROOT / relative_root
        if root.is_file() and root.suffix == ".py":
            if root.relative_to(ROOT) not in STATIC_RELEASE_SCALE_POLICY_EXCLUSIONS:
                paths.add(root)
            continue
        if root.is_dir():
            paths.update(
                path
                for path in root.rglob("*.py")
                if path.is_file()
                and path.relative_to(ROOT) not in STATIC_RELEASE_SCALE_POLICY_EXCLUSIONS
            )
    return tuple(sorted(paths))


def _parse_python_source(path: Path) -> ast.Module:
    return ast.parse(path.read_text(encoding="utf-8"), filename=str(path))


def _assigned_int_constants(tree: ast.Module) -> dict[str, int]:
    constants: dict[str, int] = {}
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign):
            value = _literal_int(node.value)
            if value is None:
                continue
            for target in node.targets:
                if isinstance(target, ast.Name):
                    constants[target.id] = value
        elif isinstance(node, ast.AnnAssign):
            value = _literal_int(node.value)
            if value is not None and isinstance(node.target, ast.Name):
                constants[node.target.id] = value
    return constants


def _literal_int(node: ast.AST | None) -> int | None:
    if isinstance(node, ast.Constant) and isinstance(node.value, int):
        return int(node.value)
    if (
        isinstance(node, ast.UnaryOp)
        and isinstance(node.op, ast.USub)
        and isinstance(node.operand, ast.Constant)
        and isinstance(node.operand.value, int)
    ):
        return -int(node.operand.value)
    return None


def _resolved_int(node: ast.AST | None, constants: dict[str, int]) -> int | None:
    value = _literal_int(node)
    if value is not None:
        return value
    if isinstance(node, ast.Name):
        return constants.get(node.id)
    return None


def _dimension_prefix(name: str, suffixes: tuple[str, ...]) -> str | None:
    for suffix in suffixes:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return None


def _fixture_constant_dimension_issues(path: Path, tree: ast.Module) -> list[str]:
    constants = _assigned_int_constants(tree)
    row_suffixes = (
        "_N_SITES",
        "_N_ROWS",
        "_SITE_COUNT",
        "_ROW_COUNT",
        "_FEATURE_COUNT",
        "_N_FEATURES",
    )
    sample_suffixes = (
        "_N_SAMPLES",
        "_SAMPLE_COUNT",
        "_COLUMN_COUNT",
        "_N_COLUMNS",
    )
    row_prefixes = {
        prefix
        for name, value in constants.items()
        if value == PROHIBITED_RELEASE_SCALE_ROWS
        for prefix in (_dimension_prefix(name, row_suffixes),)
        if prefix is not None
    }
    sample_prefixes = {
        prefix
        for name, value in constants.items()
        if value == PROHIBITED_RELEASE_SCALE_SAMPLES
        for prefix in (_dimension_prefix(name, sample_suffixes),)
        if prefix is not None
    }
    return [
        f"{path.relative_to(ROOT).as_posix()} defines prohibited 50,000x48 "
        f"fixture constants with prefix {prefix!r}"
        for prefix in sorted(row_prefixes & sample_prefixes)
    ]


def _fixture_construction_dimension_issues(path: Path, tree: ast.Module) -> list[str]:
    constants = _assigned_int_constants(tree)
    issues: list[str] = []
    row_keywords = {
        "n_sites",
        "site_count",
        "n_rows",
        "row_count",
        "n_features",
        "feature_count",
        "rows",
    }
    sample_keywords = {
        "n_samples",
        "sample_count",
        "n_columns",
        "column_count",
        "columns",
    }
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        keyword_values = {
            str(keyword.arg): _resolved_int(keyword.value, constants)
            for keyword in node.keywords
            if keyword.arg is not None
        }
        has_prohibited_rows = any(
            keyword_values.get(name) == PROHIBITED_RELEASE_SCALE_ROWS
            for name in row_keywords
        )
        has_prohibited_samples = any(
            keyword_values.get(name) == PROHIBITED_RELEASE_SCALE_SAMPLES
            for name in sample_keywords
        )
        has_prohibited_shape = any(
            _tuple_ints(argument, constants)
            == (PROHIBITED_RELEASE_SCALE_ROWS, PROHIBITED_RELEASE_SCALE_SAMPLES)
            for argument in (*node.args, *(keyword.value for keyword in node.keywords))
        )
        if has_prohibited_rows and has_prohibited_samples or has_prohibited_shape:
            issues.append(
                f"{path.relative_to(ROOT).as_posix()}:{node.lineno} constructs a "
                "prohibited required 50,000x48 fixture"
            )
    return issues


def _tuple_ints(
    node: ast.AST,
    constants: dict[str, int],
) -> tuple[int, ...] | None:
    if not isinstance(node, ast.Tuple | ast.List):
        return None
    values: list[int] = []
    for element in node.elts:
        value = _resolved_int(element, constants)
        if value is None:
            return None
        values.append(value)
    return tuple(values)


def test_default_pytest_keeps_parity_out_of_fast_local_loop() -> None:
    assert _pytest_config()["addopts"] == '-m "not parity"'


def test_release_check_command_is_maintained_project_authority() -> None:
    target_line = _make_target_line("release-check")
    target_dependencies = tuple(target_line.split(":", maxsplit=1)[1].split())

    assert AUTHORITATIVE_RELEASE_COMMAND == "make release-check"
    assert target_dependencies == RELEASE_CHECK_DEPENDENCIES
    assert _make_target_body("release-check") == ""


@pytest.mark.parametrize(
    "relative_path",
    PUBLIC_RELEASE_INSTRUCTION_DOCS,
    ids=lambda path: path.as_posix(),
)
def test_maintained_release_instructions_point_to_lightweight_release_check(
    relative_path: Path,
) -> None:
    text = _read(relative_path)
    normalized = re.sub(r"\s+", " ", text.lower())

    assert AUTHORITATIVE_RELEASE_COMMAND in text
    assert "normal ci/build confidence" in normalized
    assert "formal exact-source/exact-artifact attestation" in normalized


def test_make_release_check_covers_declared_blocking_marker_policy() -> None:
    markers = _pytest_markers()
    assert markers == {
        "unit",
        "integration",
        "parity",
        "performance",
        "slow",
        "release_gate",
        "reproducibility",
        "golden",
        "activity_parity",
        "parity_diagnostic",
    }

    makefile = _read("Makefile")
    assert re.search(r"(?ms)^\.PHONY:.*test-release-gates", makefile) is not None
    assert "make test-release-gates" in _make_target_body("help")

    assert '$(PYTEST) -m "not parity"' in _make_target_body("test-unit")
    assert (
        '$(PYTEST) tests/parity -m "parity and not parity_diagnostic" -s'
        in _make_target_body("test-parity")
    )
    assert 'tests/performance -m "performance or release_gate"' in _make_target_body(
        "test-performance"
    )
    release_gates = _make_target_body("test-release-gates")
    assert '$(MKDIR_P) "$(PYTEST_REPORT_DIR)"' in release_gates
    assert "$(PYTEST) -o addopts= tests/release tests/golden" in release_gates
    assert '-m "release_gate or golden or reproducibility"' in release_gates
    assert '--junitxml "$(PYTEST_REPORT_DIR)/release-gates.xml"' in release_gates
    assert (
        "$(PYTHON) scripts/validate_reference_bundle_index.py --repo-root ."
        in _make_target_body("validate-reference-bundles")
    )
    assert (
        "$(PYTHON) scripts/verify_installed_distributions.py --dist-dir dist "
        "--repo-root . --constraint constraints/ci.txt"
    ) in _make_target_body("verify-installed-distributions")


def test_release_scale_benchmark_is_local_optional_and_outside_pytest() -> None:
    benchmark_path = ROOT / RELEASE_SCALE_BENCHMARK_PATH

    assert benchmark_path.exists()
    assert benchmark_path.is_file()
    assert benchmark_path.relative_to(ROOT).parts[0] == "benchmarks"
    assert "tests" not in benchmark_path.relative_to(ROOT).parts
    assert not (ROOT / RELEASE_SCALE_REMOVED_PYTEST_PATH).exists()


def test_release_commands_do_not_invoke_release_scale_benchmark() -> None:
    makefile = _read("Makefile")
    release_check_line = _make_target_line("release-check")
    release_check_body = _make_target_body("release-check")
    test_performance_body = _make_target_body("test-performance")
    benchmark_body = _make_target_body(RELEASE_SCALE_BENCHMARK_TARGET)

    assert RELEASE_SCALE_BENCHMARK_TARGET in makefile
    assert "benchmarks/measure_release_scale_builder_differential.py" in benchmark_body
    assert RELEASE_SCALE_BENCHMARK_TARGET not in release_check_line
    assert release_check_body == ""
    for token in RELEASE_SCALE_BENCHMARK_TOKENS:
        assert token not in test_performance_body


def test_required_make_targets_do_not_depend_on_release_scale_benchmark() -> None:
    for target_name in ("test-performance", "test-release-gates", "release-check"):
        closure = _make_dependency_closure(target_name)
        assert RELEASE_SCALE_BENCHMARK_TARGET not in closure


def test_required_test_sources_do_not_define_50k_by_48_fixture_constants() -> None:
    issues: list[str] = []
    for path in _iter_required_test_python_sources():
        issues.extend(
            _fixture_constant_dimension_issues(path, _parse_python_source(path))
        )

    assert issues == []


def test_required_test_sources_do_not_construct_50k_by_48_fixtures() -> None:
    issues: list[str] = []
    for path in _iter_required_test_python_sources():
        issues.extend(
            _fixture_construction_dimension_issues(path, _parse_python_source(path))
        )

    assert issues == []


def test_github_actions_do_not_invoke_release_scale_benchmark() -> None:
    workflow_paths = sorted((ROOT / ".github/workflows").glob("*.yml"))
    assert workflow_paths

    for workflow_path in workflow_paths:
        workflow = workflow_path.read_text(encoding="utf-8").replace("\r\n", "\n")
        for token in RELEASE_SCALE_BENCHMARK_TOKENS:
            assert token not in workflow


def test_release_scale_policy_docs_are_local_optional_not_release_blocking() -> None:
    performance_doc = _read("docs/performance.md")
    normalized = performance_doc.lower()
    compact = re.sub(r"\s+", " ", performance_doc.lower())

    assert "optional local release-scale benchmark" in normalized
    assert "excluded from `make test-performance`" in normalized
    assert "excluded from `make test-performance`, `make release-check`" in compact
    assert "tests/performance/test_end_to_end_release_scale_contract.py" not in (
        performance_doc
    )
    assert "release-scale gate" not in normalized
    assert "tracemalloc peak < 4" not in normalized
    assert "two consecutive successful" not in normalized


def test_make_build_is_conventional_and_git_independent() -> None:
    build = _make_target_body("build")
    verifier = _make_target_body("verify-installed-distributions")

    assert "build" in _make_target_dependencies("verify-installed-distributions")
    assert "$(RM_RF) dist" in build
    assert "$(BUILD)" in build
    assert "dist/*.whl" in build
    assert "dist/*.tar.gz" in build
    assert "$(TWINE) check dist/*" in build
    assert (
        "$(PYTHON) scripts/validate_reference_bundle_distribution.py "
        "--no-git-index-compare dist/*"
    ) in build
    assert build.index("$(RM_RF) dist") < build.index("$(BUILD)")
    assert build.index("$(BUILD)") < build.index("$(TWINE) check dist/*")
    assert build.index("$(TWINE) check dist/*") < build.index(
        "scripts/validate_reference_bundle_distribution.py"
    )
    assert "scripts/verify_installed_distributions.py" not in build
    assert "scripts/validate_reference_bundle_distribution.py" not in verifier
    assert "scripts/verify_installed_distributions.py" in verifier
    assert "--dist-dir dist" in verifier
    assert "--repo-root ." in verifier
    assert "--constraint constraints/ci.txt" in verifier


def test_installed_distribution_verifier_is_standalone_release_tooling() -> None:
    verifier_path = ROOT / "scripts/verify_installed_distributions.py"
    source = verifier_path.read_text(encoding="utf-8")
    tree = ast.parse(source, filename=str(verifier_path))
    imported_modules: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported_modules.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module is not None:
            imported_modules.add(node.module)

    assert not any(
        module == "tests" or module.startswith("tests.") for module in imported_modules
    )
    assert "conftest" not in imported_modules
    assert "PYTHONPATH=src" not in source
    assert 'env.pop("PYTHONPATH", None)' in source
    assert '"-I"' in source
    assert "TemporaryDirectory" in source
    assert "phospy.__file__" in source
    assert "load_bundled_reference_manifest" in source
    assert "hashlib.sha256" in source
    assert "DifferentialAnalysisWorkflow" in source
    assert "KinaseWorkflow" in source


def test_publish_workflow_builds_once_and_publishes_uploaded_dist() -> None:
    workflow = _read(".github/workflows/publish.yml")
    build = _workflow_job_block(workflow, "build")
    verifier = _workflow_job_block(workflow, "installed-distribution-verification")
    testpypi = _workflow_job_block(workflow, "publish-to-testpypi")
    pypi = _workflow_job_block(workflow, "publish-to-pypi")

    assert "run: make release-check" in build
    assert "name: python-package-distributions" in build
    assert _has_needs(verifier, "build")
    _assert_supported_python_matrix(verifier)
    assert "python scripts/verify_installed_distributions.py" in verifier
    assert "--dist-dir dist" in verifier
    assert _has_needs(testpypi, "build")
    assert _has_needs(testpypi, "installed-distribution-verification")
    assert _has_needs(pypi, "build")
    assert _has_needs(pypi, "installed-distribution-verification")
    assert "packages-dir: dist/" in testpypi
    assert "packages-dir: dist/" in pypi
    assert "id-token: write" in testpypi
    assert "id-token: write" in pypi


def test_ci_keeps_supported_python_source_tests_and_single_build_smoke() -> None:
    workflow = _read(".github/workflows/ci.yml")
    clean_install = _workflow_job_block(workflow, "clean-constrained-install")
    default_suite = _workflow_job_block(workflow, "default-suite")
    activity_parity = _workflow_job_block(workflow, "activity-parity-gate")
    hard_parity = _workflow_job_block(workflow, "parity-tests")
    diagnostics = _workflow_job_block(workflow, "parity-diagnostics")
    performance = _workflow_job_block(workflow, "performance-contracts")
    reference_bundles = _workflow_job_block(workflow, "reference-bundles")
    fixture_integrity = _workflow_job_block(workflow, "fixture-integrity")
    release_gates = _workflow_job_block(workflow, "release-gates")
    build = _workflow_job_block(workflow, "build-distributions")
    installed = _workflow_job_block(workflow, "installed-distribution-verification")

    _assert_supported_python_matrix(clean_install)
    _assert_supported_python_matrix(default_suite)
    _assert_supported_python_matrix(activity_parity)
    _assert_supported_python_matrix(hard_parity)
    _assert_supported_python_matrix(performance)
    _assert_supported_python_matrix(release_gates)
    _assert_supported_python_matrix(installed)
    assert "python-version: '3.10'" in diagnostics
    assert "timeout-minutes: 90" in performance
    assert "make test-performance" in performance
    assert "make validate-reference-bundles" in reference_bundles
    assert "runs-on: ${{ matrix.os }}" in fixture_integrity
    assert "os: [ubuntu-latest, windows-latest]" in fixture_integrity
    assert "python-version: '3.12'" in fixture_integrity
    assert (
        "test_manifest_fixture_byte_reproducibility.py::"
        "test_manifest_governed_fixtures_use_canonical_lf_bytes_and_valid_hashes"
    ) in fixture_integrity
    assert (
        "test_large_limma_trend_fixture_manifest_hashes_match_checked_in_files"
        in fixture_integrity
    )
    assert "make test-release-gates" in release_gates
    assert '-m "parity and activity_parity"' in activity_parity
    assert '-m "parity and not parity_diagnostic"' in hard_parity
    assert "continue-on-error: true" not in hard_parity
    assert "continue-on-error: true" in diagnostics
    assert '-m "parity_diagnostic"' in diagnostics
    assert "make build" in build
    assert "Smoke-test installed wheel outside checkout" not in build
    assert _has_needs(installed, "build-distributions")
    assert "uses: actions/download-artifact@v6" in installed
    assert "python scripts/verify_installed_distributions.py" in installed
    assert "--dist-dir dist" in installed
    assert '--repo-root "$GITHUB_WORKSPACE"' in installed
    assert "--constraint constraints/ci.txt" in installed
