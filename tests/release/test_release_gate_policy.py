from __future__ import annotations

import ast
import re
import shlex
import textwrap
import tomllib
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pytest

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
    "test-contract",
    "test-parity",
    "test-performance",
    "docs-build",
    "validate-reference-bundles",
    "test-release-gates",
    "build",
    "verify-installed-distributions",
)
PROHIBITED_RELEASE_SCALE_ROWS = 50_000
PROHIBITED_RELEASE_SCALE_SAMPLES = 48
PROHIBITED_RELEASE_SCALE_SHAPE = (
    PROHIBITED_RELEASE_SCALE_ROWS,
    PROHIBITED_RELEASE_SCALE_SAMPLES,
)
REQUIRED_TEST_SOURCE_ROOTS = (
    Path("tests/unit"),
    Path("tests/integration"),
    Path("tests/parity"),
    Path("tests/workflows"),
    Path("tests/validation"),
    Path("tests/science"),
    Path("tests/architecture"),
    Path("tests/contract"),
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
RELEASE_REACHABILITY_ENTRY_TARGETS = (
    "test-performance",
    "test-release-gates",
    "release-check",
)
ROW_DIMENSION_KEYS = frozenset(
    {
        "features",
        "feature_count",
        "n_features",
        "n_rows",
        "n_sites",
        "row_count",
        "rows",
        "site_count",
        "sites",
    }
)
SAMPLE_DIMENSION_KEYS = frozenset(
    {
        "cols",
        "column_count",
        "columns",
        "kinase_count",
        "n_cols",
        "n_columns",
        "n_kinases",
        "n_samples",
        "sample_count",
        "samples",
    }
)
SHAPE_DIMENSION_KEYS = frozenset({"matrix_shape", "shape", "size"})
AUDITABLE_LOCAL_IMPORT_ROOTS = frozenset({"benchmarks", "scripts", "tests"})
COMMAND_PATH_TOKEN_ROOTS = AUDITABLE_LOCAL_IMPORT_ROOTS | frozenset({"src"})
PYTHON_SOURCE_EXCLUDED_PARTS = frozenset(
    {
        ".git",
        ".aiassistant",
        ".hypothesis",
        ".idea",
        ".mypy_cache",
        ".pytest_cache",
        ".ruff_cache",
        ".tox",
        ".venv",
        "__pycache__",
        "build",
        "dist",
        "site",
    }
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


def _release_reachable_workload_issues(root: Path = ROOT) -> list[str]:
    paths = _release_reachable_python_sources(root)
    project = _StaticPythonProject(root)
    paths = _python_import_closure(project, paths)
    return _workload_issues_for_paths(root=root, paths=paths, project=project)


def _release_reachable_python_sources(root: Path) -> set[Path]:
    makefile_path = root / "Makefile"
    make_targets: dict[str, _MakeTarget] = {}
    make_variables: dict[str, str] = {}
    if makefile_path.is_file():
        makefile = makefile_path.read_text(encoding="utf-8").replace("\r\n", "\n")
        make_targets = _parse_make_targets(makefile)
        make_variables = _parse_make_variables(makefile)

    paths: set[Path] = set()
    pending_targets = [
        target
        for target in RELEASE_REACHABILITY_ENTRY_TARGETS
        if target in make_targets
    ]
    seen_targets: set[str] = set()

    workflow_command_blocks = _workflow_command_blocks(root)
    for command_block in workflow_command_blocks:
        command_paths, invoked_targets = _python_sources_and_make_targets_from_commands(
            command_block,
            root=root,
            make_variables=make_variables,
        )
        paths.update(command_paths)
        pending_targets.extend(
            target for target in invoked_targets if target in make_targets
        )

    while pending_targets:
        target_name = pending_targets.pop()
        if target_name in seen_targets:
            continue
        seen_targets.add(target_name)
        target = make_targets.get(target_name)
        if target is None:
            continue
        pending_targets.extend(
            dependency
            for dependency in target.dependencies
            if dependency in make_targets and dependency not in seen_targets
        )
        command_paths, invoked_targets = _python_sources_and_make_targets_from_commands(
            target.body,
            root=root,
            make_variables=make_variables,
        )
        paths.update(command_paths)
        pending_targets.extend(
            target for target in invoked_targets if target in make_targets
        )
    return paths


def _parse_make_targets(makefile: str) -> dict[str, _MakeTarget]:
    targets: dict[str, _MakeTarget] = {}
    lines = makefile.splitlines()
    index = 0
    while index < len(lines):
        line = lines[index]
        match = re.match(r"^([A-Za-z0-9_.%/-]+)\s*:(?P<deps>.*)$", line)
        if match is None or line.startswith("\t") or line.startswith("."):
            index += 1
            continue
        name = match.group(1)
        dependencies = tuple(
            dependency
            for dependency in match.group("deps").split()
            if dependency and not dependency.startswith("$")
        )
        body: list[str] = []
        index += 1
        while index < len(lines) and lines[index].startswith("\t"):
            body.append(lines[index][1:])
            index += 1
        targets[name] = _MakeTarget(
            name=name,
            dependencies=dependencies,
            body="\n".join(body),
        )
    return targets


def _parse_make_variables(makefile: str) -> dict[str, str]:
    variables: dict[str, str] = {}
    for line in makefile.splitlines():
        match = re.match(r"^([A-Za-z0-9_]+)\s*(?:\?|:)?=\s*(?P<value>.*)$", line)
        if match is not None:
            variables[match.group(1)] = match.group("value").strip()
    return variables


def _workflow_command_blocks(root: Path) -> tuple[str, ...]:
    workflow_root = root / ".github" / "workflows"
    if not workflow_root.is_dir():
        return ()
    blocks: list[str] = []
    for workflow_path in sorted(workflow_root.glob("*.yml")):
        blocks.extend(_extract_run_blocks(workflow_path.read_text(encoding="utf-8")))
    return tuple(blocks)


def _extract_run_blocks(workflow: str) -> tuple[str, ...]:
    lines = workflow.replace("\r\n", "\n").splitlines()
    blocks: list[str] = []
    index = 0
    while index < len(lines):
        line = lines[index]
        match = re.match(r"^(?P<indent>\s*)(?:-\s*)?run:\s*(?P<value>.*)$", line)
        if match is None:
            index += 1
            continue
        value = match.group("value").strip()
        indent = len(match.group("indent"))
        if value not in {"|", ">"}:
            blocks.append(value)
            index += 1
            continue
        block_lines: list[str] = []
        index += 1
        while index < len(lines):
            next_line = lines[index]
            if next_line.strip() and len(next_line) - len(next_line.lstrip()) <= indent:
                break
            block_lines.append(next_line.strip())
            index += 1
        blocks.append("\n".join(block_lines))
    return tuple(blocks)


def _python_sources_and_make_targets_from_commands(
    command_text: str,
    *,
    root: Path,
    make_variables: Mapping[str, str],
) -> tuple[set[Path], set[str]]:
    expanded = _expand_make_variables(command_text, make_variables)
    paths = _python_sources_from_command_text(expanded, root=root)
    targets = _make_targets_from_command_text(expanded)
    if _looks_like_pytest_command(expanded) and not paths:
        paths.update(_pytest_default_source_paths(root))
    return paths, targets


def _expand_make_variables(command_text: str, variables: Mapping[str, str]) -> str:
    expanded = command_text
    for _ in range(4):
        previous = expanded
        for name, value in variables.items():
            expanded = expanded.replace(f"$({name})", value)
        if expanded == previous:
            break
    return expanded


def _python_sources_from_command_text(command_text: str, *, root: Path) -> set[Path]:
    paths: set[Path] = set()
    for token in _command_tokens(command_text):
        cleaned = _clean_command_path_token(token)
        if not _is_command_path_candidate(cleaned):
            continue
        paths.update(_resolve_command_path_sources(cleaned, root=root))
    for match in re.finditer(
        r"(?P<path>(?:^|\s)(?:tests|scripts|benchmarks)[/\\][^\s'\";|&]+)",
        command_text,
    ):
        cleaned = _clean_command_path_token(match.group("path").strip())
        paths.update(_resolve_command_path_sources(cleaned, root=root))
    return paths


def _command_tokens(command_text: str) -> tuple[str, ...]:
    tokens: list[str] = []
    for line in command_text.splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        try:
            line_tokens = shlex.split(stripped, posix=True)
        except ValueError:
            line_tokens = stripped.split()
        tokens.extend(_without_inline_python_command_arguments(line_tokens))
    return tuple(tokens)


def _without_inline_python_command_arguments(tokens: Iterable[str]) -> tuple[str, ...]:
    filtered: list[str] = []
    token_list = tuple(tokens)
    index = 0
    python_command_seen = False
    while index < len(token_list):
        token = token_list[index]
        if token in {"&&", "||", ";", "|", "|&"}:
            filtered.append(token)
            python_command_seen = False
            index += 1
            continue
        if python_command_seen and token in {"-c", "--command"}:
            filtered.append(token)
            index += 2
            while index < len(token_list) and token_list[index] not in {
                "&&",
                "||",
                ";",
                "|",
                "|&",
            }:
                index += 1
            python_command_seen = False
            continue
        filtered.append(token)
        if _is_python_command_token(token):
            python_command_seen = True
        index += 1
    return tuple(filtered)


def _is_python_command_token(token: str) -> bool:
    cleaned = token.strip().strip("'\"").replace("\\", "/")
    if cleaned in {"$PYTHON", "${PYTHON}", "$(PYTHON)", "python", "python.exe", "py"}:
        return True
    return re.search(r"(?:^|/)python(?:\d+(?:\.\d+)?)?(?:\.exe)?$", cleaned) is not None


def _clean_command_path_token(token: str) -> str:
    cleaned = token.strip().strip("'\"")
    if "::" in cleaned:
        cleaned = cleaned.split("::", maxsplit=1)[0]
    cleaned = cleaned.replace("\\", "/")
    cleaned = cleaned.removeprefix("./")
    return cleaned


def _is_command_path_candidate(token: str) -> bool:
    if not token or token in {".", "./"} or token.startswith("-"):
        return False
    if "\x00" in token or re.search(r"[\r\n;&|`$<>{}()[\]*?]", token):
        return False
    if token.startswith(("http://", "https://", "git@")):
        return False
    if token.startswith("../"):
        return True
    if token.endswith(".py"):
        return True
    first_part = token.split("/", maxsplit=1)[0]
    return first_part in COMMAND_PATH_TOKEN_ROOTS


def _resolve_command_path_sources(token: str, *, root: Path) -> set[Path]:
    if not token or token in {".", "./"} or token.startswith("-"):
        return set()
    path = root / token
    try:
        resolved_root = root.resolve()
        resolved_path = path.resolve()
        resolved_path.relative_to(resolved_root)
    except (OSError, ValueError):
        return set()
    try:
        path_is_file = path.is_file()
    except OSError:
        return set()
    if path_is_file and path.suffix == ".py":
        return {resolved_path}
    try:
        path_is_dir = path.is_dir()
    except OSError:
        return set()
    if path_is_dir:
        sources: set[Path] = set()
        try:
            candidate_paths = tuple(path.rglob("*.py"))
        except OSError:
            return set()
        for item in candidate_paths:
            try:
                if item.is_file() and _is_source_path_under_root(item, root):
                    sources.add(item.resolve())
            except OSError:
                continue
        return sources
    return set()


def _is_source_path_under_root(path: Path, root: Path) -> bool:
    try:
        relative = path.resolve().relative_to(root.resolve())
    except ValueError:
        return False
    return not any(part in PYTHON_SOURCE_EXCLUDED_PARTS for part in relative.parts)


def _make_targets_from_command_text(command_text: str) -> set[str]:
    targets: set[str] = set()
    for line in command_text.splitlines():
        match = re.search(
            r"(?:^|\s)(?:make|\$\(MAKE\))\s+(?P<target>[A-Za-z0-9_.%/-]+)",
            line,
        )
        if match is not None:
            targets.add(match.group("target"))
    return targets


def _looks_like_pytest_command(command_text: str) -> bool:
    return (
        re.search(r"(?:^|\s)(?:pytest|python\s+-m\s+pytest|\$\(PYTEST\))", command_text)
        is not None
    )


def _pytest_default_source_paths(root: Path) -> set[Path]:
    pyproject_path = root / "pyproject.toml"
    if not pyproject_path.is_file():
        return set()
    with pyproject_path.open("rb") as handle:
        pyproject = tomllib.load(handle)
    configured = (
        pyproject.get("tool", {})
        .get("pytest", {})
        .get("ini_options", {})
        .get("testpaths", [])
    )
    paths: set[Path] = set()
    for relative in configured:
        paths.update(_resolve_command_path_sources(str(relative), root=root))
    return paths


def _assert_supported_python_matrix(job_block: str) -> None:
    assert "python-version: ['3.11', '3.12']" in job_block


_UNKNOWN = object()


@dataclass(frozen=True)
class _ModuleAnalysis:
    path: Path
    values: Mapping[str, object]
    module_aliases: Mapping[str, Mapping[str, object]]
    function_returns: Mapping[str, object]
    imported_paths: tuple[Path, ...]


@dataclass(frozen=True)
class _MakeTarget:
    name: str
    dependencies: tuple[str, ...]
    body: str


@dataclass(frozen=True)
class _StaticEvalContext:
    values: Mapping[str, object]
    module_aliases: Mapping[str, Mapping[str, object]]
    function_returns: Mapping[str, object]


class _StaticPythonProject:
    def __init__(self, root: Path) -> None:
        self.root = root.resolve()
        self.module_index = _build_python_module_index(self.root)
        self._analysis_by_path: dict[Path, _ModuleAnalysis] = {}
        self._active_paths: set[Path] = set()

    def analysis_for_path(self, path: Path) -> _ModuleAnalysis:
        resolved = path.resolve()
        cached = self._analysis_by_path.get(resolved)
        if cached is not None:
            return cached
        if resolved in self._active_paths:
            return _ModuleAnalysis(
                path=resolved,
                values={},
                module_aliases={},
                function_returns={},
                imported_paths=(),
            )
        self._active_paths.add(resolved)
        try:
            analysis = self._analyze_path(resolved)
            self._analysis_by_path[resolved] = analysis
            return analysis
        finally:
            self._active_paths.remove(resolved)

    def analysis_for_module(self, module_name: str) -> _ModuleAnalysis | None:
        path = self.module_index.get(module_name)
        if path is None:
            return None
        return self.analysis_for_path(path)

    def _analyze_path(self, path: Path) -> _ModuleAnalysis:
        try:
            tree = _parse_python_source(path)
        except SyntaxError:
            return _ModuleAnalysis(
                path=path,
                values={},
                module_aliases={},
                function_returns={},
                imported_paths=(),
            )
        module_name = _module_name_for_path(path, self.root)
        imported_values, module_aliases, imported_paths = _collect_import_bindings(
            tree,
            current_module=module_name,
            project=self,
        )
        values: dict[str, object] = dict(imported_values)
        class_defaults: dict[str, dict[str, object]] = {}
        for node in tree.body:
            if isinstance(node, ast.ClassDef):
                defaults = _class_default_mapping(node, values, module_aliases, {})
                if defaults:
                    values[node.name] = defaults
                    class_defaults[node.name] = defaults

        context = _StaticEvalContext(
            values=values,
            module_aliases=module_aliases,
            function_returns={},
        )
        for node in tree.body:
            if isinstance(node, ast.Assign):
                _bind_assignment_targets(
                    values, node.targets, _resolve_value(node.value, context)
                )
                context = _StaticEvalContext(
                    values=values,
                    module_aliases=module_aliases,
                    function_returns={},
                )
            elif isinstance(node, ast.AnnAssign):
                _bind_assignment_target(
                    values, node.target, _resolve_value(node.value, context)
                )
                context = _StaticEvalContext(
                    values=values,
                    module_aliases=module_aliases,
                    function_returns={},
                )

        function_returns: dict[str, object] = {}
        for node in tree.body:
            if isinstance(node, ast.FunctionDef):
                returned = _static_function_return_value(
                    node,
                    _StaticEvalContext(
                        values=values,
                        module_aliases=module_aliases,
                        function_returns=function_returns,
                    ),
                )
                if returned is not _UNKNOWN:
                    function_returns[node.name] = returned

        values.update(function_returns)
        values.update(class_defaults)
        return _ModuleAnalysis(
            path=path,
            values=values,
            module_aliases=module_aliases,
            function_returns=function_returns,
            imported_paths=tuple(sorted(imported_paths)),
        )


def _iter_all_python_sources(root: Path) -> tuple[Path, ...]:
    paths: list[Path] = []
    for path in root.rglob("*.py"):
        try:
            relative = path.relative_to(root)
        except ValueError:
            continue
        if any(part in PYTHON_SOURCE_EXCLUDED_PARTS for part in relative.parts):
            continue
        if path.is_file():
            paths.append(path.resolve())
    return tuple(sorted(paths))


def _build_python_module_index(root: Path) -> dict[str, Path]:
    module_index: dict[str, Path] = {}
    for path in _iter_all_python_sources(root):
        module_name = _module_name_for_path(path, root)
        if module_name:
            module_index.setdefault(module_name, path)
            module_index.setdefault(path.stem, path)
    return module_index


def _module_name_for_path(path: Path, root: Path) -> str:
    try:
        relative = path.resolve().relative_to(root.resolve())
    except ValueError:
        return path.stem
    parts = list(relative.with_suffix("").parts)
    if parts and parts[-1] == "__init__":
        parts = parts[:-1]
    return ".".join(parts)


def _collect_import_bindings(
    tree: ast.Module,
    *,
    current_module: str,
    project: _StaticPythonProject,
) -> tuple[dict[str, object], dict[str, Mapping[str, object]], set[Path]]:
    values: dict[str, object] = {}
    module_aliases: dict[str, Mapping[str, object]] = {}
    imported_paths: set[Path] = set()
    for node in tree.body:
        if isinstance(node, ast.ImportFrom):
            module_name = _resolve_imported_module_name(
                current_module=current_module,
                module=node.module,
                level=node.level,
            )
            if module_name is None:
                continue
            analysis = project.analysis_for_module(module_name)
            if analysis is None:
                continue
            if _is_auditable_local_import(analysis.path, project.root):
                imported_paths.add(analysis.path)
            for alias in node.names:
                if alias.name == "*":
                    values.update(analysis.values)
                    continue
                value = analysis.values.get(alias.name, _UNKNOWN)
                if value is not _UNKNOWN:
                    values[alias.asname or alias.name] = value
        elif isinstance(node, ast.Import):
            for alias in node.names:
                analysis = project.analysis_for_module(alias.name)
                if analysis is None:
                    continue
                if _is_auditable_local_import(analysis.path, project.root):
                    imported_paths.add(analysis.path)
                bound_name = alias.asname or alias.name.rsplit(".", maxsplit=1)[-1]
                module_aliases[bound_name] = analysis.values
    return values, module_aliases, imported_paths


def _resolve_imported_module_name(
    *,
    current_module: str,
    module: str | None,
    level: int,
) -> str | None:
    if level == 0:
        return module
    current_parts = current_module.split(".")
    package_parts = current_parts[:-1]
    if level > len(package_parts) + 1:
        return None
    base_parts = package_parts[: len(package_parts) - level + 1]
    if module:
        base_parts.extend(module.split("."))
    return ".".join(part for part in base_parts if part)


def _is_auditable_local_import(path: Path, root: Path) -> bool:
    try:
        relative = path.resolve().relative_to(root.resolve())
    except ValueError:
        return False
    if any(part in PYTHON_SOURCE_EXCLUDED_PARTS for part in relative.parts):
        return False
    if not relative.parts:
        return False
    return relative.parts[0] in AUDITABLE_LOCAL_IMPORT_ROOTS or len(relative.parts) == 1


def _iter_required_test_python_sources(
    root: Path = ROOT,
    required_roots: Iterable[Path] = REQUIRED_TEST_SOURCE_ROOTS,
) -> tuple[Path, ...]:
    paths: set[Path] = set()
    for relative_root in required_roots:
        source_root = root / relative_root
        if source_root.is_file() and source_root.suffix == ".py":
            if (
                source_root.relative_to(root)
                not in STATIC_RELEASE_SCALE_POLICY_EXCLUSIONS
            ):
                paths.add(source_root.resolve())
            continue
        if source_root.is_dir():
            paths.update(
                path
                for path in source_root.rglob("*.py")
                if path.is_file()
                and path.relative_to(root) not in STATIC_RELEASE_SCALE_POLICY_EXCLUSIONS
                and not any(
                    part in PYTHON_SOURCE_EXCLUDED_PARTS
                    for part in path.relative_to(root).parts
                )
            )
    return tuple(sorted(paths))


def _parse_python_source(path: Path) -> ast.Module:
    return ast.parse(path.read_text(encoding="utf-8"), filename=str(path))


def _required_source_workload_issues(
    root: Path = ROOT,
    required_roots: Iterable[Path] = REQUIRED_TEST_SOURCE_ROOTS,
) -> list[str]:
    project = _StaticPythonProject(root)
    initial_paths = set(_iter_required_test_python_sources(root, required_roots))
    paths = _python_import_closure(project, initial_paths)
    return _workload_issues_for_paths(root=root, paths=paths, project=project)


def _workload_issues_for_paths(
    *,
    root: Path,
    paths: Iterable[Path],
    project: _StaticPythonProject | None = None,
) -> list[str]:
    resolved_project = _StaticPythonProject(root) if project is None else project
    issues: list[str] = []
    for path in sorted({item.resolve() for item in paths if item.is_file()}):
        try:
            relative = path.relative_to(root.resolve())
        except ValueError:
            continue
        if relative in STATIC_RELEASE_SCALE_POLICY_EXCLUSIONS:
            continue
        issues.extend(_source_workload_issues(path, root, resolved_project))
    return issues


def _python_import_closure(
    project: _StaticPythonProject,
    initial_paths: Iterable[Path],
) -> set[Path]:
    seen: set[Path] = set()
    pending = [path.resolve() for path in initial_paths if path.is_file()]
    while pending:
        path = pending.pop()
        if path in seen:
            continue
        seen.add(path)
        analysis = project.analysis_for_path(path)
        for imported_path in analysis.imported_paths:
            if imported_path not in seen:
                pending.append(imported_path)
    return seen


def _source_workload_issues(
    path: Path,
    root: Path,
    project: _StaticPythonProject,
) -> list[str]:
    try:
        tree = _parse_python_source(path)
    except SyntaxError as exc:
        return [
            f"{_relative_path(path, root)}:{exc.lineno or 1} cannot be parsed "
            "for release-scale workload audit"
        ]
    analysis = project.analysis_for_path(path)
    context = _StaticEvalContext(
        values=analysis.values,
        module_aliases=analysis.module_aliases,
        function_returns=analysis.function_returns,
    )
    return _scan_python_statements_for_workload(
        tree.body,
        context=context,
        path=path,
        root=root,
    )


def _relative_path(path: Path, root: Path) -> str:
    try:
        return path.resolve().relative_to(root.resolve()).as_posix()
    except ValueError:
        return path.as_posix()


def _scan_python_statements_for_workload(
    statements: Iterable[ast.stmt],
    *,
    context: _StaticEvalContext,
    path: Path,
    root: Path,
) -> list[str]:
    issues: list[str] = []
    local_values: dict[str, object] = dict(context.values)
    local_context = _StaticEvalContext(
        values=local_values,
        module_aliases=context.module_aliases,
        function_returns=context.function_returns,
    )
    for statement in statements:
        if isinstance(statement, ast.FunctionDef):
            issues.extend(
                _workload_issues_in_expression_nodes(
                    [*statement.decorator_list, *statement.args.defaults],
                    context=local_context,
                    path=path,
                    root=root,
                )
            )
            function_values = dict(local_values)
            function_values.update(_function_default_bindings(statement, local_context))
            issues.extend(
                _scan_python_statements_for_workload(
                    statement.body,
                    context=_StaticEvalContext(
                        values=function_values,
                        module_aliases=context.module_aliases,
                        function_returns=context.function_returns,
                    ),
                    path=path,
                    root=root,
                )
            )
            continue
        if isinstance(statement, ast.ClassDef):
            issues.extend(
                _workload_issues_in_expression_nodes(
                    [*statement.decorator_list, *statement.bases],
                    context=local_context,
                    path=path,
                    root=root,
                )
            )
            class_defaults = _class_default_mapping(
                statement,
                local_values,
                context.module_aliases,
                context.function_returns,
            )
            if _mapping_has_prohibited_shape(class_defaults):
                issues.append(
                    f"{_relative_path(path, root)}:{statement.lineno} defines a "
                    "configuration object with prohibited effective "
                    "50,000x48 dimensions"
                )
            issues.extend(
                _scan_python_statements_for_workload(
                    statement.body,
                    context=_StaticEvalContext(
                        values={**local_values, **class_defaults},
                        module_aliases=context.module_aliases,
                        function_returns=context.function_returns,
                    ),
                    path=path,
                    root=root,
                )
            )
            local_values[statement.name] = class_defaults
            local_context = _StaticEvalContext(
                values=local_values,
                module_aliases=context.module_aliases,
                function_returns=context.function_returns,
            )
            continue

        issues.extend(
            _workload_issues_in_expression_nodes(
                [statement],
                context=local_context,
                path=path,
                root=root,
            )
        )
        if isinstance(statement, ast.Assign):
            resolved = _resolve_value(statement.value, local_context)
            _bind_assignment_targets(local_values, statement.targets, resolved)
            local_context = _StaticEvalContext(
                values=local_values,
                module_aliases=context.module_aliases,
                function_returns=context.function_returns,
            )
        elif isinstance(statement, ast.AnnAssign):
            resolved = _resolve_value(statement.value, local_context)
            _bind_assignment_target(local_values, statement.target, resolved)
            local_context = _StaticEvalContext(
                values=local_values,
                module_aliases=context.module_aliases,
                function_returns=context.function_returns,
            )
    return issues


def _workload_issues_in_expression_nodes(
    nodes: Iterable[ast.AST | None],
    *,
    context: _StaticEvalContext,
    path: Path,
    root: Path,
) -> list[str]:
    issues: list[str] = []
    for node in nodes:
        if node is None:
            continue
        for descendant in ast.walk(node):
            if not isinstance(descendant, ast.Call):
                continue
            if _call_has_prohibited_shape(descendant, context):
                issues.append(
                    f"{_relative_path(path, root)}:{descendant.lineno} reaches "
                    "prohibited effective 50,000x48 workload dimensions"
                )
    return issues


def _class_default_mapping(
    node: ast.ClassDef,
    values: Mapping[str, object],
    module_aliases: Mapping[str, Mapping[str, object]],
    function_returns: Mapping[str, object],
) -> dict[str, object]:
    context = _StaticEvalContext(
        values=values,
        module_aliases=module_aliases,
        function_returns=function_returns,
    )
    defaults: dict[str, object] = {}
    for statement in node.body:
        if isinstance(statement, ast.Assign):
            resolved = _resolve_value(statement.value, context)
            _bind_assignment_targets(defaults, statement.targets, resolved)
        elif isinstance(statement, ast.AnnAssign):
            resolved = _resolve_value(statement.value, context)
            _bind_assignment_target(defaults, statement.target, resolved)
    return defaults


def _function_default_bindings(
    node: ast.FunctionDef,
    context: _StaticEvalContext,
) -> dict[str, object]:
    bindings: dict[str, object] = {}
    positional_args = [arg.arg for arg in node.args.posonlyargs + node.args.args]
    default_offset = len(positional_args) - len(node.args.defaults)
    for index, default in enumerate(node.args.defaults, start=default_offset):
        if index >= 0:
            value = _resolve_value(default, context)
            if value is not _UNKNOWN:
                bindings[positional_args[index]] = value
    for arg, default in zip(
        node.args.kwonlyargs,
        node.args.kw_defaults,
        strict=True,
    ):
        value = _resolve_value(default, context)
        if value is not _UNKNOWN:
            bindings[arg.arg] = value
    return bindings


def _static_function_return_value(
    node: ast.FunctionDef,
    context: _StaticEvalContext,
) -> object:
    local_values = dict(context.values)
    local_values.update(_function_default_bindings(node, context))
    local_context = _StaticEvalContext(
        values=local_values,
        module_aliases=context.module_aliases,
        function_returns=context.function_returns,
    )
    returned: object = _UNKNOWN
    for statement in node.body:
        if isinstance(statement, ast.Assign):
            resolved = _resolve_value(statement.value, local_context)
            _bind_assignment_targets(local_values, statement.targets, resolved)
            local_context = _StaticEvalContext(
                values=local_values,
                module_aliases=context.module_aliases,
                function_returns=context.function_returns,
            )
        elif isinstance(statement, ast.AnnAssign):
            resolved = _resolve_value(statement.value, local_context)
            _bind_assignment_target(local_values, statement.target, resolved)
            local_context = _StaticEvalContext(
                values=local_values,
                module_aliases=context.module_aliases,
                function_returns=context.function_returns,
            )
        elif isinstance(statement, ast.Return):
            value = _resolve_value(statement.value, local_context)
            if value is _UNKNOWN:
                return _UNKNOWN
            if returned is _UNKNOWN:
                returned = value
            elif returned != value:
                return _UNKNOWN
    return returned


def _bind_assignment_targets(
    values: dict[str, object],
    targets: Iterable[ast.expr],
    resolved: object,
) -> None:
    for target in targets:
        _bind_assignment_target(values, target, resolved)


def _bind_assignment_target(
    values: dict[str, object],
    target: ast.expr,
    resolved: object,
) -> None:
    if resolved is _UNKNOWN:
        return
    if isinstance(target, ast.Name):
        values[target.id] = resolved
    elif isinstance(target, ast.Tuple | ast.List) and isinstance(
        resolved, tuple | list
    ):
        for element, item in zip(target.elts, resolved, strict=False):
            _bind_assignment_target(values, element, item)


def _call_has_prohibited_shape(
    node: ast.Call,
    context: _StaticEvalContext,
) -> bool:
    positional_values = _resolved_positional_values(node, context)
    if _sequence_has_prohibited_shape(positional_values):
        return True
    if any(_value_has_prohibited_shape(value) for value in positional_values):
        return True
    keyword_values = _resolved_keyword_mapping(node, context)
    if _mapping_has_prohibited_shape(keyword_values):
        return True
    return any(_value_has_prohibited_shape(value) for value in keyword_values.values())


def _resolved_positional_values(
    node: ast.Call,
    context: _StaticEvalContext,
) -> list[object]:
    values: list[object] = []
    for argument in node.args:
        if isinstance(argument, ast.Starred):
            expanded = _resolve_value(argument.value, context)
            if isinstance(expanded, tuple | list):
                values.extend(expanded)
            else:
                values.append(expanded)
            continue
        values.append(_resolve_value(argument, context))
    return values


def _resolved_keyword_mapping(
    node: ast.Call,
    context: _StaticEvalContext,
) -> dict[str, object]:
    values: dict[str, object] = {}
    for keyword in node.keywords:
        value = _resolve_value(keyword.value, context)
        if keyword.arg is None:
            if isinstance(value, Mapping):
                values.update({str(key): item for key, item in value.items()})
            continue
        values[str(keyword.arg)] = value
    return values


def _value_has_prohibited_shape(value: object) -> bool:
    if _is_int_sequence(value) and tuple(value) == PROHIBITED_RELEASE_SCALE_SHAPE:
        return True
    if isinstance(value, Mapping):
        return _mapping_has_prohibited_shape(value)
    if isinstance(value, tuple | list):
        return any(_value_has_prohibited_shape(item) for item in value)
    return False


def _sequence_has_prohibited_shape(values: Iterable[object]) -> bool:
    int_values = [value for value in values if _is_plain_int(value)]
    return any(
        tuple(int_values[index : index + 2]) == PROHIBITED_RELEASE_SCALE_SHAPE
        for index in range(max(len(int_values) - 1, 0))
    )


def _mapping_has_prohibited_shape(value: Mapping[object, object]) -> bool:
    normalized = {str(key): item for key, item in value.items()}
    row_value = _first_dimension_value(normalized, ROW_DIMENSION_KEYS)
    sample_value = _first_dimension_value(normalized, SAMPLE_DIMENSION_KEYS)
    if (row_value, sample_value) == PROHIBITED_RELEASE_SCALE_SHAPE:
        return True
    for key in SHAPE_DIMENSION_KEYS:
        if _value_has_prohibited_shape(normalized.get(key, _UNKNOWN)):
            return True
    index_length = _length_value(normalized.get("index", _UNKNOWN))
    columns_length = _length_value(normalized.get("columns", _UNKNOWN))
    return (index_length, columns_length) == PROHIBITED_RELEASE_SCALE_SHAPE


def _first_dimension_value(
    values: Mapping[str, object],
    aliases: frozenset[str],
) -> int | None:
    for key, value in values.items():
        if _normalize_dimension_key(key) in aliases:
            resolved = _plain_int_or_length(value)
            if resolved is not None:
                return resolved
    return None


def _normalize_dimension_key(key: str) -> str:
    return key.strip().lower().replace("-", "_")


def _plain_int_or_length(value: object) -> int | None:
    if _is_plain_int(value):
        return int(value)
    return _length_value(value)


def _length_value(value: object) -> int | None:
    if isinstance(value, Mapping) and set(value) == {"__length__"}:
        length = value["__length__"]
        if _is_plain_int(length):
            return int(length)
    return None


def _is_int_sequence(value: object) -> bool:
    return (
        isinstance(value, tuple | list)
        and len(value) == 2
        and all(_is_plain_int(item) for item in value)
    )


def _literal_int(node: ast.AST | None) -> int | None:
    if (
        isinstance(node, ast.Constant)
        and isinstance(node.value, int)
        and not isinstance(node.value, bool)
    ):
        return int(node.value)
    if (
        isinstance(node, ast.UnaryOp)
        and isinstance(node.op, ast.USub)
        and isinstance(node.operand, ast.Constant)
        and isinstance(node.operand.value, int)
    ):
        return -int(node.operand.value)
    return None


def _is_plain_int(value: object) -> bool:
    return isinstance(value, int) and not isinstance(value, bool)


def _resolve_value(
    node: ast.AST | None,
    context: _StaticEvalContext,
) -> object:
    if node is None:
        return _UNKNOWN
    literal = _literal_int(node)
    if literal is not None:
        return literal
    if isinstance(node, ast.Constant):
        return node.value
    if isinstance(node, ast.Name):
        return context.values.get(node.id, _UNKNOWN)
    if isinstance(node, ast.Attribute):
        owner = _resolve_value(node.value, context)
        if isinstance(owner, Mapping) and node.attr in owner:
            return owner[node.attr]
        if isinstance(node.value, ast.Name):
            module_values = context.module_aliases.get(node.value.id)
            if module_values is not None:
                return module_values.get(node.attr, _UNKNOWN)
        return _UNKNOWN
    if isinstance(node, ast.Subscript):
        owner = _resolve_value(node.value, context)
        key = _resolve_subscript_key(node.slice, context)
        if isinstance(owner, Mapping) and key in owner:
            return owner[key]
        return _UNKNOWN
    if isinstance(node, ast.BinOp):
        left = _resolve_value(node.left, context)
        right = _resolve_value(node.right, context)
        if _is_plain_int(left) and _is_plain_int(right):
            return _apply_int_binop(int(left), int(right), node.op)
        return _UNKNOWN
    if isinstance(node, ast.Tuple):
        return tuple(_resolve_value(element, context) for element in node.elts)
    if isinstance(node, ast.List):
        return [_resolve_value(element, context) for element in node.elts]
    if isinstance(node, ast.Dict):
        resolved: dict[str, object] = {}
        for key_node, value_node in zip(node.keys, node.values, strict=True):
            key = _resolve_value(key_node, context)
            if isinstance(key, str):
                resolved[key] = _resolve_value(value_node, context)
        return resolved
    if isinstance(node, ast.Call):
        return _resolve_call_value(node, context)
    return _UNKNOWN


def _resolve_subscript_key(
    node: ast.AST,
    context: _StaticEvalContext,
) -> object:
    if isinstance(node, ast.Constant):
        return node.value
    return _resolve_value(node, context)


def _apply_int_binop(left: int, right: int, op: ast.operator) -> int | object:
    if isinstance(op, ast.Add):
        return left + right
    if isinstance(op, ast.Sub):
        return left - right
    if isinstance(op, ast.Mult):
        return left * right
    if isinstance(op, ast.FloorDiv) and right != 0:
        return left // right
    if isinstance(op, ast.Mod) and right != 0:
        return left % right
    return _UNKNOWN


def _resolve_call_value(node: ast.Call, context: _StaticEvalContext) -> object:
    function_name = _call_function_name(node)
    positional_values = _resolved_positional_values(node, context)
    keyword_values = _resolved_keyword_mapping(node, context)
    callable_value = _resolve_callable_value(node.func, context)
    if function_name == "int" and len(positional_values) == 1:
        value = positional_values[0]
        if _is_plain_int(value):
            return int(value)
    if function_name == "tuple" and len(positional_values) == 1:
        value = positional_values[0]
        if isinstance(value, tuple | list):
            return tuple(value)
    if function_name == "list" and len(positional_values) == 1:
        value = positional_values[0]
        if isinstance(value, tuple | list):
            return list(value)
    if function_name == "dict":
        resolved: dict[str, object] = {}
        for value in positional_values:
            if isinstance(value, Mapping):
                resolved.update({str(key): item for key, item in value.items()})
        resolved.update(keyword_values)
        return resolved
    if function_name == "range" and positional_values:
        length = _range_length(positional_values)
        if length is not None:
            return {"__length__": length}
    if function_name in context.function_returns and not positional_values:
        if not keyword_values:
            return context.function_returns[function_name]
    if (
        callable_value is not _UNKNOWN
        and not positional_values
        and not keyword_values
        and isinstance(callable_value, int | tuple | list | Mapping)
    ):
        return callable_value
    constructor_defaults = (
        callable_value
        if callable_value is not _UNKNOWN
        else context.values.get(function_name, _UNKNOWN)
    )
    if isinstance(constructor_defaults, Mapping):
        return {**constructor_defaults, **keyword_values}
    if keyword_values:
        return keyword_values
    return _UNKNOWN


def _resolve_callable_value(node: ast.expr, context: _StaticEvalContext) -> object:
    if isinstance(node, ast.Name):
        return context.values.get(node.id, _UNKNOWN)
    if isinstance(node, ast.Attribute):
        owner = _resolve_value(node.value, context)
        if isinstance(owner, Mapping) and node.attr in owner:
            return owner[node.attr]
        if isinstance(node.value, ast.Name):
            module_values = context.module_aliases.get(node.value.id)
            if module_values is not None:
                return module_values.get(node.attr, _UNKNOWN)
    return _UNKNOWN


def _call_function_name(node: ast.Call) -> str:
    if isinstance(node.func, ast.Name):
        return node.func.id
    if isinstance(node.func, ast.Attribute):
        return node.func.attr
    return ""


def _range_length(values: list[object]) -> int | None:
    if len(values) == 1 and _is_plain_int(values[0]):
        return max(int(values[0]), 0)
    if len(values) >= 2 and _is_plain_int(values[0]) and _is_plain_int(values[1]):
        start = int(values[0])
        stop = int(values[1])
        step = int(values[2]) if len(values) >= 3 and _is_plain_int(values[2]) else 1
        if step == 0:
            return None
        if step > 0:
            return max((stop - start + step - 1) // step, 0)
        return max((start - stop + abs(step) - 1) // abs(step), 0)
    return None


def _dimension_prefix(name: str, suffixes: tuple[str, ...]) -> str | None:
    for suffix in suffixes:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return None


def _fixture_constant_dimension_issues(path: Path, tree: ast.Module) -> list[str]:
    del tree
    root = _audit_root_for_path(path)
    project = _StaticPythonProject(root)
    return _source_workload_issues(path, root, project)


def _audit_root_for_path(path: Path) -> Path:
    resolved = path.resolve()
    try:
        resolved.relative_to(ROOT.resolve())
    except ValueError:
        pass
    else:
        return ROOT
    for parent in resolved.parents:
        if (parent / "tests").is_dir() or (parent / "Makefile").is_file():
            return parent
    return resolved.parent


def _assigned_int_constants(tree: ast.Module) -> dict[str, int]:
    constants: dict[str, int] = {}
    context = _StaticEvalContext(values={}, module_aliases={}, function_returns={})
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign):
            value = _resolve_value(node.value, context)
            if not _is_plain_int(value):
                continue
            for target in node.targets:
                if isinstance(target, ast.Name):
                    constants[target.id] = int(value)
        elif isinstance(node, ast.AnnAssign):
            value = _resolve_value(node.value, context)
            if _is_plain_int(value) and isinstance(node.target, ast.Name):
                constants[node.target.id] = int(value)
    return constants


def _fixture_shape_like_constant_dimension_issues(
    path: Path,
    tree: ast.Module,
) -> list[str]:
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
        f"{_relative_path(path, _audit_root_for_path(path))} defines prohibited 50,000x48 "
        f"fixture constants with prefix {prefix!r}"
        for prefix in sorted(row_prefixes & sample_prefixes)
    ]


def _fixture_construction_dimension_issues(path: Path, tree: ast.Module) -> list[str]:
    del tree
    root = _audit_root_for_path(path)
    project = _StaticPythonProject(root)
    return _source_workload_issues(path, root, project)


def _tuple_ints(
    node: ast.AST,
    constants: dict[str, int],
) -> tuple[int, ...] | None:
    if not isinstance(node, ast.Tuple | ast.List):
        return None
    context = _StaticEvalContext(
        values=constants,
        module_aliases={},
        function_returns={},
    )
    values: list[int] = []
    for element in node.elts:
        value = _resolve_value(element, context)
        if not _is_plain_int(value):
            return None
        values.append(int(value))
    return tuple(values)


def _write_policy_fixture_file(root: Path, relative_path: str, source: str) -> Path:
    path = root / relative_path
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(textwrap.dedent(source).strip() + "\n", encoding="utf-8")
    return path


def _audit_tmp_required_sources(root: Path) -> list[str]:
    return _required_source_workload_issues(
        root=root,
        required_roots=(Path("tests/performance"), Path("tests/support")),
    )


def _write_minimal_release_makefile(root: Path, body: str) -> None:
    _write_policy_fixture_file(
        root,
        "Makefile",
        body,
    )


def _inline_python_environment_report_command() -> str:
    long_report_label = "release_policy_inline_python_environment_report_" + ("x" * 320)
    return (
        'python -c "import sys, importlib.metadata as md; '
        "print('python', sys.version); "
        "print('executable', sys.executable); "
        "print('pyright', md.version('pyright')); "
        f"print('label', '{long_report_label}')\""
    )


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
        "contract",
        "slow",
        "release_gate",
        "reproducibility",
        "golden",
        "activity_parity",
        "parity_diagnostic",
    }

    makefile = _read("Makefile")
    assert re.search(r"(?ms)^\.PHONY:.*test-contract", makefile) is not None
    assert re.search(r"(?ms)^\.PHONY:.*test-release-gates", makefile) is not None
    assert re.search(r"(?ms)^\.PHONY:.*docs-build", makefile) is not None
    assert "make test-contract" in _make_target_body("help")
    assert "make test-release-gates" in _make_target_body("help")
    assert "make docs-build" in _make_target_body("help")

    assert '$(PYTEST) -m "not parity"' in _make_target_body("test-unit")
    contract = _make_target_body("test-contract")
    assert '$(MKDIR_P) "$(PYTEST_REPORT_DIR)"' in contract
    assert (
        "$(PYTEST) -o addopts= tests/contract --junitxml "
        '"$(PYTEST_REPORT_DIR)/contract.xml"'
    ) in contract
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
    assert "$(MKDOCS) build --strict" in _make_target_body("docs-build")
    assert "strict: true" in _read("mkdocs.yml")
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


@pytest.mark.parametrize(
    ("case_id", "source"),
    (
        (
            "generic_constants_positional",
            """
            ROWS = 50000
            COLS = 48

            def test_required_workload():
                build_matrix(ROWS, COLS)
            """,
        ),
        (
            "computed_dimensions_keyword_aliases",
            """
            ROWS = 25000 * 2
            COLS = 6 * 8

            def test_required_workload():
                build_matrix(rows=ROWS, cols=COLS)
            """,
        ),
        (
            "dictionary_kwargs",
            """
            CFG = {"n_sites": 50000, "n_samples": 48}

            def test_required_workload():
                build_matrix(**CFG)
            """,
        ),
        (
            "configuration_object",
            """
            class RequiredScale:
                n_sites: int = 50000
                n_samples: int = 48

            def test_required_workload():
                build_matrix(config=RequiredScale())
            """,
        ),
        (
            "returned_shape_starargs",
            """
            def required_shape():
                return (25000 * 2, 6 * 8)

            def test_required_workload():
                build_matrix(*required_shape())
            """,
        ),
        (
            "other_required_marker_parametrize",
            """
            import pytest

            pytestmark = pytest.mark.slow

            @pytest.mark.parametrize(("rows", "cols"), [(50000, 48)])
            def test_required_workload(rows, cols):
                assert (rows, cols)
            """,
        ),
        (
            "dataframe_index_column_lengths",
            """
            def test_required_workload():
                build_frame(index=range(50000), columns=range(48))
            """,
        ),
    ),
    ids=lambda item: item,
)
def test_required_source_workload_audit_detects_demonstrated_bypass_patterns(
    tmp_path: Path,
    case_id: str,
    source: str,
) -> None:
    _write_policy_fixture_file(
        tmp_path,
        f"tests/performance/test_{case_id}.py",
        source,
    )

    issues = _audit_tmp_required_sources(tmp_path)

    assert issues, case_id


def test_required_source_workload_audit_resolves_imported_dimensions(
    tmp_path: Path,
) -> None:
    _write_policy_fixture_file(
        tmp_path,
        "tests/support/release_dimensions.py",
        """
        ROWS = 25000 * 2
        COLS = 6 * 8
        """,
    )
    _write_policy_fixture_file(
        tmp_path,
        "tests/performance/test_imported_dimensions.py",
        """
        from tests.support.release_dimensions import COLS, ROWS
        import tests.support.release_dimensions as dims

        def test_required_workload():
            build_matrix(n_sites=ROWS, n_samples=COLS)
            build_matrix(dims.ROWS, dims.COLS)
        """,
    )

    issues = _audit_tmp_required_sources(tmp_path)

    assert issues


def test_required_source_workload_audit_detects_forwarded_helper_dimensions(
    tmp_path: Path,
) -> None:
    _write_policy_fixture_file(
        tmp_path,
        "tests/support/release_helper.py",
        """
        ROWS = 50000
        COLS = 48

        def build_required_matrix(builder):
            return builder(n_sites=ROWS, n_samples=COLS)
        """,
    )
    _write_policy_fixture_file(
        tmp_path,
        "tests/performance/test_forwarded_helper.py",
        """
        from tests.support.release_helper import build_required_matrix

        def test_required_workload():
            build_required_matrix(build_matrix)
        """,
    )

    issues = _audit_tmp_required_sources(tmp_path)

    assert issues


def test_required_source_workload_audit_keeps_bounded_50k_contracts_allowed(
    tmp_path: Path,
) -> None:
    _write_policy_fixture_file(
        tmp_path,
        "tests/performance/test_allowed_bounded_contracts.py",
        """
        def test_allowed_bounded_workloads():
            build_matrix(n_sites=50000, n_samples=12)
            build_matrix(n_sites=50000, n_samples=24)
            build_matrix((50000, 12))
            build_matrix((50000, 24))
        """,
    )

    assert _audit_tmp_required_sources(tmp_path) == []


def test_make_reachability_audit_detects_renamed_benchmark_target_and_script(
    tmp_path: Path,
) -> None:
    _write_policy_fixture_file(
        tmp_path,
        "benchmarks/renamed_scale_probe.py",
        """
        ROWS = 25000 * 2
        COLS = 6 * 8

        def main():
            build_matrix(ROWS, COLS)
        """,
    )
    _write_minimal_release_makefile(
        tmp_path,
        """
        test-performance: renamed-scale-target
        test-release-gates:
        release-check: test-performance test-release-gates

        renamed-scale-target:
        \t$(PYTHON) benchmarks/renamed_scale_probe.py
        """,
    )

    issues = _release_reachable_workload_issues(tmp_path)

    assert any("benchmarks/renamed_scale_probe.py" in issue for issue in issues)


def test_make_reachability_audit_detects_required_transitive_invocation(
    tmp_path: Path,
) -> None:
    _write_policy_fixture_file(
        tmp_path,
        "scripts/hidden_required_scale.py",
        """
        CFG = {"n_sites": 50000, "n_samples": 48}

        def main():
            build_matrix(**CFG)
        """,
    )
    _write_minimal_release_makefile(
        tmp_path,
        """
        test-performance: intermediate-target
        test-release-gates:
        release-check: test-performance test-release-gates

        intermediate-target: hidden-target

        hidden-target:
        \t$(PYTHON) scripts/hidden_required_scale.py
        """,
    )

    issues = _release_reachable_workload_issues(tmp_path)

    assert any("scripts/hidden_required_scale.py" in issue for issue in issues)


def test_workflow_reachability_audit_detects_direct_required_script(
    tmp_path: Path,
) -> None:
    _write_policy_fixture_file(
        tmp_path,
        "scripts/workflow_scale.py",
        """
        def main():
            build_matrix(n_sites=50000, n_samples=48)
        """,
    )
    _write_policy_fixture_file(
        tmp_path,
        ".github/workflows/ci.yml",
        """
        name: CI
        on: [push]
        jobs:
          performance:
            runs-on: ubuntu-latest
            steps:
              - run: python scripts/workflow_scale.py
        """,
    )

    issues = _release_reachable_workload_issues(tmp_path)

    assert any("scripts/workflow_scale.py" in issue for issue in issues)


def test_workflow_reachability_ignores_inline_python_c_command_tokens(
    tmp_path: Path,
) -> None:
    inline_report = _inline_python_environment_report_command()
    _write_policy_fixture_file(
        tmp_path,
        ".github/workflows/ci.yml",
        f"""
        name: CI
        on: [push]
        jobs:
          type-check:
            runs-on: ubuntu-latest
            steps:
              - name: Show type-check environment
                run: |
                  source .venv/bin/activate
                  {inline_report}
        """,
    )

    issues = _release_reachable_workload_issues(tmp_path)

    assert issues == []


def test_workflow_reachability_still_detects_script_path_near_inline_python_command(
    tmp_path: Path,
) -> None:
    inline_report = _inline_python_environment_report_command()
    _write_policy_fixture_file(
        tmp_path,
        "scripts/workflow_scale_near_inline_python.py",
        """
        def main():
            build_matrix(n_sites=50000, n_samples=48)
        """,
    )
    _write_policy_fixture_file(
        tmp_path,
        ".github/workflows/ci.yml",
        f"""
        name: CI
        on: [push]
        jobs:
          performance:
            runs-on: ubuntu-latest
            steps:
              - name: Run scale check
                run: |
                  {inline_report}
                  python scripts/workflow_scale_near_inline_python.py
        """,
    )

    issues = _release_reachable_workload_issues(tmp_path)

    assert any(
        "scripts/workflow_scale_near_inline_python.py" in issue for issue in issues
    )


def test_workflow_reachability_audit_detects_make_target_invoked_from_workflow(
    tmp_path: Path,
) -> None:
    _write_policy_fixture_file(
        tmp_path,
        "benchmarks/workflow_renamed.py",
        """
        class Scale:
            n_sites = 50000
            n_samples = 48

        def main():
            build_matrix(config=Scale())
        """,
    )
    _write_minimal_release_makefile(
        tmp_path,
        """
        test-performance:
        test-release-gates:
        release-check: test-performance test-release-gates

        renamed-local-benchmark:
        \t$(PYTHON) benchmarks/workflow_renamed.py
        """,
    )
    _write_policy_fixture_file(
        tmp_path,
        ".github/workflows/ci.yml",
        """
        name: CI
        on: [workflow_dispatch]
        jobs:
          manual-benchmark:
            runs-on: ubuntu-latest
            steps:
              - run: make renamed-local-benchmark
        """,
    )

    issues = _release_reachable_workload_issues(tmp_path)

    assert any("benchmarks/workflow_renamed.py" in issue for issue in issues)


def test_required_test_sources_do_not_reach_50k_by_48_workload() -> None:
    issues = _required_source_workload_issues()

    assert issues == []


def test_release_commands_and_workflows_do_not_reach_50k_by_48_workload() -> None:
    issues = _release_reachable_workload_issues()

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
    assert "ticket_1_posthoc_peptide_to_site_boundary" in source
    assert "withdrawn_asserted" in source


def test_publish_workflow_builds_once_and_publishes_uploaded_dist() -> None:
    workflow = _read(".github/workflows/publish.yml")
    build = _workflow_job_block(workflow, "build")
    verifier = _workflow_job_block(workflow, "installed-distribution-verification")
    testpypi = _workflow_job_block(workflow, "publish-to-testpypi")
    pypi = _workflow_job_block(workflow, "publish-to-pypi")

    assert "run: make release-check" in build
    assert ".[dev,test,parquet,docs]" in build
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


def test_ci_testing_audit_freshness_commands_are_not_duplicated() -> None:
    testing_audit = _workflow_job_block(
        _read(".github/workflows/ci.yml"),
        "testing-audit-freshness",
    )
    expected_commands = (
        "python tools/testing/generate_test_inventory.py",
        "python tools/testing/find_validator_test_patterns.py",
        "python tools/testing/find_dataframe_ownership_tests.py",
        "python tools/testing/find_diagnostic_assertion_clusters.py",
        "python tools/testing/find_orchestration_test_candidates.py",
    )

    for command in expected_commands:
        assert testing_audit.count(command) == 1


def test_ci_keeps_supported_python_source_tests_and_single_build_smoke() -> None:
    workflow = _read(".github/workflows/ci.yml")
    unsupported_python_310 = "3." + "10"
    lint = _workflow_job_block(workflow, "lint")
    clean_install = _workflow_job_block(workflow, "clean-constrained-install")
    minimum = _workflow_job_block(workflow, "minimum-dependency-suite")
    benchmark = _workflow_job_block(workflow, "benchmark-smoke")
    testing_audit = _workflow_job_block(workflow, "testing-audit-freshness")
    default_suite = _workflow_job_block(workflow, "default-suite")
    contract = _workflow_job_block(workflow, "public-consumer-contracts")
    documentation = _workflow_job_block(workflow, "documentation")
    activity_parity = _workflow_job_block(workflow, "activity-parity-gate")
    hard_parity = _workflow_job_block(workflow, "parity-tests")
    diagnostics = _workflow_job_block(workflow, "parity-diagnostics")
    performance = _workflow_job_block(workflow, "performance-contracts")
    reference_bundles = _workflow_job_block(workflow, "reference-bundles")
    fixture_integrity = _workflow_job_block(workflow, "fixture-integrity")
    release_gates = _workflow_job_block(workflow, "release-gates")
    adaptive = _workflow_job_block(workflow, "adaptive-mode-standard-install")
    build = _workflow_job_block(workflow, "build-distributions")
    installed = _workflow_job_block(workflow, "installed-distribution-verification")

    assert f"python-version: '{unsupported_python_310}'" not in workflow
    _assert_supported_python_matrix(clean_install)
    _assert_supported_python_matrix(default_suite)
    _assert_supported_python_matrix(contract)
    _assert_supported_python_matrix(activity_parity)
    _assert_supported_python_matrix(hard_parity)
    _assert_supported_python_matrix(performance)
    _assert_supported_python_matrix(release_gates)
    _assert_supported_python_matrix(installed)
    for lowest_supported_job in (
        lint,
        minimum,
        benchmark,
        testing_audit,
        documentation,
        reference_bundles,
        adaptive,
        diagnostics,
    ):
        assert "python-version: '3.11'" in lowest_supported_job
    assert "make test-contract" in contract
    assert "public-consumer-contracts-py${{ matrix.python-version }}" in contract
    assert 'pip install -e ".[docs]"' in documentation
    assert "make docs-build" in documentation
    assert "$(MKDOCS) build --strict" in _make_target_body("docs-build")
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
