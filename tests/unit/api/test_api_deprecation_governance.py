from __future__ import annotations

import ast
import importlib
import os
import re
import subprocess
import sys
import tomllib
from pathlib import Path

from phospy.api._compat import compatibility_export

ROOT = Path(__file__).resolve().parents[3]
COMPATIBILITY_TEST = (
    ROOT / "tests" / "unit" / "api" / ("test_api_compatibility_deprecations.py")
)


def test_ordinary_tests_do_not_import_deprecated_api_compatibility_routes() -> None:
    violations: list[str] = []
    for path in sorted((ROOT / "tests").rglob("*.py")):
        if path == COMPATIBILITY_TEST:
            continue
        source = path.read_text(encoding="utf-8")
        violations.extend(_deprecated_api_imports(path, source))

    assert violations == []


def test_documentation_examples_use_canonical_api_import_routes() -> None:
    violations: list[str] = []
    paths = (
        ROOT / "README.md",
        ROOT / "CHANGELOG.md",
        *sorted((ROOT / "docs").rglob("*.md")),
        *sorted((ROOT / "examples").rglob("*.md")),
        *sorted((ROOT / "examples").rglob("*.py")),
    )

    for path in paths:
        source = path.read_text(encoding="utf-8")
        if path.suffix == ".py":
            violations.extend(_deprecated_api_imports(path, source))
            continue
        for line_number, block in _python_code_blocks(source):
            block_path = Path(f"{path}:{line_number}")
            violations.extend(_deprecated_api_imports(block_path, block))

    assert violations == []


def test_benchmarks_use_supported_api_import_routes() -> None:
    violations: list[str] = []
    for path in sorted((ROOT / "benchmarks").rglob("*.py")):
        violations.extend(
            _deprecated_api_imports(path, path.read_text(encoding="utf-8"))
        )

    assert violations == []


def test_pytest_errors_on_unexpected_phospy_deprecation_warnings(
    tmp_path: Path,
) -> None:
    config = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    filterwarnings = config["tool"]["pytest"]["ini_options"]["filterwarnings"]

    assert "error::phospy._deprecations.PhosPyDeprecationWarning" in filterwarnings

    test_file = tmp_path / "test_uncaptured_phospy_deprecation.py"
    test_file.write_text(
        "\n".join(
            [
                "from phospy._deprecations import warn_deprecated",
                "",
                "def test_uncaptured_phospy_deprecation_fails():",
                "    warn_deprecated(",
                "        'science.differential.DifferentialAnalysis',",
                "        stacklevel=1,",
                "    )",
                "",
            ]
        ),
        encoding="utf-8",
    )
    environment = dict(os.environ)
    environment["PYTHONPATH"] = os.pathsep.join(
        [str(ROOT / "src"), environment.get("PYTHONPATH", "")]
    )
    completed = subprocess.run(
        [
            sys.executable,
            "-m",
            "pytest",
            "-c",
            str(ROOT / "pyproject.toml"),
            str(test_file),
            "-q",
        ],
        cwd=ROOT,
        env=environment,
        text=True,
        capture_output=True,
        check=False,
    )

    output = f"{completed.stdout}\n{completed.stderr}"
    assert completed.returncode != 0, output
    assert "PhosPyDeprecationWarning" in output
    assert "DifferentialAnalysis" in output


def _python_code_blocks(source: str) -> tuple[tuple[int, str], ...]:
    blocks: list[tuple[int, str]] = []
    for match in re.finditer(
        r"```python[^\n]*\n(?P<code>.*?)\n```",
        source,
        flags=re.DOTALL,
    ):
        line_number = source[: match.start("code")].count("\n") + 1
        blocks.append((line_number, match.group("code")))
    return tuple(blocks)


def _deprecated_api_imports(path: Path, source: str) -> list[str]:
    try:
        tree = ast.parse(source)
    except SyntaxError:
        return []

    violations: list[str] = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.ImportFrom):
            continue
        if node.module is None or not (
            node.module == "phospy.api" or node.module.startswith("phospy.api.")
        ):
            continue
        if any(alias.name == "*" for alias in node.names):
            continue

        module = importlib.import_module(node.module)
        for alias in node.names:
            if alias.name in module.__dict__:
                continue
            compat = compatibility_export(old_module=node.module, name=alias.name)
            if compat is None:
                continue
            violations.append(
                f"{path}:{node.lineno}: {node.module}.{alias.name} is deprecated; "
                f"use {compat.replacement_module}.{alias.name}"
            )
    return violations
