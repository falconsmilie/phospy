"""Static supported-import checks for current user-facing Markdown."""

from __future__ import annotations

import ast
import re
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True, slots=True)
class DocumentedImport:
    path: str
    line: int
    module: str
    names: tuple[str, ...]


SUPPORTED_AGGREGATES = frozenset({"phospy", "phospy.api", "phospy.advanced"})
_NAV_MARKDOWN = re.compile(r":\s+([^#\s]+\.md)\s*$", re.MULTILINE)


def user_facing_markdown_paths(root: Path) -> tuple[Path, ...]:
    """Return the explicit current documentation set, excluding history/evidence."""

    mkdocs = (root / "mkdocs.yml").read_text(encoding="utf-8")
    nav_paths = {
        root / "docs" / match.group(1) for match in _NAV_MARKDOWN.finditer(mkdocs)
    }
    current_docs = set((root / "docs").glob("*.md"))
    current_docs.update((root / "docs" / "api").rglob("*.md"))
    return tuple(
        sorted(
            {
                root / "README.md",
                root / "examples/README.md",
                root / "benchmarks/README.md",
                *current_docs,
                *nav_paths,
            }
        )
    )


def documented_import_problems(
    *,
    path: str,
    markdown: str,
    inventories: Mapping[str, frozenset[str]],
) -> tuple[list[str], list[DocumentedImport]]:
    """Validate every PhosPy import in fenced Markdown examples."""

    problems: list[str] = []
    imports: list[DocumentedImport] = []
    for start_line, code in _fenced_blocks(markdown):
        direct_import = re.search(r"\b(?:from|import)\s+phospy\b", code)
        dynamic_import = "phospy" in code and re.search(
            r"\b(?:import_module|__import__|getattr)\b", code
        )
        if not direct_import and not dynamic_import:
            continue
        try:
            tree = ast.parse(code, filename=path)
        except SyntaxError as exc:
            problems.append(
                f"{path}:{start_line + (exc.lineno or 1) - 1}: malformed PhosPy import example: {exc.msg}"
            )
            continue
        for node in ast.walk(tree):
            line = start_line + getattr(node, "lineno", 1) - 1
            if isinstance(node, ast.Import):
                for alias in node.names:
                    if not alias.name.startswith("phospy"):
                        continue
                    if alias.name not in SUPPORTED_AGGREGATES:
                        problems.append(
                            f"{path}:{line}: owner-internal or retired import module: {alias.name}"
                        )
                    if alias.asname is not None:
                        problems.append(
                            f"{path}:{line}: alias-based PhosPy import is forbidden: {alias.asname}"
                        )
                    imports.append(DocumentedImport(path, line, alias.name, ()))
            elif isinstance(node, ast.ImportFrom):
                module = node.module or ""
                if not module.startswith("phospy"):
                    continue
                names = tuple(alias.name for alias in node.names)
                if node.level != 0 or module not in SUPPORTED_AGGREGATES:
                    problems.append(
                        f"{path}:{line}: owner-internal, retired, or relative import module: {module}"
                    )
                for alias in node.names:
                    if alias.name == "*":
                        problems.append(
                            f"{path}:{line}: wildcard PhosPy import is forbidden"
                        )
                    elif alias.asname is not None:
                        problems.append(
                            f"{path}:{line}: alias-based PhosPy import is forbidden: {alias.asname}"
                        )
                    elif (
                        module in inventories and alias.name not in inventories[module]
                    ):
                        problems.append(
                            f"{path}:{line}: {alias.name!r} is not exported by {module}"
                        )
                imports.append(DocumentedImport(path, line, module, names))
            elif isinstance(node, ast.Call) and _is_dynamic_or_reflective_import(node):
                problems.append(
                    f"{path}:{line}: dynamic or reflective PhosPy import is forbidden"
                )
    return problems, imports


def _fenced_blocks(markdown: str) -> tuple[tuple[int, str], ...]:
    blocks: list[tuple[int, str]] = []
    fence: str | None = None
    start_line = 0
    lines: list[str] = []
    for line_number, line in enumerate(markdown.splitlines(), start=1):
        stripped = line.lstrip()
        marker = (
            "```"
            if stripped.startswith("```")
            else "~~~"
            if stripped.startswith("~~~")
            else None
        )
        if fence is None and marker is not None:
            fence = marker
            start_line = line_number + 1
            lines = []
        elif fence is not None and stripped.startswith(fence):
            blocks.append((start_line, "\n".join(lines)))
            fence = None
            lines = []
        elif fence is not None:
            lines.append(line)
    if fence is not None:
        blocks.append((start_line, "\n".join(lines)))
    return tuple(blocks)


def _is_dynamic_or_reflective_import(node: ast.Call) -> bool:
    if isinstance(node.func, ast.Name) and node.func.id in {"__import__", "getattr"}:
        return any(_mentions_phospy(argument) for argument in node.args)
    if isinstance(node.func, ast.Attribute) and node.func.attr == "import_module":
        return any(_mentions_phospy(argument) for argument in node.args)
    return False


def _mentions_phospy(node: ast.AST) -> bool:
    return any(
        (isinstance(child, ast.Name) and child.id == "phospy")
        or (
            isinstance(child, ast.Constant)
            and isinstance(child.value, str)
            and "phospy" in child.value
        )
        for child in ast.walk(node)
    )
