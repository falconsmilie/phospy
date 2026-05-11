"""Documentation terminology guardrails."""

from __future__ import annotations

from pathlib import Path

_ROOT = Path(__file__).resolve().parents[2]
_DOCS_ROOT = _ROOT / "docs"
_SRC_ROOT = _ROOT / "src" / "phospy"


def _public_module_paths() -> tuple[Path, ...]:
    explicit_paths = (
        _SRC_ROOT / "__init__.py",
        _SRC_ROOT / "datasets" / "models.py",
        _SRC_ROOT / "activities" / "threshold_membership.py",
    )
    api_paths = tuple((_SRC_ROOT / "api").rglob("*.py"))
    workflow_public_paths = tuple((_SRC_ROOT / "workflows").glob("*/public.py"))
    return explicit_paths + api_paths + workflow_public_paths


def _read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="ignore")


def test_docs_markdown_avoids_canonical_term() -> None:
    for markdown_path in _DOCS_ROOT.rglob("*.md"):
        assert "canonical" not in _read_text(markdown_path).lower(), (
            f"documentation must avoid the term 'canonical': {markdown_path.as_posix()}"
        )


def test_selected_public_modules_avoid_canonical_term() -> None:
    for module_path in _public_module_paths():
        assert "canonical" not in _read_text(module_path).lower(), (
            "selected public-facing module text must avoid the term "
            f"'canonical': {module_path.as_posix()}"
        )
