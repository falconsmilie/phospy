"""Documentation terminology guardrails."""

from __future__ import annotations

from pathlib import Path

_ROOT = Path(__file__).resolve().parents[2]
_DOCS_ROOT = _ROOT / "docs"
_SRC_ROOT = _ROOT / "src" / "phospy"
_README = _ROOT / "README.md"


def _public_module_paths() -> tuple[Path, ...]:
    explicit_paths = (
        _SRC_ROOT / "__init__.py",
        _SRC_ROOT / "science" / "datasets" / "models.py",
        _SRC_ROOT / "science" / "activities" / "threshold_membership.py",
    )
    api_paths = tuple((_SRC_ROOT / "api").rglob("*.py"))
    workflow_public_paths = tuple((_SRC_ROOT / "workflows").glob("*/public.py"))
    return explicit_paths + api_paths + workflow_public_paths


def _read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="ignore")


def test_docs_markdown_avoids_canonical_term() -> None:
    markdown_paths = tuple(_DOCS_ROOT.rglob("*.md")) + (_README,)
    for markdown_path in markdown_paths:
        assert "canonical" not in _read_text(markdown_path).lower(), (
            f"documentation must avoid the term 'canonical': {markdown_path.as_posix()}"
        )


def test_selected_public_modules_avoid_canonical_term() -> None:
    for module_path in _public_module_paths():
        assert "canonical" not in _read_text(module_path).lower(), (
            "selected public-facing module text must avoid the term "
            f"'canonical': {module_path.as_posix()}"
        )


def test_batch_correction_docs_keep_readiness_residualisation_and_native_correction_distinct() -> (
    None
):
    docs_text = "\n\n".join(
        _read_text(path)
        for path in (
            _README,
            _DOCS_ROOT / "workflow_contracts.md",
            _DOCS_ROOT / "validation.md",
            _DOCS_ROOT / "api" / "dataset-build-workflow.md",
            _DOCS_ROOT / "scientific-coverage.md",
            _DOCS_ROOT / "parity.md",
        )
    )
    normalised = " ".join(docs_text.split()).casefold()

    assert (
        "`linear_residualize_batch`, a limited fixed-effect residualisation"
        in normalised
    )
    assert "`ruv_readiness` diagnostics are report-only" in normalised
    assert (
        "native phospy sps/ruv-style preprocessing correction estimates "
        "unwanted factors from eligible control-site residuals after "
        "protected-design handling"
    ) in normalised
    assert "batch terms are resolved for validation and diagnostics" in normalised
    assert (
        "not directly residualized as fixed effects by the native correction"
        in normalised
        or "not directly residualized as fixed effects by native sps/ruv-style "
        "correction"
        in normalised
    )
