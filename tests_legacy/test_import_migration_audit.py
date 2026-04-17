from __future__ import annotations

import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
AUDITED_PATHS = (
    REPO_ROOT / "src",
    REPO_ROOT / "examples",
    REPO_ROOT / "docs" / "api.md",
    REPO_ROOT / "README.md",
)
LEGACY_IMPORT_PATTERNS = (
    r"(?<![\w.])phospy\.workflow(?![\w.])",
    r"(?<![\w.])phospy\.dataset(?![\w.])",
    r"(?<![\w.])phospy\.dataset_loader(?![\w.])",
    r"(?<![\w.])phospy\.dataset_schema(?![\w.])",
    r"(?<![\w.])phospy\.dataset_site_matrix(?![\w.])",
    r"(?<![\w.])phospy\.analysis(?![\w.])",
    r"(?<![\w.])phospy\.scoring(?![\w.])",
    r"(?<![\w.])phospy\.prediction\.models(?![\w.])",
    r"(?<![\w.])phospy\.prediction\.service(?![\w.])",
    r"(?<![\w.])phospy\.prediction\.workflows(?![\w.])",
    r"(?<![\w.])phospy\.publishing(?![\w.])",
    r"(?<![\w.])phospy\.writers(?![\w.])",
    r"(?<![\w.])phospy\.constants(?![\w.])",
    r"(?<![\w.])phospy\.types(?![\w.])",
    r"(?<![\w.])phospy\.validation\.errors(?![\w.])",
    r"(?<![\w.])phospy\.preprocessing_services(?![\w.])",
    r"(?<![\w.])phospy\.core_processing(?![\w.])",
    r"(?<![\w.])phospy\.dataset_preprocessing(?![\w.])",
    r"(?<![\w.])phospy\.site_matrix_builder(?![\w.])",
    r"(?<![\w.])phospy\._preprocessing_primitives(?![\w.])",
    r"(?<![\w.])phospy\._protein_correction(?![\w.])",
    r"(?<![\w.])phospy\.matrices(?![\w.])",
    r"(?<![\w.])phospy\.motifs(?![\w.])",
    r"(?<![\w.])phospy\.profiles(?![\w.])",
    r"(?<![\w.])phospy\.orchestration(?![\w.])",
)


def _iter_audited_files() -> list[Path]:
    files: list[Path] = []
    for path in AUDITED_PATHS:
        if path.is_file():
            files.append(path)
            continue
        for file_path in path.rglob("*"):
            if file_path.suffix not in {".py", ".md"}:
                continue
            if "__pycache__" in file_path.parts:
                continue
            files.append(file_path)
    return files


def test_no_removed_flat_modules_are_referenced_in_source_docs_or_examples() -> None:
    checked_files = _iter_audited_files()

    for file_path in checked_files:
        contents = file_path.read_text(encoding="utf-8")
        for legacy_pattern in LEGACY_IMPORT_PATTERNS:
            assert re.search(legacy_pattern, contents) is None, (
                f"Legacy import pattern {legacy_pattern!r} is still referenced in "
                f"{file_path.relative_to(REPO_ROOT)}"
            )


def test_examples_and_api_docs_prefer_domain_package_imports() -> None:
    preferred_files = [
        REPO_ROOT / "README.md",
        REPO_ROOT / "docs" / "api.md",
        *sorted((REPO_ROOT / "examples").glob("*.py")),
    ]

    for file_path in preferred_files:
        contents = file_path.read_text(encoding="utf-8")
        assert "from phospy import" not in contents, (
            f"Root package imports remain in {file_path.relative_to(REPO_ROOT)}"
        )
