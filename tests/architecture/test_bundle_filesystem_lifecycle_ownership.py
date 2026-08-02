from __future__ import annotations

from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = PROJECT_ROOT / "src" / "phospy"

_FORBIDDEN_RESULT_MODEL_TOKENS = (
    "write_bundle_atomically",
    "verify_bundle_integrity",
    "build_table_file_entry",
    "build_json_file_entry",
    "mkdtemp",
    "rmtree",
    "tempfile",
)


def test_bundle_filesystem_lifecycle_remains_in_io_not_result_models() -> None:
    model_paths = [
        *sorted((SRC_ROOT / "contracts" / "results").glob("*.py")),
        *sorted((SRC_ROOT / "science").glob("**/models.py")),
    ]

    violations: list[str] = []
    for path in model_paths:
        source = path.read_text(encoding="utf-8")
        for token in _FORBIDDEN_RESULT_MODEL_TOKENS:
            if token in source:
                relative_path = path.relative_to(PROJECT_ROOT).as_posix()
                violations.append(f"{relative_path}: {token}")

    assert violations == []
