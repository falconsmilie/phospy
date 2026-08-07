from __future__ import annotations

from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = PROJECT_ROOT / "src" / "phospy"
BUNDLE_IO_ROOT = SRC_ROOT / "io" / "bundles"

_FORBIDDEN_RESULT_MODEL_TOKENS = (
    "write_bundle_atomically",
    "verify_bundle_integrity",
    "build_table_file_entry",
    "build_json_file_entry",
    "mkdtemp",
    "rmtree",
    "tempfile",
)
_BUNDLE_TRANSACTION_LIFECYCLE_TOKENS = (
    "write_bundle_atomically",
    "_BundleTransactionState",
    "_promote_staged_bundle",
    "_restore_backup_after_failed_promotion",
    "bundle writer returned path outside staged root",
    "Pass overwrite=True to replace it transactionally",
    "failed to promote staged bundle",
    "recovery backup",
    "rollback restored",
    ".previous-",
    ".tmp-",
    "mkdtemp",
    "shutil.rmtree",
    "replace(backup_root",
    "replace(final_root",
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


def test_bundle_filesystem_lifecycle_remains_under_io_bundles() -> None:
    source_paths = [
        path
        for path in sorted(SRC_ROOT.rglob("*.py"))
        if not _is_relative_to(path, BUNDLE_IO_ROOT)
    ]

    violations: list[str] = []
    for path in source_paths:
        source = path.read_text(encoding="utf-8")
        for token in _BUNDLE_TRANSACTION_LIFECYCLE_TOKENS:
            if token in source:
                relative_path = path.relative_to(PROJECT_ROOT).as_posix()
                violations.append(f"{relative_path}: {token}")

    assert violations == []


def _is_relative_to(path: Path, root: Path) -> bool:
    try:
        path.relative_to(root)
    except ValueError:
        return False
    return True
