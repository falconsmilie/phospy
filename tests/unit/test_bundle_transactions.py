from __future__ import annotations

import os
from collections.abc import Callable
from pathlib import Path
from typing import Any

import pytest

from phospy.errors.input import PhosPyInputError
from phospy.io.bundles._shared import transactions
from phospy.io.bundles._shared.transactions import write_bundle_atomically

MANIFEST_FILENAME = "manifest.json"


def test_new_write_promotes_staged_bundle_and_returns_final_paths(
    tmp_path: Path,
) -> None:
    bundle_root = tmp_path / "bundle"

    written = write_bundle_atomically(
        output_root=bundle_root,
        manifest_filename=MANIFEST_FILENAME,
        overwrite=False,
        write_staged_bundle=lambda staged_root: _write_minimal_bundle(
            staged_root,
            label="new",
        ),
    )

    assert _payload_text(bundle_root) == "new"
    assert written == {
        "manifest": bundle_root / MANIFEST_FILENAME,
        "payload": bundle_root / "data" / "payload.txt",
    }
    assert all(path.exists() for path in written.values())
    assert not list(tmp_path.glob(".bundle.tmp-*"))
    assert not list(tmp_path.glob(".bundle.previous-*"))


def test_successful_overwrite_removes_stale_files_and_backup(
    tmp_path: Path,
) -> None:
    bundle_root = tmp_path / "bundle"
    _create_existing_bundle(bundle_root, label="original")
    stale_path = bundle_root / "stale.txt"
    stale_path.write_text("stale", encoding="utf-8")

    written = write_bundle_atomically(
        output_root=bundle_root,
        manifest_filename=MANIFEST_FILENAME,
        overwrite=True,
        write_staged_bundle=lambda staged_root: _write_minimal_bundle(
            staged_root,
            label="replacement",
        ),
    )

    assert _payload_text(bundle_root) == "replacement"
    assert not stale_path.exists()
    assert written["payload"] == bundle_root / "data" / "payload.txt"
    assert not list(tmp_path.glob(".bundle.tmp-*"))
    assert not list(tmp_path.glob(".bundle.previous-*"))


def test_staged_serialization_failure_leaves_original_untouched(
    tmp_path: Path,
) -> None:
    bundle_root = tmp_path / "bundle"
    _create_existing_bundle(bundle_root, label="original")

    def fail_serialization(staged_root: Path) -> dict[str, Path]:
        _write_minimal_bundle(staged_root, label="partial")
        raise PhosPyInputError("simulated staged serialization failure")

    with pytest.raises(
        PhosPyInputError,
        match="simulated staged serialization failure",
    ):
        write_bundle_atomically(
            output_root=bundle_root,
            manifest_filename=MANIFEST_FILENAME,
            overwrite=True,
            write_staged_bundle=fail_serialization,
        )

    assert _payload_text(bundle_root) == "original"
    assert not list(tmp_path.glob(".bundle.tmp-*"))
    assert not list(tmp_path.glob(".bundle.previous-*"))


def test_new_write_replace_failure_publishes_no_target_and_cleans_stage(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    bundle_root = tmp_path / "bundle"
    staged_root: Path | None = None
    _fail_path_replace_when(
        monkeypatch,
        lambda source, target: (
            staged_root is not None and source == staged_root and target == bundle_root
        ),
        PermissionError("simulated new-write promotion failure"),
    )

    def writer(root: Path) -> dict[str, Path]:
        nonlocal staged_root
        staged_root = root
        return _write_minimal_bundle(root, label="new")

    with pytest.raises(
        PhosPyInputError,
        match="failed to promote staged bundle.*transaction_state='staged'",
    ):
        write_bundle_atomically(
            output_root=bundle_root,
            manifest_filename=MANIFEST_FILENAME,
            overwrite=False,
            write_staged_bundle=writer,
        )

    assert not bundle_root.exists()
    assert not list(tmp_path.glob(".bundle.tmp-*"))


def test_overwrite_original_to_backup_replace_failure_keeps_original(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    bundle_root = tmp_path / "bundle"
    _create_existing_bundle(bundle_root, label="original")
    _fail_path_replace_when(
        monkeypatch,
        lambda source, target: (
            source == bundle_root and target.name.startswith(".bundle.previous-")
        ),
        PermissionError("simulated backup move failure"),
    )

    with pytest.raises(
        PhosPyInputError,
        match="failed to move existing bundle.*recovery backup",
    ):
        write_bundle_atomically(
            output_root=bundle_root,
            manifest_filename=MANIFEST_FILENAME,
            overwrite=True,
            write_staged_bundle=lambda staged_root: _write_minimal_bundle(
                staged_root,
                label="replacement",
            ),
        )

    assert _payload_text(bundle_root) == "original"
    assert not list(tmp_path.glob(".bundle.tmp-*"))
    assert not list(tmp_path.glob(".bundle.previous-*"))


def test_overwrite_staged_to_final_replace_failure_restores_original(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    bundle_root = tmp_path / "bundle"
    _create_existing_bundle(bundle_root, label="original")
    staged_root: Path | None = None
    _fail_path_replace_when(
        monkeypatch,
        lambda source, target: (
            staged_root is not None and source == staged_root and target == bundle_root
        ),
        PermissionError("simulated staged promotion failure"),
    )

    def writer(root: Path) -> dict[str, Path]:
        nonlocal staged_root
        staged_root = root
        return _write_minimal_bundle(root, label="replacement")

    with pytest.raises(
        PhosPyInputError,
        match="rollback restored the original bundle",
    ):
        write_bundle_atomically(
            output_root=bundle_root,
            manifest_filename=MANIFEST_FILENAME,
            overwrite=True,
            write_staged_bundle=writer,
        )

    assert _payload_text(bundle_root) == "original"
    assert not list(tmp_path.glob(".bundle.tmp-*"))
    assert not list(tmp_path.glob(".bundle.previous-*"))


def test_overwrite_promotion_and_rollback_failure_retains_recovery_backup(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    bundle_root = tmp_path / "bundle"
    _create_existing_bundle(bundle_root, label="original")
    staged_root: Path | None = None
    _fail_path_replace_when(
        monkeypatch,
        lambda source, target: (
            (
                staged_root is not None
                and source == staged_root
                and target == bundle_root
            )
            or (source.name.startswith(".bundle.previous-") and target == bundle_root)
        ),
        PermissionError("simulated promotion and rollback failure"),
    )

    def writer(root: Path) -> dict[str, Path]:
        nonlocal staged_root
        staged_root = root
        return _write_minimal_bundle(root, label="replacement")

    with pytest.raises(
        PhosPyInputError,
        match="rollback failed.*Recovery backup retained",
    ) as exc_info:
        write_bundle_atomically(
            output_root=bundle_root,
            manifest_filename=MANIFEST_FILENAME,
            overwrite=True,
            write_staged_bundle=writer,
        )

    backup_root = _single_recovery_backup(tmp_path, bundle_name="bundle")
    assert str(backup_root) in str(exc_info.value)
    assert not bundle_root.exists()
    assert _payload_text(backup_root) == "original"
    assert not list(tmp_path.glob(".bundle.tmp-*"))


def test_backup_cleanup_failure_preserves_promoted_bundle_and_retains_backup(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    bundle_root = tmp_path / "bundle"
    _create_existing_bundle(bundle_root, label="original")
    _fail_rmtree_when(
        monkeypatch,
        lambda path: path.name.startswith(".bundle.previous-"),
        PermissionError("simulated backup cleanup failure"),
    )

    with pytest.raises(
        PhosPyInputError,
        match="promoted successfully.*backup cleanup failed",
    ) as exc_info:
        write_bundle_atomically(
            output_root=bundle_root,
            manifest_filename=MANIFEST_FILENAME,
            overwrite=True,
            write_staged_bundle=lambda staged_root: _write_minimal_bundle(
                staged_root,
                label="replacement",
            ),
        )

    backup_root = _single_recovery_backup(tmp_path, bundle_name="bundle")
    assert str(backup_root) in str(exc_info.value)
    assert _payload_text(bundle_root) == "replacement"
    assert _payload_text(backup_root) == "original"
    assert not list(tmp_path.glob(".bundle.tmp-*"))


def test_staged_directory_cleanup_failure_reports_retained_stage_and_keeps_original(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    bundle_root = tmp_path / "bundle"
    _create_existing_bundle(bundle_root, label="original")
    _fail_rmtree_when(
        monkeypatch,
        lambda path: path.name.startswith(".bundle.tmp-"),
        PermissionError("simulated staged cleanup failure"),
    )

    def fail_serialization(staged_root: Path) -> dict[str, Path]:
        _write_minimal_bundle(staged_root, label="partial")
        raise PhosPyInputError("simulated serialization failure")

    with pytest.raises(
        PhosPyInputError,
        match="Failed to clean staged bundle directory",
    ) as exc_info:
        write_bundle_atomically(
            output_root=bundle_root,
            manifest_filename=MANIFEST_FILENAME,
            overwrite=True,
            write_staged_bundle=fail_serialization,
        )

    staged_roots = list(tmp_path.glob(".bundle.tmp-*"))
    assert len(staged_roots) == 1
    assert str(staged_roots[0]) in str(exc_info.value)
    assert _payload_text(bundle_root) == "original"
    assert _payload_text(staged_roots[0]) == "partial"
    assert not list(tmp_path.glob(".bundle.previous-*"))


def test_existing_target_file_is_rejected_before_staging(tmp_path: Path) -> None:
    bundle_root = tmp_path / "bundle"
    bundle_root.write_text("not a directory", encoding="utf-8")

    with pytest.raises(
        PhosPyInputError,
        match="bundle output path exists and is not a directory",
    ):
        write_bundle_atomically(
            output_root=bundle_root,
            manifest_filename=MANIFEST_FILENAME,
            overwrite=True,
            write_staged_bundle=lambda staged_root: _write_minimal_bundle(
                staged_root,
                label="replacement",
            ),
        )

    assert bundle_root.read_text(encoding="utf-8") == "not a directory"
    assert not list(tmp_path.glob(".bundle.tmp-*"))
    assert not list(tmp_path.glob(".bundle.previous-*"))


@pytest.mark.skipif(
    os.name != "nt",
    reason="Windows file-lock behavior is platform-specific",
)
def test_windows_locked_existing_bundle_blocks_backup_move_and_keeps_original(
    tmp_path: Path,
) -> None:
    import msvcrt

    bundle_root = tmp_path / "bundle"
    _create_existing_bundle(bundle_root, label="original")
    locked_path = bundle_root / "data" / "payload.txt"

    with locked_path.open("r+", encoding="utf-8") as handle:
        handle.seek(0)
        msvcrt.locking(handle.fileno(), msvcrt.LK_NBLCK, 1)
        try:
            with pytest.raises(
                PhosPyInputError,
                match="failed to move existing bundle.*recovery backup",
            ):
                write_bundle_atomically(
                    output_root=bundle_root,
                    manifest_filename=MANIFEST_FILENAME,
                    overwrite=True,
                    write_staged_bundle=lambda staged_root: _write_minimal_bundle(
                        staged_root,
                        label="replacement",
                    ),
                )
        finally:
            handle.seek(0)
            msvcrt.locking(handle.fileno(), msvcrt.LK_UNLCK, 1)

    assert _payload_text(bundle_root) == "original"
    assert not list(tmp_path.glob(".bundle.tmp-*"))
    assert not list(tmp_path.glob(".bundle.previous-*"))


def _create_existing_bundle(bundle_root: Path, *, label: str) -> None:
    bundle_root.mkdir(parents=True)
    _write_minimal_bundle(bundle_root, label=label)


def _write_minimal_bundle(bundle_root: Path, *, label: str) -> dict[str, Path]:
    data_root = bundle_root / "data"
    data_root.mkdir()
    payload_path = data_root / "payload.txt"
    payload_path.write_text(label, encoding="utf-8")
    manifest_path = bundle_root / MANIFEST_FILENAME
    manifest_path.write_text('{"manifest": true}', encoding="utf-8")
    return {"manifest": manifest_path, "payload": payload_path}


def _payload_text(bundle_root: Path) -> str:
    return (bundle_root / "data" / "payload.txt").read_text(encoding="utf-8")


def _single_recovery_backup(tmp_path: Path, *, bundle_name: str) -> Path:
    backup_roots = list(tmp_path.glob(f".{bundle_name}.previous-*"))
    assert len(backup_roots) == 1
    return backup_roots[0]


def _fail_path_replace_when(
    monkeypatch: pytest.MonkeyPatch,
    predicate: Callable[[Path, Path], bool],
    error: OSError,
) -> None:
    original_replace = Path.replace

    def replace(self: Path, target: str | Path) -> Path:
        source = Path(self)
        target_path = Path(target)
        if predicate(source, target_path):
            raise error
        return original_replace(self, target)

    monkeypatch.setattr(Path, "replace", replace)


def _fail_rmtree_when(
    monkeypatch: pytest.MonkeyPatch,
    predicate: Callable[[Path], bool],
    error: OSError,
) -> None:
    original_rmtree = transactions.shutil.rmtree

    def rmtree(path: str | Path, *args: Any, **kwargs: Any) -> None:
        root = Path(path)
        if predicate(root):
            raise error
        original_rmtree(path, *args, **kwargs)

    monkeypatch.setattr(transactions.shutil, "rmtree", rmtree)
