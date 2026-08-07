"""Transactional directory lifecycle helpers for result bundles."""

from __future__ import annotations

import shutil
import tempfile
import uuid
from collections.abc import Callable
from enum import Enum
from pathlib import Path

from phospy.errors.input import PhosPyInputError


class _BundleTransactionState(str, Enum):
    """Internal transaction lifecycle checkpoints for directory promotion."""

    STAGED = "staged"
    ORIGINAL_MOVED_TO_BACKUP = "original moved to backup"
    STAGED_PROMOTED = "staged promoted"
    ROLLBACK_RESTORED = "rollback restored"
    RECOVERY_BACKUP_RETAINED = "recovery backup retained"


def write_bundle_atomically(
    *,
    output_root: Path,
    manifest_filename: str,
    overwrite: bool,
    write_staged_bundle: Callable[[Path], dict[str, Path]],
) -> dict[str, Path]:
    """Write a bundle in a sibling temporary directory, then promote it.

    The caller writes all bundle files, including the manifest, into the staged
    root. This helper promotes only after the manifest exists and never mutates
    the target directory during staged serialization.
    """

    final_root = Path(output_root)
    _validate_destination(final_root, overwrite=overwrite)
    parent = final_root.parent
    try:
        parent.mkdir(parents=True, exist_ok=True)
        staged_root = Path(
            tempfile.mkdtemp(prefix=f".{final_root.name}.tmp-", dir=parent)
        )
    except OSError as exc:
        raise PhosPyInputError(
            f"failed to create temporary bundle directory for '{final_root}': {exc}"
        ) from exc

    try:
        staged_written = write_staged_bundle(staged_root)
        staged_manifest = staged_root / manifest_filename
        if not staged_manifest.is_file():
            raise PhosPyInputError(
                f"bundle writer did not complete staged manifest: {staged_manifest}"
            )
        final_written = _rebase_written_paths(
            staged_written,
            staged_root=staged_root,
            final_root=final_root,
        )
        _promote_staged_bundle(
            staged_root=staged_root,
            final_root=final_root,
            overwrite=overwrite,
        )
        return final_written
    except Exception as exc:
        _cleanup_staged_after_failure(staged_root, original_error=exc)
        raise


def _validate_destination(final_root: Path, *, overwrite: bool) -> None:
    if not final_root.exists():
        return
    if not final_root.is_dir():
        raise PhosPyInputError(
            f"bundle output path exists and is not a directory: {final_root}"
        )
    if not overwrite:
        raise PhosPyInputError(
            f"bundle output directory already exists: {final_root}. "
            "Pass overwrite=True to replace it transactionally."
        )


def _promote_staged_bundle(
    *,
    staged_root: Path,
    final_root: Path,
    overwrite: bool,
) -> None:
    state = _BundleTransactionState.STAGED
    if not final_root.exists():
        try:
            staged_root.replace(final_root)
            state = _BundleTransactionState.STAGED_PROMOTED
            return
        except OSError as exc:
            raise PhosPyInputError(
                "failed to promote staged bundle "
                f"'{staged_root}' to '{final_root}' "
                f"(transaction_state='{state.value}'): {exc}"
            ) from exc

    if not overwrite:
        raise PhosPyInputError(
            f"bundle output directory already exists: {final_root}. "
            "Pass overwrite=True to replace it transactionally."
        )

    backup_root = final_root.parent / f".{final_root.name}.previous-{uuid.uuid4().hex}"
    try:
        final_root.replace(backup_root)
        state = _BundleTransactionState.ORIGINAL_MOVED_TO_BACKUP
    except OSError as exc:
        raise PhosPyInputError(
            "failed to move existing bundle "
            f"'{final_root}' to recovery backup '{backup_root}' "
            f"(transaction_state='{state.value}'): {exc}"
        ) from exc

    try:
        staged_root.replace(final_root)
        state = _BundleTransactionState.STAGED_PROMOTED
    except OSError as exc:
        if not final_root.exists() and backup_root.exists():
            _restore_backup_after_failed_promotion(
                backup_root=backup_root,
                final_root=final_root,
                staged_root=staged_root,
                promotion_error=exc,
                state=state,
            )
            state = _BundleTransactionState.ROLLBACK_RESTORED
            raise PhosPyInputError(
                "failed to promote staged bundle "
                f"'{staged_root}' to '{final_root}', but rollback restored "
                f"the original bundle from recovery backup '{backup_root}' "
                f"(transaction_state='{state.value}'): {exc}"
            ) from exc
        raise PhosPyInputError(
            "failed to promote staged bundle "
            f"'{staged_root}' to '{final_root}' after moving the original to "
            f"recovery backup '{backup_root}'. The recovery backup was retained "
            f"(transaction_state='"
            f"{_BundleTransactionState.RECOVERY_BACKUP_RETAINED.value}'): {exc}"
        ) from exc

    _remove_backup_after_success(
        backup_root=backup_root,
        final_root=final_root,
        state=state,
    )


def _restore_backup_after_failed_promotion(
    *,
    backup_root: Path,
    final_root: Path,
    staged_root: Path,
    promotion_error: OSError,
    state: _BundleTransactionState,
) -> None:
    try:
        backup_root.replace(final_root)
    except OSError as rollback_error:
        retained_state = _BundleTransactionState.RECOVERY_BACKUP_RETAINED
        raise PhosPyInputError(
            "failed to promote staged bundle "
            f"'{staged_root}' to '{final_root}' after moving the original to "
            f"recovery backup '{backup_root}', and rollback failed. Recovery "
            f"backup retained at exact path '{backup_root}' "
            f"(transaction_state='{retained_state.value}'; "
            f"previous_state='{state.value}'; "
            f"promotion_error='{promotion_error}'; "
            f"rollback_error='{rollback_error}')"
        ) from rollback_error


def _remove_backup_after_success(
    *,
    backup_root: Path,
    final_root: Path,
    state: _BundleTransactionState,
) -> None:
    if not backup_root.exists():
        return
    try:
        shutil.rmtree(backup_root)
    except OSError as exc:
        raise PhosPyInputError(
            "staged bundle was promoted successfully, but obsolete recovery "
            f"backup cleanup failed. Promoted bundle remains at '{final_root}', "
            f"and the retained backup is at exact path '{backup_root}' "
            f"(transaction_state='{state.value}'): {exc}"
        ) from exc


def _cleanup_staged_after_failure(
    staged_root: Path,
    *,
    original_error: Exception,
) -> None:
    if not staged_root.exists():
        return
    try:
        shutil.rmtree(staged_root)
    except OSError as cleanup_error:
        raise PhosPyInputError(
            f"{original_error} Failed to clean staged bundle directory "
            f"'{staged_root}'. The failed staged directory was retained for "
            f"inspection (cleanup_error='{cleanup_error}')."
        ) from original_error


def _rebase_written_paths(
    written: dict[str, Path],
    *,
    staged_root: Path,
    final_root: Path,
) -> dict[str, Path]:
    staged_resolved = staged_root.resolve()
    rebased: dict[str, Path] = {}
    for key, path in written.items():
        try:
            relative_path = Path(path).resolve().relative_to(staged_resolved)
        except ValueError as exc:
            raise PhosPyInputError(
                f"bundle writer returned path outside staged root for key '{key}': {path}"
            ) from exc
        rebased[key] = final_root / relative_path
    return rebased
