"""Transactional directory lifecycle helpers for generated output trees."""

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

    return write_output_tree_atomically(
        output_root=output_root,
        overwrite=overwrite,
        artifact_label="bundle",
        required_relative_paths=(Path(manifest_filename),),
        write_staged_tree=write_staged_bundle,
    )


def write_output_tree_atomically(
    *,
    output_root: Path,
    overwrite: bool,
    artifact_label: str,
    required_relative_paths: tuple[Path, ...],
    write_staged_tree: Callable[[Path], dict[str, Path]],
    allow_non_directory_overwrite: bool = False,
) -> dict[str, Path]:
    """Write a generated output tree transactionally and return final paths.

    The caller writes every output into ``staged_root``. This helper validates
    required staged files before any final-path mutation, then promotes the
    complete staged tree with backup/rollback semantics for overwrite.
    """

    final_root = Path(output_root)
    _validate_destination(
        final_root,
        overwrite=overwrite,
        artifact_label=artifact_label,
        allow_non_directory_overwrite=allow_non_directory_overwrite,
    )
    parent = final_root.parent
    try:
        parent.mkdir(parents=True, exist_ok=True)
        staged_root = Path(
            tempfile.mkdtemp(prefix=f".{final_root.name}.tmp-", dir=parent)
        )
    except OSError as exc:
        raise PhosPyInputError(
            "failed to create temporary "
            f"{artifact_label} directory for '{final_root}': {exc}"
        ) from exc

    try:
        staged_written = write_staged_tree(staged_root)
        for relative_path in required_relative_paths:
            staged_required_file = staged_root / relative_path
            if staged_required_file.is_file():
                continue
            raise PhosPyInputError(
                f"{artifact_label} writer did not complete staged required file: "
                f"{staged_required_file}"
            )
        final_written = _rebase_written_paths(
            staged_written,
            staged_root=staged_root,
            final_root=final_root,
            artifact_label=artifact_label,
        )
        _promote_staged_tree(
            staged_root=staged_root,
            final_root=final_root,
            overwrite=overwrite,
            artifact_label=artifact_label,
        )
        return final_written
    except Exception as exc:
        _cleanup_staged_after_failure(
            staged_root,
            original_error=exc,
            artifact_label=artifact_label,
        )
        raise


def _validate_destination(
    final_root: Path,
    *,
    overwrite: bool,
    artifact_label: str,
    allow_non_directory_overwrite: bool,
) -> None:
    if not _path_lexists(final_root):
        return
    if final_root.is_symlink():
        if not overwrite:
            raise PhosPyInputError(
                f"{artifact_label} output path already exists: {final_root}. "
                "Pass overwrite=True to replace it transactionally."
            )
        return
    if final_root.is_dir():
        if not overwrite:
            raise PhosPyInputError(
                f"{artifact_label} output directory already exists: {final_root}. "
                "Pass overwrite=True to replace it transactionally."
            )
        return
    if allow_non_directory_overwrite and overwrite:
        return
    if overwrite:
        raise PhosPyInputError(
            f"{artifact_label} output path exists and is not a directory: {final_root}"
        )
    raise PhosPyInputError(
        f"{artifact_label} output path already exists: {final_root}. "
        "Pass overwrite=True to replace it transactionally."
    )


def _promote_staged_tree(
    *,
    staged_root: Path,
    final_root: Path,
    overwrite: bool,
    artifact_label: str,
) -> None:
    state = _BundleTransactionState.STAGED
    if not _path_lexists(final_root):
        try:
            staged_root.replace(final_root)
            state = _BundleTransactionState.STAGED_PROMOTED
            return
        except OSError as exc:
            raise PhosPyInputError(
                f"failed to promote staged {artifact_label} "
                f"'{staged_root}' to '{final_root}' "
                f"(transaction_state='{state.value}'): {exc}"
            ) from exc

    if not overwrite:
        raise PhosPyInputError(
            f"{artifact_label} output path already exists: {final_root}. "
            "Pass overwrite=True to replace it transactionally."
        )

    backup_root = final_root.parent / f".{final_root.name}.previous-{uuid.uuid4().hex}"
    try:
        final_root.replace(backup_root)
        state = _BundleTransactionState.ORIGINAL_MOVED_TO_BACKUP
    except OSError as exc:
        raise PhosPyInputError(
            f"failed to move existing {artifact_label} "
            f"'{final_root}' to recovery backup '{backup_root}' "
            f"(transaction_state='{state.value}'): {exc}"
        ) from exc

    try:
        staged_root.replace(final_root)
        state = _BundleTransactionState.STAGED_PROMOTED
    except OSError as exc:
        if not _path_lexists(final_root) and _path_lexists(backup_root):
            _restore_backup_after_failed_promotion(
                backup_root=backup_root,
                final_root=final_root,
                staged_root=staged_root,
                promotion_error=exc,
                state=state,
                artifact_label=artifact_label,
            )
            state = _BundleTransactionState.ROLLBACK_RESTORED
            raise PhosPyInputError(
                f"failed to promote staged {artifact_label} "
                f"'{staged_root}' to '{final_root}', but rollback restored "
                f"the original {artifact_label} from recovery backup '{backup_root}' "
                f"(transaction_state='{state.value}'): {exc}"
            ) from exc
        raise PhosPyInputError(
            f"failed to promote staged {artifact_label} "
            f"'{staged_root}' to '{final_root}' after moving the original to "
            f"recovery backup '{backup_root}'. The recovery backup was retained "
            f"(transaction_state='"
            f"{_BundleTransactionState.RECOVERY_BACKUP_RETAINED.value}'): {exc}"
        ) from exc

    _remove_backup_after_success(
        backup_root=backup_root,
        final_root=final_root,
        state=state,
        artifact_label=artifact_label,
    )


def _restore_backup_after_failed_promotion(
    *,
    backup_root: Path,
    final_root: Path,
    staged_root: Path,
    promotion_error: OSError,
    state: _BundleTransactionState,
    artifact_label: str,
) -> None:
    try:
        backup_root.replace(final_root)
    except OSError as rollback_error:
        retained_state = _BundleTransactionState.RECOVERY_BACKUP_RETAINED
        raise PhosPyInputError(
            f"failed to promote staged {artifact_label} "
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
    artifact_label: str,
) -> None:
    if not _path_lexists(backup_root):
        return
    try:
        _remove_recovery_backup(backup_root)
    except OSError as exc:
        raise PhosPyInputError(
            f"staged {artifact_label} was promoted successfully, but obsolete recovery "
            f"backup cleanup failed. Promoted bundle remains at '{final_root}', "
            f"and the retained backup is at exact path '{backup_root}' "
            f"(transaction_state='{state.value}'): {exc}"
        ) from exc


def _remove_recovery_backup(backup_root: Path) -> None:
    if backup_root.is_symlink() or backup_root.is_file():
        backup_root.unlink()
        return
    shutil.rmtree(backup_root)


def _cleanup_staged_after_failure(
    staged_root: Path,
    *,
    original_error: Exception,
    artifact_label: str = "bundle",
) -> None:
    if not _path_lexists(staged_root):
        return
    try:
        shutil.rmtree(staged_root)
    except OSError as cleanup_error:
        raise PhosPyInputError(
            f"{original_error} Failed to clean staged {artifact_label} directory "
            f"'{staged_root}'. The failed staged directory was retained for "
            f"inspection (cleanup_error='{cleanup_error}')."
        ) from original_error


def _rebase_written_paths(
    written: dict[str, Path],
    *,
    staged_root: Path,
    final_root: Path,
    artifact_label: str = "bundle",
) -> dict[str, Path]:
    staged_resolved = staged_root.resolve()
    rebased: dict[str, Path] = {}
    for key, path in written.items():
        try:
            relative_path = Path(path).resolve().relative_to(staged_resolved)
        except ValueError as exc:
            raise PhosPyInputError(
                f"{artifact_label} writer returned path outside staged root "
                f"for key '{key}': {path}"
            ) from exc
        rebased[key] = final_root / relative_path
    return rebased


def _path_lexists(path: Path) -> bool:
    return path.exists() or path.is_symlink()
