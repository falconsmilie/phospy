"""Transactional directory lifecycle helpers for result bundles."""

from __future__ import annotations

import shutil
import tempfile
import uuid
from collections.abc import Callable
from pathlib import Path

from phospy.errors.input import PhosPyInputError


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
    except Exception:
        if staged_root.exists():
            shutil.rmtree(staged_root, ignore_errors=True)
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
    if not final_root.exists():
        try:
            staged_root.replace(final_root)
            return
        except OSError as exc:
            raise PhosPyInputError(
                f"failed to promote staged bundle '{staged_root}' to '{final_root}': {exc}"
            ) from exc

    if not overwrite:
        raise PhosPyInputError(
            f"bundle output directory already exists: {final_root}. "
            "Pass overwrite=True to replace it transactionally."
        )

    backup_root = final_root.parent / f".{final_root.name}.previous-{uuid.uuid4().hex}"
    try:
        final_root.replace(backup_root)
        staged_root.replace(final_root)
    except OSError as exc:
        if not final_root.exists() and backup_root.exists():
            try:
                backup_root.replace(final_root)
            except OSError:
                pass
        raise PhosPyInputError(
            f"failed to promote staged bundle '{staged_root}' to '{final_root}': {exc}"
        ) from exc
    finally:
        if backup_root.exists():
            shutil.rmtree(backup_root, ignore_errors=True)


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
