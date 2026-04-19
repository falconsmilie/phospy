"""Filesystem path constraints for bundle-relative references."""

from __future__ import annotations

from pathlib import Path

from phospy.errors.input import PhosPyInputError


def resolve_bundle_relative_path(
    bundle_root: Path,
    relative_path: str,
    *,
    field_name: str,
) -> Path:
    """Resolve a manifest relative path and enforce bundle-root containment."""

    candidate = Path(relative_path)
    if candidate.is_absolute():
        raise PhosPyInputError(f"{field_name} must be a relative bundle path")
    resolved_root = bundle_root.resolve()
    resolved_candidate = (bundle_root / candidate).resolve()
    try:
        resolved_candidate.relative_to(resolved_root)
    except ValueError as exc:
        raise PhosPyInputError(
            f"{field_name} points outside the bundle root: {relative_path}"
        ) from exc
    return resolved_candidate
