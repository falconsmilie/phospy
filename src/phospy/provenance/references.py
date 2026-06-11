"""Reference-source provenance helpers."""

from __future__ import annotations

from hashlib import sha256
from pathlib import Path

from phospy.errors.input import PhosPyInputError
from phospy.provenance.models import JsonValue


def fingerprint_local_reference_source_file(
    path: Path,
    *,
    role: str,
) -> dict[str, JsonValue]:
    """Return JSON-safe identity metadata for one local reference source file."""

    normalized_path = Path(path)
    try:
        raw_bytes = normalized_path.read_bytes()
    except FileNotFoundError as exc:
        raise PhosPyInputError(
            f"reference source file does not exist: {normalized_path}"
        ) from exc
    except PermissionError as exc:
        raise PhosPyInputError(
            f"permission denied while reading reference source file: {normalized_path}"
        ) from exc
    except OSError as exc:
        raise PhosPyInputError(
            f"failed to read reference source file '{normalized_path}': {exc}"
        ) from exc
    return {
        "role": role,
        "path": str(normalized_path),
        "resolved_path": str(normalized_path.resolve()),
        "sha256": sha256(raw_bytes).hexdigest(),
        "bytes": len(raw_bytes),
    }


__all__ = ["fingerprint_local_reference_source_file"]
