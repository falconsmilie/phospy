"""JSON file read/write utilities with bundle-level error translation."""

from __future__ import annotations

import json
from pathlib import Path

from phospy.errors.input import PhosPyInputError


def read_json(path: Path, *, label: str) -> object:
    """Read and decode JSON payload from disk."""

    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError as exc:
        raise PhosPyInputError(f"{label} does not exist: {path}") from exc
    except PermissionError as exc:
        raise PhosPyInputError(
            f"permission denied while reading {label}: {path}"
        ) from exc
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise PhosPyInputError(f"failed to parse {label} '{path}': {exc}") from exc


def write_json(path: Path, payload: object, *, label: str) -> None:
    """Encode and write JSON payload to disk."""

    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        path.write_text(
            json.dumps(payload, indent=2, sort_keys=True),
            encoding="utf-8",
        )
    except (OSError, TypeError, ValueError) as exc:
        raise PhosPyInputError(f"failed to write {label} '{path}': {exc}") from exc
