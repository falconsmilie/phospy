"""Policy-free ordered mechanics for DataFrame-or-path input boundaries."""

from __future__ import annotations

from os import PathLike
from pathlib import Path
from typing import TypeVar, cast

import pandas as pd

from phospy.errors.input import PhosPyInputError

_InputBoundaryError = TypeVar("_InputBoundaryError", bound=PhosPyInputError)


def validate_dataframe_or_path(
    value: object | None,
    *,
    field_name: str,
    allow_none: bool,
    error_type: type[_InputBoundaryError],
) -> pd.DataFrame | str | Path | None:
    """Apply the shared ordered source checks with a caller-owned error policy."""

    if value is None:
        if allow_none:
            return None
        raise error_type(f"dataset build request {field_name} is required")
    if isinstance(value, pd.DataFrame):
        return value
    if isinstance(value, (str, Path, PathLike)):
        if isinstance(value, str) and not value.strip():
            raise error_type(f"dataset build request {field_name} path cannot be empty")
        if isinstance(value, (str, Path)):
            return value
        return Path(cast(PathLike[str], value))
    raise error_type(
        f"dataset build request {field_name} must be a pandas DataFrame or a file "
        "path (str/pathlib.Path)"
    )


__all__ = ["validate_dataframe_or_path"]
