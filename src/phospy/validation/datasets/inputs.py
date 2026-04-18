"""Dataset-build input source validation."""

from __future__ import annotations

from os import PathLike
from pathlib import Path

import pandas as pd

from phospy.errors.input import UnsupportedInputFormatError


class DatasetInputSourceValidator:
    """Validate supported source types for dataset build inputs."""

    def run(
        self,
        value: object | None,
        *,
        field_name: str,
        allow_none: bool = False,
    ) -> object | None:
        if value is None:
            if allow_none:
                return None
            raise UnsupportedInputFormatError(
                f"dataset build request {field_name} is required"
            )
        if isinstance(value, pd.DataFrame):
            return value
        if isinstance(value, (str, Path, PathLike)):
            if isinstance(value, str) and not value.strip():
                raise UnsupportedInputFormatError(
                    f"dataset build request {field_name} path cannot be empty"
                )
            return value
        raise UnsupportedInputFormatError(
            f"dataset build request {field_name} must be a pandas DataFrame or a file "
            "path (str/pathlib.Path)"
        )
