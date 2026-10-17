"""Dataset-build request input source validation."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from phospy.errors.input import UnsupportedInputFormatError
from phospy.validation.common.input_sources import validate_dataframe_or_path


class DatasetInputSourceValidator:
    """Validate supported source types for dataset build inputs."""

    def run(
        self,
        value: object | None,
        *,
        field_name: str,
        allow_none: bool = False,
    ) -> pd.DataFrame | str | Path | None:
        return validate_dataframe_or_path(
            value,
            field_name=field_name,
            allow_none=allow_none,
            error_type=UnsupportedInputFormatError,
        )
