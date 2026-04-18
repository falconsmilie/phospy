"""Internal validator for dataset build requests."""

from __future__ import annotations

from os import PathLike
from pathlib import Path

import pandas as pd

from phospy.api.requests import DatasetBuildRequest
from phospy.errors.input import PhosPyInputError, UnsupportedInputFormatError


class DatasetBuildRequestValidator:
    """Validate the supported subset of `DatasetBuildRequest`."""

    def run(self, request: DatasetBuildRequest) -> DatasetBuildRequest:
        if not isinstance(request, DatasetBuildRequest):
            raise PhosPyInputError("builder input must be a DatasetBuildRequest")
        self._require_supported_source(request.phospho, field_name="phospho")
        self._require_supported_source(
            request.site_metadata, field_name="site_metadata"
        )
        self._require_optional_supported_source(
            request.sample_metadata, field_name="sample_metadata"
        )
        self._require_optional_supported_source(request.total, field_name="total")
        return request

    @staticmethod
    def _require_supported_source(value: object, *, field_name: str) -> None:
        if isinstance(value, str | Path | PathLike) and not str(value).strip():
            raise UnsupportedInputFormatError(
                f"dataset build request {field_name} path cannot be empty"
            )
        if isinstance(value, Path | PathLike | str):
            return
        if isinstance(value, pd.DataFrame):
            return
        raise UnsupportedInputFormatError(
            f"dataset build request {field_name} must be a pandas DataFrame or a file "
            "path (str/pathlib.Path)"
        )

    @classmethod
    def _require_optional_supported_source(
        cls, value: object | None, *, field_name: str
    ) -> None:
        if value is None:
            return
        cls._require_supported_source(value, field_name=field_name)
