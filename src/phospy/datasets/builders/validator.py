"""Internal validator for dataset build requests."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from phospy.api.requests import DatasetBuildRequest
from phospy.errors.input import PhosPyInputError, UnsupportedInputFormatError


class DatasetBuildRequestValidator:
    """Validate the supported subset of `DatasetBuildRequest`."""

    def run(self, request: DatasetBuildRequest) -> DatasetBuildRequest:
        if not isinstance(request, DatasetBuildRequest):
            raise PhosPyInputError("builder input must be a DatasetBuildRequest")
        self._require_dataframe(request.phospho, field_name="phospho")
        if request.site_metadata is None:
            raise UnsupportedInputFormatError(
                "site_metadata is required in the current dataset builder path"
            )
        self._require_dataframe(request.site_metadata, field_name="site_metadata")
        self._require_optional_dataframe(
            request.sample_metadata, field_name="sample_metadata"
        )
        self._require_optional_dataframe(request.total, field_name="total")
        return request

    @staticmethod
    def _require_dataframe(value: object, *, field_name: str) -> None:
        if isinstance(value, pd.DataFrame):
            return
        if isinstance(value, (str, Path)):
            raise UnsupportedInputFormatError(
                f"{field_name} file-path inputs are not supported yet in this rewrite phase"
            )
        raise UnsupportedInputFormatError(
            f"{field_name} must be provided as a pandas DataFrame in this rewrite phase"
        )

    @classmethod
    def _require_optional_dataframe(
        cls, value: object | None, *, field_name: str
    ) -> None:
        if value is None:
            return
        cls._require_dataframe(value, field_name=field_name)
