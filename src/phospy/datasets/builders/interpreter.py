"""Internal interpreter for dataset build requests."""

from __future__ import annotations

import pandas as pd

from phospy.api.requests import DatasetBuildRequest
from phospy.datasets.builders.contracts import InterpretedDatasetBuildRequest
from phospy.errors.input import UnsupportedInputFormatError


class DatasetBuildRequestInterpreter:
    """Resolve validated builder request data into execution inputs."""

    def run(self, request: DatasetBuildRequest) -> InterpretedDatasetBuildRequest:
        if not isinstance(request.phospho, pd.DataFrame):
            raise UnsupportedInputFormatError(
                "phospho must be provided as a pandas DataFrame after validation"
            )
        if not isinstance(request.site_metadata, pd.DataFrame):
            raise UnsupportedInputFormatError(
                "site_metadata must be provided as a pandas DataFrame after validation"
            )
        if request.sample_metadata is not None and not isinstance(
            request.sample_metadata, pd.DataFrame
        ):
            raise UnsupportedInputFormatError(
                "sample_metadata must be provided as a pandas DataFrame after validation"
            )
        if request.total is not None and not isinstance(request.total, pd.DataFrame):
            raise UnsupportedInputFormatError(
                "total must be provided as a pandas DataFrame after validation"
            )
        return InterpretedDatasetBuildRequest(
            phospho=request.phospho,
            site_metadata=request.site_metadata,
            sample_metadata=request.sample_metadata,
            total=request.total,
            organism=request.organism,
            transformation_state=request.transformation_state,
        )
