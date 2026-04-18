"""Internal validator for dataset build requests."""

from __future__ import annotations

from phospy.api.requests import DatasetBuildRequest
from phospy.errors.input import PhosPyInputError
from phospy.references.models import Organism
from phospy.transformations.models import TransformationState
from phospy.validation.datasets.inputs import DatasetInputSourceValidator


class DatasetBuildRequestValidator:
    """Validate the supported subset of `DatasetBuildRequest`."""

    def __init__(
        self, *, source_validator: DatasetInputSourceValidator | None = None
    ) -> None:
        self._source_validator = source_validator or DatasetInputSourceValidator()

    def run(self, request: DatasetBuildRequest) -> DatasetBuildRequest:
        if not isinstance(request, DatasetBuildRequest):
            raise PhosPyInputError("builder input must be a DatasetBuildRequest")
        self._source_validator.run(request.phospho, field_name="phospho")
        self._source_validator.run(request.site_metadata, field_name="site_metadata")
        self._source_validator.run(
            request.sample_metadata, field_name="sample_metadata", allow_none=True
        )
        self._source_validator.run(request.total, field_name="total", allow_none=True)
        if request.organism is not None and not isinstance(request.organism, Organism):
            raise PhosPyInputError("dataset build request organism must be an Organism")
        if request.transformation_state is not None and not isinstance(
            request.transformation_state, TransformationState
        ):
            raise PhosPyInputError(
                "dataset build request transformation_state must be a TransformationState"
            )
        return request
