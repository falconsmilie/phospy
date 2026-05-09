"""Internal validator for dataset build requests."""

from __future__ import annotations

from phospy.api.configs import (
    DATASET_COMPARISON_BUILDING_POLICY_SAMPLE_METADATA_PAIRS,
    DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_NONE,
)
from phospy.api.requests import DatasetBuildRequest
from phospy.errors.input import PhosPyInputError
from phospy.references.models import Organism
from phospy.transformations.models import QuantitativeMeaning
from phospy.validation.datasets.inputs import DatasetInputSourceValidator
from phospy.validation.datasets.preprocessing import DatasetPreprocessingConfigValidator


class DatasetBuildRequestValidator:
    """Validate the supported subset of `DatasetBuildRequest`."""

    def __init__(
        self,
        *,
        source_validator: DatasetInputSourceValidator | None = None,
        preprocessing_validator: DatasetPreprocessingConfigValidator | None = None,
    ) -> None:
        self._source_validator = source_validator or DatasetInputSourceValidator()
        self._preprocessing_validator = (
            preprocessing_validator or DatasetPreprocessingConfigValidator()
        )

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
        _validate_quantitative_meaning(request.quantitative_meaning)
        self._preprocessing_validator.run(request.preprocessing_config)
        requested_total_policy = (
            request.preprocessing_config.total_protein_correction.policy
        )
        resolved_total_policy = requested_total_policy
        if (
            resolved_total_policy != DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_NONE
            and request.total is None
        ):
            raise PhosPyInputError(
                "dataset build request preprocessing_config.total_protein_correction."
                f"policy={requested_total_policy!r} requires total input data"
            )
        if (
            request.preprocessing_config.comparisons.policy
            == DATASET_COMPARISON_BUILDING_POLICY_SAMPLE_METADATA_PAIRS
            and request.sample_metadata is None
        ):
            raise PhosPyInputError(
                "dataset build request preprocessing_config.comparisons."
                "policy='sample_metadata_pairs' requires sample_metadata input data"
            )
        return request


def _validate_quantitative_meaning(
    quantitative_meaning: QuantitativeMeaning | str | None,
) -> None:
    if quantitative_meaning is None or isinstance(
        quantitative_meaning, QuantitativeMeaning
    ):
        return
    try:
        QuantitativeMeaning(str(quantitative_meaning))
    except ValueError as exc:
        supported = ", ".join(member.value for member in QuantitativeMeaning)
        raise PhosPyInputError(
            f"dataset build request quantitative_meaning must be one of: {supported}"
        ) from exc
