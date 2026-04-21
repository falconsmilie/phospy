"""Internal interpreter for dataset build requests."""

from __future__ import annotations

from phospy.api.configs import DATASET_SITE_MATRIX_POLICY_BUILD_FROM_METADATA
from phospy.api.requests import DatasetBuildRequest
from phospy.datasets.builders.contracts import InterpretedDatasetBuildRequest
from phospy.datasets.builders.normalizer import DatasetConventionNormalizer
from phospy.datasets.builders.reader import DatasetInputReader
from phospy.datasets.builders.sequence_derivation import SiteSequenceDeriver
from phospy.datasets.preprocessing.models import PreprocessingPlan


class DatasetBuildRequestInterpreter:
    """Resolve validated builder request data into execution inputs."""

    def __init__(
        self,
        *,
        reader: DatasetInputReader | None = None,
        normalizer: DatasetConventionNormalizer | None = None,
        site_sequence_deriver: SiteSequenceDeriver | None = None,
    ) -> None:
        self._reader = reader or DatasetInputReader()
        self._normalizer = normalizer or DatasetConventionNormalizer()
        self._site_sequence_deriver = site_sequence_deriver or SiteSequenceDeriver()

    def run(self, request: DatasetBuildRequest) -> InterpretedDatasetBuildRequest:
        phospho = self._reader.run(request.phospho, field_name="phospho")
        site_metadata = self._reader.run(
            request.site_metadata,
            field_name="site_metadata",
        )
        sample_metadata = (
            None
            if request.sample_metadata is None
            else self._reader.run(
                request.sample_metadata,
                field_name="sample_metadata",
            )
        )
        total = (
            None
            if request.total is None
            else self._reader.run(
                request.total,
                field_name="total",
            )
        )
        normalized = self._normalizer.run(
            phospho=phospho,
            site_metadata=site_metadata,
            sample_metadata=sample_metadata,
            total=total,
        )
        preprocessing_plan = PreprocessingPlan.from_config(request.preprocessing_config)
        enriched_site_metadata = self._site_sequence_deriver.run(
            normalized.site_metadata,
            organism=request.organism,
            allow_partial=(
                preprocessing_plan.site_matrix_policy
                == DATASET_SITE_MATRIX_POLICY_BUILD_FROM_METADATA
            ),
        )
        return InterpretedDatasetBuildRequest(
            phospho=normalized.phospho,
            site_metadata=enriched_site_metadata,
            sample_metadata=normalized.sample_metadata,
            total=normalized.total,
            organism=request.organism,
            preprocessing_plan=preprocessing_plan,
        )
