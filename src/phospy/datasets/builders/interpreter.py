"""Internal interpreter for dataset build requests."""

from __future__ import annotations

from typing import NoReturn

from phospy.api.configs import DATASET_SITE_MATRIX_POLICY_BUILD_FROM_METADATA
from phospy.api.requests import DatasetBuildRequest
from phospy.datasets.builders.contracts import InterpretedDatasetBuildRequest
from phospy.datasets.builders.normalizer import DatasetConventionNormalizer
from phospy.datasets.builders.reader import DatasetInputReader
from phospy.datasets.builders.sequence_derivation import SiteSequenceDeriver
from phospy.datasets.preprocessing.models import PreprocessingPlan
from phospy.errors.input import PhosPyInputError


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
        try:
            normalized = self._normalizer.run(
                phospho=phospho,
                site_metadata=site_metadata,
                sample_metadata=sample_metadata,
                total=total,
            )
        except (TypeError, ValueError, KeyError) as exc:
            self._raise_wrapped_input_error(
                stage_name="dataset_builder.normalization",
                field_name="dataset build request input tables",
                operation="normalizing input indices and metadata column conventions",
                next_action=(
                    "ensure phospho/site_metadata/sample_metadata/total tables use "
                    "supported rectangular DataFrame shapes, canonical site labels, "
                    "and non-conflicting metadata column conventions"
                ),
                original_error=exc,
            )
        preprocessing_plan = PreprocessingPlan.from_config(request.preprocessing_config)
        defer_missing_site_sequence_fill = (
            self._should_defer_missing_site_sequence_fill(preprocessing_plan)
        )
        try:
            enriched_site_metadata = self._site_sequence_deriver.run(
                normalized.site_metadata,
                organism=request.organism,
                allow_partial=(
                    preprocessing_plan.site_matrix_policy
                    == DATASET_SITE_MATRIX_POLICY_BUILD_FROM_METADATA
                    or preprocessing_plan.site_sequence_resolution_enabled
                ),
                derive_missing_from_reference=not defer_missing_site_sequence_fill,
            )
        except (TypeError, ValueError, KeyError) as exc:
            self._raise_wrapped_input_error(
                stage_name="dataset_builder.site_sequence_derivation",
                field_name="dataset build request site_metadata",
                operation="deriving and validating site_sequence values",
                next_action=(
                    "ensure site_metadata contains canonical non-empty gene/site "
                    "fields and supported site identifiers for sequence derivation"
                ),
                original_error=exc,
            )
        return InterpretedDatasetBuildRequest(
            phospho=normalized.phospho,
            site_metadata=enriched_site_metadata,
            sample_metadata=normalized.sample_metadata,
            total=normalized.total,
            organism=request.organism,
            preprocessing_plan=preprocessing_plan,
        )

    @staticmethod
    def _should_defer_missing_site_sequence_fill(plan: PreprocessingPlan) -> bool:
        return bool(plan.site_sequence_resolution_enabled)

    @staticmethod
    def _raise_wrapped_input_error(
        *,
        stage_name: str,
        field_name: str,
        operation: str,
        next_action: str,
        original_error: Exception,
    ) -> NoReturn:
        original_message = " ".join(str(original_error).split())
        raise PhosPyInputError(
            f"{stage_name} failed while {operation} for {field_name}. "
            f"Original error: {type(original_error).__name__}: {original_message}. "
            f"Next action: {next_action}"
        ) from original_error
