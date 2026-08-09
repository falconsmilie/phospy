"""Internal interpreter for dataset build requests."""

from __future__ import annotations

from typing import cast

from phospy.errors.input import PhosPyInputError
from phospy.science.datasets.builders.contracts import (
    DatasetBuildRequestProtocol,
    InterpretedDatasetBuildRequest,
)
from phospy.science.datasets.builders.interpreter_collaborators import (
    DatasetBuildPreprocessingPlanner,
    DatasetBuildSiteSequenceResolver,
    DatasetBuildSourceResolver,
)
from phospy.science.datasets.builders.normalizer import DatasetConventionNormalizer
from phospy.science.datasets.builders.reader import DatasetInputReader
from phospy.science.datasets.builders.sequence_derivation import SiteSequenceDeriver
from phospy.science.datasets.preprocessing.correction_output import (
    CorrectedPreprocessingOutput,
)
from phospy.science.evidence.dataset_resolution.resolver import (
    PeptideEvidenceDatasetResolver,
)
from phospy.science.transformations.models import (
    IntensityScaleKind,
    QuantitativeMeaning,
    caller_declarable_quantitative_meaning_values,
)


class DatasetBuildRequestInterpreter:
    """Resolve validated builder request data into execution inputs."""

    def __init__(
        self,
        *,
        source_resolver: DatasetBuildSourceResolver | None = None,
        preprocessing_planner: DatasetBuildPreprocessingPlanner | None = None,
        site_sequence_resolver: DatasetBuildSiteSequenceResolver | None = None,
        reader: DatasetInputReader | None = None,
        normalizer: DatasetConventionNormalizer | None = None,
        site_sequence_deriver: SiteSequenceDeriver | None = None,
        peptide_evidence_resolver: PeptideEvidenceDatasetResolver | None = None,
    ) -> None:
        self._source_resolver = source_resolver or DatasetBuildSourceResolver(
            reader=reader,
            normalizer=normalizer,
            peptide_evidence_resolver=peptide_evidence_resolver,
        )
        self._preprocessing_planner = (
            preprocessing_planner or DatasetBuildPreprocessingPlanner()
        )
        self._site_sequence_resolver = site_sequence_resolver or (
            DatasetBuildSiteSequenceResolver(
                site_sequence_deriver=site_sequence_deriver
            )
        )

    def run(
        self,
        request: DatasetBuildRequestProtocol,
    ) -> InterpretedDatasetBuildRequest:
        declared_input_intensity_scale_kind = _resolve_input_intensity_scale_kind(
            request.input_intensity_scale
        )
        resolved_sources = self._source_resolver.run(
            request,
            declared_input_intensity_scale_kind=declared_input_intensity_scale_kind,
        )
        preprocessing_plan = self._preprocessing_planner.run(
            request,
            phospho=resolved_sources.phospho,
            sample_metadata=resolved_sources.sample_metadata,
            declared_input_intensity_scale_kind=(declared_input_intensity_scale_kind),
        )
        sequence_resolution = self._site_sequence_resolver.run(
            site_metadata=resolved_sources.site_metadata,
            organism=request.organism,
            preprocessing_plan=preprocessing_plan,
        )
        return InterpretedDatasetBuildRequest(
            phospho=resolved_sources.phospho,
            site_metadata=sequence_resolution.site_metadata,
            sample_metadata=resolved_sources.sample_metadata,
            total=resolved_sources.total,
            organism=request.organism,
            preprocessing_plan=preprocessing_plan,
            declared_input_intensity_scale_kind=declared_input_intensity_scale_kind,
            declared_input_intensity_scale_source=(
                "dataset_build_request.input_intensity_scale"
                if request.input_intensity_scale is not None
                else None
            ),
            site_identifier_normalisation=resolved_sources.site_identifier_normalisation,
            site_sequence_derivation=sequence_resolution.site_sequence_derivation,
            quantitative_meaning=_resolve_quantitative_meaning(
                request.quantitative_meaning
            ),
            site_resolution_mode=resolved_sources.site_resolution_mode,
            multi_site_policy=resolved_sources.multi_site_policy,
            allow_opaque_site_values=bool(request.allow_opaque_site_values),
            allow_suspicious_declared_input_intensity_scale=(
                request.allow_suspicious_declared_input_intensity_scale
            ),
            peptide_evidence_resolution=resolved_sources.peptide_evidence_resolution,
            corrected_preprocessing_output=cast(
                CorrectedPreprocessingOutput | None,
                request.corrected_preprocessing_output,
            ),
        )


def _resolve_quantitative_meaning(
    quantitative_meaning: QuantitativeMeaning | str | None,
) -> QuantitativeMeaning | None:
    if quantitative_meaning is None:
        return None
    if isinstance(quantitative_meaning, QuantitativeMeaning):
        resolved = quantitative_meaning
    else:
        try:
            resolved = QuantitativeMeaning(str(quantitative_meaning))
        except ValueError as exc:
            supported = ", ".join(member.value for member in QuantitativeMeaning)
            raise PhosPyInputError(
                "dataset build request quantitative_meaning must be one of: "
                f"{supported}"
            ) from exc
    allowed = caller_declarable_quantitative_meaning_values()
    if resolved.value not in allowed:
        raise PhosPyInputError(
            "dataset build request quantitative_meaning may only declare direct "
            "input meanings: " + ", ".join(allowed) + f"; got {resolved.value!r}"
        )
    return resolved


def _resolve_input_intensity_scale_kind(
    input_intensity_scale: IntensityScaleKind | str | None,
) -> IntensityScaleKind | None:
    if input_intensity_scale is None:
        return None
    if isinstance(input_intensity_scale, IntensityScaleKind):
        return input_intensity_scale
    try:
        return IntensityScaleKind(str(input_intensity_scale))
    except ValueError as exc:
        supported = ", ".join(member.value for member in IntensityScaleKind)
        raise PhosPyInputError(
            f"dataset build request input_intensity_scale must be one of: {supported}"
        ) from exc
