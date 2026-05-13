"""Internal interpreter for dataset build requests."""

from __future__ import annotations

from typing import NoReturn

from phospy.api.requests import DatasetBuildRequest
from phospy.datasets.builders.contracts import (
    DatasetInput,
    InterpretedDatasetBuildRequest,
)
from phospy.datasets.builders.normalizer import DatasetConventionNormalizer
from phospy.datasets.builders.reader import DatasetInputReader
from phospy.datasets.builders.sequence_derivation import SiteSequenceDeriver
from phospy.datasets.preprocessing.models import PreprocessingPlan
from phospy.datasets.preprocessing.policy_models import SiteMatrixPolicy
from phospy.errors.input import PhosPyInputError
from phospy.evidence.dataset_resolution import (
    DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE,
    DATASET_SITE_RESOLUTION_MODE_SITE_LEVEL_RESOLVED,
    PeptideEvidenceDatasetResolver,
    build_multi_site_handling_config_for_dataset_policy,
)
from phospy.evidence.models import PeptideEvidenceTable
from phospy.transformations.models import IntensityScaleKind, QuantitativeMeaning


class DatasetBuildRequestInterpreter:
    """Resolve validated builder request data into execution inputs."""

    def __init__(
        self,
        *,
        reader: DatasetInputReader | None = None,
        normalizer: DatasetConventionNormalizer | None = None,
        site_sequence_deriver: SiteSequenceDeriver | None = None,
        peptide_evidence_resolver: PeptideEvidenceDatasetResolver | None = None,
    ) -> None:
        self._reader = reader or DatasetInputReader()
        self._normalizer = normalizer or DatasetConventionNormalizer()
        self._site_sequence_deriver = site_sequence_deriver or SiteSequenceDeriver()
        self._peptide_evidence_resolver = (
            peptide_evidence_resolver or PeptideEvidenceDatasetResolver()
        )

    def run(self, request: DatasetBuildRequest) -> InterpretedDatasetBuildRequest:
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
        site_resolution_mode = str(request.site_resolution_mode).strip()
        peptide_evidence_resolution_payload: dict[str, object] | None = None
        multi_site_policy: str | None = None
        if site_resolution_mode == DATASET_SITE_RESOLUTION_MODE_SITE_LEVEL_RESOLVED:
            phospho = self._reader.run(
                _require_dataset_input(
                    request.phospho,
                    field_name="dataset build request phospho",
                ),
                field_name="phospho",
            )
            site_metadata = self._reader.run(
                _require_dataset_input(
                    request.site_metadata,
                    field_name="dataset build request site_metadata",
                ),
                field_name="site_metadata",
            )
        elif site_resolution_mode == DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE:
            peptide_evidence = self._reader.run(
                _require_dataset_input(
                    request.peptide_evidence,
                    field_name="dataset build request peptide_evidence",
                ),
                field_name="peptide_evidence",
            )
            peptide_site_mapping = (
                None
                if request.peptide_site_mapping is None
                else self._reader.run(
                    request.peptide_site_mapping,
                    field_name="peptide_site_mapping",
                )
            )
            try:
                evidence = PeptideEvidenceTable(
                    frame=peptide_evidence,
                    sample_intensity_columns=(
                        ()
                        if request.peptide_evidence_sample_intensity_columns is None
                        else request.peptide_evidence_sample_intensity_columns
                    ),
                    site_mapping=peptide_site_mapping,
                    multi_site_handling_config=build_multi_site_handling_config_for_dataset_policy(
                        multi_site_policy=str(request.multi_site_policy)
                    ),
                )
                resolved = self._peptide_evidence_resolver.run(
                    evidence=evidence,
                    multi_site_policy=str(request.multi_site_policy),
                )
            except (TypeError, ValueError, KeyError, PhosPyInputError) as exc:
                self._raise_wrapped_input_error(
                    stage_name="dataset_builder.peptide_evidence_resolution",
                    field_name="dataset build request peptide_evidence",
                    operation=(
                        "resolving peptide evidence into site-level phospho and "
                        "site_metadata tables"
                    ),
                    next_action=(
                        "provide peptide_evidence with required evidence columns, "
                        "valid peptide_evidence_sample_intensity_columns, and a "
                        "compatible multi_site_policy/site mapping"
                    ),
                    original_error=exc,
                )
            phospho = resolved.phospho
            site_metadata = resolved.site_metadata
            multi_site_policy = str(request.multi_site_policy)
            peptide_evidence_resolution_payload = resolved.summary.to_payload()
        else:  # pragma: no cover - validator owns this branch; keep defensive.
            raise PhosPyInputError(
                "dataset build request site_resolution_mode is unsupported after "
                "validation"
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
                    is SiteMatrixPolicy.BUILD_FROM_METADATA
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
        site_sequence_derivation_payload = _resolve_site_sequence_derivation_payload(
            self._site_sequence_deriver
        )
        return InterpretedDatasetBuildRequest(
            phospho=normalized.phospho,
            site_metadata=enriched_site_metadata,
            sample_metadata=normalized.sample_metadata,
            total=normalized.total,
            organism=request.organism,
            preprocessing_plan=preprocessing_plan,
            declared_input_intensity_scale_kind=_resolve_input_intensity_scale_kind(
                request.input_intensity_scale
            ),
            declared_input_intensity_scale_source=(
                "dataset_build_request.input_intensity_scale"
                if request.input_intensity_scale is not None
                else None
            ),
            site_identifier_normalisation=normalized.site_identifier_normalisation,
            site_sequence_derivation=site_sequence_derivation_payload,
            quantitative_meaning=_resolve_quantitative_meaning(
                request.quantitative_meaning
            ),
            site_resolution_mode=site_resolution_mode,
            multi_site_policy=multi_site_policy,
            peptide_evidence_resolution=peptide_evidence_resolution_payload,
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


def _resolve_quantitative_meaning(
    quantitative_meaning: QuantitativeMeaning | str | None,
) -> QuantitativeMeaning | None:
    if quantitative_meaning is None:
        return None
    if isinstance(quantitative_meaning, QuantitativeMeaning):
        return quantitative_meaning
    try:
        return QuantitativeMeaning(str(quantitative_meaning))
    except ValueError as exc:
        supported = ", ".join(member.value for member in QuantitativeMeaning)
        raise PhosPyInputError(
            f"dataset build request quantitative_meaning must be one of: {supported}"
        ) from exc


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


def _resolve_site_sequence_derivation_payload(
    site_sequence_deriver: object,
) -> dict[str, object] | None:
    report = getattr(site_sequence_deriver, "last_report", None)
    if report is None:
        return None
    to_payload = getattr(report, "to_payload", None)
    if callable(to_payload):
        payload = to_payload()
        if isinstance(payload, dict):
            return payload
    return None


def _require_dataset_input(
    value: DatasetInput | None, *, field_name: str
) -> DatasetInput:
    if value is None:
        raise PhosPyInputError(f"{field_name} is required")
    return value
