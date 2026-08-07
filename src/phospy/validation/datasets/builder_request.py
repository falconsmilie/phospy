"""Internal validator for dataset build requests."""

from __future__ import annotations

from collections.abc import Sequence
from typing import cast

from phospy.contracts.configs import (
    DATASET_COMPARISON_BUILDING_POLICY_SAMPLE_METADATA_PAIRS,
    DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_NONE,
)
from phospy.contracts.requests import DatasetBuildRequest
from phospy.errors.input import PhosPyInputError
from phospy.science.datasets.preprocessing.correction_output import (
    validate_corrected_preprocessing_output,
)
from phospy.science.evidence.dataset_resolution import (
    DATASET_SITE_RESOLUTION_MODE_SITE_LEVEL_RESOLVED,
    SUPPORTED_DATASET_MULTI_SITE_POLICIES,
    SUPPORTED_DATASET_SITE_RESOLUTION_MODES,
)
from phospy.science.references.models import Organism
from phospy.science.transformations.models import (
    IntensityScaleKind,
    QuantitativeMeaning,
    caller_declarable_quantitative_meaning_values,
)
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

    def run(self, request: object) -> DatasetBuildRequest:
        if not isinstance(request, DatasetBuildRequest):
            raise PhosPyInputError("builder input must be a DatasetBuildRequest")
        site_resolution_mode = _validate_site_resolution_mode(
            request.site_resolution_mode
        )
        if site_resolution_mode == DATASET_SITE_RESOLUTION_MODE_SITE_LEVEL_RESOLVED:
            self._source_validator.run(request.phospho, field_name="phospho")
            self._source_validator.run(
                request.site_metadata, field_name="site_metadata"
            )
            if request.peptide_evidence is not None:
                raise PhosPyInputError(
                    "dataset build request peptide_evidence requires "
                    "site_resolution_mode='peptide_evidence'"
                )
            if request.multi_site_policy is not None:
                raise PhosPyInputError(
                    "dataset build request multi_site_policy is only supported when "
                    "site_resolution_mode='peptide_evidence'"
                )
            if request.peptide_evidence_sample_intensity_columns is not None:
                raise PhosPyInputError(
                    "dataset build request peptide_evidence_sample_intensity_columns "
                    "is only supported when site_resolution_mode='peptide_evidence'"
                )
            if request.peptide_site_mapping is not None:
                raise PhosPyInputError(
                    "dataset build request peptide_site_mapping is only supported "
                    "when site_resolution_mode='peptide_evidence'"
                )
        else:
            if request.phospho is not None or request.site_metadata is not None:
                raise PhosPyInputError(
                    "dataset build request with "
                    "site_resolution_mode='peptide_evidence' must not provide "
                    "phospho/site_metadata inputs"
                )
            self._source_validator.run(
                request.peptide_evidence,
                field_name="peptide_evidence",
            )
            self._source_validator.run(
                request.peptide_site_mapping,
                field_name="peptide_site_mapping",
                allow_none=True,
            )
            if request.multi_site_policy is None:
                raise PhosPyInputError(
                    "dataset build request peptide_evidence input requires "
                    "multi_site_policy"
                )
            _validate_dataset_multi_site_policy(request.multi_site_policy)
            _validate_peptide_sample_columns(
                request.peptide_evidence_sample_intensity_columns
            )
        self._source_validator.run(
            request.sample_metadata, field_name="sample_metadata", allow_none=True
        )
        self._source_validator.run(request.total, field_name="total", allow_none=True)
        organism = cast(object, request.organism)
        if organism is not None and not isinstance(organism, Organism):
            raise PhosPyInputError("dataset build request organism must be an Organism")
        _validate_input_intensity_scale(request.input_intensity_scale)
        allow_suspicious_declared_scale = cast(
            object,
            request.allow_suspicious_declared_input_intensity_scale,
        )
        if not isinstance(allow_suspicious_declared_scale, bool):
            raise PhosPyInputError(
                "dataset build request "
                "allow_suspicious_declared_input_intensity_scale must be a bool"
            )
        _validate_quantitative_meaning(request.quantitative_meaning)
        validate_corrected_preprocessing_output(request.corrected_preprocessing_output)
        allow_opaque_site_values = cast(object, request.allow_opaque_site_values)
        if not isinstance(allow_opaque_site_values, bool):
            raise PhosPyInputError(
                "dataset build request allow_opaque_site_values must be a bool"
            )
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
        if (
            request.preprocessing_config.group_coverage_filter.enabled
            and request.sample_metadata is None
        ):
            raise PhosPyInputError(
                "dataset build request preprocessing_config.group_coverage_filter "
                "requires sample_metadata input data"
            )
        return request


def _validate_site_resolution_mode(site_resolution_mode: object) -> str:
    if (
        isinstance(site_resolution_mode, str)
        and site_resolution_mode in SUPPORTED_DATASET_SITE_RESOLUTION_MODES
    ):
        return site_resolution_mode
    supported = ", ".join(
        repr(value) for value in SUPPORTED_DATASET_SITE_RESOLUTION_MODES
    )
    raise PhosPyInputError(
        f"dataset build request site_resolution_mode must be one of: {supported}"
    )


def _validate_dataset_multi_site_policy(multi_site_policy: object) -> None:
    if multi_site_policy == "keep_joint":
        raise PhosPyInputError(
            "dataset build request multi_site_policy='keep_joint' is no longer "
            "supported for AnalysisReadyDatasetBuilder peptide-evidence requests "
            "because unresolved joint evidence cannot satisfy the strict "
            "site-level identity contract. Use multi_site_policy='split' to "
            "allocate ambiguous evidence to strict site rows, 'reject' to fail "
            "on ambiguous peptide rows, or 'exclude_from_sequence_scoring' to "
            "remove ambiguous rows from the analysis-ready build."
        )
    if (
        isinstance(multi_site_policy, str)
        and multi_site_policy in SUPPORTED_DATASET_MULTI_SITE_POLICIES
    ):
        return
    supported = ", ".join(
        repr(value) for value in SUPPORTED_DATASET_MULTI_SITE_POLICIES
    )
    raise PhosPyInputError(
        f"dataset build request multi_site_policy must be one of: {supported}"
    )


def _validate_peptide_sample_columns(value: object) -> None:
    if value is None:
        raise PhosPyInputError(
            "dataset build request peptide_evidence_sample_intensity_columns "
            "is required when site_resolution_mode='peptide_evidence'"
        )
    if isinstance(value, str) or not isinstance(value, Sequence):
        raise PhosPyInputError(
            "dataset build request peptide_evidence_sample_intensity_columns must be "
            "a sequence of non-empty string column names"
        )
    columns = tuple(cast(Sequence[object], value))
    if not columns:
        raise PhosPyInputError(
            "dataset build request peptide_evidence_sample_intensity_columns must "
            "contain at least one column name"
        )
    try:
        has_duplicates = len(set(columns)) != len(columns)
    except TypeError as exc:
        raise PhosPyInputError(
            "dataset build request peptide_evidence_sample_intensity_columns must "
            "contain hashable string column names"
        ) from exc
    if has_duplicates:
        raise PhosPyInputError(
            "dataset build request peptide_evidence_sample_intensity_columns must "
            "not contain duplicate column names"
        )
    for column in columns:
        if not isinstance(column, str) or not column.strip():
            raise PhosPyInputError(
                "dataset build request peptide_evidence_sample_intensity_columns "
                "must contain non-empty string column names"
            )


def _validate_quantitative_meaning(
    quantitative_meaning: QuantitativeMeaning | str | None,
) -> None:
    if quantitative_meaning is None:
        return
    try:
        resolved = (
            quantitative_meaning
            if isinstance(quantitative_meaning, QuantitativeMeaning)
            else QuantitativeMeaning(str(quantitative_meaning))
        )
    except ValueError as exc:
        supported = ", ".join(member.value for member in QuantitativeMeaning)
        raise PhosPyInputError(
            f"dataset build request quantitative_meaning must be one of: {supported}"
        ) from exc
    allowed = caller_declarable_quantitative_meaning_values()
    if resolved.value not in allowed:
        raise PhosPyInputError(
            "dataset build request quantitative_meaning may only declare direct "
            "input meanings: " + ", ".join(allowed) + f"; got {resolved.value!r}"
        )


def _validate_input_intensity_scale(
    input_intensity_scale: IntensityScaleKind | str | None,
) -> None:
    if input_intensity_scale is None:
        return
    if isinstance(input_intensity_scale, IntensityScaleKind):
        return
    try:
        IntensityScaleKind(str(input_intensity_scale))
    except ValueError as exc:
        supported = ", ".join(member.value for member in IntensityScaleKind)
        raise PhosPyInputError(
            f"dataset build request input_intensity_scale must be one of: {supported}"
        ) from exc
