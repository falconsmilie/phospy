"""Validation for enrichment workflow requests."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from numbers import Integral
from typing import Literal, cast

from phospy.contracts.configs import EnrichmentConfig
from phospy.contracts.enrichment_identifier_sets import (
    EnrichmentIdentifierSetProvenance,
    EnrichmentIdentifierSetSourceType,
)
from phospy.contracts.requests import EnrichmentWorkflowRequest
from phospy.errors.validation import WorkflowValidationError
from phospy.provenance.models import InputIntensityScaleEvidence
from phospy.science.enrichment.models import (
    ENRICHMENT_COLLECTION_KIND_GENE_SET,
    ENRICHMENT_COLLECTION_KIND_PTM_SET,
    GENE_LEVEL_ENRICHMENT_IDENTIFIER_KINDS,
    PTM_LEVEL_ENRICHMENT_IDENTIFIER_KINDS,
    SUPPORTED_ENRICHMENT_IDENTIFIER_KINDS,
    SUPPORTED_ENRICHMENT_METHODS,
    SUPPORTED_MULTIPLE_TESTING_CORRECTIONS,
    EnrichmentCollectionKind,
    EnrichmentIdentifierKind,
    EnrichmentSet,
    EnrichmentSetCollection,
)
from phospy.science.transformations.models import (
    IntensityScaleEvidenceLevel,
    IntensityScaleKind,
)
from phospy.validation.common.dataframes import require_columns, require_dataframe

EnrichmentSelectedIdentifierSource = Literal["selected_identifiers", "input_table"]


@dataclass(frozen=True, slots=True)
class ValidatedEnrichmentWorkflowRequest:
    """Validated enrichment request passed to interpretation/execution."""

    request: EnrichmentWorkflowRequest
    identifier_column: str
    identifier_kind: EnrichmentIdentifierKind
    set_collection: EnrichmentSetCollection
    background_universe: tuple[str, ...]
    selected_identifiers: tuple[str, ...]
    config: EnrichmentConfig
    selected_identifier_source: EnrichmentSelectedIdentifierSource
    selected_identifier_provenance: EnrichmentIdentifierSetProvenance | None
    background_identifier_provenance: EnrichmentIdentifierSetProvenance | None


class EnrichmentWorkflowValidator:
    """Validate `EnrichmentWorkflowRequest` before interpretation or execution."""

    def run(self, request: object) -> ValidatedEnrichmentWorkflowRequest:
        if not isinstance(request, EnrichmentWorkflowRequest):
            raise WorkflowValidationError(
                "enrichment workflow input must be an EnrichmentWorkflowRequest"
            )

        identifier_column = _require_non_empty_string(
            request.identifier_column,
            field_name="enrichment workflow request identifier_column",
        )
        identifier_kind = _require_identifier_kind(
            request.identifier_kind,
            field_name="enrichment workflow request identifier_kind",
        )
        set_collection = _validate_set_collection(
            request.set_collection,
            identifier_kind=identifier_kind,
        )
        if (
            request.background_identifier_provenance is not None
            and request.background_universe is None
        ):
            raise WorkflowValidationError(
                "Background identifier-set provenance is invalid: "
                "background_identifier_provenance requires explicit "
                "background_universe"
            )
        background_universe = _normalise_identifier_sequence(
            request.background_universe,
            field_name="enrichment workflow request background_universe",
            allow_empty=False,
        )
        config = _validate_config(request.config)
        selected_identifiers, selected_identifier_source = (
            _resolve_selected_identifiers(
                request=request,
                identifier_column=identifier_column,
            )
        )
        selected_identifier_provenance = _validate_identifier_set_provenance(
            request.selected_identifier_provenance,
            role="Selected",
            field_name="selected_identifier_provenance",
            normalized_identifier_count=len(selected_identifiers),
        )
        background_identifier_provenance = _validate_identifier_set_provenance(
            request.background_identifier_provenance,
            role="Background",
            field_name="background_identifier_provenance",
            normalized_identifier_count=len(background_universe),
        )

        return ValidatedEnrichmentWorkflowRequest(
            request=request,
            identifier_column=identifier_column,
            identifier_kind=identifier_kind,
            set_collection=set_collection,
            background_universe=background_universe,
            selected_identifiers=selected_identifiers,
            config=config,
            selected_identifier_source=selected_identifier_source,
            selected_identifier_provenance=selected_identifier_provenance,
            background_identifier_provenance=background_identifier_provenance,
        )


def _validate_set_collection(
    value: object,
    *,
    identifier_kind: EnrichmentIdentifierKind,
) -> EnrichmentSetCollection:
    if not isinstance(value, EnrichmentSetCollection):
        raise WorkflowValidationError(
            "enrichment workflow request set_collection must be "
            "EnrichmentSetCollection, GeneSetCollection, or PtmSetCollection"
        )
    if value.identifier_kind != identifier_kind:
        raise WorkflowValidationError(
            "enrichment workflow request set_collection.identifier_kind must match "
            "request identifier_kind; "
            f"observed={value.identifier_kind!r}, expected={identifier_kind!r}"
        )

    expected_collection_kind = _expected_collection_kind(identifier_kind)
    if value.collection_kind != expected_collection_kind:
        raise WorkflowValidationError(
            "enrichment workflow request set_collection.collection_kind must match "
            "request identifier_kind; "
            f"observed={value.collection_kind!r}, "
            f"expected={expected_collection_kind!r}"
        )

    enrichment_sets = tuple(value.enrichment_sets)
    if not enrichment_sets:
        raise WorkflowValidationError(
            "enrichment workflow request set_collection must contain at least one set"
        )
    observed_set_ids: set[str] = set()
    for enrichment_set in enrichment_sets:
        _validate_enrichment_set(
            enrichment_set,
            collection_identifier_kind=value.identifier_kind,
            observed_set_ids=observed_set_ids,
        )
    return value


def _validate_enrichment_set(
    value: object,
    *,
    collection_identifier_kind: EnrichmentIdentifierKind,
    observed_set_ids: set[str],
) -> None:
    if not isinstance(value, EnrichmentSet):
        raise WorkflowValidationError(
            "enrichment workflow request set_collection must contain "
            "EnrichmentSet values"
        )
    set_id = _require_non_empty_string(
        value.set_id,
        field_name="enrichment workflow request set_collection.set_id",
    )
    if set_id in observed_set_ids:
        raise WorkflowValidationError(
            "enrichment workflow request set_collection set_id values must be unique"
        )
    observed_set_ids.add(set_id)
    _require_non_empty_string(
        value.name,
        field_name=f"enrichment workflow request set_collection[{set_id!r}].name",
    )
    if value.identifier_kind != collection_identifier_kind:
        raise WorkflowValidationError(
            "enrichment workflow request set_collection cannot mix "
            "identifier_kind values; "
            f"set_id={set_id!r}, observed={value.identifier_kind!r}, "
            f"expected={collection_identifier_kind!r}"
        )
    _normalise_identifier_sequence(
        value.identifiers,
        field_name=(
            f"enrichment workflow request set_collection[{set_id!r}].identifiers"
        ),
        allow_empty=False,
    )


def _resolve_selected_identifiers(
    *,
    request: EnrichmentWorkflowRequest,
    identifier_column: str,
) -> tuple[tuple[str, ...], EnrichmentSelectedIdentifierSource]:
    has_input_table = request.input_table is not None
    has_selected_identifiers = request.selected_identifiers is not None
    if has_input_table == has_selected_identifiers:
        raise WorkflowValidationError(
            "enrichment workflow request requires exactly one of input_table or "
            "selected_identifiers"
        )
    if has_selected_identifiers:
        return (
            _normalise_identifier_sequence(
                request.selected_identifiers,
                field_name="enrichment workflow request selected_identifiers",
                allow_empty=False,
            ),
            "selected_identifiers",
        )

    input_table = require_dataframe(
        request.input_table,
        field_name="enrichment workflow request input_table",
        allow_empty=False,
        error_type=WorkflowValidationError,
    )
    require_columns(
        input_table,
        field_name="enrichment workflow request input_table",
        required_columns=(identifier_column,),
        error_type=WorkflowValidationError,
    )
    return (
        _normalise_identifier_sequence(
            tuple(input_table.loc[:, identifier_column].tolist()),
            field_name=(
                f"enrichment workflow request input_table[{identifier_column!r}]"
            ),
            allow_empty=False,
        ),
        "input_table",
    )


def _validate_config(value: object) -> EnrichmentConfig:
    if not isinstance(value, EnrichmentConfig):
        raise WorkflowValidationError(
            "enrichment workflow request config must be EnrichmentConfig"
        )
    if value.method not in SUPPORTED_ENRICHMENT_METHODS:
        supported = ", ".join(repr(method) for method in SUPPORTED_ENRICHMENT_METHODS)
        raise WorkflowValidationError(
            "enrichment workflow request config.method must be one of: " + supported
        )
    if value.multiple_testing_correction not in SUPPORTED_MULTIPLE_TESTING_CORRECTIONS:
        supported = ", ".join(
            repr(method) for method in SUPPORTED_MULTIPLE_TESTING_CORRECTIONS
        )
        raise WorkflowValidationError(
            "enrichment workflow request config.multiple_testing_correction must be "
            "one of: " + supported
        )
    min_set_size = _validate_optional_set_size_threshold(
        value.min_set_size,
        field_name="enrichment workflow request config.min_set_size",
    )
    max_set_size = _validate_optional_set_size_threshold(
        value.max_set_size,
        field_name="enrichment workflow request config.max_set_size",
    )
    if (
        min_set_size is not None
        and max_set_size is not None
        and min_set_size > max_set_size
    ):
        raise WorkflowValidationError(
            "enrichment workflow request config.min_set_size must be less than "
            "or equal to config.max_set_size"
        )
    return value


def _validate_identifier_set_provenance(
    value: object | None,
    *,
    role: Literal["Selected", "Background"],
    field_name: str,
    normalized_identifier_count: int,
) -> EnrichmentIdentifierSetProvenance | None:
    if value is None:
        return None
    if not isinstance(value, EnrichmentIdentifierSetProvenance):
        raise WorkflowValidationError(
            f"{role} identifier-set provenance is invalid: {field_name} must be "
            "EnrichmentIdentifierSetProvenance"
        )
    source_type = _validate_identifier_set_source_type(
        value.source_type,
        role=role,
    )
    source_label = _require_provenance_text(
        value.source_label,
        role=role,
        field_name="source_label",
    )
    identifier_count = _validate_identifier_count(
        value.identifier_count,
        role=role,
    )
    if identifier_count != normalized_identifier_count:
        raise WorkflowValidationError(
            f"{role} identifier-set provenance count mismatch: "
            f"declared identifier_count={identifier_count}, "
            "normalized request identifier count="
            f"{normalized_identifier_count}."
        )
    upstream_workflow_id = _validate_optional_provenance_text(
        value.upstream_workflow_id,
        role=role,
        field_name="upstream_workflow_id",
    )
    upstream_result_id = _validate_optional_provenance_text(
        value.upstream_result_id,
        role=role,
        field_name="upstream_result_id",
    )
    evidence = value.input_intensity_scale_evidence
    if source_type is EnrichmentIdentifierSetSourceType.PHOSPY_DERIVED_QUANTITATIVE:
        evidence = _validate_input_intensity_scale_evidence(
            evidence,
            role=role,
        )
    elif evidence is not None:
        raise WorkflowValidationError(
            f"{role} identifier-set provenance is invalid: "
            "input_intensity_scale_evidence is only valid when "
            "source_type='phospy_derived_quantitative'; "
            f"observed source_type={source_type.value!r}."
        )
    return EnrichmentIdentifierSetProvenance(
        source_type=source_type,
        source_label=source_label,
        identifier_count=identifier_count,
        upstream_workflow_id=upstream_workflow_id,
        upstream_result_id=upstream_result_id,
        input_intensity_scale_evidence=evidence,
    )


def _validate_identifier_set_source_type(
    value: object,
    *,
    role: Literal["Selected", "Background"],
) -> EnrichmentIdentifierSetSourceType:
    if isinstance(value, EnrichmentIdentifierSetSourceType):
        return value
    try:
        return EnrichmentIdentifierSetSourceType(str(value).strip())
    except ValueError as exc:
        supported = ", ".join(item.value for item in EnrichmentIdentifierSetSourceType)
        raise WorkflowValidationError(
            f"{role} identifier-set provenance is invalid: source_type must be "
            f"one of: {supported}; observed={value!r}"
        ) from exc


def _validate_identifier_count(
    value: object,
    *,
    role: Literal["Selected", "Background"],
) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise WorkflowValidationError(
            f"{role} identifier-set provenance is invalid: identifier_count must "
            "be a non-negative integer"
        )
    count = int(value)
    if count < 0:
        raise WorkflowValidationError(
            f"{role} identifier-set provenance is invalid: identifier_count must "
            "be nonnegative"
        )
    return count


def _validate_input_intensity_scale_evidence(
    value: object | None,
    *,
    role: Literal["Selected", "Background"],
) -> InputIntensityScaleEvidence:
    if value is None:
        raise WorkflowValidationError(
            f"{role} identifier-set provenance is invalid: "
            "source_type='phospy_derived_quantitative' requires "
            "input_intensity_scale_evidence."
        )
    if not isinstance(value, InputIntensityScaleEvidence):
        raise WorkflowValidationError(
            f"{role} identifier-set provenance is invalid: "
            "input_intensity_scale_evidence must be InputIntensityScaleEvidence"
        )
    input_intensity_scale = _validate_intensity_scale(
        value.input_intensity_scale,
        role=role,
    )
    evidence_level = _validate_intensity_scale_evidence_level(
        value.input_intensity_scale_evidence_level,
        role=role,
    )
    evidence_source = _require_provenance_text(
        value.input_intensity_scale_source,
        role=role,
        field_name="input_intensity_scale_evidence.input_intensity_scale_source",
    )
    source_detail = _validate_optional_provenance_text(
        value.input_intensity_scale_source_detail,
        role=role,
        field_name="input_intensity_scale_evidence.input_intensity_scale_source_detail",
    )
    return InputIntensityScaleEvidence(
        input_intensity_scale=input_intensity_scale.value,
        input_intensity_scale_evidence_level=evidence_level.value,
        input_intensity_scale_source=evidence_source,
        input_intensity_scale_source_detail=source_detail,
    )


def _validate_intensity_scale(
    value: object,
    *,
    role: Literal["Selected", "Background"],
) -> IntensityScaleKind:
    try:
        return IntensityScaleKind(str(value).strip())
    except ValueError as exc:
        supported = ", ".join(item.value for item in IntensityScaleKind)
        raise WorkflowValidationError(
            f"{role} identifier-set provenance is invalid: "
            "input_intensity_scale_evidence.input_intensity_scale must be one "
            f"of: {supported}; observed={value!r}"
        ) from exc


def _validate_intensity_scale_evidence_level(
    value: object,
    *,
    role: Literal["Selected", "Background"],
) -> IntensityScaleEvidenceLevel:
    try:
        return IntensityScaleEvidenceLevel(str(value).strip())
    except ValueError as exc:
        supported = ", ".join(item.value for item in IntensityScaleEvidenceLevel)
        raise WorkflowValidationError(
            f"{role} identifier-set provenance is invalid: "
            "input_intensity_scale_evidence."
            "input_intensity_scale_evidence_level must be one of: "
            f"{supported}; observed={value!r}"
        ) from exc


def _require_provenance_text(
    value: object,
    *,
    role: Literal["Selected", "Background"],
    field_name: str,
) -> str:
    if not isinstance(value, str) or not value.strip():
        raise WorkflowValidationError(
            f"{role} identifier-set provenance is invalid: {field_name} must be "
            "nonblank"
        )
    return value.strip()


def _validate_optional_provenance_text(
    value: object | None,
    *,
    role: Literal["Selected", "Background"],
    field_name: str,
) -> str | None:
    if value is None:
        return None
    return _require_provenance_text(value, role=role, field_name=field_name)


def _require_identifier_kind(
    value: object,
    *,
    field_name: str,
) -> EnrichmentIdentifierKind:
    if value not in SUPPORTED_ENRICHMENT_IDENTIFIER_KINDS:
        supported = ", ".join(
            repr(identifier_kind)
            for identifier_kind in SUPPORTED_ENRICHMENT_IDENTIFIER_KINDS
        )
        raise WorkflowValidationError(f"{field_name} must be one of: {supported}")
    return cast(EnrichmentIdentifierKind, value)


def _expected_collection_kind(
    identifier_kind: EnrichmentIdentifierKind,
) -> EnrichmentCollectionKind:
    if identifier_kind in GENE_LEVEL_ENRICHMENT_IDENTIFIER_KINDS:
        return ENRICHMENT_COLLECTION_KIND_GENE_SET
    if identifier_kind in PTM_LEVEL_ENRICHMENT_IDENTIFIER_KINDS:
        return ENRICHMENT_COLLECTION_KIND_PTM_SET
    raise WorkflowValidationError(
        f"enrichment workflow request identifier_kind is unsupported: {identifier_kind!r}"
    )


def _normalise_identifier_sequence(
    value: object,
    *,
    field_name: str,
    allow_empty: bool,
) -> tuple[str, ...]:
    if isinstance(value, str | bytes | bytearray) or not isinstance(value, Sequence):
        raise WorkflowValidationError(
            f"{field_name} must be a sequence of non-empty strings"
        )
    seen: set[str] = set()
    normalised: list[str] = []
    for raw_identifier in value:
        identifier = _require_non_empty_string(
            raw_identifier,
            field_name=f"{field_name}[]",
        )
        if identifier in seen:
            continue
        seen.add(identifier)
        normalised.append(identifier)
    if not allow_empty and not normalised:
        raise WorkflowValidationError(f"{field_name} must not be empty")
    return tuple(normalised)


def _require_non_empty_string(value: object, *, field_name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise WorkflowValidationError(f"{field_name} must be a non-empty string")
    return value.strip()


def _validate_optional_set_size_threshold(
    value: object | None,
    *,
    field_name: str,
) -> int | None:
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, int):
        raise WorkflowValidationError(f"{field_name} must be an int or None")
    if value < 1:
        raise WorkflowValidationError(
            f"{field_name} must be greater than or equal to 1"
        )
    return value


__all__ = [
    "EnrichmentSelectedIdentifierSource",
    "EnrichmentWorkflowValidator",
    "ValidatedEnrichmentWorkflowRequest",
]
