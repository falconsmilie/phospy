from __future__ import annotations

import pandas.testing as pdt
import pytest

from phospy.api import (
    ContractValidationError,
    EnrichmentConfig,
    EnrichmentIdentifierSetProvenance,
    EnrichmentIdentifierSetSourceType,
    EnrichmentWorkflow,
    EnrichmentWorkflowRequest,
    GeneSetCollection,
    WorkflowValidationError,
)
from phospy.api.configs import ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL
from phospy.provenance.models import InputIntensityScaleEvidence
from phospy.validation.workflows.enrichment import EnrichmentWorkflowValidator
from phospy.workflows.enrichment.interpreter import EnrichmentWorkflowInterpreter
from phospy.workflows.enrichment.validator import (
    EnrichmentWorkflowValidator as InternalEnrichmentWorkflowValidator,
)
from phospy.workflows.intensity_scale_evidence import (
    INPUT_INTENSITY_SCALE_DECLARED_CAVEAT_CODE,
)


def _gene_collection() -> GeneSetCollection:
    return GeneSetCollection(
        sets={
            "MAPK_PATHWAY": ("AKT1", "MAPK1"),
            "MTOR_SIGNALING": ("MTOR",),
        },
        identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
        term_names={
            "MAPK_PATHWAY": "MAPK pathway",
            "MTOR_SIGNALING": "MTOR signaling",
        },
        source_name="unit_test",
    )


def _observed_evidence() -> InputIntensityScaleEvidence:
    return InputIntensityScaleEvidence(
        input_intensity_scale="log2",
        input_intensity_scale_evidence_level="observed_transformation",
        input_intensity_scale_source="phospy.science.transformations.transformers.log2",
        input_intensity_scale_source_detail="log2_transform",
    )


def _declared_evidence() -> InputIntensityScaleEvidence:
    return InputIntensityScaleEvidence(
        input_intensity_scale="linear",
        input_intensity_scale_evidence_level="declared_by_user",
        input_intensity_scale_source="declared_by_user",
        input_intensity_scale_source_detail="unit-test declaration",
    )


def _provenance(
    source_type: EnrichmentIdentifierSetSourceType,
    *,
    count: int,
    label: str = "unit-test identifiers",
    evidence: InputIntensityScaleEvidence | None = None,
) -> EnrichmentIdentifierSetProvenance:
    return EnrichmentIdentifierSetProvenance(
        source_type=source_type,
        source_label=label,
        identifier_count=count,
        upstream_workflow_id="workflow-1",
        upstream_result_id="result-1",
        input_intensity_scale_evidence=evidence,
    )


def _request(
    *,
    selected_provenance: object | None = None,
    background_provenance: object | None = None,
    background_universe: object = ("AKT1", "MAPK1", "MTOR", "MTOR"),
) -> EnrichmentWorkflowRequest:
    return EnrichmentWorkflowRequest(
        identifier_column="gene_symbol",
        identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
        set_collection=_gene_collection(),
        selected_identifiers=(" AKT1 ", "MAPK1", "AKT1"),
        background_universe=background_universe,  # type: ignore[arg-type]
        config=EnrichmentConfig(),
        selected_identifier_provenance=selected_provenance,  # type: ignore[arg-type]
        background_identifier_provenance=background_provenance,  # type: ignore[arg-type]
    )


def _declared_role_caveats(result) -> dict[str, object]:
    return {
        caveat.details["identifier_set_role"]: caveat
        for caveat in result.caveats
        if caveat.code == INPUT_INTENSITY_SCALE_DECLARED_CAVEAT_CODE
        and "identifier_set_role" in caveat.details
    }


def test_legacy_request_without_provenance_still_runs() -> None:
    request = _request()

    validated = EnrichmentWorkflowValidator().run(request)
    result = EnrichmentWorkflow().run(request)

    assert validated.selected_identifier_provenance is None
    assert validated.background_identifier_provenance is None
    assert result.selected_identifier_provenance is None
    assert result.background_identifier_provenance is None
    assert result.provenance is not None
    assert (
        result.provenance.workflow_parameters["selected_identifier_provenance"] is None
    )
    assert (
        result.provenance.workflow_parameters["background_identifier_provenance"]
        is None
    )
    assert result.records


@pytest.mark.parametrize(
    "source_type",
    (
        EnrichmentIdentifierSetSourceType.MANUAL,
        EnrichmentIdentifierSetSourceType.RAW_IDENTIFIER_LIST,
    ),
)
def test_manual_and_raw_selected_provenance_validate_without_evidence(
    source_type: EnrichmentIdentifierSetSourceType,
) -> None:
    provenance = _provenance(source_type, count=2)

    validated = EnrichmentWorkflowValidator().run(
        _request(selected_provenance=provenance)
    )

    assert validated.selected_identifier_provenance == provenance
    assert (
        validated.selected_identifier_provenance.input_intensity_scale_evidence is None
    )


def test_phospy_derived_selected_observed_evidence_propagates_without_caveat() -> None:
    provenance = _provenance(
        EnrichmentIdentifierSetSourceType.PHOSPY_DERIVED_QUANTITATIVE,
        count=2,
        evidence=_observed_evidence(),
    )

    result = EnrichmentWorkflow().run(_request(selected_provenance=provenance))

    assert result.selected_identifier_provenance == provenance
    assert result.selected_identifier_provenance.input_intensity_scale_evidence == (
        _observed_evidence()
    )
    assert _declared_role_caveats(result) == {}


def test_phospy_derived_selected_declared_evidence_emits_selected_caveat() -> None:
    provenance = _provenance(
        EnrichmentIdentifierSetSourceType.PHOSPY_DERIVED_QUANTITATIVE,
        count=2,
        evidence=_declared_evidence(),
    )

    result = EnrichmentWorkflow().run(_request(selected_provenance=provenance))
    caveat = _declared_role_caveats(result)["selected"]

    assert caveat.details["workflow_scope"] == "enrichment"
    assert caveat.details["identifier_set_role"] == "selected"
    assert caveat.details["input_intensity_scale"] == "linear"
    assert caveat.details["input_intensity_scale_source"] == "declared_by_user"
    assert caveat.details["input_intensity_scale_source_detail"] == (
        "unit-test declaration"
    )
    assert caveat.details["scale_declared_not_observed"] is True


def test_phospy_derived_background_observed_evidence_propagates_without_caveat() -> (
    None
):
    provenance = _provenance(
        EnrichmentIdentifierSetSourceType.PHOSPY_DERIVED_QUANTITATIVE,
        count=3,
        label="background from quantitative workflow",
        evidence=_observed_evidence(),
    )

    result = EnrichmentWorkflow().run(_request(background_provenance=provenance))

    assert result.background_identifier_provenance == provenance
    assert _declared_role_caveats(result) == {}


def test_mixed_selected_manual_and_background_derived_inputs_survive_interpretation() -> (
    None
):
    selected = _provenance(EnrichmentIdentifierSetSourceType.MANUAL, count=2)
    background = _provenance(
        EnrichmentIdentifierSetSourceType.PHOSPY_DERIVED_QUANTITATIVE,
        count=3,
        label="quantitative background",
        evidence=_observed_evidence(),
    )

    validated = InternalEnrichmentWorkflowValidator().run(
        _request(selected_provenance=selected, background_provenance=background)
    )
    interpreted = EnrichmentWorkflowInterpreter().run(validated)

    assert interpreted.selected_identifier_provenance == selected
    assert interpreted.background_identifier_provenance == background


@pytest.mark.parametrize(
    ("role", "kwargs", "pattern"),
    (
        (
            "selected",
            {
                "selected_provenance": _provenance(
                    EnrichmentIdentifierSetSourceType.PHOSPY_DERIVED_QUANTITATIVE,
                    count=2,
                )
            },
            "Selected identifier-set provenance.*requires input_intensity_scale_evidence",
        ),
        (
            "background",
            {
                "background_provenance": _provenance(
                    EnrichmentIdentifierSetSourceType.PHOSPY_DERIVED_QUANTITATIVE,
                    count=3,
                    label="background",
                )
            },
            "Background identifier-set provenance.*requires input_intensity_scale_evidence",
        ),
    ),
)
def test_phospy_derived_provenance_requires_intensity_evidence(
    role: str,
    kwargs: dict[str, object],
    pattern: str,
) -> None:
    _ = role
    with pytest.raises(WorkflowValidationError, match=pattern):
        EnrichmentWorkflowValidator().run(_request(**kwargs))


@pytest.mark.parametrize(
    ("kwargs", "pattern"),
    (
        (
            {
                "selected_provenance": _provenance(
                    EnrichmentIdentifierSetSourceType.MANUAL, count=1
                )
            },
            "Selected identifier-set provenance count mismatch.*declared identifier_count=1",
        ),
        (
            {
                "background_provenance": _provenance(
                    EnrichmentIdentifierSetSourceType.MANUAL, count=2
                )
            },
            "Background identifier-set provenance count mismatch.*declared identifier_count=2",
        ),
    ),
)
def test_identifier_count_mismatches_fail(
    kwargs: dict[str, object],
    pattern: str,
) -> None:
    with pytest.raises(WorkflowValidationError, match=pattern):
        EnrichmentWorkflowValidator().run(_request(**kwargs))


def test_background_provenance_without_corresponding_background_fails() -> None:
    with pytest.raises(
        WorkflowValidationError,
        match="Background identifier-set provenance.*background_universe",
    ):
        EnrichmentWorkflowValidator().run(
            _request(
                background_universe=None,
                background_provenance=_provenance(
                    EnrichmentIdentifierSetSourceType.MANUAL,
                    count=0,
                ),
            )
        )


def test_blank_source_label_fails_with_role_and_field() -> None:
    with pytest.raises(
        ContractValidationError,
        match="enrichment identifier-set provenance source_label",
    ):
        _provenance(
            EnrichmentIdentifierSetSourceType.MANUAL,
            count=2,
            label=" ",
        )


def test_invalid_source_type_during_provenance_deserialization_fails() -> None:
    with pytest.raises(ContractValidationError, match="source_type must be one of"):
        EnrichmentIdentifierSetProvenance(
            source_type="spreadsheet_guess",  # type: ignore[arg-type]
            source_label="bad",
            identifier_count=2,
        )


def test_untyped_request_provenance_mapping_fails_workflow_validation() -> None:
    request = _request(
        selected_provenance={
            "source_type": "manual",
            "source_label": "typed model required",
            "identifier_count": 2,
        }
    )

    with pytest.raises(
        WorkflowValidationError,
        match="selected_identifier_provenance must be EnrichmentIdentifierSetProvenance",
    ):
        EnrichmentWorkflowValidator().run(request)


def test_mutated_invalid_source_type_fails_workflow_validation() -> None:
    provenance = _provenance(EnrichmentIdentifierSetSourceType.MANUAL, count=2)
    object.__setattr__(provenance, "source_type", "spreadsheet_guess")

    with pytest.raises(
        WorkflowValidationError,
        match="Selected identifier-set provenance.*source_type",
    ):
        EnrichmentWorkflowValidator().run(_request(selected_provenance=provenance))


@pytest.mark.parametrize(
    ("evidence", "pattern"),
    (
        (
            InputIntensityScaleEvidence(
                input_intensity_scale="natural_log",
                input_intensity_scale_evidence_level="observed_transformation",
                input_intensity_scale_source="transformer",
            ),
            "input_intensity_scale.*one of",
        ),
        (
            InputIntensityScaleEvidence(
                input_intensity_scale="log2",
                input_intensity_scale_evidence_level="visual_guess",
                input_intensity_scale_source="transformer",
            ),
            "input_intensity_scale_evidence_level.*one of",
        ),
    ),
)
def test_invalid_intensity_evidence_fields_fail(
    evidence: InputIntensityScaleEvidence,
    pattern: str,
) -> None:
    provenance = _provenance(
        EnrichmentIdentifierSetSourceType.PHOSPY_DERIVED_QUANTITATIVE,
        count=2,
        evidence=evidence,
    )

    with pytest.raises(WorkflowValidationError, match=pattern):
        EnrichmentWorkflowValidator().run(_request(selected_provenance=provenance))


def test_blank_intensity_evidence_source_fails() -> None:
    evidence = _observed_evidence()
    object.__setattr__(evidence, "input_intensity_scale_source", " ")
    provenance = _provenance(
        EnrichmentIdentifierSetSourceType.PHOSPY_DERIVED_QUANTITATIVE,
        count=2,
        evidence=evidence,
    )

    with pytest.raises(
        WorkflowValidationError,
        match="input_intensity_scale_source",
    ):
        EnrichmentWorkflowValidator().run(_request(selected_provenance=provenance))


@pytest.mark.parametrize(
    "source_type",
    (
        EnrichmentIdentifierSetSourceType.MANUAL,
        EnrichmentIdentifierSetSourceType.RAW_IDENTIFIER_LIST,
    ),
)
def test_manual_and_raw_provenance_reject_intensity_evidence(
    source_type: EnrichmentIdentifierSetSourceType,
) -> None:
    provenance = _provenance(
        source_type,
        count=2,
        evidence=_observed_evidence(),
    )

    with pytest.raises(
        WorkflowValidationError,
        match="input_intensity_scale_evidence is only valid",
    ):
        EnrichmentWorkflowValidator().run(_request(selected_provenance=provenance))


def test_both_provenance_objects_survive_result_and_run_provenance_serialization() -> (
    None
):
    selected = _provenance(
        EnrichmentIdentifierSetSourceType.PHOSPY_DERIVED_QUANTITATIVE,
        count=2,
        label="selected quantitative hits",
        evidence=_observed_evidence(),
    )
    background = _provenance(
        EnrichmentIdentifierSetSourceType.PHOSPY_DERIVED_QUANTITATIVE,
        count=3,
        label="background quantitative universe",
        evidence=_observed_evidence(),
    )

    result = EnrichmentWorkflow().run(
        _request(selected_provenance=selected, background_provenance=background)
    )

    assert result.selected_identifier_provenance == selected
    assert result.background_identifier_provenance == background
    assert result.provenance is not None
    workflow_parameters = result.provenance.workflow_parameters
    selected_payload = workflow_parameters["selected_identifier_provenance"]
    background_payload = workflow_parameters["background_identifier_provenance"]
    assert selected_payload["source_type"] == "phospy_derived_quantitative"
    assert selected_payload["source_label"] == "selected quantitative hits"
    assert selected_payload["identifier_count"] == 2
    assert selected_payload["upstream_workflow_id"] == "workflow-1"
    assert selected_payload["upstream_result_id"] == "result-1"
    assert selected_payload["input_intensity_scale_evidence"] == (
        _observed_evidence().to_payload()
    )
    assert background_payload["source_label"] == "background quantitative universe"
    assert background_payload["identifier_count"] == 3


def test_declared_selected_and_background_evidence_emit_role_specific_caveats() -> None:
    selected = _provenance(
        EnrichmentIdentifierSetSourceType.PHOSPY_DERIVED_QUANTITATIVE,
        count=2,
        label="declared selected",
        evidence=_declared_evidence(),
    )
    background = _provenance(
        EnrichmentIdentifierSetSourceType.PHOSPY_DERIVED_QUANTITATIVE,
        count=3,
        label="declared background",
        evidence=_declared_evidence(),
    )

    result = EnrichmentWorkflow().run(
        _request(selected_provenance=selected, background_provenance=background)
    )

    by_role = _declared_role_caveats(result)
    assert set(by_role) == {"selected", "background"}
    assert by_role["selected"].details["identifier_set_source_label"] == (
        "declared selected"
    )
    assert by_role["background"].details["identifier_set_source_label"] == (
        "declared background"
    )


def test_enrichment_records_and_numeric_statistics_are_identical_with_provenance() -> (
    None
):
    baseline = EnrichmentWorkflow().run(_request())
    with_provenance = EnrichmentWorkflow().run(
        _request(
            selected_provenance=_provenance(
                EnrichmentIdentifierSetSourceType.PHOSPY_DERIVED_QUANTITATIVE,
                count=2,
                evidence=_observed_evidence(),
            ),
            background_provenance=_provenance(
                EnrichmentIdentifierSetSourceType.PHOSPY_DERIVED_QUANTITATIVE,
                count=3,
                evidence=_observed_evidence(),
            ),
        )
    )

    assert with_provenance.records == baseline.records
    pdt.assert_frame_equal(with_provenance.table, baseline.table)
