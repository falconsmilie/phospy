from __future__ import annotations

import json
from dataclasses import fields

import phospy.advanced
import phospy.api
from phospy.quantitative.evidence import InputIntensityScaleEvidence
from phospy.quantitative.models import IntensityScaleKind
from phospy.quantitative.serialization import intensity_scale_state_to_payload
from phospy.quantitative.transformers import IdentityTransformer
from phospy.statistics.multiple_testing import MultipleTestingCorrection
from tests.support.intensity_scale_states import (
    supported_log2_intensity_scale_state,
)


def test_input_intensity_scale_evidence_shape_payload_order_and_facade_identities() -> (
    None
):
    assert tuple(field.name for field in fields(InputIntensityScaleEvidence)) == (
        "input_intensity_scale",
        "input_intensity_scale_evidence_level",
        "input_intensity_scale_source",
        "input_intensity_scale_source_detail",
    )
    without_detail = InputIntensityScaleEvidence(
        input_intensity_scale="log2",
        input_intensity_scale_evidence_level="observed_transformation",
        input_intensity_scale_source="transformed_by_phospy",
    )
    with_detail = InputIntensityScaleEvidence(
        input_intensity_scale="log2",
        input_intensity_scale_evidence_level="declared_by_user",
        input_intensity_scale_source="declared_by_user",
        input_intensity_scale_source_detail="request.input_intensity_scale",
    )

    assert tuple(without_detail.to_payload()) == (
        "input_intensity_scale",
        "input_intensity_scale_evidence_level",
        "input_intensity_scale_source",
    )
    assert "input_intensity_scale_source_detail" not in without_detail.to_payload()
    assert tuple(with_detail.to_payload()) == (
        "input_intensity_scale",
        "input_intensity_scale_evidence_level",
        "input_intensity_scale_source",
        "input_intensity_scale_source_detail",
    )
    assert phospy.api.InputIntensityScaleEvidence is InputIntensityScaleEvidence
    assert phospy.api.IntensityScaleKind is IntensityScaleKind
    assert phospy.advanced.MultipleTestingCorrection is MultipleTestingCorrection


def test_intensity_scale_serialization_order_and_json_bytes_are_exact() -> None:
    payload = intensity_scale_state_to_payload(
        supported_log2_intensity_scale_state(has_total_matrix=False)
    )

    assert tuple(payload) == (
        "phospho",
        "total",
        "quantity",
        "quantitative_meaning_provenance",
    )
    assert tuple(payload["phospho"]) == ("kind", "transformed", "established_by")
    assert json.dumps(payload, separators=(",", ":"), ensure_ascii=False) == (
        '{"phospho":{"kind":"log2","transformed":true,'
        '"established_by":"trusted.input"},"total":null,'
        '"quantity":"phosphosite_log_abundance",'
        '"quantitative_meaning_provenance":{"schema_version":1,'
        '"source_quantity":null,"target_quantity":"phosphosite_log_abundance",'
        '"operation_id":"phospy.dataset_builder.quantitative_meaning.'
        'infer_from_scale_contract","producer_id":"phospy.science.datasets.'
        'preprocessing.state_builder","evidence_mode":"inferred_from_scale_contract",'
        '"parameters":{"scale_kind":"log2","scale_label":"log2"},'
        '"input_table_fingerprints":[],"output_table_fingerprint":null,'
        '"trace_id":null,"diagnostic_caveat_codes":[]}}'
    )


def test_identity_transformer_identifier_remains_historical() -> None:
    assert IdentityTransformer.stable_identifier == (
        "phospy.science.transformations.transformers.identity.IdentityTransformer"
    )
