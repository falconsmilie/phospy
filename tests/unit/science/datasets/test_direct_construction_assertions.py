from __future__ import annotations

import json
from collections.abc import Mapping
from dataclasses import replace

import pandas as pd
import pytest

from phospy.api import Organism
from phospy.errors import DatasetValidationError, PhosPyInputError
from phospy.provenance import (
    RunProvenance,
    TrustedDatasetConstructionAssertions,
    TrustedDatasetConstructionEvidence,
    from_payload,
    to_payload,
)
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
)
from tests.support.site_keys import protein_site_key_index, site_key_context_columns


def _site_index() -> pd.Index:
    return protein_site_key_index(
        protein_identifiers=["MAPK14", "AKT1"],
        sites=["Y182", "T308"],
    )


def _phospho(index: pd.Index) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample_a": [1.0, 2.0],
            "sample_b": [1.5, 2.5],
        },
        index=index.copy(),
    )


def _site_metadata(index: pd.Index) -> pd.DataFrame:
    display_ids = ["MAPK14;Y182;", "AKT1;T308;"]
    return pd.DataFrame(
        {
            "site_key": index.astype(str).tolist(),
            "display_id": display_ids,
            **site_key_context_columns(index),
            "gene_symbol": ["MAPK14", "AKT1"],
            "protein_id": ["MAPK14", "AKT1"],
            "site": ["Y182", "T308"],
            "site_sequence": [
                "AAAAAAAAAAAAAAAYAAAAAAAAAAAAAAA",
                "AAAAAAAAAAAAAAATAAAAAAAAAAAAAAA",
            ],
            "localisation_confidence": [0.95, 0.9],
        },
        index=index.copy(),
    )


def _complete_assertions(
    *,
    localisation: TrustedDatasetConstructionEvidence | None = None,
    reference_context: TrustedDatasetConstructionEvidence | None = None,
    identity: TrustedDatasetConstructionEvidence | None = None,
    aligned_structure: TrustedDatasetConstructionEvidence | None = None,
) -> TrustedDatasetConstructionAssertions:
    return TrustedDatasetConstructionAssertions(
        identity=identity
        or TrustedDatasetConstructionEvidence.evidence(
            source="protein-scoped site_key export",
            details={"site_key_schema": "protein-scoped-v1"},
        ),
        intensity_scale=TrustedDatasetConstructionEvidence.evidence(
            source="IntensityScaleState established before trusted construction",
            policy="require_established_intensity_scale_state",
            details={"scale": "linear"},
        ),
        quantitative_meaning=TrustedDatasetConstructionEvidence.evidence(
            source="curated log2 intensity export",
            policy="analysis-ready quantitative matrix",
        ),
        aligned_structure=aligned_structure
        or TrustedDatasetConstructionEvidence.evidence(
            source="pre-export table alignment audit",
            policy="require_identical_site_indexes_and_sample_axes",
            details={
                "phospho_index": "site_key",
                "site_metadata_index": "site_key",
                "sample_axis": "phospho.columns",
            },
        ),
        localisation=localisation
        or TrustedDatasetConstructionEvidence.evidence(
            source="localisation_confidence column",
            policy="require_threshold",
            threshold=0.75,
        ),
        sequence=TrustedDatasetConstructionEvidence.evidence(
            source="site_sequence column curated before PhosPy import",
        ),
        reference_context=reference_context
        or TrustedDatasetConstructionEvidence.waiver(
            reason="source export did not retain reference context metadata",
        ),
        asserted_by="unit-test",
        assertion_source="curated analysis-ready export",
    )


def _trusted_dataset(
    assertions: TrustedDatasetConstructionAssertions | None = None,
    provenance: RunProvenance | None = None,
) -> AnalysisReadyPhosphoDataset:
    index = _site_index()
    return AnalysisReadyPhosphoDataset.from_trusted_tables(
        phospho=_phospho(index),
        site_metadata=_site_metadata(index),
        organism=Organism.RAT,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
        trusted_construction_assertions=assertions or _complete_assertions(),
        provenance=provenance,
    )


def test_from_trusted_tables_records_typed_construction_assertions() -> None:
    assertions = _complete_assertions()

    dataset = _trusted_dataset(assertions)

    assert dataset.trusted_construction_assertions == assertions
    assert dataset.provenance is not None
    construction = dataset.provenance.workflow_parameters["construction"]
    assert isinstance(construction, Mapping)
    payload = construction["trusted_construction_assertions"]
    assert isinstance(payload, Mapping)
    assert payload["schema_version"] == 3
    assert payload["assertion_metadata_provided"] is True
    assert payload["identity"]["source"] == "protein-scoped site_key export"
    assert payload["intensity_scale"]["policy"] == (
        "require_established_intensity_scale_state"
    )
    assert payload["quantitative_meaning"]["policy"] == (
        "analysis-ready quantitative matrix"
    )
    assert payload["aligned_structure"]["policy"] == (
        "require_identical_site_indexes_and_sample_axes"
    )
    assert payload["localisation"]["source"] == "localisation_confidence column"
    assert payload["localisation"]["policy"] == "require_threshold"
    assert payload["localisation"]["threshold"] == 0.75
    assert payload["sequence"]["source"] == (
        "site_sequence column curated before PhosPy import"
    )
    assert payload["reference_context"]["kind"] == "waiver"
    assert payload["sequence_user_asserted"] is True
    assert payload["identity_user_asserted"] is True
    assert payload["intensity_scale_user_asserted"] is True
    assert payload["quantitative_meaning_user_asserted"] is True
    assert payload["aligned_structure_user_asserted"] is True
    assert payload["localisation_user_asserted"] is True
    assert payload["reference_context_user_asserted"] is True
    assert payload["asserted_by"] == "unit-test"
    assert payload["assertion_source"] == "curated analysis-ready export"
    assert payload["waived_assertions"] == ["reference_context"]
    assert payload["missing_assertions"] == []
    assert construction["missing_trusted_assertions"] == []
    assert construction["trusted_construction_assertion_fingerprint"] == (
        assertions.assertion_fingerprint
    )


def test_from_trusted_tables_rejects_aligned_structure_waiver() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="aligned_structure cannot be waived",
    ):
        _complete_assertions(
            aligned_structure=TrustedDatasetConstructionEvidence.waiver(
                reason="caller wants to skip alignment audit"
            )
        )


def test_from_trusted_tables_rejects_missing_localisation_evidence() -> None:
    index = _site_index()

    with pytest.raises(
        DatasetValidationError,
        match="from_trusted_tables requires.*localisation",
    ):
        AnalysisReadyPhosphoDataset.from_trusted_tables(
            phospho=_phospho(index),
            site_metadata=_site_metadata(index),
            organism=Organism.RAT,
            intensity_scale_state=supported_linear_intensity_scale_state(
                has_total_matrix=False
            ),
            processing_state=supported_linear_processing_state(has_total_matrix=False),
        )


def test_from_trusted_tables_accepts_and_serializes_localisation_waiver() -> None:
    assertions = _complete_assertions(
        localisation=TrustedDatasetConstructionEvidence.waiver(
            reason="historical source lacks localisation confidence export",
            policy="trusted_curation_waiver",
        )
    )

    dataset = _trusted_dataset(assertions)

    assert dataset.provenance is not None
    construction = dataset.provenance.workflow_parameters["construction"]
    payload = construction["trusted_construction_assertions"]
    assert isinstance(payload, Mapping)
    assert payload["localisation"]["kind"] == "waiver"
    assert payload["localisation"]["waiver_reason"] == (
        "historical source lacks localisation confidence export"
    )
    assert "localisation" in payload["waived_assertions"]
    assert payload["missing_assertions"] == []


def test_from_trusted_tables_rejects_false_table_fingerprint() -> None:
    assertions = _complete_assertions()
    trusted = _trusted_dataset(assertions)
    assert trusted.provenance is not None
    bad_output_tables = tuple(
        replace(fingerprint, exact_hash_value="0" * 64)
        if fingerprint.name == "dataset.phospho"
        else fingerprint
        for fingerprint in trusted.provenance.output_tables
    )
    bad_provenance = replace(trusted.provenance, output_tables=bad_output_tables)

    with pytest.raises(
        DatasetValidationError,
        match=r"run_provenance\.output_tables\.dataset\.phospho.*exact_hash_value",
    ):
        _trusted_dataset(assertions=assertions, provenance=bad_provenance)


def test_from_trusted_tables_serializes_seven_assertion_dimensions() -> None:
    dataset = _trusted_dataset()
    assert dataset.provenance is not None

    payload = to_payload(dataset.provenance)
    json.dumps(payload)
    restored = from_payload(payload)

    assert to_payload(restored) == payload
    construction = payload["workflow_parameters"]["construction"]
    trusted_assertions = construction["trusted_construction_assertions"]
    assert trusted_assertions["schema_version"] == 3
    assert trusted_assertions["missing_assertions"] == []
    for dimension in (
        "identity",
        "intensity_scale",
        "quantitative_meaning",
        "aligned_structure",
        "localisation",
        "sequence",
        "reference_context",
    ):
        assert trusted_assertions[dimension] is not None
        assert trusted_assertions[f"{dimension}_user_asserted"] is True


def test_direct_constructor_emits_deprecation_warning() -> None:
    index = _site_index()

    with pytest.warns(
        DeprecationWarning,
        match="AnalysisReadyPhosphoDataset\\(\\.\\.\\.\\) direct construction "
        "is deprecated",
    ):
        AnalysisReadyPhosphoDataset(
            phospho=_phospho(index),
            site_metadata=_site_metadata(index),
            organism=Organism.RAT,
            intensity_scale_state=supported_linear_intensity_scale_state(
                has_total_matrix=False
            ),
            processing_state=supported_linear_processing_state(has_total_matrix=False),
        )


def test_caller_mutable_assertion_details_cannot_alter_provenance() -> None:
    details = {"columns": ["site_key"], "schema": {"version": 1}}
    assertions = _complete_assertions(
        identity=TrustedDatasetConstructionEvidence.evidence(
            source="mutable caller metadata",
            details=details,
        )
    )

    dataset = _trusted_dataset(assertions)
    details["columns"].append("mutated")
    details["schema"]["version"] = 99

    assert dataset.provenance is not None
    construction = dataset.provenance.workflow_parameters["construction"]
    payload = construction["trusted_construction_assertions"]
    identity_payload = payload["identity"]
    assert identity_payload["details"] == {
        "columns": ["site_key"],
        "schema": {"version": 1},
    }
    with pytest.raises(TypeError):
        assertions.identity.details["schema"] = {"version": 2}  # type: ignore[index]


def test_from_trusted_tables_rejects_untyped_assertion_mapping() -> None:
    index = _site_index()

    with pytest.raises(
        DatasetValidationError,
        match="TrustedDatasetConstructionAssertions",
    ):
        AnalysisReadyPhosphoDataset.from_trusted_tables(
            phospho=_phospho(index),
            site_metadata=_site_metadata(index),
            organism=Organism.RAT,
            intensity_scale_state=supported_linear_intensity_scale_state(
                has_total_matrix=False
            ),
            processing_state=supported_linear_processing_state(has_total_matrix=False),
            trusted_construction_assertions={
                "sequence_user_asserted": True,
                "identity_user_asserted": True,
                "quantitative_meaning_user_asserted": True,
                "reference_context_user_asserted": True,
            },
        )
