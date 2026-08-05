from __future__ import annotations

import inspect
import json
import warnings
from collections.abc import Mapping
from dataclasses import replace
from typing import cast

import pandas as pd
import pytest

from phospy.api import Organism
from phospy.errors import DatasetValidationError, PhosPyInputError
from phospy.provenance import (
    ReferenceProvenance,
    RunProvenance,
    TableFingerprint,
    TrustedDatasetConstructionAssertions,
    TrustedDatasetConstructionEvidence,
    from_payload,
    to_payload,
)
from phospy.provenance.reference_context import ReferenceContext
from phospy.science.datasets.construction import service as construction_service
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.science.transformations.state_coherence import observe_numeric_domain
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
        numeric_semantic_domain=TrustedDatasetConstructionEvidence.evidence(
            source="pre-export numeric-semantic domain audit",
            policy="analysis_ready_numeric_semantic_domain_checked",
            details=_numeric_domain_details(_phospho(_site_index())),
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


def _trusted_dataset_with_all_tables() -> AnalysisReadyPhosphoDataset:
    index = _site_index()
    phospho = _phospho(index)
    return AnalysisReadyPhosphoDataset.from_trusted_tables(
        phospho=phospho,
        site_metadata=_site_metadata(index),
        sample_metadata=pd.DataFrame(
            {"condition": ["control", "treated"]},
            index=phospho.columns.copy(),
        ),
        total=pd.DataFrame(
            {
                "sample_a": [10.0, 20.0],
                "sample_b": [11.0, 21.0],
            },
            index=pd.Index(["MAPK14", "AKT1"], name="protein_id"),
        ),
        comparisons=pd.DataFrame(
            {"treated_vs_control": [0.5, 0.7]},
            index=index.copy(),
        ),
        imputation_observation_mask=pd.DataFrame(
            {
                "sample_a": [True, False],
                "sample_b": [True, True],
            },
            index=index.copy(),
            columns=phospho.columns.copy(),
        ),
        organism=Organism.RAT,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=True
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=True),
        trusted_construction_assertions=_complete_assertions(),
    )


def _public_constructor_payload_from_dataset(
    dataset: AnalysisReadyPhosphoDataset,
) -> dict[str, object]:
    return {
        "phospho": dataset.phospho,
        "site_metadata": dataset.site_metadata,
        "sample_metadata": dataset.sample_metadata,
        "total": dataset.total,
        "comparisons": dataset.comparisons,
        "imputation_observation_mask": dataset.imputation_observed_mask_dataframe(),
        "organism": dataset.organism,
        "intensity_scale_state": dataset.intensity_scale_state,
        "processing_state": dataset.processing_state,
        "provenance": dataset.provenance,
        "trusted_construction_assertions": dataset.trusted_construction_assertions,
    }


def _provenance_with_construction_payload(
    assertions: TrustedDatasetConstructionAssertions,
    construction_update: Mapping[str, object],
) -> RunProvenance:
    trusted = _trusted_dataset(assertions)
    assert trusted.provenance is not None
    construction = _construction_payload(trusted.provenance)
    construction.update(construction_update)
    return replace(
        trusted.provenance,
        workflow_parameters={"construction": construction},
    )


def _construction_payload(provenance: RunProvenance) -> dict[str, object]:
    construction = provenance.workflow_parameters["construction"]
    assert isinstance(construction, Mapping)
    return dict(cast(Mapping[str, object], construction))


def _assertion_payload(value: object) -> Mapping[str, object]:
    assert isinstance(value, Mapping)
    return cast(Mapping[str, object], value)


def _reference_context(organism: object) -> ReferenceContext:
    return ReferenceContext(
        organism=organism,
        protein_namespace="gene_symbol",
        source_name="unit-reference",
        source_version="v1",
        proteome_version=None,
        reference_table_sha256="a" * 64,
    )


def _numeric_domain_details(phospho: pd.DataFrame) -> Mapping[str, object]:
    observation = observe_numeric_domain(phospho, table_name="dataset.phospho")
    return {
        "table": observation.table_name,
        "observed_numeric_domain": observation.observed_domain.value,
        "value_count": observation.value_count,
        "negative_count": observation.negative_count,
        "zero_count": observation.zero_count,
        "positive_count": observation.positive_count,
        "min": observation.minimum,
        "max": observation.maximum,
    }


def _mutate_payload_table(payload: dict[str, object], table_name: str) -> None:
    table_by_name = {
        "dataset.phospho": "phospho",
        "dataset.site_metadata": "site_metadata",
        "dataset.sample_metadata": "sample_metadata",
        "dataset.total": "total",
        "dataset.comparisons": "comparisons",
        "dataset.imputation_observation_mask": "imputation_observation_mask",
    }
    frame = payload[table_by_name[table_name]]
    assert isinstance(frame, pd.DataFrame)
    if table_name == "dataset.site_metadata":
        frame["curation_score"] = [0.1, 0.2]
    elif table_name == "dataset.sample_metadata":
        frame.loc[:, "condition"] = ["changed", "treated"]
    elif table_name == "dataset.imputation_observation_mask":
        frame.iloc[0, 0] = not bool(frame.iloc[0, 0])
    else:
        frame.iloc[0, 0] = float(frame.iloc[0, 0]) + 10.0


def _fingerprint_named(
    fingerprints: tuple[TableFingerprint, ...],
    table_name: str,
) -> TableFingerprint:
    for fingerprint in fingerprints:
        if fingerprint.name == table_name:
            return fingerprint
    raise AssertionError(f"missing fingerprint for {table_name}")


def _assert_stale_provenance_error(
    exc: DatasetValidationError,
    *,
    provenance: RunProvenance,
    table_name: str,
) -> None:
    expected = _fingerprint_named(provenance.output_tables, table_name)
    message = str(exc)
    assert table_name in message
    assert f"expected exact digest {expected.exact_hash_value}" in message
    assert "actual exact digest " in message
    assert f"expected tolerance digest {expected.tolerance_hash_value}" in message
    assert "actual tolerance digest " in message


def _assert_from_trusted_tables_rejects_stale_table(table_name: str) -> None:
    trusted = _trusted_dataset_with_all_tables()
    assert trusted.provenance is not None
    assert trusted.trusted_construction_assertions is not None
    payload = _public_constructor_payload_from_dataset(trusted)
    _mutate_payload_table(payload, table_name)

    with pytest.raises(DatasetValidationError) as exc_info:
        AnalysisReadyPhosphoDataset.from_trusted_tables(
            phospho=cast(pd.DataFrame, payload["phospho"]),
            site_metadata=cast(pd.DataFrame, payload["site_metadata"]),
            sample_metadata=cast(pd.DataFrame | None, payload["sample_metadata"]),
            total=cast(pd.DataFrame | None, payload["total"]),
            comparisons=cast(pd.DataFrame | None, payload["comparisons"]),
            imputation_observation_mask=cast(
                pd.DataFrame | None,
                payload["imputation_observation_mask"],
            ),
            organism=Organism.RAT,
            intensity_scale_state=trusted.intensity_scale_state,
            processing_state=trusted.processing_state,
            provenance=trusted.provenance,
            trusted_construction_assertions=trusted.trusted_construction_assertions,
        )
    _assert_stale_provenance_error(
        exc_info.value,
        provenance=trusted.provenance,
        table_name=table_name,
    )


def test_from_trusted_tables_rejects_stale_phospho_provenance() -> None:
    _assert_from_trusted_tables_rejects_stale_table("dataset.phospho")


def test_from_trusted_tables_rejects_supplied_builder_claimed_provenance() -> None:
    trusted = _trusted_dataset()
    assert trusted.provenance is not None
    assert trusted.trusted_construction_assertions is not None
    construction = dict(
        cast(
            Mapping[str, object], trusted.provenance.workflow_parameters["construction"]
        )
    )
    construction["method"] = "AnalysisReadyDatasetBuilder.run"
    construction["builder_used"] = True
    builder_claiming_provenance = replace(
        trusted.provenance,
        workflow_name="dataset_builder",
        workflow_parameters={"construction": construction},
    )

    with pytest.raises(DatasetValidationError, match="builder-observed construction"):
        _trusted_dataset(
            assertions=trusted.trusted_construction_assertions,
            provenance=builder_claiming_provenance,
        )


def test_from_trusted_tables_rejects_stale_site_metadata_provenance() -> None:
    _assert_from_trusted_tables_rejects_stale_table("dataset.site_metadata")


@pytest.mark.parametrize(
    "table_name",
    [
        pytest.param("dataset.sample_metadata", id="sample_metadata"),
        pytest.param("dataset.total", id="total"),
        pytest.param("dataset.comparisons", id="comparisons"),
        pytest.param(
            "dataset.imputation_observation_mask",
            id="imputation_observation_mask",
        ),
    ],
)
def test_from_trusted_tables_rejects_stale_optional_table_provenance(
    table_name: str,
) -> None:
    _assert_from_trusted_tables_rejects_stale_table(table_name)


def test_direct_constructor_cannot_be_unsealed_by_public_argument() -> None:
    trusted = _trusted_dataset_with_all_tables()
    payload = _public_constructor_payload_from_dataset(trusted)

    with pytest.raises(
        TypeError,
        match="AnalysisReadyDatasetBuilder.*from_trusted_tables",
    ):
        AnalysisReadyPhosphoDataset(
            **payload,
            _emit_direct_constructor_deprecation=False,
        )

    with pytest.raises(
        TypeError,
        match="AnalysisReadyDatasetBuilder.*from_trusted_tables",
    ):
        AnalysisReadyPhosphoDataset(**payload)


def test_from_trusted_tables_remains_warning_free_and_fingerprint_strict() -> None:
    trusted = _trusted_dataset_with_all_tables()
    assert trusted.provenance is not None
    assert trusted.trusted_construction_assertions is not None
    payload = _public_constructor_payload_from_dataset(trusted)
    _mutate_payload_table(payload, "dataset.phospho")

    with warnings.catch_warnings(record=True) as recorded:
        warnings.simplefilter("always")
        with pytest.raises(DatasetValidationError) as exc_info:
            AnalysisReadyPhosphoDataset.from_trusted_tables(
                phospho=cast(pd.DataFrame, payload["phospho"]),
                site_metadata=cast(pd.DataFrame, payload["site_metadata"]),
                sample_metadata=cast(pd.DataFrame, payload["sample_metadata"]),
                total=cast(pd.DataFrame, payload["total"]),
                comparisons=cast(pd.DataFrame, payload["comparisons"]),
                imputation_observation_mask=cast(
                    pd.DataFrame,
                    payload["imputation_observation_mask"],
                ),
                organism=Organism.RAT,
                intensity_scale_state=trusted.intensity_scale_state,
                processing_state=trusted.processing_state,
                provenance=trusted.provenance,
                trusted_construction_assertions=(
                    trusted.trusted_construction_assertions
                ),
            )

    assert recorded == []
    _assert_stale_provenance_error(
        exc_info.value,
        provenance=trusted.provenance,
        table_name="dataset.phospho",
    )


def test_exported_dataset_signature_has_no_private_validation_controls() -> None:
    parameters = inspect.signature(AnalysisReadyPhosphoDataset).parameters

    assert tuple(parameters) == ("args", "kwargs")
    assert "_emit_direct_constructor_deprecation" not in parameters
    assert "_enforce_trusted_table_fingerprints" not in parameters


def test_validated_table_aggregate_rejects_direct_construction() -> None:
    with pytest.raises(TypeError, match="private analysis-ready dataset construction"):
        construction_service._ValidatedAnalysisReadyTables()


def test_private_construction_service_is_not_publicly_exported() -> None:
    assert "_AnalysisReadyDatasetConstructionService" not in getattr(
        construction_service,
        "__all__",
        (),
    )
    assert "_ValidatedAnalysisReadyTables" not in getattr(
        construction_service,
        "__all__",
        (),
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
    assert payload["schema_version"] == 4
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
    assert payload["numeric_semantic_domain"]["policy"] == (
        "analysis_ready_numeric_semantic_domain_checked"
    )
    assert payload["sequence_user_asserted"] is True
    assert payload["identity_user_asserted"] is True
    assert payload["intensity_scale_user_asserted"] is True
    assert payload["quantitative_meaning_user_asserted"] is True
    assert payload["aligned_structure_user_asserted"] is True
    assert payload["localisation_user_asserted"] is True
    assert payload["reference_context_user_asserted"] is True
    assert payload["numeric_semantic_domain_user_asserted"] is True
    assert payload["asserted_by"] == "unit-test"
    assert payload["assertion_source"] == "curated analysis-ready export"
    assert payload["waived_assertions"] == ["reference_context"]
    assert payload["missing_assertions"] == []
    assert construction["missing_trusted_assertions"] == []
    assert construction["trusted_construction_assertion_fingerprint"] == (
        assertions.assertion_fingerprint
    )


def test_from_trusted_tables_rejects_digest_only_supplied_assertion_provenance() -> (
    None
):
    assertions = _complete_assertions()
    provenance = _provenance_with_construction_payload(
        assertions,
        {
            "trusted_construction_assertions": None,
            "trusted_construction_assertion_fingerprint": (
                assertions.assertion_fingerprint
            ),
        },
    )
    construction = _construction_payload(provenance)
    digest_only_provenance = replace(
        provenance,
        workflow_parameters={
            "construction": {
                "trusted_construction_assertion_fingerprint": (
                    construction["trusted_construction_assertion_fingerprint"]
                )
            }
        },
    )

    with pytest.raises(
        DatasetValidationError,
        match="trusted_construction_assertions",
    ):
        _trusted_dataset(assertions=assertions, provenance=digest_only_provenance)


def test_from_trusted_tables_rejects_missing_serialized_assertion_dimension() -> None:
    assertions = _complete_assertions()
    payload = assertions.to_payload()
    payload.pop("sequence")
    payload["sequence_user_asserted"] = False
    provenance = _provenance_with_construction_payload(
        assertions,
        {
            "trusted_construction_assertions": payload,
            "trusted_construction_assertion_fingerprint": (
                assertions.assertion_fingerprint
            ),
        },
    )

    with pytest.raises(DatasetValidationError, match="sequence"):
        _trusted_dataset(assertions=assertions, provenance=provenance)


def test_from_trusted_tables_rejects_altered_evidence_with_stale_fingerprint() -> None:
    assertions = _complete_assertions()
    payload = assertions.to_payload()
    identity = dict(_assertion_payload(payload["identity"]))
    identity["source"] = "tampered identity evidence"
    payload["identity"] = identity
    provenance = _provenance_with_construction_payload(
        assertions,
        {
            "trusted_construction_assertions": payload,
            "trusted_construction_assertion_fingerprint": (
                assertions.assertion_fingerprint
            ),
        },
    )

    with pytest.raises(DatasetValidationError, match="assertion_fingerprint"):
        _trusted_dataset(assertions=assertions, provenance=provenance)


def test_from_trusted_tables_rejects_altered_evidence_with_forged_fingerprint() -> None:
    assertions = _complete_assertions()
    altered_assertions = _complete_assertions(
        identity=TrustedDatasetConstructionEvidence.evidence(
            source="forged identity evidence",
            details={"site_key_schema": "forged"},
        )
    )
    provenance = _provenance_with_construction_payload(
        assertions,
        {
            "trusted_construction_assertions": altered_assertions.to_payload(),
            "trusted_construction_assertion_fingerprint": (
                altered_assertions.assertion_fingerprint
            ),
        },
    )

    with pytest.raises(DatasetValidationError, match="does not match"):
        _trusted_dataset(assertions=assertions, provenance=provenance)


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


def test_from_trusted_tables_rejects_incomplete_missing_assertion_bundle() -> None:
    with pytest.raises(
        DatasetValidationError,
        match="from_trusted_tables requires.*trusted_construction_assertions",
    ):
        _trusted_dataset(assertions=TrustedDatasetConstructionAssertions.missing())


def test_from_trusted_tables_accepts_and_serializes_localisation_waiver() -> None:
    assertions = _complete_assertions(
        localisation=TrustedDatasetConstructionEvidence.waiver(
            reason="historical source lacks localisation confidence export",
            policy="trusted_curation_waiver",
        )
    )

    dataset = _trusted_dataset(assertions)

    assert dataset.provenance is not None
    construction = _construction_payload(dataset.provenance)
    payload = construction["trusted_construction_assertions"]
    assert isinstance(payload, Mapping)
    localisation = _assertion_payload(payload["localisation"])
    assert localisation["kind"] == "waiver"
    assert localisation["waiver_reason"] == (
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


def test_from_trusted_tables_serializes_trusted_assertion_dimensions() -> None:
    dataset = _trusted_dataset()
    assert dataset.provenance is not None

    payload = to_payload(dataset.provenance)
    json.dumps(payload)
    restored = from_payload(payload)

    assert to_payload(restored) == payload
    workflow_parameters = _assertion_payload(payload["workflow_parameters"])
    construction = _assertion_payload(workflow_parameters["construction"])
    trusted_assertions = _assertion_payload(
        construction["trusted_construction_assertions"]
    )
    assert trusted_assertions["schema_version"] == 4
    assert trusted_assertions["missing_assertions"] == []
    for dimension in (
        "identity",
        "intensity_scale",
        "quantitative_meaning",
        "aligned_structure",
        "localisation",
        "sequence",
        "reference_context",
        "numeric_semantic_domain",
    ):
        assert trusted_assertions[dimension] is not None
        assert trusted_assertions[f"{dimension}_user_asserted"] is True


def test_trusted_construction_assertions_deserialize_canonical_round_trip() -> None:
    assertions = _complete_assertions(
        localisation=TrustedDatasetConstructionEvidence.waiver(
            reason="historical source lacks localisation confidence export",
            policy="trusted_curation_waiver",
            details={"approved_by": "unit-test"},
        )
    )

    restored = TrustedDatasetConstructionAssertions.from_payload(
        assertions.to_payload()
    )

    assert restored == assertions
    assert restored.to_payload() == assertions.to_payload()
    assert restored.assertion_fingerprint == assertions.assertion_fingerprint


def test_replayed_dataset_rejects_table_fingerprint_mismatch() -> None:
    assertions = _complete_assertions()
    trusted = _trusted_dataset(assertions)
    assert trusted.provenance is not None
    payload = to_payload(trusted.provenance)
    output_tables = payload["output_tables"]
    assert isinstance(output_tables, list)
    phospho_fingerprint = dict(_assertion_payload(output_tables[0]))
    phospho_fingerprint["exact_hash_value"] = "0" * 64
    output_tables[0] = phospho_fingerprint
    replayed_provenance = from_payload(payload)

    with pytest.raises(
        DatasetValidationError,
        match=r"run_provenance\.output_tables\.dataset\.phospho.*exact_hash_value",
    ):
        _trusted_dataset(assertions=assertions, provenance=replayed_provenance)


def test_replayed_trusted_construction_rejects_reference_organism_contradiction() -> (
    None
):
    assertions = _complete_assertions()
    trusted = _trusted_dataset(assertions)
    assert trusted.provenance is not None
    contradictory_provenance = replace(
        trusted.provenance,
        reference_context=_reference_context("human"),
        reference=ReferenceProvenance(
            source_type="explicit",
            organism=Organism.HUMAN,
            bundle_id=None,
            table_fingerprints=(),
            reference_context=_reference_context("human"),
        ),
    )
    replayed_provenance = from_payload(to_payload(contradictory_provenance))

    with pytest.raises(DatasetValidationError, match="dataset organism identity"):
        _trusted_dataset(assertions=assertions, provenance=replayed_provenance)


def test_direct_constructor_fails_immediately_and_names_supported_paths() -> None:
    index = _site_index()

    with pytest.raises(
        TypeError,
        match="AnalysisReadyDatasetBuilder.*from_trusted_tables",
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
    construction = _construction_payload(dataset.provenance)
    payload = _assertion_payload(construction["trusted_construction_assertions"])
    identity_payload = _assertion_payload(payload["identity"])
    assert identity_payload["details"] == {
        "columns": ["site_key"],
        "schema": {"version": 1},
    }
    assert assertions.identity is not None
    with pytest.raises(TypeError):
        cast(dict[str, object], assertions.identity.details)["schema"] = {"version": 2}
    assert dataset.trusted_construction_assertions is not None
    assert dataset.trusted_construction_assertions.identity is not None
    with pytest.raises(TypeError):
        cast(
            dict[str, object],
            dataset.trusted_construction_assertions.identity.details,
        )["schema"] = {"version": 3}


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
            trusted_construction_assertions=cast(
                TrustedDatasetConstructionAssertions,
                {
                    "sequence_user_asserted": True,
                    "identity_user_asserted": True,
                    "quantitative_meaning_user_asserted": True,
                    "reference_context_user_asserted": True,
                },
            ),
        )
