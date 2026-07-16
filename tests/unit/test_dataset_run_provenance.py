from __future__ import annotations

import json
from collections.abc import Mapping, Sequence

import pandas as pd
import pytest

from phospy import AnalysisReadyDatasetBuilder, AnalysisReadyPhosphoDataset
from phospy.api import (
    DatasetBuildRequest,
    DatasetIntensityTransformConfig,
    DatasetMissingDataConfig,
    DatasetPreprocessingConfig,
    Organism,
)
from phospy.errors.input import PhosPyInputError
from phospy.provenance import environment as provenance_environment
from phospy.provenance.environment import (
    collect_batch_correction_environment_provenance,
    collect_environment_provenance,
)
from phospy.provenance.scientific_policy_models import ScientificPolicyId
from phospy.provenance.serialization import from_payload, to_payload
from phospy.science.datasets.preprocessing.models import (
    PREPROCESSING_STAGE_ORDER_RATIONALE_MINPROB_MISSING_DATA,
    PREPROCESSING_STAGE_ORDER_RATIONALE_NON_MINPROB_MISSING_DATA,
)
from phospy.science.sites.site_keys import (
    ProteinScopedPhosphositeKey,
    encode_site_key,
)
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
)


def _gene_site_key(
    *,
    organism: Organism,
    gene_symbol: str,
    residue: str,
    position: int,
) -> str:
    return encode_site_key(
        ProteinScopedPhosphositeKey(
            organism=organism.value,
            protein_namespace="protein_id",
            protein_identifier=gene_symbol,
            residue=residue,
            position=position,
        )
    )


def _direct_analysis_ready_dataset() -> AnalysisReadyPhosphoDataset:
    display_ids = ["MAPK14;Y182;", "AKT1;T308;"]
    site_keys = [
        _gene_site_key(
            organism=Organism.RAT,
            gene_symbol="MAPK14",
            residue="Y",
            position=182,
        ),
        _gene_site_key(
            organism=Organism.RAT,
            gene_symbol="AKT1",
            residue="T",
            position=308,
        ),
    ]
    site_index = pd.Index(site_keys, name="site_key")
    sample_index = pd.Index(["sample_a", "sample_b"], name="sample_id")
    return AnalysisReadyPhosphoDataset(
        phospho=pd.DataFrame(
            {
                "sample_a": [1.0, 2.0],
                "sample_b": [1.5, 2.5],
            },
            index=site_index.copy(),
        ),
        site_metadata=pd.DataFrame(
            {
                "site_key": site_keys,
                "display_id": display_ids,
                "organism": ["rat", "rat"],
                "protein_namespace": ["protein_id", "protein_id"],
                "protein_identifier": ["MAPK14", "AKT1"],
                "gene_symbol": ["MAPK14", "AKT1"],
                "site": ["Y182", "T308"],
                "site_sequence": [
                    "AAAAAAAAAAAAAAAYAAAAAAAAAAAAAAA",
                    "AAAAAAAAAAAAAAATAAAAAAAAAAAAAAA",
                ],
                "localisation_confidence": [0.95, 0.9],
            },
            index=site_index.copy(),
        ),
        sample_metadata=pd.DataFrame(
            {"condition": ["control", "treated"]},
            index=sample_index.copy(),
        ),
        comparisons=pd.DataFrame(
            {"treated_vs_control": [0.5, 0.7]},
            index=site_index.copy(),
        ),
        imputation_observation_mask=pd.DataFrame(
            {
                "sample_a": [True, False],
                "sample_b": [True, True],
            },
            index=site_index.copy(),
            columns=sample_index.copy(),
        ),
        organism=Organism.RAT,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
    )


def test_collect_environment_provenance_reports_expected_keys() -> None:
    provenance_environment.clear_environment_provenance_cache()
    environment = collect_environment_provenance(use_cache=False)
    dependency_names = set(environment.dependency_versions)
    assert environment.package_name == "phospy"
    assert environment.python_version
    assert {"numpy", "pandas", "scipy", "scikit-learn"}.issubset(dependency_names)
    assert {"pyarrow", "openpyxl"}.issubset(dependency_names)
    assert {"platform", "system", "machine"}.issubset(set(environment.platform))
    assert environment.schema_version >= 2
    assert {
        "blas_name",
        "lapack_name",
        "blas_version",
        "lapack_version",
    }.issubset(set(environment.blas_lapack))
    assert set(environment.thread_environment) == set(
        provenance_environment.THREAD_ENVIRONMENT_VARIABLES
    )
    assert {"language_code", "encoding", "lc_all", "preferred_encoding"}.issubset(
        set(environment.locale)
    )
    assert {"algorithm", "value", "sources"}.issubset(
        set(environment.constraints_fingerprint)
    )


def test_collect_environment_provenance_falls_back_to_project_metadata_version(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    provenance_environment.clear_environment_provenance_cache()

    def _version_lookup(distribution_name: str) -> str | None:
        if distribution_name == "phospy":
            return None
        return "test-version"

    monkeypatch.setattr(
        provenance_environment, "_distribution_version", _version_lookup
    )

    environment = collect_environment_provenance(
        dependency_names=(),
        use_cache=False,
    )

    assert environment.package_version == "1.6.0"
    assert environment.package_version != "unknown"


def test_collect_batch_correction_environment_tolerates_unavailable_dependency(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    provenance_environment.clear_environment_provenance_cache()
    real_lookup = provenance_environment._distribution_version

    def _version_lookup(distribution_name: str) -> str | None:
        if distribution_name == "scikit-learn":
            return None
        return real_lookup(distribution_name)

    monkeypatch.setattr(
        provenance_environment, "_distribution_version", _version_lookup
    )

    environment = collect_batch_correction_environment_provenance(use_cache=False)

    assert set(environment.dependency_versions) == {
        "numpy",
        "pandas",
        "scipy",
        "scikit-learn",
    }
    assert environment.dependency_versions["scikit-learn"] is None


def test_collect_environment_provenance_tolerates_missing_optional_engines(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    provenance_environment.clear_environment_provenance_cache()
    optional_dependencies = {"pyarrow", "openpyxl"}
    real_lookup = provenance_environment._distribution_version

    def _version_lookup(distribution_name: str) -> str | None:
        if distribution_name in optional_dependencies:
            return None
        return real_lookup(distribution_name)

    monkeypatch.setattr(
        provenance_environment, "_distribution_version", _version_lookup
    )
    environment = collect_environment_provenance(use_cache=False)

    for dependency_name in optional_dependencies:
        assert dependency_name in environment.dependency_versions
        assert environment.dependency_versions[dependency_name] is None


def test_collect_environment_provenance_tolerates_missing_optional_backend_metadata(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    provenance_environment.clear_environment_provenance_cache()
    monkeypatch.setattr(
        provenance_environment,
        "_blas_lapack_provenance",
        lambda: {
            "source": "numpy_backend_unavailable",
            "blas_name": None,
            "blas_version": None,
            "blas_detection_method": None,
            "blas_openblas_configuration": None,
            "lapack_name": None,
            "lapack_version": None,
            "lapack_detection_method": None,
            "lapack_openblas_configuration": None,
        },
    )
    environment = collect_environment_provenance(use_cache=False)
    assert environment.blas_lapack["source"] == "numpy_backend_unavailable"
    assert environment.blas_lapack["blas_name"] is None
    assert environment.blas_lapack["lapack_name"] is None


def test_collect_environment_provenance_captures_thread_environment_variables(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    provenance_environment.clear_environment_provenance_cache()
    monkeypatch.setenv("OMP_NUM_THREADS", "2")
    monkeypatch.setenv("OPENBLAS_NUM_THREADS", "3")
    monkeypatch.setenv("MKL_NUM_THREADS", "4")
    monkeypatch.setenv("NUMEXPR_NUM_THREADS", "5")

    environment = collect_environment_provenance(use_cache=False)

    assert environment.thread_environment == {
        "OMP_NUM_THREADS": "2",
        "OPENBLAS_NUM_THREADS": "3",
        "MKL_NUM_THREADS": "4",
        "NUMEXPR_NUM_THREADS": "5",
    }


def test_environment_provenance_payload_serializes_extended_metadata() -> None:
    provenance_environment.clear_environment_provenance_cache()
    environment = collect_environment_provenance(use_cache=False)
    payload = to_payload(
        from_payload(
            {
                "environment": {
                    "schema_version": environment.schema_version,
                    "package_name": environment.package_name,
                    "package_version": environment.package_version,
                    "python_version": environment.python_version,
                    "dependency_versions": environment.dependency_versions,
                    "platform": environment.platform,
                    "blas_lapack": environment.blas_lapack,
                    "thread_environment": environment.thread_environment,
                    "timezone": environment.timezone,
                    "locale": environment.locale,
                    "constraints_fingerprint": environment.constraints_fingerprint,
                },
                "input_tables": [],
                "preprocessing_stages": [],
                "reference": None,
                "workflow_name": None,
                "workflow_parameters": {},
                "random_state": None,
                "random_seed_policy": None,
                "output_tables": [],
                "scientific_policies": [],
            }
        )
    )
    environment_payload = payload["environment"]
    assert isinstance(environment_payload, dict)
    assert environment_payload["schema_version"] >= 2
    assert "blas_lapack" in environment_payload
    assert "thread_environment" in environment_payload
    assert "timezone" in environment_payload
    assert "locale" in environment_payload
    assert "constraints_fingerprint" in environment_payload


def test_collect_environment_provenance_uses_process_cache(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    provenance_environment.clear_environment_provenance_cache()
    real_collector = provenance_environment._collect_environment_provenance_uncached
    invocation_count = 0

    def _counting_collector(*, package_name: str, dependency_names: tuple[str, ...]):
        nonlocal invocation_count
        invocation_count += 1
        return real_collector(
            package_name=package_name,
            dependency_names=dependency_names,
        )

    monkeypatch.setattr(
        provenance_environment,
        "_collect_environment_provenance_uncached",
        _counting_collector,
    )
    first = collect_environment_provenance()
    second = collect_environment_provenance()
    assert invocation_count == 1
    assert first == second


def test_direct_dataset_construction_without_provenance_records_marker() -> None:
    dataset = _direct_analysis_ready_dataset()

    provenance = dataset.provenance

    assert provenance is not None
    assert provenance.workflow_name == "analysis_ready_dataset_direct_construction"
    assert provenance.preprocessing_stages == ()
    assert provenance.reference is None
    assert provenance.random_state is None
    assert provenance.random_seed_policy is None
    assert provenance.scientific_policies == ()
    construction = provenance.workflow_parameters["construction"]
    assert isinstance(construction, Mapping)
    assert construction["method"] == "AnalysisReadyPhosphoDataset.__init__"
    assert construction["source"] == "direct_trusted_construction"
    assert construction["builder_used"] is False
    assert construction["warning"] == (
        "Direct construction cannot prove biological correctness of "
        "caller-provided analysis-ready state."
    )
    assert construction["trusted_assertion_metadata_provided"] is False
    assert construction["missing_trusted_assertions"] == (
        "identity_user_asserted",
        "quantitative_meaning_user_asserted",
        "localisation_user_asserted",
        "sequence_user_asserted",
        "reference_context_user_asserted",
    )
    assertions = construction["trusted_construction_assertions"]
    assert isinstance(assertions, Mapping)
    assert assertions["assertion_metadata_provided"] is False
    assert {item.name for item in provenance.input_tables} == {
        "dataset.phospho",
        "dataset.site_metadata",
        "dataset.sample_metadata",
        "dataset.comparisons",
        "dataset.imputation_observation_mask",
    }
    assert provenance.input_tables == provenance.output_tables


def test_direct_dataset_construction_provenance_serializes_round_trip() -> None:
    dataset = _direct_analysis_ready_dataset()
    assert dataset.provenance is not None

    payload = to_payload(dataset.provenance)
    json.dumps(payload)
    restored = from_payload(payload)

    assert to_payload(restored) == payload
    assert payload["workflow_name"] == "analysis_ready_dataset_direct_construction"
    assert payload["preprocessing_stages"] == []
    workflow_parameters = payload["workflow_parameters"]
    assert isinstance(workflow_parameters, dict)
    construction = workflow_parameters["construction"]
    assert isinstance(construction, dict)
    trusted_assertions = construction["trusted_construction_assertions"]
    assert isinstance(trusted_assertions, dict)
    assert trusted_assertions["assertion_metadata_provided"] is False
    assert trusted_assertions["missing_assertions"] == [
        "identity_user_asserted",
        "quantitative_meaning_user_asserted",
        "localisation_user_asserted",
        "sequence_user_asserted",
        "reference_context_user_asserted",
    ]


def test_dataset_builder_emits_run_provenance_and_stage_details() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, float("nan"), float("nan")],
            "sample_b": [2.0, 2.0, float("nan")],
            "sample_c": [3.0, 3.0, float("nan")],
        },
        index=pd.Index(["MAPK14;Y182;", "AKT1;T308;", "GSK3B;S9;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "AKT1", "GSK3B"],
            "site": ["Y182", "T308", "S9"],
            "site_sequence": ["SEQ_A", "SEQ_R", "SEQ_C"],
            "protein_id": ["MAPK14", "AKT1", "GSK3B"],
            "localisation_confidence": [0.95, 0.9, 0.92],
        },
        index=phospho.index.copy(),
    )

    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
            input_intensity_scale="linear",
            preprocessing_config=DatasetPreprocessingConfig(
                missing_data=DatasetMissingDataConfig(
                    policy="impute_row_median",
                    min_observed_values=2,
                )
            ),
        )
    )

    provenance = built.provenance
    assert provenance is not None
    assert provenance.workflow_name == "dataset_builder"
    assert len(provenance.input_tables) >= 2
    assert len(provenance.output_tables) >= 2
    first_input_fingerprint = provenance.input_tables[0]
    assert first_input_fingerprint.exact_hash_algorithm == "sha256-stable-json-v1"
    assert isinstance(first_input_fingerprint.exact_hash_value, str)
    assert len(first_input_fingerprint.exact_hash_value) == 64
    assert (
        first_input_fingerprint.tolerance_hash_algorithm == "sha256-float-round-8dp-v1"
    )
    assert isinstance(first_input_fingerprint.tolerance_hash_value, str)
    assert len(first_input_fingerprint.tolerance_hash_value) == 64
    assert provenance.preprocessing_stages
    missing_stage = next(
        stage
        for stage in provenance.preprocessing_stages
        if stage.stage == "missing_data"
    )
    assert missing_stage.operation == "impute_row_median"
    assert missing_stage.schema_version >= 2
    consumed_names = {item.name for item in missing_stage.consumed_input_tables}
    produced_names = {item.name for item in missing_stage.produced_output_tables}
    assert "dataset.phospho" in consumed_names
    assert "dataset.site_metadata" in consumed_names
    assert "dataset.phospho" in produced_names
    assert missing_stage.backend in {"pandas", "numpy", None}
    assert missing_stage.determinism == "deterministic"
    assert missing_stage.random_seed is None
    assert isinstance(missing_stage.phospho_input_hash, str)
    assert isinstance(missing_stage.phospho_output_hash, str)
    assert missing_stage.input_hash != missing_stage.phospho_input_hash
    assert missing_stage.output_hash != missing_stage.phospho_output_hash
    assert _gene_site_key(
        organism=Organism.RAT,
        gene_symbol="GSK3B",
        residue="S",
        position=9,
    ) in set(missing_stage.dropped_row_ids)
    assert missing_stage.imputed_cell_count >= 1
    assert _gene_site_key(
        organism=Organism.RAT,
        gene_symbol="AKT1",
        residue="T",
        position=308,
    ) in set(missing_stage.imputed_row_ids)
    policy_ids = {policy.id for policy in provenance.scientific_policies}
    assert ScientificPolicyId.PREPROCESSING_STAGE_ORDER in policy_ids


def test_dataset_builder_records_reference_context_from_bundled_sequence_derivation() -> (
    None
):
    phospho = pd.DataFrame(
        {"sample_a": [1.0], "sample_b": [2.0]},
        index=pd.Index(["AAK1;S677;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["AAK1"],
            "site": ["S677"],
            "protein_id": ["AAK1"],
            "localisation_confidence": [0.95],
        },
        index=phospho.index.copy(),
    )

    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
            input_intensity_scale="linear",
        )
    )

    provenance = built.provenance
    assert provenance is not None
    assert provenance.reference is None
    assert provenance.reference_context is not None
    assert built.reference_context is provenance.reference_context
    assert provenance.reference_context.organism == Organism.RAT.value
    assert provenance.reference_context.reference_context_id.startswith(
        "reference-context-v1:"
    )
    restored = from_payload(to_payload(provenance))
    assert restored.reference_context is not None
    assert (
        restored.reference_context.reference_context_id
        == provenance.reference_context.reference_context_id
    )


def test_dataset_builder_records_unknown_reference_context_when_not_reference_derived() -> (
    None
):
    phospho = pd.DataFrame(
        {"sample_a": [1.0], "sample_b": [2.0]},
        index=pd.Index(["MAPK14;Y182;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14"],
            "site": ["Y182"],
            "site_sequence": ["AAAAAAAAAAAAAAAYAAAAAAAAAAAAAAA"],
            "protein_id": ["MAPK14"],
            "localisation_confidence": [0.95],
        },
        index=phospho.index.copy(),
    )

    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
            input_intensity_scale="linear",
        )
    )

    assert built.provenance is not None
    assert built.provenance.reference_context is None
    assert built.reference_context is None
    payload = to_payload(built.provenance)
    assert "reference_context" in payload
    assert payload["reference_context"] is None


def test_dataset_builder_provenance_records_site_identifier_normalisation_changes() -> (
    None
):
    phospho = pd.DataFrame(
        {" sample_a ": [1.0], " sample_b ": [2.0]},
        index=pd.Index([" mapk14 ; y182 "], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["mapk14"],
            "site": ["y182"],
            "site_sequence": ["SEQ_A"],
            "protein_id": ["mapk14"],
            "localisation_confidence": [0.95],
        },
        index=pd.Index([" mapk14 ; y182 "], name="site_id"),
    )

    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
            input_intensity_scale="linear",
        )
    )

    assert built.provenance is not None
    normalisation_payload = built.provenance.workflow_parameters.get(
        "site_identifier_normalisation"
    )
    assert isinstance(normalisation_payload, Mapping)
    assert normalisation_payload["changed_identifier_count"] >= 2
    records = normalisation_payload["records"]
    assert isinstance(records, Sequence)
    assert not isinstance(records, (str, bytes, bytearray))
    fields = {record["field_name"] for record in records}
    assert "dataset build request phospho.index" in fields
    assert "dataset build request site_metadata.index" in fields
    assert all(record["normalised_value"] == "MAPK14;Y182;" for record in records)


def test_dataset_stage_order_policy_changes_with_preprocessing_plan() -> None:
    phospho = pd.DataFrame(
        {"sample_a": [1.0], "sample_b": [2.0]},
        index=pd.Index(["MAPK14;Y182;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14"],
            "site": ["Y182"],
            "site_sequence": ["SEQ_A"],
            "protein_id": ["MAPK14"],
            "localisation_confidence": [0.95],
        },
        index=phospho.index.copy(),
    )
    default_result = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
            input_intensity_scale="linear",
        )
    )
    transformed_result = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
            preprocessing_config=DatasetPreprocessingConfig(
                intensity_transform=DatasetIntensityTransformConfig(
                    policy="log2",
                    pseudocount=1.0,
                )
            ),
        )
    )

    assert default_result.provenance is not None
    assert transformed_result.provenance is not None
    default_policy = next(
        policy
        for policy in default_result.provenance.scientific_policies
        if policy.id == ScientificPolicyId.PREPROCESSING_STAGE_ORDER
    )
    transformed_policy = next(
        policy
        for policy in transformed_result.provenance.scientific_policies
        if policy.id == ScientificPolicyId.PREPROCESSING_STAGE_ORDER
    )
    assert (
        default_policy.parameters["configured_stage_order"]
        == "localisation_confidence -> missing_data"
    )
    assert transformed_policy.parameters["configured_stage_order"] == (
        "localisation_confidence -> missing_data -> intensity_transform"
    )


def test_run_provenance_serializes_resolved_stage_order_for_minprob_with_log2() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, float("nan"), 2.0],
            "sample_b": [2.0, 3.0, float("nan")],
            "sample_c": [4.0, 5.0, 6.0],
        },
        index=pd.Index(["MAPK14;Y182;", "AKT1;T308;", "GSK3B;S9;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "AKT1", "GSK3B"],
            "site": ["Y182", "T308", "S9"],
            "site_sequence": ["SEQ_A", "SEQ_R", "SEQ_C"],
            "protein_id": ["MAPK14", "AKT1", "GSK3B"],
            "localisation_confidence": [0.95, 0.9, 0.92],
        },
        index=phospho.index.copy(),
    )
    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
            preprocessing_config=DatasetPreprocessingConfig(
                intensity_transform=DatasetIntensityTransformConfig(
                    policy="log2",
                    pseudocount=1.0,
                ),
                missing_data=DatasetMissingDataConfig(
                    policy="impute_minprob",
                    q=0.01,
                    width=0.3,
                    seed=123,
                    max_missing_fraction_per_row=0.5,
                ),
            ),
        )
    )
    assert built.provenance is not None
    preprocessing_plan = built.provenance.workflow_parameters["preprocessing_plan"]
    assert isinstance(preprocessing_plan, Mapping)
    resolved_stage_order = preprocessing_plan["resolved_stage_order"]
    assert isinstance(resolved_stage_order, Sequence)
    assert not isinstance(resolved_stage_order, (str, bytes, bytearray))
    assert [item["stage"] for item in resolved_stage_order[:3]] == [
        "localisation_confidence",
        "intensity_transform",
        "missing_data",
    ]
    assert resolved_stage_order[2]["rationale"] == (
        PREPROCESSING_STAGE_ORDER_RATIONALE_MINPROB_MISSING_DATA
    )
    missing_stage = next(
        stage
        for stage in built.provenance.preprocessing_stages
        if stage.stage == "missing_data"
    )
    assert missing_stage.determinism == "seeded_stochastic"
    assert missing_stage.random_seed == 123


def test_run_provenance_serializes_resolved_stage_order_for_non_minprob_with_log2() -> (
    None
):
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, float("nan"), 2.0],
            "sample_b": [2.0, 3.0, 4.0],
            "sample_c": [4.0, 5.0, 6.0],
        },
        index=pd.Index(["MAPK14;Y182;", "AKT1;T308;", "GSK3B;S9;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "AKT1", "GSK3B"],
            "site": ["Y182", "T308", "S9"],
            "site_sequence": ["SEQ_A", "SEQ_R", "SEQ_C"],
            "protein_id": ["MAPK14", "AKT1", "GSK3B"],
            "localisation_confidence": [0.95, 0.9, 0.92],
        },
        index=phospho.index.copy(),
    )
    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
            preprocessing_config=DatasetPreprocessingConfig(
                intensity_transform=DatasetIntensityTransformConfig(
                    policy="log2",
                    pseudocount=1.0,
                ),
                missing_data=DatasetMissingDataConfig(
                    policy="impute_row_median",
                    min_observed_values=1,
                ),
            ),
        )
    )
    assert built.provenance is not None
    preprocessing_plan = built.provenance.workflow_parameters["preprocessing_plan"]
    assert isinstance(preprocessing_plan, Mapping)
    resolved_stage_order = preprocessing_plan["resolved_stage_order"]
    assert isinstance(resolved_stage_order, Sequence)
    assert not isinstance(resolved_stage_order, (str, bytes, bytearray))
    assert [item["stage"] for item in resolved_stage_order[:3]] == [
        "localisation_confidence",
        "missing_data",
        "intensity_transform",
    ]
    assert resolved_stage_order[1]["rationale"] == (
        PREPROCESSING_STAGE_ORDER_RATIONALE_NON_MINPROB_MISSING_DATA
    )


def test_run_provenance_serialization_round_trip_preserves_payload() -> None:
    phospho = pd.DataFrame(
        {"sample_a": [1.0], "sample_b": [2.0]},
        index=pd.Index(["MAPK14;Y182;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14"],
            "site": ["Y182"],
            "site_sequence": ["SEQ_A"],
            "protein_id": ["MAPK14"],
            "localisation_confidence": [0.95],
        },
        index=phospho.index.copy(),
    )
    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
            input_intensity_scale="linear",
        )
    )
    assert built.provenance is not None
    payload = to_payload(built.provenance)
    input_tables_payload = payload["input_tables"]
    assert isinstance(input_tables_payload, list)
    assert input_tables_payload
    first_input = input_tables_payload[0]
    assert isinstance(first_input, dict)
    assert first_input["exact_hash_algorithm"] == "sha256-stable-json-v1"
    assert first_input["tolerance_hash_algorithm"] == "sha256-float-round-8dp-v1"
    assert isinstance(first_input["exact_hash_value"], str)
    assert isinstance(first_input["tolerance_hash_value"], str)
    json.dumps(payload)
    restored = from_payload(payload)
    assert to_payload(restored) == payload


def test_run_provenance_from_payload_rejects_legacy_stage_shape() -> None:
    phospho = pd.DataFrame(
        {"sample_a": [1.0], "sample_b": [2.0]},
        index=pd.Index(["MAPK14;Y182;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14"],
            "site": ["Y182"],
            "site_sequence": ["SEQ_A"],
            "protein_id": ["MAPK14"],
            "localisation_confidence": [0.95],
        },
        index=phospho.index.copy(),
    )
    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
            input_intensity_scale="linear",
            preprocessing_config=DatasetPreprocessingConfig(
                missing_data=DatasetMissingDataConfig(
                    policy="impute_row_median",
                    min_observed_values=1,
                )
            ),
        )
    )
    assert built.provenance is not None
    payload = to_payload(built.provenance)
    environment_payload = payload["environment"]
    assert isinstance(environment_payload, dict)
    environment_payload.pop("schema_version", None)
    environment_payload.pop("platform", None)
    environment_payload.pop("blas_lapack", None)
    environment_payload.pop("thread_environment", None)
    environment_payload.pop("timezone", None)
    environment_payload.pop("locale", None)
    environment_payload.pop("constraints_fingerprint", None)
    stages = payload["preprocessing_stages"]
    assert isinstance(stages, list)
    for stage in stages:
        assert isinstance(stage, dict)
        stage.pop("schema_version", None)
        stage.pop("consumed_input_tables", None)
        stage.pop("produced_output_tables", None)
        stage.pop("backend", None)
        stage.pop("random_seed", None)
        stage.pop("determinism", None)
        stage.pop("phospho_input_hash", None)
        stage.pop("phospho_output_hash", None)
    payload.pop("scientific_policies", None)

    with pytest.raises(
        PhosPyInputError,
        match=(
            "Legacy provenance schemas are no longer supported. "
            "Regenerate the result with the current PhosPy version."
        ),
    ):
        from_payload(payload)
