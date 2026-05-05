from __future__ import annotations

import json

import pandas as pd
import pytest

from phospy import AnalysisReadyDatasetBuilder
from phospy.api import (
    DatasetBuildRequest,
    DatasetIntensityTransformConfig,
    DatasetMissingDataConfig,
    DatasetPreprocessingConfig,
    Organism,
)
from phospy.datasets.preprocessing.models import (
    PREPROCESSING_STAGE_ORDER_RATIONALE_MINPROB_MISSING_DATA,
    PREPROCESSING_STAGE_ORDER_RATIONALE_NON_MINPROB_MISSING_DATA,
)
from phospy.provenance import environment as provenance_environment
from phospy.provenance.environment import collect_environment_provenance
from phospy.provenance.serialization import from_payload, to_payload
from phospy.scientific_policies import ScientificPolicyId


def test_collect_environment_provenance_reports_expected_keys() -> None:
    environment = collect_environment_provenance()
    dependency_names = set(environment.dependency_versions)
    assert environment.package_name == "phospy"
    assert environment.python_version
    assert {"numpy", "pandas", "scipy", "scikit-learn"}.issubset(dependency_names)
    assert {"pyarrow", "openpyxl"}.issubset(dependency_names)
    assert {"platform", "system", "machine"}.issubset(set(environment.platform))


def test_collect_environment_provenance_tolerates_missing_optional_engines(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    optional_dependencies = {"pyarrow", "openpyxl"}
    real_lookup = provenance_environment._distribution_version

    def _version_lookup(distribution_name: str) -> str | None:
        if distribution_name in optional_dependencies:
            return None
        return real_lookup(distribution_name)

    monkeypatch.setattr(
        provenance_environment, "_distribution_version", _version_lookup
    )
    environment = collect_environment_provenance()

    for dependency_name in optional_dependencies:
        assert dependency_name in environment.dependency_versions
        assert environment.dependency_versions[dependency_name] is None


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
            "site_sequence": ["SEQ_A", "SEQ_B", "SEQ_C"],
        },
        index=phospho.index.copy(),
    )

    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
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
    assert missing_stage.is_deterministic is True
    assert missing_stage.random_seed is None
    assert "GSK3B;S9;" in set(missing_stage.dropped_row_ids)
    assert missing_stage.imputed_cell_count >= 1
    assert "AKT1;T308;" in set(missing_stage.imputed_row_ids)
    policy_ids = {policy.id for policy in provenance.scientific_policies}
    assert ScientificPolicyId.PREPROCESSING_STAGE_ORDER in policy_ids


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
        },
        index=phospho.index.copy(),
    )
    default_result = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
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
    assert default_policy.parameters["configured_stage_order"] == "missing_data"
    assert transformed_policy.parameters["configured_stage_order"] == (
        "missing_data -> intensity_transform"
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
            "site_sequence": ["SEQ_A", "SEQ_B", "SEQ_C"],
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
    assert isinstance(preprocessing_plan, dict)
    resolved_stage_order = preprocessing_plan["resolved_stage_order"]
    assert isinstance(resolved_stage_order, list)
    assert [item["stage"] for item in resolved_stage_order[:2]] == [
        "intensity_transform",
        "missing_data",
    ]
    assert resolved_stage_order[1]["rationale"] == (
        PREPROCESSING_STAGE_ORDER_RATIONALE_MINPROB_MISSING_DATA
    )


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
            "site_sequence": ["SEQ_A", "SEQ_B", "SEQ_C"],
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
    assert isinstance(preprocessing_plan, dict)
    resolved_stage_order = preprocessing_plan["resolved_stage_order"]
    assert isinstance(resolved_stage_order, list)
    assert [item["stage"] for item in resolved_stage_order[:2]] == [
        "missing_data",
        "intensity_transform",
    ]
    assert resolved_stage_order[0]["rationale"] == (
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
        },
        index=phospho.index.copy(),
    )
    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
        )
    )
    assert built.provenance is not None
    payload = to_payload(built.provenance)
    json.dumps(payload)
    restored = from_payload(payload)
    assert to_payload(restored) == payload


def test_run_provenance_from_payload_accepts_legacy_stage_shape() -> None:
    phospho = pd.DataFrame(
        {"sample_a": [1.0], "sample_b": [2.0]},
        index=pd.Index(["MAPK14;Y182;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14"],
            "site": ["Y182"],
            "site_sequence": ["SEQ_A"],
        },
        index=phospho.index.copy(),
    )
    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
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
    environment_payload.pop("platform", None)
    stages = payload["preprocessing_stages"]
    assert isinstance(stages, list)
    for stage in stages:
        assert isinstance(stage, dict)
        stage.pop("schema_version", None)
        stage.pop("consumed_input_tables", None)
        stage.pop("produced_output_tables", None)
        stage.pop("backend", None)
        stage.pop("random_seed", None)
        stage.pop("is_deterministic", None)
    payload.pop("scientific_policies", None)

    restored = from_payload(payload)
    stage = next(item for item in restored.preprocessing_stages if item.stage)
    assert stage.schema_version == 1
    assert stage.consumed_input_tables == ()
    assert stage.produced_output_tables == ()
    assert stage.backend is None
    assert stage.random_seed is None
    assert stage.is_deterministic is True
    assert restored.scientific_policies == ()
