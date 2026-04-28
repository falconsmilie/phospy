from __future__ import annotations

import pandas as pd

from phospy import AnalysisReadyDatasetBuilder
from phospy.api import (
    DatasetBuildRequest,
    DatasetMissingDataConfig,
    DatasetPreprocessingConfig,
    Organism,
)
from phospy.provenance.environment import collect_environment_provenance
from phospy.provenance.serialization import from_payload, to_payload


def test_collect_environment_provenance_reports_expected_keys() -> None:
    environment = collect_environment_provenance()
    assert environment.package_name == "phospy"
    assert environment.python_version
    assert {"numpy", "pandas", "scikit-learn"}.issubset(
        set(environment.dependency_versions.keys())
    )


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

    restored = from_payload(payload)
    stage = next(item for item in restored.preprocessing_stages if item.stage)
    assert stage.schema_version == 1
    assert stage.consumed_input_tables == ()
    assert stage.produced_output_tables == ()
    assert stage.backend is None
    assert stage.random_seed is None
    assert stage.is_deterministic is True
