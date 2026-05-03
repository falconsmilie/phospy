from __future__ import annotations

import pandas as pd

from phospy import AnalysisReadyDatasetBuilder
from phospy.api import (
    DatasetBuildRequest,
    DatasetComparisonBuildingConfig,
    DatasetIntensityTransformConfig,
    DatasetMissingDataConfig,
    DatasetPreprocessingConfig,
    DatasetSiteMatrixConfig,
    DatasetTotalProteinCorrectionConfig,
    Organism,
)
from phospy.provenance.models import PreprocessingStageProvenance, TableFingerprint


def _base_phospho() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample_a": [15.0, 7.0],
            "sample_b": [31.0, 15.0],
        },
        index=pd.Index(["MAPK14;Y182;", "AKT1;T308;"], name="site_id"),
    )


def _base_site_metadata(index: pd.Index) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "AKT1"],
            "site": ["Y182", "T308"],
            "site_sequence": ["SEQ_A", "SEQ_B"],
        },
        index=index.copy(),
    )


def _base_sample_metadata(columns: pd.Index) -> pd.DataFrame:
    return pd.DataFrame(
        {"comparison_group": ["treated", "control"]},
        index=columns.copy(),
    )


def _base_total(columns: pd.Index) -> pd.DataFrame:
    return pd.DataFrame(
        {
            str(columns[0]): [3.0, 1.0],
            str(columns[1]): [7.0, 3.0],
        },
        index=pd.Index(["MAPK14", "AKT1"], name="protein_id"),
    )


def _build_dataset(
    *,
    phospho: pd.DataFrame | None = None,
    site_metadata: pd.DataFrame | None = None,
    sample_metadata: pd.DataFrame | None = None,
    total: pd.DataFrame | None = None,
    preprocessing_config: DatasetPreprocessingConfig,
):
    resolved_phospho = _base_phospho() if phospho is None else phospho
    resolved_site_metadata = (
        _base_site_metadata(resolved_phospho.index)
        if site_metadata is None
        else site_metadata
    )
    return AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=resolved_phospho,
            site_metadata=resolved_site_metadata,
            sample_metadata=sample_metadata,
            total=total,
            organism=Organism.RAT,
            preprocessing_config=preprocessing_config,
        )
    )


def _stage(
    dataset,
    stage_name: str,
) -> PreprocessingStageProvenance:
    assert dataset.provenance is not None
    return next(
        stage
        for stage in dataset.provenance.preprocessing_stages
        if stage.stage == stage_name
    )


def _hash_by_name(
    fingerprints: tuple[TableFingerprint, ...],
    table_name: str,
) -> str:
    return next(item.hash_value for item in fingerprints if item.name == table_name)


def test_site_metadata_change_updates_site_matrix_stage_provenance() -> None:
    phospho = _base_phospho()
    site_metadata_a = _base_site_metadata(phospho.index)
    site_metadata_b = site_metadata_a.copy(deep=True)
    site_metadata_b.loc["AKT1;T308;", "site"] = "S473"

    config = DatasetPreprocessingConfig(
        site_matrix=DatasetSiteMatrixConfig(policy="build_from_metadata")
    )
    built_a = _build_dataset(
        phospho=phospho,
        site_metadata=site_metadata_a,
        preprocessing_config=config,
    )
    built_b = _build_dataset(
        phospho=phospho,
        site_metadata=site_metadata_b,
        preprocessing_config=config,
    )

    stage_a = _stage(built_a, "site_matrix")
    stage_b = _stage(built_b, "site_matrix")
    assert _hash_by_name(
        stage_a.consumed_input_tables, "dataset.site_metadata"
    ) != _hash_by_name(stage_b.consumed_input_tables, "dataset.site_metadata")


def test_sample_metadata_change_updates_comparisons_stage_provenance() -> None:
    phospho = _base_phospho()
    sample_metadata_a = _base_sample_metadata(phospho.columns)
    sample_metadata_b = sample_metadata_a.copy(deep=True)
    sample_metadata_b.loc["sample_b", "comparison_group"] = "vehicle"

    config = DatasetPreprocessingConfig(
        comparisons=DatasetComparisonBuildingConfig(policy="sample_metadata_pairs")
    )
    built_a = _build_dataset(
        phospho=phospho,
        sample_metadata=sample_metadata_a,
        preprocessing_config=config,
    )
    built_b = _build_dataset(
        phospho=phospho,
        sample_metadata=sample_metadata_b,
        preprocessing_config=config,
    )

    stage_a = _stage(built_a, "comparisons")
    stage_b = _stage(built_b, "comparisons")
    assert _hash_by_name(
        stage_a.consumed_input_tables, "dataset.sample_metadata"
    ) != _hash_by_name(stage_b.consumed_input_tables, "dataset.sample_metadata")
    assert _hash_by_name(
        stage_a.produced_output_tables, "dataset.comparisons"
    ) != _hash_by_name(stage_b.produced_output_tables, "dataset.comparisons")


def test_total_input_change_updates_total_correction_stage_provenance() -> None:
    phospho = _base_phospho()
    total_a = _base_total(phospho.columns)
    total_b = total_a.copy(deep=True)
    total_b.loc["AKT1", "sample_b"] = 4.0
    config = DatasetPreprocessingConfig(
        intensity_transform=DatasetIntensityTransformConfig(
            policy="log2",
            pseudocount=1.0,
        ),
        total_protein_correction=DatasetTotalProteinCorrectionConfig(
            policy="subtract_log_total"
        ),
    )
    built_a = _build_dataset(
        phospho=phospho,
        total=total_a,
        preprocessing_config=config,
    )
    built_b = _build_dataset(
        phospho=phospho,
        total=total_b,
        preprocessing_config=config,
    )

    stage_a = _stage(built_a, "total_protein_correction")
    stage_b = _stage(built_b, "total_protein_correction")
    assert _hash_by_name(
        stage_a.consumed_input_tables, "dataset.total"
    ) != _hash_by_name(stage_b.consumed_input_tables, "dataset.total")


def test_stage_configuration_change_updates_stage_provenance() -> None:
    phospho = _base_phospho()
    config_a = DatasetPreprocessingConfig(
        intensity_transform=DatasetIntensityTransformConfig(
            policy="log2",
            pseudocount=1.0,
        )
    )
    config_b = DatasetPreprocessingConfig(
        intensity_transform=DatasetIntensityTransformConfig(
            policy="log2",
            pseudocount=2.0,
        )
    )
    built_a = _build_dataset(phospho=phospho, preprocessing_config=config_a)
    built_b = _build_dataset(phospho=phospho, preprocessing_config=config_b)

    stage_a = _stage(built_a, "intensity_transform")
    stage_b = _stage(built_b, "intensity_transform")
    assert stage_a.parameters != stage_b.parameters
    assert _hash_by_name(
        stage_a.produced_output_tables, "dataset.phospho"
    ) != _hash_by_name(stage_b.produced_output_tables, "dataset.phospho")


def test_missing_data_missingness_mask_hash_is_stable() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [15.0, float("nan")],
            "sample_b": [float("nan"), 15.0],
        },
        index=pd.Index(["MAPK14;Y182;", "AKT1;T308;"], name="site_id"),
    )
    config = DatasetPreprocessingConfig(
        missing_data=DatasetMissingDataConfig(
            policy="impute_row_median",
            min_observed_values=1,
        )
    )

    built_a = _build_dataset(phospho=phospho, preprocessing_config=config)
    built_b = _build_dataset(phospho=phospho, preprocessing_config=config)

    stage_a = _stage(built_a, "missing_data")
    stage_b = _stage(built_b, "missing_data")
    diagnostics_a = stage_a.diagnostics or {}
    diagnostics_b = stage_b.diagnostics or {}

    assert diagnostics_a.get("missingness_mask_hash") == diagnostics_b.get(
        "missingness_mask_hash"
    )
