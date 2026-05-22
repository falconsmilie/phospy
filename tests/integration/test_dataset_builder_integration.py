from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pandas.testing as pdt
import pytest

from phospy import (
    AnalysisReadyDatasetBuilder,
)
from phospy.api import (
    DATASET_TOTAL_PROTEIN_CORRECTION_UNMATCHED_POLICY_ALLOW_UNCORRECTED,
    DatasetBuildRequest,
    DatasetComparisonBuildingConfig,
    DatasetIntensityTransformConfig,
    DatasetMissingDataConfig,
    DatasetNormalisationConfig,
    DatasetPreprocessingConfig,
    DatasetRuvReadinessConfig,
    DatasetSiteMatrixConfig,
    DatasetSiteSequenceResolutionConfig,
    DatasetTotalProteinCorrectionConfig,
    DatasetTotalProteinCorrectionIdentityConfig,
    Organism,
    ReferencePreset,
)
from phospy.errors import PhosPyInputError
from phospy.io.publishers.workflows import publish_dataset
from phospy.science.references.resolution import ReferenceResolver
from phospy.science.transformations.models import (
    IntensityScaleState,
    MatrixIntensityScaleState,
)
from tests.support.rewrite_fixture_data import load_rat_l6_phospho, site_metadata_for

pytestmark = pytest.mark.integration


def test_dataset_builder_populates_preprocessing_report_for_successful_build() -> None:
    phospho = load_rat_l6_phospho().head(16).copy(deep=True)
    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata_for(phospho),
            organism=Organism.RAT,
            input_intensity_scale="linear",
        )
    )

    pdt.assert_frame_equal(built.phospho, phospho)
    report = built.preprocessing_report
    assert report is not None
    assert report.row_counts.shape[0] >= 1
    assert report.operations.shape[0] >= 1
    assert report.row_audit is not None
    assert {
        "stage",
        "input_rows",
        "output_rows",
        "dropped_rows",
    }.issubset(set(report.row_counts.columns))
    assert {
        "step_order",
        "stage",
        "operation",
        "parameters",
        "input_rows",
        "output_rows",
        "notes",
    }.issubset(set(report.operations.columns))
    assert {
        "stage",
        "action",
        "reason",
        "source_row_id",
        "site_id",
        "retained",
        "retained_row_id",
        "source_rows",
        "retained_row",
        "parameter_snapshot",
    }.issubset(set(report.row_audit.columns))
    assert report.duplicate_site_resolution is not None
    assert report.metadata_conflicts is not None
    assert report.comparison_group_stats is not None
    assert report.comparison_pair_stats is not None
    assert {
        "site_id",
        "source_row_id",
        "retained",
        "resolution_policy",
        "retained_reason",
        "dropped_reason",
        "observed_values",
        "mean_signal",
    }.issubset(set(report.duplicate_site_resolution.columns))
    assert {
        "site_id",
        "field",
        "values",
        "n_distinct_values",
        "source_row_ids",
    }.issubset(set(report.metadata_conflicts.columns))
    assert {"site_id", "group", "n", "mean", "sd", "sem"}.issubset(
        set(report.comparison_group_stats.columns)
    )
    assert {
        "site_id",
        "comparison",
        "left_n",
        "right_n",
        "left_mean",
        "right_mean",
        "left_sd",
        "right_sd",
        "left_sem",
        "right_sem",
        "effect_size",
    }.issubset(set(report.comparison_pair_stats.columns))
    assert report.comparison_group_stats.empty
    assert report.comparison_pair_stats.empty
    assert "missing_data" in set(report.row_counts.loc[:, "stage"])
    assert "site_matrix" in set(report.row_counts.loc[:, "stage"])
    assert "final_dataset_construction" in set(report.row_counts.loc[:, "stage"])


def test_dataset_builder_builds_analysis_ready_dataset_from_fixture() -> None:
    phospho = load_rat_l6_phospho().head(32).copy(deep=True)
    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata_for(phospho),
            organism=Organism.RAT,
            input_intensity_scale="linear",
        )
    )
    pdt.assert_frame_equal(built.phospho, phospho)
    assert list(built.site_metadata.columns) == [
        "gene_symbol",
        "site",
        "site_sequence",
        "localisation_confidence",
    ]
    assert built.intensity_scale_state == IntensityScaleState(
        phospho=MatrixIntensityScaleState.linear(
            established_by=(
                "phospy.science.datasets.builders.executor.input_intensity_scale"
            )
        ),
        total=None,
    )
    assert built.processing_state.intensity_scale == built.intensity_scale_state
    assert built.processing_state.missing_data.policy == "forbid"
    assert built.processing_state.missing_data.imputed is False
    assert built.processing_state.missing_data.complete_matrix is True
    assert built.processing_state.normalisation.policy == "none"
    assert built.processing_state.total_protein_correction.policy == "none"
    assert built.processing_state.total_protein_correction.applied is False
    assert built.processing_state.total_protein_correction.formula is None
    assert built.processing_state.total_protein_correction.requires_log_scale is False
    assert built.processing_state.total_protein_correction.input_scale is None
    assert built.processing_state.total_protein_correction.output_scale is None
    assert (
        built.processing_state.total_protein_correction.quantitative_meaning
        == "phosphosite_abundance"
    )
    correction_diagnostics = built.processing_state.total_protein_correction.diagnostics
    assert correction_diagnostics is not None
    assert correction_diagnostics.get("diagnostics_schema_version") == 1
    assert correction_diagnostics.get("quantitative_meaning") == "phosphosite_abundance"
    assert built.processing_state.site_matrix.policy == "as_input"
    assert built.processing_state.site_matrix.constructed is False
    assert built.processing_state.comparisons.policy == "none"
    assert built.processing_state.comparisons.pairs is None


def test_dataset_builder_reports_ruv_readiness_without_mutating_phospho_values() -> (
    None
):
    phospho = load_rat_l6_phospho().head(8).copy(deep=True)
    site_metadata = site_metadata_for(phospho).copy(deep=True)

    baseline = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho.copy(deep=True),
            site_metadata=site_metadata.copy(deep=True),
            organism=Organism.RAT,
            input_intensity_scale="linear",
        )
    )
    with_readiness = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho.copy(deep=True),
            site_metadata=site_metadata.copy(deep=True),
            organism=Organism.RAT,
            input_intensity_scale="linear",
            preprocessing_config=DatasetPreprocessingConfig(
                ruv_readiness=DatasetRuvReadinessConfig(enabled=True),
            ),
        )
    )

    pdt.assert_frame_equal(with_readiness.phospho, baseline.phospho)
    assert with_readiness.processing_state.ruv_readiness.enabled is True
    assert with_readiness.processing_state.ruv_readiness.ready is False
    assert "control feature column missing" in set(
        with_readiness.processing_state.ruv_readiness.reasons
    )


def test_dataset_builder_establishes_intensity_scale_state_via_supported_path() -> None:
    phospho = load_rat_l6_phospho().head(8).copy(deep=True)
    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata_for(phospho),
            organism=Organism.RAT,
            input_intensity_scale="linear",
        )
    )
    assert built.intensity_scale_state.label == "linear"
    assert built.intensity_scale_state.is_established
    assert built.intensity_scale_state.established_via is not None
    assert (
        built.intensity_scale_state.phospho.established_by
        == "phospy.science.datasets.builders.executor.input_intensity_scale"
    )


def test_dataset_builder_preserves_total_matrix_and_establishes_linear_state() -> None:
    phospho = load_rat_l6_phospho().head(8).copy(deep=True)
    total = pd.DataFrame(
        {
            sample_name: [float(i + 1), float(i + 2)]
            for i, sample_name in enumerate(phospho.columns.astype(str))
        },
        index=pd.Index(["MAPK14", "AKT1"], name="protein_id"),
    )
    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata_for(phospho),
            total=total,
            organism=Organism.RAT,
            input_intensity_scale="linear",
        )
    )
    assert built.total is not None
    pdt.assert_frame_equal(built.total, total)
    assert built.intensity_scale_state == IntensityScaleState(
        phospho=MatrixIntensityScaleState.linear(
            established_by=(
                "phospy.science.datasets.builders.executor.input_intensity_scale"
            )
        ),
        total=MatrixIntensityScaleState.linear(
            established_by=(
                "phospy.science.datasets.builders.executor.input_intensity_scale"
            )
        ),
    )
    assert built.intensity_scale_state.label == "linear"


def test_dataset_builder_applies_subtract_log_total_after_log2_transform() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [15.0, 7.0, 3.0],
            "sample_b": [31.0, 15.0, 7.0],
        },
        index=pd.Index(
            ["MAPK14;Y182;", "AKT1;T308;", "GSK3B;S9;"],
            name="site_id",
        ),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "AKT1", "GSK3B"],
            "site": ["Y182", "T308", "S9"],
            "site_sequence": ["SEQ_A", "SEQ_R", "SEQ_C"],
        },
        index=phospho.index.copy(),
    )
    site_metadata.loc[:, "localisation_confidence"] = [0.95] * site_metadata.shape[0]
    gene_symbols = site_metadata.loc[:, "gene_symbol"].astype(str)
    unique_genes = pd.Index(["MAPK14", "AKT1", "GSK3B"], name="protein_id")
    total = pd.DataFrame(
        {
            "sample_a": [3.0, 1.0, 1.0],
            "sample_b": [7.0, 3.0, 1.0],
        },
        index=unique_genes,
    )

    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            total=total,
            organism=Organism.RAT,
            preprocessing_config=DatasetPreprocessingConfig(
                intensity_transform=DatasetIntensityTransformConfig(
                    policy="log2",
                    pseudocount=1.0,
                ),
                total_protein_correction=DatasetTotalProteinCorrectionConfig(
                    policy="subtract_log_total"
                ),
            ),
        )
    )

    total_by_site = total.reindex(gene_symbols.tolist())
    total_by_site.index = phospho.index
    expected = np.log2(phospho + 1.0) - np.log2(total_by_site + 1.0)
    pdt.assert_frame_equal(built.phospho, expected)
    assert built.total is not None
    expected_total = np.log2(total + 1.0)
    pdt.assert_frame_equal(built.total, expected_total)
    assert built.intensity_scale_state.label == "log2"
    assert built.intensity_scale_state.kind.value == "log2"
    assert built.intensity_scale_state.quantity.value == "phospho_total_log_ratio"
    assert (
        built.processing_state.total_protein_correction.policy == "subtract_log_total"
    )
    assert built.processing_state.total_protein_correction.applied is True
    assert (
        built.processing_state.total_protein_correction.formula
        == "log2_phospho - log2_total"
    )
    assert built.processing_state.total_protein_correction.requires_log_scale is True
    assert built.processing_state.total_protein_correction.input_scale == "log2"
    assert built.processing_state.total_protein_correction.output_scale == "log2_ratio"
    assert (
        built.processing_state.total_protein_correction.quantitative_meaning
        == "phospho_total_log_ratio"
    )
    assert (
        built.processing_state.intensity_scale.quantity.value
        == "phospho_total_log_ratio"
    )
    correction_diagnostics = (
        built.processing_state.total_protein_correction.diagnostics or {}
    )
    assert correction_diagnostics.get("quantitative_meaning") == (
        "phospho_total_log_ratio"
    )
    assert correction_diagnostics.get("diagnostics_schema_version") == 1
    assert correction_diagnostics.get("matched_rows") == 3
    assert isinstance(correction_diagnostics.get("input_phospho_hash"), str)
    assert isinstance(correction_diagnostics.get("output_phospho_hash"), str)
    assert built.preprocessing_report is not None
    operations = built.preprocessing_report.operations
    intensity_stage_order = int(
        operations.loc[
            operations.loc[:, "stage"] == "intensity_transform",
            "step_order",
        ].iloc[0]
    )
    total_stage_order = int(
        operations.loc[
            operations.loc[:, "stage"] == "total_protein_correction",
            "step_order",
        ].iloc[0]
    )
    assert intensity_stage_order < total_stage_order
    assert built.provenance is not None
    total_stage = next(
        stage
        for stage in built.provenance.preprocessing_stages
        if stage.stage == "total_protein_correction"
    )
    diagnostics = total_stage.diagnostics or {}
    assert diagnostics["policy"] == "subtract_log_total"
    assert diagnostics["requested_policy"] == "subtract_log_total"
    assert diagnostics["resolved_policy"] == "subtract_log_total"
    assert diagnostics["formula"] == "log2_phospho - log2_total"
    assert diagnostics["requires_log_scale"] is True
    assert diagnostics["input_scale"] == "log2"
    assert diagnostics["output_scale"] == "log2_ratio"
    assert diagnostics["diagnostics_schema_version"] == 1
    assert diagnostics["quantitative_meaning"] == "phospho_total_log_ratio"
    assert isinstance(diagnostics.get("total_table_hash"), str)
    assert isinstance(diagnostics.get("input_phospho_hash"), str)
    assert isinstance(diagnostics.get("output_phospho_hash"), str)
    final_operation = operations.loc[
        operations.loc[:, "stage"] == "final_dataset_construction"
    ].iloc[0]
    assert final_operation["parameters"]["quantitative_meaning"] == (
        "phospho_total_log_ratio"
    )
    assert built.provenance is not None
    assert (
        built.provenance.workflow_parameters.get("quantitative_meaning")
        == "phospho_total_log_ratio"
    )


def test_dataset_builder_marks_mixed_quantitative_meaning_when_uncorrected_rows_are_retained(
    tmp_path: Path,
) -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [15.0, 7.0],
            "sample_b": [31.0, 15.0],
        },
        index=pd.Index(["MAPK14;Y182;", "AKT1;T308;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "AKT1"],
            "site": ["Y182", "T308"],
            "site_sequence": ["SEQ_A", "SEQ_R"],
            "protein_id": ["MAPK14", "AKT1"],
        },
        index=phospho.index.copy(),
    )
    site_metadata.loc[:, "localisation_confidence"] = [0.95] * site_metadata.shape[0]
    total = pd.DataFrame(
        {
            "sample_a": [3.0],
            "sample_b": [7.0],
        },
        index=pd.Index(["MAPK14"], name="protein_id"),
    )

    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            total=total,
            preprocessing_config=DatasetPreprocessingConfig(
                intensity_transform=DatasetIntensityTransformConfig(
                    policy="log2",
                    pseudocount=1.0,
                ),
                total_protein_correction=DatasetTotalProteinCorrectionConfig(
                    policy="subtract_log_total",
                    identity=DatasetTotalProteinCorrectionIdentityConfig(
                        unmatched_policy=DATASET_TOTAL_PROTEIN_CORRECTION_UNMATCHED_POLICY_ALLOW_UNCORRECTED
                    ),
                ),
            ),
        )
    )
    mixed_meaning = "mixed_phospho_total_log_ratio_and_phosphosite_log_abundance"
    assert built.intensity_scale_state.quantity is not None
    assert built.intensity_scale_state.quantity.value == mixed_meaning
    correction = built.processing_state.total_protein_correction
    assert correction.quantitative_meaning == mixed_meaning
    assert correction.diagnostics is not None
    assert correction.diagnostics.get("uncorrected_row_count") == 1
    assert correction.diagnostics.get("corrected_row_count") == 1
    assert correction.diagnostics.get("unmatched_policy") == "allow_uncorrected"
    assert correction.diagnostics.get("corrected_phosphosite_row_ids") == [
        "MAPK14;Y182;"
    ]
    assert correction.diagnostics.get(
        "corrected_phosphosite_to_total_protein_row_id"
    ) == {"MAPK14;Y182;": "MAPK14"}
    assert correction.diagnostics.get("uncorrected_phosphosite_row_reasons") == {
        "AKT1;T308;": (
            "no_matching_total_protein_row_retained_by_"
            "unmatched_policy_allow_uncorrected"
        )
    }

    written = publish_dataset(
        built,
        tmp_path / "published_mixed",
        output_format="csv",
    )
    manifest = json.loads(written["dataset.manifest"].read_text(encoding="utf-8"))
    assert manifest["quantitative_meaning"] == mixed_meaning
    assert manifest["processing_state"]["intensity_scale"]["quantity"] == mixed_meaning
    correction_payload = manifest["processing_state"]["total_protein_correction"]
    assert correction_payload["quantitative_meaning"] == mixed_meaning
    assert correction_payload["diagnostics"]["quantitative_meaning"] == mixed_meaning


def test_dataset_builder_requires_total_when_subtract_log_total_is_requested() -> None:
    phospho = load_rat_l6_phospho().head(4).copy(deep=True)
    with pytest.raises(
        PhosPyInputError,
        match="policy='subtract_log_total' requires total input data",
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=site_metadata_for(phospho),
                organism=Organism.RAT,
                preprocessing_config=DatasetPreprocessingConfig(
                    intensity_transform=DatasetIntensityTransformConfig(
                        policy="log2",
                        pseudocount=1.0,
                    ),
                    total_protein_correction=DatasetTotalProteinCorrectionConfig(
                        policy="subtract_log_total"
                    ),
                ),
            )
        )


def test_dataset_builder_rejects_removed_ratio_to_total_alias() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="preprocessing_config.total_protein_correction.policy must be one of",
    ):
        DatasetTotalProteinCorrectionConfig(policy="ratio_to_total")  # type: ignore[arg-type]


def test_dataset_builder_supports_documented_alias_and_index_derivation_conventions() -> (
    None
):
    phospho = load_rat_l6_phospho().head(8).copy(deep=True)
    canonical_site_metadata = site_metadata_for(phospho)
    site_metadata = pd.DataFrame(
        {
            "centralized_sequence": canonical_site_metadata.loc[
                :, "site_sequence"
            ].tolist(),
        },
        index=phospho.index.copy(),
    )
    site_metadata.loc[:, "localisation_confidence"] = [0.95] * site_metadata.shape[0]
    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
            input_intensity_scale="linear",
        )
    )
    assert list(built.phospho.index) == list(phospho.index)
    assert list(built.site_metadata.columns) == [
        "site_sequence",
        "localisation_confidence",
        "gene_symbol",
        "site",
    ]
    assert (
        built.site_metadata.loc[:, "site_sequence"].tolist()
        == canonical_site_metadata.loc[:, "site_sequence"].tolist()
    )


def test_dataset_builder_preserves_explicit_protein_identity_column() -> None:
    phospho = load_rat_l6_phospho().head(4).copy(deep=True)
    site_metadata = site_metadata_for(phospho).copy(deep=True)
    site_metadata.loc[:, "protein_id"] = [
        f"PROT_{position:03d}" for position in range(site_metadata.shape[0])
    ]
    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
            input_intensity_scale="linear",
        )
    )
    assert "protein_id" in built.site_metadata.columns
    assert (
        built.site_metadata.loc[:, "protein_id"].tolist()
        == site_metadata.loc[:, "protein_id"].tolist()
    )


def test_dataset_builder_supports_row_median_missing_data_preprocessing_policy() -> (
    None
):
    phospho = load_rat_l6_phospho().head(4).copy(deep=True)
    phospho.iloc[0, 0] = float("nan")
    phospho.iloc[1, :] = float("nan")
    original_index = phospho.index.copy()
    expected_imputed = phospho.loc[original_index[0]].median(skipna=True)

    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata_for(phospho),
            organism=Organism.RAT,
            input_intensity_scale="linear",
            preprocessing_config=DatasetPreprocessingConfig(
                missing_data=DatasetMissingDataConfig(
                    policy="impute_row_median",
                    min_observed_values=2,
                ),
            ),
        )
    )

    assert built.phospho.index.tolist() == [
        site_id for site_id in original_index.tolist() if site_id != original_index[1]
    ]
    assert built.phospho.isna().to_numpy().sum() == 0
    assert built.phospho.loc[original_index[0], phospho.columns[0]] == expected_imputed
    assert built.processing_state.missing_data.policy == "impute_row_median"
    assert built.processing_state.missing_data.imputed is True


def test_dataset_builder_runs_site_sequence_resolution_before_minprob_diagnostics(
    tmp_path: Path,
) -> None:
    fasta_path = tmp_path / "proteins.fasta"
    fasta_path.write_text(
        ">P1 protein_1\nAAAASAAAA\n>P2 protein_2\nCCCCCTCCCC\n>P3 protein_3\nAAAASAAAA\n",
        encoding="utf-8",
    )
    phospho = pd.DataFrame(
        {
            "sample_a": [10.0, float("nan"), 8.0],
            "sample_b": [9.0, 7.0, float("nan")],
            "sample_c": [11.0, 6.0, 5.0],
        },
        index=pd.Index(["MAPK14;S5;", "GSK3B;T6;", "AKT1;S5;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "GSK3B", "AKT1"],
            "site": ["S5", "T6", "S5"],
            "protein_accession": ["P1", "P2", "P3"],
            "site_sequence": [pd.NA, pd.NA, pd.NA],
        },
        index=phospho.index.copy(),
    )
    site_metadata.loc[:, "localisation_confidence"] = [0.95] * site_metadata.shape[0]

    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            preprocessing_config=DatasetPreprocessingConfig(
                site_sequence_resolution=DatasetSiteSequenceResolutionConfig(
                    fasta_path=str(fasta_path),
                    flank_size=2,
                ),
                intensity_transform=DatasetIntensityTransformConfig(
                    policy="log2",
                    pseudocount=1.0,
                ),
                missing_data=DatasetMissingDataConfig(
                    policy="impute_minprob",
                    q=0.01,
                    width=0.3,
                    seed=12345,
                    max_missing_fraction_per_row=0.5,
                ),
            ),
        )
    )

    assert built.provenance is not None
    stage_order = [stage.stage for stage in built.provenance.preprocessing_stages]
    assert stage_order == [
        "site_sequence_resolution",
        "localisation_confidence",
        "intensity_transform",
        "missing_data",
    ]

    sequence_stage = built.provenance.preprocessing_stages[0]
    missing_data_stage = next(
        stage
        for stage in built.provenance.preprocessing_stages
        if stage.stage == "missing_data"
    )
    sequence_diagnostics = sequence_stage.diagnostics or {}
    missing_data_diagnostics = missing_data_stage.diagnostics or {}
    assert sequence_diagnostics.get("resolved_site_count") == 3
    assert sequence_diagnostics.get("unresolved_site_count") == 0
    assert missing_data_diagnostics.get("imputation_method_id") == "minprob"


def test_dataset_builder_treats_fasta_resolution_as_authoritative_for_missing_sequences(
    tmp_path: Path,
) -> None:
    fasta_path = tmp_path / "proteins.fasta"
    fasta_path.write_text(
        f">P1\n{'A' * 181}Y{'C' * 24}\n",
        encoding="utf-8",
    )
    phospho = pd.DataFrame(
        {"sample_a": [10.0], "sample_b": [9.0]},
        index=pd.Index(["MAPK14;Y182;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14"],
            "site": ["Y182"],
            "protein_accession": ["P1"],
            "site_sequence": [pd.NA],
        },
        index=phospho.index.copy(),
    )
    site_metadata.loc[:, "localisation_confidence"] = [0.95] * site_metadata.shape[0]

    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
            input_intensity_scale="linear",
            preprocessing_config=DatasetPreprocessingConfig(
                site_sequence_resolution=DatasetSiteSequenceResolutionConfig(
                    fasta_path=str(fasta_path),
                    flank_size=2,
                ),
            ),
        )
    )

    assert built.site_metadata.loc["MAPK14;Y182;", "site_sequence"] == "AAYCC"
    assert built.provenance is not None
    sequence_stage = next(
        stage
        for stage in built.provenance.preprocessing_stages
        if stage.stage == "site_sequence_resolution"
    )
    diagnostics = sequence_stage.diagnostics or {}
    assert diagnostics.get("filled_missing_count") == 1
    assert diagnostics.get("existing_sequence_conflict_count") == 0


def test_dataset_builder_site_sequence_resolution_processing_state_tracks_diagnostics(
    tmp_path: Path,
) -> None:
    fasta_path = tmp_path / "proteins.fasta"
    fasta_path.write_text(
        ">P1 protein_1\nAAAASAAAA\n>P2 protein_2\nCCCCCTCCCC\n",
        encoding="utf-8",
    )
    phospho = pd.DataFrame(
        {
            "sample_a": [10.0, 8.0, 7.0, 6.0],
            "sample_b": [9.0, 7.0, 6.0, 5.0],
        },
        index=pd.Index(
            ["MAPK14;S5;", "GSK3B;T6;", "BAD;T6;", "MAPK1;S5;"],
            name="site_id",
        ),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "GSK3B", "BAD", "MAPK1"],
            "site": ["S5", "T6", "T6", "S5"],
            "protein_accession": ["P1", "P2", "P404", "P1"],
            "site_sequence": [pd.NA, "XXXXX", "KEEP", "AASAA"],
        },
        index=phospho.index.copy(),
    )
    site_metadata.loc[:, "localisation_confidence"] = [0.95] * site_metadata.shape[0]

    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            preprocessing_config=DatasetPreprocessingConfig(
                site_sequence_resolution=DatasetSiteSequenceResolutionConfig(
                    fasta_path=str(fasta_path),
                    mode="replace_existing",
                    flank_size=2,
                ),
            ),
            input_intensity_scale="linear",
        )
    )

    assert built.provenance is not None
    sequence_stages = [
        stage
        for stage in built.provenance.preprocessing_stages
        if stage.stage == "site_sequence_resolution"
    ]
    assert len(sequence_stages) == 1
    sequence_diagnostics = dict(sequence_stages[0].diagnostics or {})
    state = built.processing_state.site_sequence_resolution

    assert state.configured is True
    assert state.mode == "replace_existing"
    assert state.flank_size == 2
    assert state.fasta_source_path == sequence_diagnostics.get("fasta_source_path")
    assert state.fasta_source_label == sequence_diagnostics.get("fasta_source_label")
    assert state.fasta_sha256 == sequence_diagnostics.get("fasta_sha256")
    assert state.resolver_version == sequence_diagnostics.get("resolver_version")
    assert state.resolved_site_count == sequence_diagnostics.get("resolved_site_count")
    assert state.unresolved_site_count == sequence_diagnostics.get(
        "unresolved_site_count"
    )
    assert state.unresolved_counts_by_reason == sequence_diagnostics.get(
        "unresolved_counts_by_reason"
    )
    assert state.filled_missing_count == sequence_diagnostics.get(
        "filled_missing_count"
    )
    assert state.replaced_existing_count == sequence_diagnostics.get(
        "replaced_existing_count"
    )
    assert state.preserved_existing_count == sequence_diagnostics.get(
        "preserved_existing_count"
    )
    assert state.existing_sequence_conflict_count == sequence_diagnostics.get(
        "existing_sequence_conflict_count"
    )


def test_dataset_builder_supports_site_matrix_build_from_metadata_policy() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 2.0, 3.0],
            "sample_b": [1.5, 2.5, float("nan")],
        },
        index=pd.Index(["row_a", "row_b", "row_c"], name="input_row"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "MAPK14", "AKT1"],
            "site": ["Y182", "Y182", "T308"],
            "site_sequence": ["SEQ_A", "SEQ_R", "SEQ_C"],
        },
        index=phospho.index.copy(),
    )
    site_metadata.loc[:, "localisation_confidence"] = [0.95] * site_metadata.shape[0]

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
                ),
                site_matrix=DatasetSiteMatrixConfig(policy="build_from_metadata"),
            ),
        )
    )

    assert built.phospho.index.tolist() == ["MAPK14;Y182;"]
    assert built.phospho.iloc[0, 0] == pytest.approx(2.0)
    assert built.phospho.iloc[0, 1] == pytest.approx(2.5)
    assert built.site_metadata.index.tolist() == ["MAPK14;Y182;"]
    assert built.site_metadata.loc["MAPK14;Y182;", "site_sequence"] == "SEQ_R"
    assert built.processing_state.site_matrix.policy == "build_from_metadata"
    assert built.processing_state.site_matrix.constructed is True


def test_dataset_builder_supports_site_matrix_duplicate_aggregation_policy() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 3.0],
            "sample_b": [2.0, 4.0],
        },
        index=pd.Index(["row_a", "row_b"], name="input_row"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "MAPK14"],
            "site": ["Y182", "Y182"],
            "site_sequence": ["SEQ_A", "SEQ_A"],
            "source_uid": ["UID_A", "UID_B"],
        },
        index=phospho.index.copy(),
    )
    site_metadata.loc[:, "localisation_confidence"] = [0.95] * site_metadata.shape[0]

    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
            input_intensity_scale="linear",
            preprocessing_config=DatasetPreprocessingConfig(
                site_matrix=DatasetSiteMatrixConfig(
                    policy="build_from_metadata",
                    duplicate_site_policy="aggregate_mean",
                )
            ),
        )
    )

    assert built.phospho.index.tolist() == ["MAPK14;Y182;"]
    assert built.phospho.loc["MAPK14;Y182;", "sample_a"] == pytest.approx(2.0)
    assert built.phospho.loc["MAPK14;Y182;", "sample_b"] == pytest.approx(3.0)
    assert pd.isna(built.site_metadata.loc["MAPK14;Y182;", "source_uid"])
    assert built.preprocessing_report is not None
    assert built.preprocessing_report.row_audit is not None
    assert built.preprocessing_report.duplicate_site_resolution is not None
    assert built.preprocessing_report.metadata_conflicts is not None
    duplicate_rows = built.preprocessing_report.duplicate_site_resolution
    assert duplicate_rows.shape[0] == 2
    assert duplicate_rows["source_row_id"].tolist() == ["row_a", "row_b"]
    assert duplicate_rows["retained"].tolist() == [True, True]
    assert duplicate_rows["n_aggregated_rows"].tolist() == [2, 2]
    conflicts = built.preprocessing_report.metadata_conflicts
    assert not conflicts.empty
    assert "source_uid" in set(conflicts.loc[:, "field"])
    aggregated = built.preprocessing_report.row_audit.loc[
        (built.preprocessing_report.row_audit.loc[:, "stage"] == "site_matrix")
        & (built.preprocessing_report.row_audit.loc[:, "action"] == "aggregated")
    ]
    assert set(aggregated.loc[:, "source_row_id"].astype(str)) == {"row_a", "row_b"}


def test_dataset_builder_site_matrix_derivation_keeps_all_fully_resolvable_rows() -> (
    None
):
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 2.0],
            "sample_b": [1.5, 2.5],
        },
        index=pd.Index(["MAPK14;Y182;", "GSK3B;S9;"], name="input_row"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "GSK3B"],
            "site": ["Y182", "S9"],
        },
        index=phospho.index.copy(),
    )
    site_metadata.loc[:, "localisation_confidence"] = [0.95] * site_metadata.shape[0]

    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
            input_intensity_scale="linear",
            preprocessing_config=DatasetPreprocessingConfig(
                site_matrix=DatasetSiteMatrixConfig(policy="build_from_metadata")
            ),
        )
    )

    assert built.phospho.index.tolist() == ["GSK3B;S9;", "MAPK14;Y182;"]
    assert built.site_metadata.index.tolist() == ["GSK3B;S9;", "MAPK14;Y182;"]
    assert built.site_metadata.loc[:, "site_sequence"].isna().sum() == 0
    assert built.preprocessing_report is not None
    dropped = built.preprocessing_report.row_audit.loc[
        (built.preprocessing_report.row_audit.loc[:, "stage"] == "site_matrix")
        & (built.preprocessing_report.row_audit.loc[:, "action"] == "dropped")
    ]
    assert dropped.empty


def test_dataset_builder_site_matrix_derivation_excludes_only_unresolved_rows() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 2.0, 3.0],
            "sample_b": [1.5, 2.5, 3.5],
        },
        index=pd.Index(["MAPK14;Y182;", "FAKE1;S1;", "GSK3B;S9;"], name="input_row"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "FAKE1", "GSK3B"],
            "site": ["Y182", "S1", "S9"],
        },
        index=phospho.index.copy(),
    )
    site_metadata.loc[:, "localisation_confidence"] = [0.95] * site_metadata.shape[0]

    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
            input_intensity_scale="linear",
            preprocessing_config=DatasetPreprocessingConfig(
                site_matrix=DatasetSiteMatrixConfig(policy="build_from_metadata")
            ),
        )
    )

    assert built.phospho.index.tolist() == ["GSK3B;S9;", "MAPK14;Y182;"]
    assert built.site_metadata.index.tolist() == ["GSK3B;S9;", "MAPK14;Y182;"]
    assert built.preprocessing_report is not None
    dropped = built.preprocessing_report.row_audit.loc[
        (built.preprocessing_report.row_audit.loc[:, "stage"] == "site_matrix")
        & (built.preprocessing_report.row_audit.loc[:, "action"] == "dropped")
    ]
    assert set(dropped.loc[:, "source_row_id"].astype(str)) == {"FAKE1;S1;"}
    assert "missing or blank" in str(dropped.iloc[0]["reason"])


def test_dataset_builder_site_matrix_derivation_uses_row_metadata_site_identity() -> (
    None
):
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 2.0, 3.0],
            "sample_b": [1.5, 2.5, 3.5],
        },
        index=pd.Index(["row_a", "row_b", "row_c"], name="input_row"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "FAKE1", "GSK3B"],
            "site": ["Y182", "S1", "S9"],
        },
        index=phospho.index.copy(),
    )
    site_metadata.loc[:, "localisation_confidence"] = [0.95] * site_metadata.shape[0]

    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
            input_intensity_scale="linear",
            preprocessing_config=DatasetPreprocessingConfig(
                site_matrix=DatasetSiteMatrixConfig(policy="build_from_metadata")
            ),
        )
    )

    assert built.phospho.index.tolist() == ["GSK3B;S9;", "MAPK14;Y182;"]
    assert built.site_metadata.index.tolist() == ["GSK3B;S9;", "MAPK14;Y182;"]
    assert built.site_metadata.loc[:, "site_sequence"].isna().sum() == 0
    assert built.preprocessing_report is not None
    dropped = built.preprocessing_report.row_audit.loc[
        (built.preprocessing_report.row_audit.loc[:, "stage"] == "site_matrix")
        & (built.preprocessing_report.row_audit.loc[:, "action"] == "dropped")
    ]
    assert set(dropped.loc[:, "source_row_id"].astype(str)) == {"row_b"}
    assert "missing or blank" in str(dropped.iloc[0]["reason"])


def test_dataset_builder_site_matrix_excludes_unusable_supplied_sequence_rows() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 2.0],
            "sample_b": [1.5, 2.5],
        },
        index=pd.Index(["MAPK14;Y182;", "FAKE1;S1;"], name="input_row"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "FAKE1"],
            "site": ["Y182", "S1"],
            "site_sequence": [("A" * 15) + "Y" + ("A" * 15), "  "],
        },
        index=phospho.index.copy(),
    )
    site_metadata.loc[:, "localisation_confidence"] = [0.95] * site_metadata.shape[0]

    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            input_intensity_scale="linear",
            preprocessing_config=DatasetPreprocessingConfig(
                site_matrix=DatasetSiteMatrixConfig(policy="build_from_metadata")
            ),
        )
    )

    assert built.phospho.index.tolist() == ["MAPK14;Y182;"]
    assert built.site_metadata.index.tolist() == ["MAPK14;Y182;"]
    assert built.preprocessing_report is not None
    dropped = built.preprocessing_report.row_audit.loc[
        (built.preprocessing_report.row_audit.loc[:, "stage"] == "site_matrix")
        & (built.preprocessing_report.row_audit.loc[:, "action"] == "dropped")
    ]
    assert set(dropped.loc[:, "source_row_id"].astype(str)) == {"FAKE1;S1;"}
    assert "missing or blank" in str(dropped.iloc[0]["reason"])


def test_dataset_builder_site_matrix_derivation_reports_no_rows_when_fully_unresolvable() -> (
    None
):
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 2.0],
            "sample_b": [1.5, 2.5],
        },
        index=pd.Index(["FAKE1;S1;", "FAKE2;T2;"], name="input_row"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["FAKE1", "FAKE2"],
            "site": ["S1", "T2"],
        },
        index=phospho.index.copy(),
    )
    site_metadata.loc[:, "localisation_confidence"] = [0.95] * site_metadata.shape[0]

    with pytest.raises(
        PhosPyInputError,
        match=(
            "site-matrix construction produced no retained rows after filtering; "
            "input_rows=2, dropped_missing_sequence=2"
        ),
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=site_metadata,
                organism=Organism.RAT,
                preprocessing_config=DatasetPreprocessingConfig(
                    site_matrix=DatasetSiteMatrixConfig(policy="build_from_metadata")
                ),
            )
        )


@pytest.mark.parametrize(
    "missing_data_policy",
    ("retain_missing", "require_min_observed_values"),
)
def test_dataset_builder_rejects_incompatible_site_matrix_missing_data_modes_early(
    missing_data_policy: str,
) -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0],
            "sample_b": [2.0],
        },
        index=pd.Index(["row_a"], name="input_row"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14"],
            "site": ["Y182"],
            "site_sequence": ["SEQ_A"],
        },
        index=phospho.index.copy(),
    )
    site_metadata.loc[:, "localisation_confidence"] = [0.95] * site_metadata.shape[0]

    with pytest.raises(
        PhosPyInputError,
        match=(
            f"missing_data_policy='{missing_data_policy}' is not supported for strict "
            "AnalysisReadyPhosphoDataset construction"
        ),
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=site_metadata,
                organism=Organism.RAT,
                preprocessing_config=DatasetPreprocessingConfig(
                    site_matrix=DatasetSiteMatrixConfig(
                        policy="build_from_metadata",
                        missing_data_policy=missing_data_policy,  # type: ignore[arg-type]
                    )
                ),
            )
        )


@pytest.mark.parametrize(
    "missing_data_policy",
    ("retain_missing", "require_min_observed_values"),
)
def test_dataset_builder_rejects_dead_end_site_matrix_missing_modes_before_dataset_boundary(
    missing_data_policy: str,
) -> None:
    phospho = pd.DataFrame(
        {"sample_a": [1.0], "sample_b": [float("nan")]},
        index=pd.Index(["row_a"], name="input_row"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14"],
            "site": ["Y182"],
            "site_sequence": ["SEQ_A"],
        },
        index=phospho.index.copy(),
    )
    site_metadata.loc[:, "localisation_confidence"] = [0.95] * site_metadata.shape[0]

    with pytest.raises(
        PhosPyInputError,
        match=(
            f"missing_data_policy='{missing_data_policy}' is not supported for strict "
            "AnalysisReadyPhosphoDataset construction"
        ),
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=site_metadata,
                organism=Organism.RAT,
                preprocessing_config=DatasetPreprocessingConfig(
                    site_matrix=DatasetSiteMatrixConfig(
                        policy="build_from_metadata",
                        missing_data_policy=missing_data_policy,  # type: ignore[arg-type]
                    )
                ),
            )
        )


def test_dataset_builder_builds_inferred_comparisons_from_sample_metadata() -> None:
    phospho = pd.DataFrame(
        {
            "sample_1": [8.0, 2.0],
            "sample_2": [8.0, 4.0],
            "sample_3": [5.0, 1.0],
            "sample_4": [5.0, 1.0],
        },
        index=pd.Index(["MAPK14;Y182;", "AKT1;T308;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "AKT1"],
            "site": ["Y182", "T308"],
            "site_sequence": ["SEQ_A", "SEQ_R"],
        },
        index=phospho.index.copy(),
    )
    site_metadata.loc[:, "localisation_confidence"] = [0.95] * site_metadata.shape[0]
    sample_metadata = pd.DataFrame(
        {"comparison_group": ["group1", "group1", "group4", "group4"]},
        index=phospho.columns.copy(),
    )

    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            sample_metadata=sample_metadata,
            organism=Organism.RAT,
            input_intensity_scale="linear",
            preprocessing_config=DatasetPreprocessingConfig(
                comparisons=DatasetComparisonBuildingConfig(
                    policy="sample_metadata_pairs"
                )
            ),
        )
    )

    assert built.comparisons is not None
    expected = pd.DataFrame(
        {"p_group1_group4": [3.0, 2.0]},
        index=phospho.index.copy(),
    )
    pdt.assert_frame_equal(built.comparisons, expected)
    assert built.preprocessing_report is not None
    assert built.preprocessing_report.row_audit is not None
    assert built.preprocessing_report.comparison_group_stats is not None
    assert built.preprocessing_report.comparison_pair_stats is not None
    group_stats = built.preprocessing_report.comparison_group_stats
    pair_stats = built.preprocessing_report.comparison_pair_stats
    assert not group_stats.empty
    assert not pair_stats.empty
    assert {"site_id", "group", "n", "mean", "sd", "sem"}.issubset(
        set(group_stats.columns)
    )
    assert {
        "site_id",
        "comparison",
        "left_n",
        "right_n",
        "left_mean",
        "right_mean",
        "left_sd",
        "right_sd",
        "left_sem",
        "right_sem",
        "effect_size",
    }.issubset(set(pair_stats.columns))
    comparison_long = built.comparisons.reset_index().melt(
        id_vars=["site_id"],
        var_name="comparison",
        value_name="expected_effect_size",
    )
    merged = pair_stats.merge(
        comparison_long,
        how="inner",
        on=["site_id", "comparison"],
    )
    assert merged.shape[0] == built.comparisons.shape[0]
    assert (merged.loc[:, "effect_size"] == merged.loc[:, "expected_effect_size"]).all()
    assert built.processing_state.comparisons.policy == "sample_metadata_pairs"
    assert built.processing_state.comparisons.sample_group_column == "comparison_group"
    assert built.processing_state.comparisons.pairs is None


def test_dataset_builder_rejects_comparison_groups_missing_from_metadata() -> None:
    phospho = pd.DataFrame(
        {"sample_a": [7.0], "sample_b": [4.0]},
        index=pd.Index(["PRKACA;S339;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["PRKACA"],
            "site": ["S339"],
            "site_sequence": ["AAAAAA"],
        },
        index=phospho.index.copy(),
    )
    site_metadata.loc[:, "localisation_confidence"] = [0.95] * site_metadata.shape[0]
    sample_metadata = pd.DataFrame(
        {"comparison_group": ["sample_a", "sample_b"]},
        index=phospho.columns.copy(),
    )

    with pytest.raises(
        PhosPyInputError,
        match="references unknown sample groups",
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=site_metadata,
                sample_metadata=sample_metadata,
                organism=Organism.RAT,
                preprocessing_config=DatasetPreprocessingConfig(
                    comparisons=DatasetComparisonBuildingConfig(
                        policy="sample_metadata_pairs",
                        pairs=(("sample_a", "missing_group"),),
                    )
                ),
            )
        )


def test_dataset_builder_rejects_site_matrix_build_without_site_sequence_column() -> (
    None
):
    phospho = load_rat_l6_phospho().head(4).copy(deep=True)
    site_metadata = site_metadata_for(phospho).drop(columns=["site_sequence"])

    with pytest.raises(
        PhosPyInputError,
        match=(
            "site-matrix construction produced no retained rows after filtering; "
            "input_rows=4, dropped_missing_sequence=4"
        ),
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=site_metadata,
                preprocessing_config=DatasetPreprocessingConfig(
                    site_matrix=DatasetSiteMatrixConfig(policy="build_from_metadata")
                ),
            )
        )


def test_dataset_builder_default_forbid_policy_keeps_missingness_strict() -> None:
    phospho = load_rat_l6_phospho().head(4).copy(deep=True)
    phospho.iloc[0, 0] = float("nan")

    with pytest.raises(
        PhosPyInputError,
        match="missing_data.policy='forbid'",
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=site_metadata_for(phospho),
                organism=Organism.RAT,
            )
        )


def test_dataset_builder_log2_preprocessing_records_operation_and_state() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 3.0],
            "sample_b": [2.0, 7.0],
        },
        index=pd.Index(["MAPK14;Y182;", "AKT1;T308;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "AKT1"],
            "site": ["Y182", "T308"],
            "site_sequence": ["SEQ_A", "SEQ_R"],
        },
        index=phospho.index.copy(),
    )
    site_metadata.loc[:, "localisation_confidence"] = [0.95] * site_metadata.shape[0]

    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            preprocessing_config=DatasetPreprocessingConfig(
                intensity_transform=DatasetIntensityTransformConfig(
                    policy="log2",
                    pseudocount=1.0,
                )
            ),
        )
    )

    assert built.intensity_scale_state.label == "log2"
    assert built.intensity_scale_state.quantity.value == "phosphosite_log_abundance"
    assert built.preprocessing_report is not None
    operations = built.preprocessing_report.operations
    log2_operation = operations.loc[
        (operations.loc[:, "stage"] == "intensity_transform")
        & (operations.loc[:, "notes"] == "stage executed")
    ]
    assert log2_operation.shape[0] == 1
    assert log2_operation.iloc[0]["operation"] == "log2"
    assert log2_operation.iloc[0]["parameters"] == {"pseudocount": 1.0}
    assert built.processing_state.intensity_scale == built.intensity_scale_state


def test_dataset_builder_distinguishes_corrected_vs_uncorrected_log2_quantity() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [15.0],
            "sample_b": [31.0],
        },
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
    site_metadata.loc[:, "localisation_confidence"] = [0.95] * site_metadata.shape[0]
    total = pd.DataFrame(
        {
            "sample_a": [3.0],
            "sample_b": [7.0],
        },
        index=pd.Index(["MAPK14"], name="protein_id"),
    )

    uncorrected = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            preprocessing_config=DatasetPreprocessingConfig(
                intensity_transform=DatasetIntensityTransformConfig(
                    policy="log2",
                    pseudocount=1.0,
                )
            ),
        )
    )
    corrected = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            total=total,
            preprocessing_config=DatasetPreprocessingConfig(
                intensity_transform=DatasetIntensityTransformConfig(
                    policy="log2",
                    pseudocount=1.0,
                ),
                total_protein_correction=DatasetTotalProteinCorrectionConfig(
                    policy="subtract_log_total"
                ),
            ),
        )
    )

    assert uncorrected.intensity_scale_state.label == "log2"
    assert corrected.intensity_scale_state.label == "log2"
    assert (
        uncorrected.intensity_scale_state.quantity.value == "phosphosite_log_abundance"
    )
    assert corrected.intensity_scale_state.quantity.value == "phospho_total_log_ratio"
    assert (
        uncorrected.intensity_scale_state.quantity
        != corrected.intensity_scale_state.quantity
    )
    assert (
        uncorrected.processing_state.total_protein_correction.quantitative_meaning
        == "phosphosite_log_abundance"
    )
    assert corrected.processing_state.total_protein_correction.quantitative_meaning == (
        "phospho_total_log_ratio"
    )


def test_dataset_builder_public_payloads_pair_scale_with_quantitative_meaning(
    tmp_path: Path,
) -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [15.0],
            "sample_b": [31.0],
        },
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
    site_metadata.loc[:, "localisation_confidence"] = [0.95] * site_metadata.shape[0]
    total = pd.DataFrame(
        {
            "sample_a": [3.0],
            "sample_b": [7.0],
        },
        index=pd.Index(["MAPK14"], name="protein_id"),
    )
    datasets = (
        (
            "uncorrected",
            AnalysisReadyDatasetBuilder().run(
                DatasetBuildRequest(
                    phospho=phospho,
                    site_metadata=site_metadata,
                    preprocessing_config=DatasetPreprocessingConfig(
                        intensity_transform=DatasetIntensityTransformConfig(
                            policy="log2",
                            pseudocount=1.0,
                        )
                    ),
                )
            ),
            "phosphosite_log_abundance",
        ),
        (
            "corrected",
            AnalysisReadyDatasetBuilder().run(
                DatasetBuildRequest(
                    phospho=phospho,
                    site_metadata=site_metadata,
                    total=total,
                    preprocessing_config=DatasetPreprocessingConfig(
                        intensity_transform=DatasetIntensityTransformConfig(
                            policy="log2",
                            pseudocount=1.0,
                        ),
                        total_protein_correction=DatasetTotalProteinCorrectionConfig(
                            policy="subtract_log_total"
                        ),
                    ),
                )
            ),
            "phospho_total_log_ratio",
        ),
    )

    for label, dataset, expected_quantitative_meaning in datasets:
        assert dataset.intensity_scale_state.label == "log2"
        assert (
            dataset.intensity_scale_state.quantity.value
            == expected_quantitative_meaning
        )
        written = publish_dataset(
            dataset,
            tmp_path / f"published_{label}",
            output_format="csv",
        )
        manifest = json.loads(written["dataset.manifest"].read_text(encoding="utf-8"))
        assert manifest["intensity_scale"] == "log2"
        assert manifest["quantitative_meaning"] == expected_quantitative_meaning
        processing_state_payload = manifest["processing_state"]
        assert (
            processing_state_payload["intensity_scale"]["quantity"]
            == expected_quantitative_meaning
        )
        correction_payload = processing_state_payload["total_protein_correction"]
        assert (
            correction_payload["quantitative_meaning"] == expected_quantitative_meaning
        )
        correction_diagnostics = correction_payload["diagnostics"]
        assert correction_diagnostics["diagnostics_schema_version"] == 1
        assert (
            correction_diagnostics["quantitative_meaning"]
            == expected_quantitative_meaning
        )
        assert manifest["provenance"]["workflow_parameters"][
            "intensity_scale_label"
        ] == ("log2")
        assert manifest["provenance"]["workflow_parameters"][
            "quantitative_meaning"
        ] == (expected_quantitative_meaning)
        correction_stage = next(
            (
                stage
                for stage in manifest["provenance"]["preprocessing_stages"]
                if stage["stage"] == "total_protein_correction"
            ),
            None,
        )
        if correction_stage is not None:
            diagnostics = correction_stage["diagnostics"] or {}
            if "output_scale" in diagnostics:
                assert diagnostics["quantitative_meaning"] == (
                    expected_quantitative_meaning
                )


def test_dataset_builder_median_center_preprocessing_records_operation() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 3.0, 5.0],
            "sample_b": [2.0, 3.0, 4.0],
        },
        index=pd.Index(["MAPK14;Y182;", "AKT1;T308;", "GSK3B;S9;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "AKT1", "GSK3B"],
            "site": ["Y182", "T308", "S9"],
            "site_sequence": ["SEQ_A", "SEQ_R", "SEQ_C"],
        },
        index=phospho.index.copy(),
    )
    site_metadata.loc[:, "localisation_confidence"] = [0.95] * site_metadata.shape[0]

    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            input_intensity_scale="linear",
            preprocessing_config=DatasetPreprocessingConfig(
                normalisation=DatasetNormalisationConfig(policy="median_center")
            ),
        )
    )

    assert built.preprocessing_report is not None
    operations = built.preprocessing_report.operations
    normalisation_operation = operations.loc[
        (operations.loc[:, "stage"] == "normalisation")
        & (operations.loc[:, "notes"] == "stage executed")
    ]
    assert normalisation_operation.shape[0] == 1
    assert normalisation_operation.iloc[0]["operation"] == "median_center"
    assert normalisation_operation.iloc[0]["parameters"] == {
        "applied": True,
        "centering_statistic": "median",
        "axis": "columns",
        "skipna": True,
    }
    assert built.processing_state.normalisation.policy == "median_center"
    assert built.provenance is not None
    normalisation_stage = next(
        stage
        for stage in built.provenance.preprocessing_stages
        if stage.stage == "normalisation"
    )
    assert normalisation_stage.operation == "median_center"
    assert dict(normalisation_stage.parameters) == {
        "applied": True,
        "centering_statistic": "median",
        "axis": "columns",
        "skipna": True,
    }
    diagnostics = normalisation_stage.diagnostics or {}
    assert diagnostics["method"] == "median_center"
    assert diagnostics["parameters"] == {
        "applied": True,
        "centering_statistic": "median",
        "axis": "columns",
        "skipna": True,
    }


def test_dataset_builder_quantile_preprocessing_records_operation() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [4.0, 1.0, 3.0],
            "sample_b": [5.0, 2.0, 1.0],
        },
        index=pd.Index(["MAPK14;Y182;", "AKT1;T308;", "GSK3B;S9;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "AKT1", "GSK3B"],
            "site": ["Y182", "T308", "S9"],
            "site_sequence": ["SEQ_A", "SEQ_R", "SEQ_C"],
        },
        index=phospho.index.copy(),
    )
    site_metadata.loc[:, "localisation_confidence"] = [0.95] * site_metadata.shape[0]

    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            input_intensity_scale="linear",
            preprocessing_config=DatasetPreprocessingConfig(
                normalisation=DatasetNormalisationConfig(policy="quantile")
            ),
        )
    )

    assert built.preprocessing_report is not None
    operations = built.preprocessing_report.operations
    normalisation_operation = operations.loc[
        (operations.loc[:, "stage"] == "normalisation")
        & (operations.loc[:, "notes"] == "stage executed")
    ]
    assert normalisation_operation.shape[0] == 1
    assert normalisation_operation.iloc[0]["operation"] == "quantile"
    assert normalisation_operation.iloc[0]["parameters"] == {
        "applied": True,
        "target_distribution": "mean_rank_distribution",
        "tie_strategy": "deterministic_rank_average",
        "dtype": "float64",
    }
    assert built.processing_state.normalisation.policy == "quantile"


def test_reference_bundle_rat_tables_are_structurally_coherent() -> None:
    bundle = ReferenceResolver().run(
        ReferencePreset.RAT,
        dataset_organism=Organism.RAT,
    )
    substrate_sites = set(
        bundle.kinase_substrate_map.loc[:, "substrate_site"].astype(str)
    )
    known_sites = set(bundle.site_sequences.index.astype(str))
    assert substrate_sites.issubset(known_sites)


def test_dataset_builder_emits_machine_readable_run_provenance() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, float("nan"), 2.0, float("nan")],
            "sample_b": [2.0, 2.0, 3.0, float("nan")],
            "sample_c": [3.0, 3.0, 4.0, float("nan")],
        },
        index=pd.Index(["row_a", "row_b", "row_c", "row_d"], name="input_row"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "MAPK14", "AKT1", "GSK3B"],
            "site": ["Y182", "Y182", "T308", "S9"],
            "site_sequence": ["SEQ_A", "SEQ_R", "", "SEQ_D"],
        },
        index=phospho.index.copy(),
    )
    site_metadata.loc[:, "localisation_confidence"] = [0.95] * site_metadata.shape[0]

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
                ),
                site_matrix=DatasetSiteMatrixConfig(policy="build_from_metadata"),
            ),
        )
    )

    assert built.preprocessing_report is not None
    assert built.preprocessing_report.row_audit is not None
    provenance = built.provenance
    assert provenance is not None
    assert provenance.environment.package_name == "phospy"
    assert provenance.workflow_name == "dataset_builder"
    assert provenance.reference is None
    assert provenance.random_state is None
    input_names = {item.name for item in provenance.input_tables}
    output_names = {item.name for item in provenance.output_tables}
    assert "dataset.phospho" in input_names
    assert "dataset.site_metadata" in input_names
    assert "dataset.phospho" in output_names
    assert "dataset.site_metadata" in output_names

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
    assert missing_stage.determinism == "pure"
    assert missing_stage.is_deterministic is True
    assert missing_stage.random_seed is None
    assert isinstance(missing_stage.phospho_input_hash, str)
    assert isinstance(missing_stage.phospho_output_hash, str)
    assert missing_stage.input_hash != missing_stage.phospho_input_hash
    assert missing_stage.output_hash != missing_stage.phospho_output_hash
    assert set(missing_stage.dropped_row_ids) == {"row_d"}
    assert missing_stage.imputed_cell_count == 1
    assert set(missing_stage.imputed_row_ids) == {"row_b"}

    site_matrix_stage = next(
        stage
        for stage in provenance.preprocessing_stages
        if stage.stage == "site_matrix"
    )
    assert site_matrix_stage.operation == "build_from_metadata"
    diagnostics = site_matrix_stage.diagnostics or {}
    assert "row_c" in set(diagnostics["dropped_missing_sequence_row_ids"])
    assert diagnostics["duplicate_site_policy"] == "max_mean_signal"
    assert diagnostics["missing_data_policy"] == "drop_any_missing"
    assert "MAPK14;Y182;" in set(diagnostics["final_constructed_site_ids"])


def test_dataset_builder_marks_minprob_stage_as_seeded_stochastic() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [10.0, float("nan"), 4.0],
            "sample_b": [9.0, 8.0, float("nan")],
            "sample_c": [11.0, 7.0, 5.0],
        },
        index=pd.Index(
            ["MAPK14;Y182;", "AKT1;T308;", "PRKACA;S339;"],
            name="site_id",
        ),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "AKT1", "PRKACA"],
            "site": ["Y182", "T308", "S339"],
            "site_sequence": ["SEQ_A", "SEQ_R", "SEQ_C"],
        },
        index=phospho.index.copy(),
    )
    site_metadata.loc[:, "localisation_confidence"] = [0.95] * site_metadata.shape[0]
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
                    max_missing_fraction_per_row=0.75,
                ),
            ),
        )
    )
    assert built.provenance is not None
    missing_stage = next(
        stage
        for stage in built.provenance.preprocessing_stages
        if stage.stage == "missing_data"
    )
    assert missing_stage.determinism == "seeded_stochastic"
    assert missing_stage.is_deterministic is False
    assert missing_stage.random_seed == 123
