from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from phospy.advanced import (
    DatasetIntensityTransformConfig,
    DatasetSiteMatrixConfig,
    DifferentialAnalysisConfig,
    EmpiricalBayesConfig,
)
from phospy.api import (
    AnalysisReadyDatasetBuilder,
    Contrast,
    DatasetPreprocessingConfig,
    DifferentialAnalysisRequest,
    DifferentialAnalysisWorkflow,
    ExperimentalDesign,
    Organism,
    PhosphositeImportResult,
    SampleDesignRecord,
)
from phospy.api.requests import (
    DATASET_MULTI_SITE_POLICY_REJECT,
    DATASET_MULTI_SITE_POLICY_SPLIT,
    DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE,
)
from phospy.contracts.results import (
    IMPORTER_QUALITY_STATUS_NOT_APPLICABLE,
    IMPORTER_QUALITY_STATUS_REPORTED,
)
from phospy.errors import PhosPyInputError
from phospy.io.readers import (
    MaxQuantColumnMapping,
    MaxQuantPhosphositeImporter,
    MaxQuantPhosphositeImportRequest,
)
from phospy.science.differential.models import (
    EMPIRICAL_BAYES_TREND_COVARIATE_QUANTIFICATION_DEPTH,
    QUANTIFICATION_DEPTH_KIND_PEPTIDE_COUNT,
    QUANTIFICATION_DEPTH_KIND_PSM_COUNT,
)

ROOT = Path(__file__).resolve().parents[2]
FIXTURES = ROOT / "tests" / "fixtures" / "maxquant"


def _depth_source(
    depth: list[object],
    *,
    proteins: list[str] | None = None,
    genes: list[str] | None = None,
    residues: list[str] | None = None,
    positions: list[object] | None = None,
) -> pd.DataFrame:
    genes = genes or ["MAPK1", "AKT1", "GSK3B", "MTOR", "RPS6KB1"][: len(depth)]
    residues = residues or ["S", "S", "S", "S", "T"][: len(depth)]
    positions = positions or [10, 473, 9, 2448, 389][: len(depth)]
    proteins = (
        proteins or ["P28482", "P31749", "P49841", "P42345", "P23443"][: len(depth)]
    )
    return pd.DataFrame(
        {
            "Proteins": proteins,
            "Gene names": genes,
            "Amino acid": residues,
            "Positions within proteins": [str(position) for position in positions],
            "Localization prob": ["0.95"] * len(depth),
            "Sequence": ["AAAAA" + residue + "AAAA" for residue in residues],
            "Modified sequence": [
                "AAAAA(ph)" + residue + "AAAA" for residue in residues
            ],
            "Sequence window": [
                ("A" * 15) + residue + ("A" * 15) for residue in residues
            ],
            "Intensity A_1": [10.0 + row for row in range(len(depth))],
            "Intensity A_2": [11.0 + row for row in range(len(depth))],
            "Intensity B_1": [20.0 + row for row in range(len(depth))],
            "Intensity B_2": [22.0 + row for row in range(len(depth))],
            "Depth": depth,
            "Potential contaminant": [""] * len(depth),
            "Reverse": [""] * len(depth),
        }
    )


def _depth_mapping(
    *,
    kind: str = QUANTIFICATION_DEPTH_KIND_PSM_COUNT,
) -> MaxQuantColumnMapping:
    return MaxQuantColumnMapping(
        quantification_depth="Depth",
        quantification_depth_kind=kind,  # type: ignore[arg-type]
    )


def _build_dataset(
    import_result: PhosphositeImportResult,
    *,
    site_resolution_mode: str,
    multi_site_policy: str | None = None,
    site_matrix_duplicate_policy: str | None = None,
):
    preprocessing_config = DatasetPreprocessingConfig(
        intensity_transform=DatasetIntensityTransformConfig(
            policy="log2",
            pseudocount=1.0,
        )
    )
    if site_matrix_duplicate_policy is not None:
        preprocessing_config = DatasetPreprocessingConfig(
            intensity_transform=DatasetIntensityTransformConfig(
                policy="log2",
                pseudocount=1.0,
            ),
            site_matrix=DatasetSiteMatrixConfig(
                policy="build_from_metadata",
                duplicate_site_policy=site_matrix_duplicate_policy,
            ),
        )
    return AnalysisReadyDatasetBuilder().run(
        import_result.to_dataset_build_request(
            site_resolution_mode=site_resolution_mode,
            multi_site_policy=(
                multi_site_policy
                if multi_site_policy is not None
                else (
                    DATASET_MULTI_SITE_POLICY_REJECT
                    if site_resolution_mode
                    == DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE
                    else None
                )
            ),
            organism=Organism.HUMAN,
            input_intensity_scale="linear",
            preprocessing_config=preprocessing_config,
        )
    )


def _depth_design() -> ExperimentalDesign:
    return ExperimentalDesign(
        samples=(
            SampleDesignRecord(
                sample_id="A_1",
                condition="A",
                biological_replicate_id="A_r1",
            ),
            SampleDesignRecord(
                sample_id="A_2",
                condition="A",
                biological_replicate_id="A_r2",
            ),
            SampleDesignRecord(
                sample_id="B_1",
                condition="B",
                biological_replicate_id="B_r1",
            ),
            SampleDesignRecord(
                sample_id="B_2",
                condition="B",
                biological_replicate_id="B_r2",
            ),
        )
    )


def test_maxquant_importer_reads_standard_phospho_sty_sites_columns() -> None:
    result = MaxQuantPhosphositeImporter().run(
        MaxQuantPhosphositeImportRequest(
            source=FIXTURES / "phospho_sty_sites_standard.txt"
        )
    )

    assert isinstance(result, PhosphositeImportResult)
    phospho = result.phospho_matrix_candidate
    assert list(phospho.columns) == ["Control", "Stim"]
    assert phospho.shape == (2, 2)
    assert float(phospho.iloc[0]["Control"]) == pytest.approx(10.0)
    metadata = result.site_metadata_candidate
    assert metadata.loc[:, "gene_symbol"].tolist() == ["MAPK1", "AKT1"]
    assert metadata.loc[:, "site"].tolist() == ["S10", "S473"]
    assert metadata.loc[:, "protein_accession"].tolist() == ["P28482", "P31749"]
    assert metadata.loc[:, "protein_id"].tolist() == ["P28482", "P31749"]
    assert metadata.index.astype(str).str.startswith("maxquant:P").all()
    assert "MAPK1;S10;" not in metadata.index.astype(str).tolist()
    assert result.localisation_confidence_column == "localisation_confidence"
    assert metadata.loc[:, "localisation_confidence"].tolist() == pytest.approx(
        [0.95, 0.91]
    )
    assert result.diagnostics["maxquant"]["filtering"]["removed_rows"] == 2


def test_maxquant_realistic_grouped_multisite_rows_remain_candidates() -> None:
    result = MaxQuantPhosphositeImporter().run(
        MaxQuantPhosphositeImportRequest(
            source=FIXTURES / "phospho_sty_sites_realistic_variants.txt"
        )
    )

    phospho = result.phospho_matrix_candidate
    assert list(phospho.columns) == ["Control", "Stim"]
    assert phospho.shape == (2, 2)
    assert float(phospho.iloc[0]["Control"]) == pytest.approx(1000.0)
    metadata = result.site_metadata_candidate
    assert metadata.loc[:, "gene_symbol"].tolist() == ["MAPK1", "AKT1"]
    assert metadata.loc[:, "site"].tolist() == ["S10", "S473,T308"]
    assert metadata.loc[:, "protein_accession"].tolist() == ["P28482", "P31749"]
    assert metadata.index.astype(str).tolist() == [
        "maxquant:P28482:S10:row1",
        "maxquant:P31749:S473,T308:row2",
    ]
    assert "site_sequence" not in metadata.columns
    assert metadata.loc[:, "localisation_confidence"].tolist() == pytest.approx(
        [0.97, 0.78]
    )

    evidence = result.peptide_evidence
    assert evidence is not None
    assert evidence.loc[:, "site_string"].tolist() == ["S10", "S473;T308"]
    assert evidence.loc[:, "site_id"].tolist() == [
        "MAPK1;S10;",
        "AKT1;S473,T308;",
    ]
    assert evidence.loc[:, "multi_site"].tolist() == [False, True]
    assert "site_sequence" not in evidence.columns

    maxquant_diagnostics = result.diagnostics["maxquant"]
    assert maxquant_diagnostics["filtering"]["removed_rows"] == 2
    adaptation = maxquant_diagnostics["adaptation"]
    assert adaptation["protein_group_rows_collapsed_to_first_accession"] == 1
    assert adaptation["gene_group_rows_collapsed_to_first_symbol"] == 1
    assert adaptation["multi_site_rows"] == 1
    assert any("protein-group rows" in warning for warning in result.warnings)
    assert any("gene-name group rows" in warning for warning in result.warnings)
    assert any("multi-site candidates" in warning for warning in result.warnings)
    assert result.quality_report.warnings == result.warnings


def test_maxquant_importer_supports_custom_column_mapping() -> None:
    source = pd.DataFrame(
        {
            "accession": ["Q16539", "P28482"],
            "symbol": ["MAPK14", "MAPK1"],
            "site_token": ["Y182", "T185"],
            "loc_percent": ["95", "92"],
            "pep": ["AAAAAYAAAA", "BBBBBTBBBB"],
            "mod_pep": ["AAAAA(ph)YAAAA", "BBBBB(ph)TBBBB"],
            "raw control": ["100", "110"],
            "raw stim": ["120", "130"],
            "contam": ["", ""],
            "rev": ["", ""],
        }
    )

    result = MaxQuantPhosphositeImporter().run(
        MaxQuantPhosphositeImportRequest(
            source=source,
            column_mapping=MaxQuantColumnMapping(
                protein_accession="accession",
                gene_symbol="symbol",
                modified_site="site_token",
                localisation_confidence="loc_percent",
                peptide_sequence="pep",
                modified_peptide_sequence="mod_pep",
                intensity_columns={
                    "raw control": "control",
                    "raw stim": "stim",
                },
                potential_contaminant="contam",
                reverse="rev",
            ),
            localisation_confidence_scale="percent",
            source_name="custom_maxquant",
        )
    )

    assert result.source_name == "custom_maxquant"
    assert result.sample_column_mapping == {
        "raw control": "control",
        "raw stim": "stim",
    }
    assert list(result.phospho_matrix_candidate.columns) == ["control", "stim"]
    assert result.site_metadata_candidate.loc[:, "site"].tolist() == ["Y182", "T185"]
    assert result.site_metadata_candidate.loc[
        :, "localisation_confidence"
    ].tolist() == (pytest.approx([0.95, 0.92]))


def test_maxquant_explicit_depth_mapping_emits_canonical_metadata_and_provenance() -> (
    None
):
    result = MaxQuantPhosphositeImporter().run(
        MaxQuantPhosphositeImportRequest(
            source=_depth_source([3, 4, 5]),
            column_mapping=_depth_mapping(kind=QUANTIFICATION_DEPTH_KIND_PEPTIDE_COUNT),
        )
    )

    metadata = result.site_metadata_candidate
    assert metadata.loc[:, "quantification_depth"].tolist() == pytest.approx(
        [3.0, 4.0, 5.0]
    )
    evidence = result.peptide_evidence
    assert evidence is not None
    assert evidence.loc[:, "quantification_depth"].tolist() == pytest.approx(
        [3.0, 4.0, 5.0]
    )
    maxquant_diagnostics = result.diagnostics["maxquant"]
    assert maxquant_diagnostics["resolved_columns"]["quantification_depth"] == "Depth"
    assert (
        maxquant_diagnostics["resolved_columns"]["quantification_depth_kind"]
        == QUANTIFICATION_DEPTH_KIND_PEPTIDE_COUNT
    )
    depth_diagnostics = maxquant_diagnostics["adaptation"]["quantification_depth"]
    assert depth_diagnostics["status"] == "reported"
    assert depth_diagnostics["source_column"] == "Depth"
    assert depth_diagnostics["output_column"] == "quantification_depth"
    assert depth_diagnostics["quantification_depth_kind"] == (
        QUANTIFICATION_DEPTH_KIND_PEPTIDE_COUNT
    )
    report_depth = result.quality_report.format_specific["maxquant"]["adaptation"][
        "quantification_depth"
    ]
    assert report_depth["source_column"] == "Depth"
    assert report_depth["output_column"] == "quantification_depth"


def test_maxquant_depth_mapping_is_optional_and_not_inferred() -> None:
    result = MaxQuantPhosphositeImporter().run(
        MaxQuantPhosphositeImportRequest(source=_depth_source([3, 4, 5]))
    )

    assert "quantification_depth" not in result.site_metadata_candidate.columns
    evidence = result.peptide_evidence
    assert evidence is not None
    assert "quantification_depth" not in evidence.columns
    maxquant_diagnostics = result.diagnostics["maxquant"]
    assert maxquant_diagnostics["resolved_columns"]["quantification_depth"] is None
    assert (
        maxquant_diagnostics["adaptation"]["quantification_depth"]["status"]
        == "not_mapped"
    )


@pytest.mark.parametrize(
    ("depth_values", "message"),
    (
        ([3, "", 5], "missing values"),
        ([3, "bad", 5], "numeric count values"),
        ([3, 0, 5], ">= 1"),
        ([3, 4.5, 5], "integer-valued"),
    ),
)
def test_maxquant_explicit_depth_mapping_rejects_invalid_counts(
    depth_values: list[object],
    message: str,
) -> None:
    with pytest.raises(PhosPyInputError, match=message):
        MaxQuantPhosphositeImporter().run(
            MaxQuantPhosphositeImportRequest(
                source=_depth_source(depth_values),
                column_mapping=_depth_mapping(),
            )
        )


def test_maxquant_depth_mapping_requires_explicit_supported_depth_kind() -> None:
    with pytest.raises(PhosPyInputError, match="quantification_depth_kind requires"):
        MaxQuantPhosphositeImporter().run(
            MaxQuantPhosphositeImportRequest(
                source=_depth_source([3, 4]),
                column_mapping=MaxQuantColumnMapping(
                    quantification_depth_kind=QUANTIFICATION_DEPTH_KIND_PSM_COUNT
                ),
            )
        )

    with pytest.raises(PhosPyInputError, match="quantification_depth_kind"):
        MaxQuantPhosphositeImporter().run(
            MaxQuantPhosphositeImportRequest(
                source=_depth_source([3, 4]),
                column_mapping=MaxQuantColumnMapping(quantification_depth="Depth"),
            )
        )

    with pytest.raises(PhosPyInputError, match="quantification_depth_kind"):
        MaxQuantPhosphositeImporter().run(
            MaxQuantPhosphositeImportRequest(
                source=_depth_source([3, 4]),
                column_mapping=MaxQuantColumnMapping(
                    quantification_depth="Depth",
                    quantification_depth_kind="spectral_count",  # type: ignore[arg-type]
                ),
            )
        )


def test_maxquant_explicit_intensity_override_selects_lfq_columns() -> None:
    result = MaxQuantPhosphositeImporter().run(
        MaxQuantPhosphositeImportRequest(
            source=FIXTURES / "phospho_sty_sites_raw_and_lfq.txt",
            column_mapping=MaxQuantColumnMapping(
                intensity_columns={
                    "LFQ intensity Control": "control_lfq",
                    "LFQ intensity Stim": "stim_lfq",
                },
            ),
        )
    )

    assert result.sample_column_mapping == {
        "LFQ intensity Control": "control_lfq",
        "LFQ intensity Stim": "stim_lfq",
    }
    phospho = result.phospho_matrix_candidate
    assert list(phospho.columns) == ["control_lfq", "stim_lfq"]
    assert phospho.iloc[:, 0].tolist() == pytest.approx([100.0, 200.0])
    assert phospho.iloc[:, 1].tolist() == pytest.approx([120.0, 210.0])


def test_maxquant_lfq_intensity_columns_are_detected_when_raw_absent() -> None:
    result = MaxQuantPhosphositeImporter().run(
        MaxQuantPhosphositeImportRequest(
            source=FIXTURES / "phospho_sty_sites_lfq_only.txt"
        )
    )

    assert result.sample_column_mapping == {
        "LFQ intensity Control": "Control",
        "LFQ intensity Stim": "Stim",
    }
    phospho = result.phospho_matrix_candidate
    assert list(phospho.columns) == ["Control", "Stim"]
    assert phospho.iloc[:, 0].tolist() == pytest.approx([101.0, 202.0])
    assert phospho.iloc[:, 1].tolist() == pytest.approx([121.0, 222.0])


def test_maxquant_mixed_raw_and_lfq_detection_requires_explicit_mapping() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="multiple intensity columns for the same inferred sample IDs",
    ):
        MaxQuantPhosphositeImporter().run(
            MaxQuantPhosphositeImportRequest(
                source=FIXTURES / "phospho_sty_sites_raw_and_lfq.txt"
            )
        )


def test_maxquant_localisation_probability_strings_become_threshold_ready_numeric() -> (
    None
):
    source = pd.DataFrame(
        {
            "Proteins": ["P28482", "P31749"],
            "Gene names": ["MAPK1", "AKT1"],
            "Amino acid": ["S", "S"],
            "Positions within proteins": ["10", "473"],
            "Phospho (STY) Probabilities": ["S(0.80)", "S(0.93)"],
            "Sequence": ["AAAAASAAAA", "BBBBBSBBBB"],
            "Modified sequence": ["AAAAA(ph)SAAAA", "BBBBB(ph)SBBBB"],
            "Intensity A": ["1.0", "2.0"],
            "Intensity B": ["3.0", "4.0"],
            "Potential contaminant": ["", ""],
            "Reverse": ["", ""],
        }
    )

    result = MaxQuantPhosphositeImporter().run(
        MaxQuantPhosphositeImportRequest(source=source)
    )

    confidence = result.site_metadata_candidate.loc[:, "localisation_confidence"]
    assert confidence.tolist() == pytest.approx([0.80, 0.93])
    assert all(isinstance(value, float) for value in confidence.tolist())
    assert all(0.0 <= float(value) <= 1.0 for value in confidence.tolist())
    assert result.diagnostics["localisation_confidence"]["scale"] == "probability"


def test_maxquant_contaminants_and_reverse_hits_are_removed_by_default() -> None:
    result = MaxQuantPhosphositeImporter().run(
        MaxQuantPhosphositeImportRequest(
            source=FIXTURES / "phospho_sty_sites_standard.txt"
        )
    )

    metadata = result.site_metadata_candidate
    assert metadata.loc[:, "gene_symbol"].tolist() == ["MAPK1", "AKT1"]
    filtering = result.diagnostics["maxquant"]["filtering"]
    assert filtering["potential_contaminant_rows"] == 1
    assert filtering["reverse_rows"] == 1
    assert filtering["removed_rows"] == 2

    report = result.quality_report
    assert report.row_count_status == IMPORTER_QUALITY_STATUS_REPORTED
    assert report.rows_read == 4
    assert report.rows_retained == 2
    assert report.rows_dropped == 2
    assert [
        (column.source_column, column.sample_id)
        for column in report.detected_intensity_columns
    ] == [
        ("Intensity Control", "Control"),
        ("Intensity Stim", "Stim"),
    ]
    assert report.missing_intensity.total_missing_values == 0
    assert report.localisation_confidence.source_column == "Localization prob"
    assert report.localisation_confidence.row_count == 2
    assert report.flagged_rows.contaminant.status == (IMPORTER_QUALITY_STATUS_REPORTED)
    assert report.flagged_rows.contaminant.count == 1
    assert report.flagged_rows.contaminant.source_column == "Potential contaminant"
    assert report.flagged_rows.contaminant.policy == "remove"
    assert report.flagged_rows.reverse.count == 1
    assert report.flagged_rows.reverse.source_column == "Reverse"
    assert report.flagged_rows.reverse.policy == "remove"
    assert report.flagged_rows.decoy.status == IMPORTER_QUALITY_STATUS_NOT_APPLICABLE
    assert report.format_specific["maxquant"]["filtering"]["removed_rows"] == 2


def test_maxquant_contaminants_and_reverse_hits_can_be_flagged() -> None:
    result = MaxQuantPhosphositeImporter().run(
        MaxQuantPhosphositeImportRequest(
            source=FIXTURES / "phospho_sty_sites_standard.txt",
            contaminant_policy="flag",
            reverse_policy="flag",
        )
    )

    metadata = result.site_metadata_candidate
    assert metadata.shape[0] == 4
    assert metadata.loc[:, "maxquant_potential_contaminant"].tolist() == [
        False,
        False,
        True,
        False,
    ]
    assert metadata.loc[:, "maxquant_reverse"].tolist() == [
        False,
        False,
        False,
        True,
    ]
    evidence = result.peptide_evidence
    assert evidence is not None
    assert "maxquant_potential_contaminant" in evidence.columns
    assert "maxquant_reverse" in evidence.columns
    assert result.quality_report.rows_retained == 4
    assert result.quality_report.rows_dropped == 0
    assert result.quality_report.flagged_rows.contaminant.policy == "flag"
    assert result.quality_report.flagged_rows.reverse.policy == "flag"


def test_maxquant_multisite_rows_are_retained_as_peptide_evidence() -> None:
    result = MaxQuantPhosphositeImporter().run(
        MaxQuantPhosphositeImportRequest(
            source=FIXTURES / "phospho_sty_sites_multisite.txt"
        )
    )

    metadata = result.site_metadata_candidate
    assert metadata.loc[:, "site"].tolist() == ["S10,T12", "S473"]
    evidence = result.peptide_evidence
    assert evidence is not None
    assert evidence.loc[:, "site_string"].tolist() == ["S10;T12", "S473"]
    assert evidence.loc[:, "site_id"].tolist() == ["MAPK1;S10,T12;", "AKT1;S473;"]
    assert evidence.loc[:, "multi_site"].tolist() == [True, False]
    assert any("multi-site candidates" in warning for warning in result.warnings)


def test_maxquant_depth_survives_dataset_build_without_protein_group_multiplication() -> (
    None
):
    result = MaxQuantPhosphositeImporter().run(
        MaxQuantPhosphositeImportRequest(
            source=_depth_source(
                [7],
                proteins=["P28482;P28482-2"],
                genes=["MAPK1"],
                residues=["S"],
                positions=[10],
            ),
            column_mapping=_depth_mapping(),
        )
    )
    dataset = _build_dataset(result, site_resolution_mode="site_level_resolved")

    metadata = dataset.site_metadata
    assert metadata.loc[:, "quantification_depth"].tolist() == pytest.approx([7.0])
    assert (
        result.diagnostics["maxquant"]["adaptation"][
            "protein_group_rows_collapsed_to_first_accession"
        ]
        == 1
    )


def test_maxquant_site_matrix_aggregate_preserves_identical_depth_without_summing() -> (
    None
):
    result = MaxQuantPhosphositeImporter().run(
        MaxQuantPhosphositeImportRequest(
            source=_depth_source(
                [7, 7],
                proteins=["P28482", "P28482"],
                genes=["MAPK1", "MAPK1"],
                residues=["S", "S"],
                positions=[10, 10],
            ),
            column_mapping=_depth_mapping(),
        )
    )

    dataset = _build_dataset(
        result,
        site_resolution_mode="site_level_resolved",
        site_matrix_duplicate_policy="aggregate_mean",
    )

    metadata = dataset.site_metadata
    assert metadata.shape[0] == 1
    assert pd.api.types.is_numeric_dtype(metadata.loc[:, "quantification_depth"].dtype)
    assert metadata.loc[:, "quantification_depth"].tolist() == pytest.approx([7.0])


def test_maxquant_site_matrix_aggregate_marks_conflicting_depth_unavailable() -> None:
    result = MaxQuantPhosphositeImporter().run(
        MaxQuantPhosphositeImportRequest(
            source=_depth_source(
                [7, 9],
                proteins=["P28482", "P28482"],
                genes=["MAPK1", "MAPK1"],
                residues=["S", "S"],
                positions=[10, 10],
            ),
            column_mapping=_depth_mapping(),
        )
    )

    dataset = _build_dataset(
        result,
        site_resolution_mode="site_level_resolved",
        site_matrix_duplicate_policy="aggregate_mean",
    )

    metadata = dataset.site_metadata
    assert metadata.shape[0] == 1
    assert pd.api.types.is_numeric_dtype(metadata.loc[:, "quantification_depth"].dtype)
    assert pd.isna(metadata.loc[:, "quantification_depth"].iloc[0])


def test_maxquant_peptide_evidence_build_preserves_identical_collapsed_depth() -> None:
    result = MaxQuantPhosphositeImporter().run(
        MaxQuantPhosphositeImportRequest(
            source=_depth_source(
                [7, 7],
                proteins=["P28482", "P28482"],
                genes=["MAPK1", "MAPK1"],
                residues=["S", "S"],
                positions=[10, 10],
            ),
            column_mapping=_depth_mapping(),
        )
    )

    dataset = _build_dataset(
        result,
        site_resolution_mode=DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE,
    )

    metadata = dataset.site_metadata
    assert metadata.shape[0] == 1
    assert metadata.loc[:, "display_id"].tolist() == ["MAPK1;S10;"]
    assert metadata.loc[:, "quantification_depth"].tolist() == pytest.approx([7.0])


def test_maxquant_peptide_evidence_build_marks_conflicting_depth_unavailable() -> None:
    result = MaxQuantPhosphositeImporter().run(
        MaxQuantPhosphositeImportRequest(
            source=_depth_source(
                [7, 9],
                proteins=["P28482", "P28482"],
                genes=["MAPK1", "MAPK1"],
                residues=["S", "S"],
                positions=[10, 10],
            ),
            column_mapping=_depth_mapping(),
        )
    )

    dataset = _build_dataset(
        result,
        site_resolution_mode=DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE,
    )

    metadata = dataset.site_metadata
    assert metadata.shape[0] == 1
    assert pd.isna(metadata.loc[:, "quantification_depth"].iloc[0])


def test_maxquant_split_multisite_depth_is_unavailable_without_count_multiplication() -> (
    None
):
    result = MaxQuantPhosphositeImporter().run(
        MaxQuantPhosphositeImportRequest(
            source=_depth_source(
                [7],
                proteins=["P28482"],
                genes=["MAPK1"],
                residues=["S"],
                positions=["10;12"],
            ),
            column_mapping=_depth_mapping(),
        )
    )

    dataset = _build_dataset(
        result,
        site_resolution_mode=DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE,
        multi_site_policy=DATASET_MULTI_SITE_POLICY_SPLIT,
    )

    metadata = dataset.site_metadata.sort_index()
    assert metadata.loc[:, "display_id"].tolist() == ["MAPK1;S10;", "MAPK1;S12;"]
    assert metadata.loc[:, "quantification_depth"].isna().all()


def test_maxquant_imported_depth_supports_downstream_depth_aware_differential() -> None:
    import_result = MaxQuantPhosphositeImporter().run(
        MaxQuantPhosphositeImportRequest(
            source=_depth_source([3, 4, 5, 6, 7]),
            column_mapping=_depth_mapping(),
        )
    )
    dataset = _build_dataset(import_result, site_resolution_mode="site_level_resolved")

    result = DifferentialAnalysisWorkflow().run(
        DifferentialAnalysisRequest(
            dataset=dataset,
            design=_depth_design(),
            contrasts=(
                Contrast(
                    name="B_vs_A",
                    numerator_condition="B",
                    denominator_condition="A",
                ),
            ),
            config=DifferentialAnalysisConfig(
                empirical_bayes=EmpiricalBayesConfig(
                    trend=True,
                    trend_covariate=(
                        EMPIRICAL_BAYES_TREND_COVARIATE_QUANTIFICATION_DEPTH
                    ),
                    quantification_depth_kind=QUANTIFICATION_DEPTH_KIND_PSM_COUNT,
                ),
            ),
        )
    )

    diagnostics = result.quantification_depth_trend_diagnostics
    assert diagnostics is not None
    assert diagnostics.quantification_depth.tolist() == pytest.approx(
        [3.0, 4.0, 5.0, 6.0, 7.0]
    )
    assert diagnostics.quantification_depth_kind == QUANTIFICATION_DEPTH_KIND_PSM_COUNT
    assert "B_vs_A" in result.contrast_tables
