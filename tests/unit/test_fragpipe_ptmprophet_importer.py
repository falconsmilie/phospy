from __future__ import annotations

from pathlib import Path
from typing import Any, cast

import pandas as pd
import pytest

from phospy.advanced import (
    DatasetIntensityTransformConfig,
    DatasetSiteMatrixConfig,
    DatasetSiteMatrixDuplicateSitePolicy,
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
    FragPipeColumnMapping,
    FragPipePTMProphetImporter,
    FragPipePTMProphetImportRequest,
)
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.science.differential.models import (
    EMPIRICAL_BAYES_TREND_COVARIATE_QUANTIFICATION_DEPTH,
    QUANTIFICATION_DEPTH_KIND_PEPTIDE_COUNT,
    QUANTIFICATION_DEPTH_KIND_PSM_COUNT,
)

ROOT = Path(__file__).resolve().parents[2]
FIXTURES = ROOT / "tests" / "fixtures" / "fragpipe"


def test_fragpipe_reader_split_modules_preserve_legacy_import_identity() -> None:
    import importlib

    import phospy.io.readers.fragpipe as fragpipe_reader
    from phospy.io.readers.fragpipe.importer import FragPipePTMProphetImporter
    from phospy.io.readers.fragpipe.models import (
        FragPipeColumnMapping,
        FragPipePTMProphetImportRequest,
    )

    for module_name in (
        "columns",
        "conversion",
        "filtering",
        "normalization",
        "raw",
        "reporting",
    ):
        importlib.import_module(f"phospy.io.readers.fragpipe.{module_name}")

    assert fragpipe_reader.FragPipeColumnMapping is FragPipeColumnMapping
    assert fragpipe_reader.FragPipePTMProphetImporter is FragPipePTMProphetImporter
    assert (
        fragpipe_reader.FragPipePTMProphetImportRequest
        is FragPipePTMProphetImportRequest
    )


def _import_fixture(
    filename: str = "ptmprophet_sites.tsv",
    **kwargs: Any,
) -> PhosphositeImportResult:
    return FragPipePTMProphetImporter().run(
        FragPipePTMProphetImportRequest(
            source=FIXTURES / filename,
            **kwargs,
        )
    )


def _depth_source(
    depth: list[object],
    *,
    proteins: list[str] | None = None,
    genes: list[str] | None = None,
    residues: list[str] | None = None,
    positions: list[object] | None = None,
) -> pd.DataFrame:
    gene_values = genes or ["MAPK1", "AKT1", "GSK3B", "MTOR", "RPS6KB1"][: len(depth)]
    residue_values = residues or ["S", "S", "S", "S", "T"][: len(depth)]
    position_values: tuple[object, ...] = (
        tuple(positions)
        if positions is not None
        else tuple([10, 473, 9, 2448, 389][: len(depth)])
    )
    protein_values = (
        proteins or ["P28482", "P31749", "P49841", "P42345", "P23443"][: len(depth)]
    )
    return pd.DataFrame(
        {
            "Protein": protein_values,
            "Gene": gene_values,
            "Peptide": ["AAAAA" + residue + "AAAA" for residue in residue_values],
            "Modified Peptide": [
                "AAAAA[p" + residue + "]AAAA" for residue in residue_values
            ],
            "PTMProphet Probability": [
                f"{residue}{position}(0.95)"
                for residue, position in zip(
                    residue_values,
                    position_values,
                    strict=True,
                )
            ],
            "Site": [
                f"{residue}{position}"
                for residue, position in zip(
                    residue_values,
                    position_values,
                    strict=True,
                )
            ],
            "Sequence Window": [
                ("A" * 15) + residue + ("A" * 15) for residue in residue_values
            ],
            "Intensity A_1": [10.0 + row for row in range(len(depth))],
            "Intensity A_2": [11.0 + row for row in range(len(depth))],
            "Intensity B_1": [20.0 + row for row in range(len(depth))],
            "Intensity B_2": [22.0 + row for row in range(len(depth))],
            "Depth": depth,
            "Spectrum": [f"scan.{row + 1}.{row + 1}.2" for row in range(len(depth))],
            "PSM ID": [f"psm-{row + 1}" for row in range(len(depth))],
        }
    )


def _depth_mapping(
    *,
    kind: str = QUANTIFICATION_DEPTH_KIND_PSM_COUNT,
) -> FragPipeColumnMapping:
    return FragPipeColumnMapping(
        quantification_depth="Depth",
        quantification_depth_kind=kind,  # type: ignore[arg-type]
    )


def _build_dataset(
    import_result: PhosphositeImportResult,
    *,
    site_resolution_mode: str,
    multi_site_policy: str | None = None,
    site_matrix_duplicate_policy: DatasetSiteMatrixDuplicateSitePolicy | None = None,
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


def _fragpipe_diagnostics(result: PhosphositeImportResult) -> dict[str, Any]:
    return cast(dict[str, Any], result.diagnostics["fragpipe"])


def _fragpipe_report(result: PhosphositeImportResult) -> dict[str, Any]:
    return cast(
        dict[str, Any],
        result.quality_report.format_specific["fragpipe_ptmprophet"],
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


def test_fragpipe_ptmprophet_importer_reads_single_site_peptide() -> None:
    result = _import_fixture()

    assert isinstance(result, PhosphositeImportResult)
    phospho = result.phospho_matrix_candidate
    assert list(phospho.columns) == ["Control", "Stim"]
    assert phospho.shape == (3, 2)
    assert float(phospho.iloc[0]["Control"]) == pytest.approx(10.0)
    metadata = result.site_metadata_candidate
    first = metadata.iloc[0]
    assert first["gene_symbol"] == "MAPK1"
    assert first["site"] == "S10"
    assert first["protein_accession"] == "P28482"
    assert first["localisation_confidence"] == pytest.approx(0.95)
    assert first["fragpipe_ptmprophet_candidate_sites"] == "S10"
    assert first["fragpipe_ptmprophet_site_probabilities"] == "S10:0.95"
    assert bool(first["fragpipe_ptmprophet_ambiguous"]) is False
    assert result.localisation_confidence_column == "localisation_confidence"
    assert _fragpipe_diagnostics(result)["filtering"]["removed_rows"] == 1


def test_fragpipe_explicit_psm_count_depth_mapping_emits_canonical_metadata() -> None:
    result = FragPipePTMProphetImporter().run(
        FragPipePTMProphetImportRequest(
            source=_depth_source([3, 4, 5]),
            column_mapping=_depth_mapping(),
            ptmprophet_position_reference="protein",
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
    fragpipe_diagnostics = _fragpipe_diagnostics(result)
    assert fragpipe_diagnostics["resolved_columns"]["quantification_depth"] == "Depth"
    assert (
        fragpipe_diagnostics["resolved_columns"]["quantification_depth_kind"]
        == QUANTIFICATION_DEPTH_KIND_PSM_COUNT
    )
    depth_diagnostics = fragpipe_diagnostics["adaptation"]["quantification_depth"]
    assert (
        depth_diagnostics["quantification_depth_kind"]
        == QUANTIFICATION_DEPTH_KIND_PSM_COUNT
    )
    report_depth = _fragpipe_report(result)["adaptation"]["quantification_depth"]
    assert (
        report_depth["quantification_depth_kind"] == QUANTIFICATION_DEPTH_KIND_PSM_COUNT
    )


def test_fragpipe_explicit_peptide_count_depth_mapping_is_supported() -> None:
    result = FragPipePTMProphetImporter().run(
        FragPipePTMProphetImportRequest(
            source=_depth_source([2, 3]),
            column_mapping=_depth_mapping(kind=QUANTIFICATION_DEPTH_KIND_PEPTIDE_COUNT),
            ptmprophet_position_reference="protein",
        )
    )

    metadata = result.site_metadata_candidate
    assert metadata.loc[:, "quantification_depth"].tolist() == pytest.approx([2.0, 3.0])
    assert (
        _fragpipe_diagnostics(result)["resolved_columns"]["quantification_depth_kind"]
        == QUANTIFICATION_DEPTH_KIND_PEPTIDE_COUNT
    )


def test_fragpipe_depth_mapping_is_optional_and_not_inferred() -> None:
    result = FragPipePTMProphetImporter().run(
        FragPipePTMProphetImportRequest(
            source=_depth_source([3, 4, 5]),
            ptmprophet_position_reference="protein",
        )
    )

    assert "quantification_depth" not in result.site_metadata_candidate.columns
    evidence = result.peptide_evidence
    assert evidence is not None
    assert "quantification_depth" not in evidence.columns
    fragpipe_diagnostics = _fragpipe_diagnostics(result)
    assert fragpipe_diagnostics["resolved_columns"]["quantification_depth"] is None
    assert (
        fragpipe_diagnostics["adaptation"]["quantification_depth"]["status"]
        == "not_mapped"
    )


def test_fragpipe_spectrum_and_psm_identifiers_are_not_counted_as_depth() -> None:
    result = FragPipePTMProphetImporter().run(
        FragPipePTMProphetImportRequest(
            source=_depth_source([3, 4]),
            ptmprophet_position_reference="protein",
        )
    )

    assert _fragpipe_diagnostics(result)["resolved_columns"]["unique_feature_id"] == (
        "Spectrum"
    )
    assert "quantification_depth" not in result.site_metadata_candidate.columns
    assert (
        _fragpipe_diagnostics(result)["adaptation"]["quantification_depth"]["status"]
        == "not_mapped"
    )


@pytest.mark.parametrize(
    ("depth_values", "message"),
    (
        ([3, "", 5], "missing values"),
        ([3, "bad", 5], "numeric count values"),
        ([3, True, 5], "not booleans"),
        ([3, "inf", 5], "finite numeric count values"),
        ([3, 0, 5], ">= 1"),
        ([3, 4.5, 5], "integer-valued"),
    ),
)
def test_fragpipe_explicit_depth_mapping_rejects_invalid_counts(
    depth_values: list[object],
    message: str,
) -> None:
    with pytest.raises(PhosPyInputError, match=message):
        FragPipePTMProphetImporter().run(
            FragPipePTMProphetImportRequest(
                source=_depth_source(depth_values),
                column_mapping=_depth_mapping(),
                ptmprophet_position_reference="protein",
            )
        )


def test_fragpipe_depth_mapping_requires_explicit_supported_depth_kind() -> None:
    with pytest.raises(PhosPyInputError, match="quantification_depth_kind requires"):
        FragPipePTMProphetImporter().run(
            FragPipePTMProphetImportRequest(
                source=_depth_source([3, 4]),
                column_mapping=FragPipeColumnMapping(
                    quantification_depth_kind=QUANTIFICATION_DEPTH_KIND_PSM_COUNT
                ),
                ptmprophet_position_reference="protein",
            )
        )

    with pytest.raises(PhosPyInputError, match="quantification_depth_kind"):
        FragPipePTMProphetImporter().run(
            FragPipePTMProphetImportRequest(
                source=_depth_source([3, 4]),
                column_mapping=FragPipeColumnMapping(quantification_depth="Depth"),
                ptmprophet_position_reference="protein",
            )
        )

    with pytest.raises(PhosPyInputError, match="quantification_depth_kind"):
        FragPipePTMProphetImporter().run(
            FragPipePTMProphetImportRequest(
                source=_depth_source([3, 4]),
                column_mapping=FragPipeColumnMapping(
                    quantification_depth="Depth",
                    quantification_depth_kind="spectral_count",  # type: ignore[arg-type]
                ),
                ptmprophet_position_reference="protein",
            )
        )


def test_fragpipe_dataset_builder_preserves_imported_depth() -> None:
    import_result = FragPipePTMProphetImporter().run(
        FragPipePTMProphetImportRequest(
            source=_depth_source([3, 4, 5]),
            column_mapping=_depth_mapping(),
            ptmprophet_position_reference="protein",
        )
    )

    dataset = _build_dataset(import_result, site_resolution_mode="site_level_resolved")

    assert dataset.site_metadata.loc[:, "quantification_depth"].tolist() == (
        pytest.approx([3.0, 4.0, 5.0])
    )


def test_fragpipe_site_matrix_aggregate_marks_conflicting_depth_unavailable() -> None:
    import_result = FragPipePTMProphetImporter().run(
        FragPipePTMProphetImportRequest(
            source=_depth_source(
                [7, 9],
                proteins=["P28482", "P28482"],
                genes=["MAPK1", "MAPK1"],
                residues=["S", "S"],
                positions=[10, 10],
            ),
            column_mapping=_depth_mapping(),
            ptmprophet_position_reference="protein",
        )
    )

    dataset = _build_dataset(
        import_result,
        site_resolution_mode="site_level_resolved",
        site_matrix_duplicate_policy="aggregate_mean",
    )

    metadata = dataset.site_metadata
    assert metadata.shape[0] == 1
    assert pd.api.types.is_numeric_dtype(metadata.loc[:, "quantification_depth"].dtype)
    assert pd.isna(metadata.loc[:, "quantification_depth"].iloc[0])


def test_fragpipe_split_multisite_depth_is_unavailable_without_count_multiplication() -> (
    None
):
    import_result = FragPipePTMProphetImporter().run(
        FragPipePTMProphetImportRequest(
            source=pd.DataFrame(
                {
                    "Protein": ["P28482"],
                    "Gene": ["MAPK1"],
                    "Peptide": ["AAAAASASAAA"],
                    "Modified Peptide": ["AAAAA[pS]A[pS]AAA"],
                    "PTMProphet Probability": ["S10(0.95);S12(0.93)"],
                    "Site": ["S10,S12"],
                    "Sequence Window": [("A" * 15) + "S" + ("A" * 15)],
                    "Intensity A_1": [10.0],
                    "Intensity A_2": [11.0],
                    "Intensity B_1": [20.0],
                    "Intensity B_2": [22.0],
                    "Depth": [7],
                }
            ),
            column_mapping=_depth_mapping(),
            ptmprophet_position_reference="protein",
        )
    )

    dataset = _build_dataset(
        import_result,
        site_resolution_mode=DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE,
        multi_site_policy=DATASET_MULTI_SITE_POLICY_SPLIT,
    )

    metadata = dataset.site_metadata.sort_values("display_id")
    assert metadata.loc[:, "display_id"].tolist() == ["MAPK1;S10;", "MAPK1;S12;"]
    assert metadata.loc[:, "quantification_depth"].isna().all()


def test_fragpipe_imported_depth_supports_downstream_depth_aware_differential() -> None:
    import_result = FragPipePTMProphetImporter().run(
        FragPipePTMProphetImportRequest(
            source=_depth_source([3, 4, 5, 6, 7]),
            column_mapping=_depth_mapping(),
            ptmprophet_position_reference="protein",
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


def test_fragpipe_ptmprophet_importer_retains_multi_site_peptide_evidence() -> None:
    result = _import_fixture()

    metadata = result.site_metadata_candidate
    assert metadata.iloc[1]["site"] == "S473,T475"
    assert metadata.iloc[1]["localisation_confidence"] == pytest.approx(0.97)
    evidence = result.peptide_evidence
    assert evidence is not None
    assert evidence.iloc[1]["site_string"] == "S473;T475"
    assert evidence.iloc[1]["site_id"] == "AKT1;S473,T475;"
    assert bool(evidence.iloc[1]["multi_site"]) is True
    assert evidence.iloc[1]["modified_peptide_sequence"] == "AA[pS]A[pT]AA"
    assert evidence.iloc[1]["fragpipe_modified_peptide_phospho_count"] == 2
    assert any("multi-site candidates" in warning for warning in result.warnings)


def test_fragpipe_ptmprophet_ambiguous_localisation_is_joint_not_first_site() -> None:
    result = _import_fixture()

    metadata = result.site_metadata_candidate
    ambiguous = metadata.iloc[2]
    assert ambiguous["gene_symbol"] == "GSK3B"
    assert ambiguous["site"] == "S10,T11"
    assert ambiguous["localisation_confidence"] == pytest.approx(0.50)
    assert ambiguous["fragpipe_ptmprophet_candidate_sites"] == "S10;T11"
    assert ambiguous["fragpipe_ptmprophet_site_probabilities"] == ("S10:0.5;T11:0.5")
    assert bool(ambiguous["fragpipe_ptmprophet_ambiguous"]) is True
    evidence = result.peptide_evidence
    assert evidence is not None
    assert evidence.iloc[2]["site_string"] == "S10;T11"
    assert evidence.iloc[2]["site_id"] == "GSK3B;S10,T11;"
    assert bool(evidence.iloc[2]["multi_site"]) is True
    assert any("ambiguous localisation" in warning for warning in result.warnings)
    assert (
        _fragpipe_diagnostics(result)["adaptation"]["ambiguous_localisation_rows"] == 1
    )


def test_fragpipe_ptmprophet_malformed_localisation_string_fails() -> None:
    source = pd.read_csv(FIXTURES / "ptmprophet_sites.tsv", sep="\t")
    source.loc[0, "PTMProphet Probability"] = "S4=not_a_probability"

    with pytest.raises(PhosPyInputError, match="malformed PTMProphet"):
        FragPipePTMProphetImporter().run(FragPipePTMProphetImportRequest(source=source))


def test_fragpipe_ptmprophet_contaminants_can_be_flagged() -> None:
    result = _import_fixture(contaminant_policy="flag")

    metadata = result.site_metadata_candidate
    assert metadata.shape[0] == 4
    assert metadata.iloc[3]["protein_accession"] == "P99999"
    assert metadata.loc[:, "fragpipe_contaminant"].tolist() == [
        False,
        False,
        False,
        True,
    ]
    evidence = result.peptide_evidence
    assert evidence is not None
    assert "fragpipe_contaminant" in evidence.columns
    assert bool(evidence.iloc[3]["fragpipe_contaminant"]) is True


def test_fragpipe_ptmprophet_supports_protein_position_localisation_strings() -> None:
    source = pd.DataFrame(
        {
            "protein": ["sp|Q16539|MK14_HUMAN"],
            "gene": ["MAPK14"],
            "peptide": ["AAAAAYAAAA"],
            "mod": ["AAAAA[pY]AAAA"],
            "ptm": ["Y182(0.93)"],
            "raw control": ["100.0"],
            "raw stim": ["120.0"],
        }
    )

    result = FragPipePTMProphetImporter().run(
        FragPipePTMProphetImportRequest(
            source=source,
            column_mapping=FragPipeColumnMapping(
                protein_accession="protein",
                gene_symbol="gene",
                peptide_sequence="peptide",
                modified_peptide_sequence="mod",
                ptmprophet_probabilities="ptm",
                intensity_columns={
                    "raw control": "control",
                    "raw stim": "stim",
                },
            ),
            ptmprophet_position_reference="protein",
            source_name="custom_fragpipe",
        )
    )

    assert result.source_name == "custom_fragpipe"
    assert result.sample_column_mapping == {
        "raw control": "control",
        "raw stim": "stim",
    }
    assert result.site_metadata_candidate.iloc[0]["site"] == "Y182"
    assert result.site_metadata_candidate.iloc[0][
        "localisation_confidence"
    ] == pytest.approx(0.93)


def test_fragpipe_explicit_protein_site_fixture_remains_import_candidate() -> None:
    result = _import_fixture(
        "ptmprophet_explicit_site_edge_cases.tsv",
        ptmprophet_position_reference="protein",
    )

    assert isinstance(result, PhosphositeImportResult)
    assert not isinstance(result, AnalysisReadyPhosphoDataset)
    metadata = result.site_metadata_candidate
    assert list(result.phospho_matrix_candidate.columns) == ["Control", "Stim"]
    assert metadata.loc[:, "protein_accession"].tolist() == [
        "P28482",
        "Q9Y243",
        "P49841",
        "Q13554",
    ]
    assert metadata.loc[:, "gene_symbol"].tolist() == [
        "MAPK1",
        "AKT3",
        "GSK3B",
        "CAMK2B",
    ]
    assert metadata.loc[:, "site"].tolist() == ["S10", "S472,T474", "S9", "S17"]
    assert metadata.iloc[0]["localisation_confidence"] == pytest.approx(0.95)
    assert metadata.iloc[0]["fragpipe_ptmprophet_site_probabilities"] == "S10:0.95"

    request = result.to_dataset_build_request(input_intensity_scale="linear")
    assert request.phospho is not None
    assert request.site_metadata is not None
    assert request.peptide_evidence is None


def test_fragpipe_peptide_position_fixture_maps_ptmprophet_string_variants() -> None:
    result = _import_fixture("ptmprophet_peptide_position_edge_cases.tsv")

    metadata = result.site_metadata_candidate
    assert metadata.loc[:, "site"].tolist() == [
        "Y182",
        "S473,T475",
        "S10",
        "S203",
    ]
    assert metadata.loc[:, "fragpipe_ptmprophet_site_probabilities"].tolist() == [
        "Y182:0.93",
        "S473:0.98;T475:0.97",
        "S10:0.88",
        "S203:0.91",
    ]
    assert metadata.loc[:, "localisation_confidence"].tolist() == pytest.approx(
        [0.93, 0.97, 0.88, 0.91]
    )

    evidence = result.peptide_evidence
    assert evidence is not None
    assert evidence.loc[:, "site_string"].tolist() == [
        "Y182",
        "S473;T475",
        "S10",
        "S203",
    ]
    assert evidence.loc[:, "multi_site"].tolist() == [False, True, False, False]
    assert _fragpipe_diagnostics(result)["adaptation"]["multi_site_rows"] == 1


def test_fragpipe_explicit_site_ambiguous_localisation_is_diagnostic() -> None:
    result = _import_fixture(
        "ptmprophet_explicit_site_edge_cases.tsv",
        ptmprophet_position_reference="protein",
    )

    metadata = result.site_metadata_candidate
    ambiguous = metadata.loc[metadata.loc[:, "gene_symbol"] == "GSK3B"].iloc[0]
    assert ambiguous["site"] == "S9"
    assert ambiguous["fragpipe_ptmprophet_candidate_sites"] == "S9;T10"
    assert ambiguous["fragpipe_ptmprophet_site_probabilities"] == "S9:0.5;T10:0.5"
    assert bool(ambiguous["fragpipe_ptmprophet_ambiguous"]) is True
    assert (
        _fragpipe_diagnostics(result)["adaptation"]["ambiguous_localisation_rows"] == 1
    )
    assert any("ambiguous localisation" in warning for warning in result.warnings)


def test_fragpipe_reports_protein_group_collapse_and_sequence_mismatch() -> None:
    result = _import_fixture(
        "ptmprophet_explicit_site_edge_cases.tsv",
        ptmprophet_position_reference="protein",
    )

    diagnostics = _fragpipe_diagnostics(result)["adaptation"]
    assert diagnostics["protein_group_rows_collapsed_to_first_accession"] == 1
    assert diagnostics["peptide_sequence_mismatch_rows"] == 1
    assert any("protein-group rows" in warning for warning in result.warnings)
    assert any("peptide sequence" in warning for warning in result.warnings)


def test_fragpipe_excludes_decoys_and_contaminants_by_default() -> None:
    result = _import_fixture(
        "ptmprophet_explicit_site_edge_cases.tsv",
        ptmprophet_position_reference="protein",
    )

    filtering = _fragpipe_diagnostics(result)["filtering"]
    assert filtering["input_row_count"] == 6
    assert filtering["contaminant_rows"] == 1
    assert filtering["decoy_rows"] == 1
    assert filtering["removed_rows"] == 2
    assert filtering["retained_row_count"] == 4
    retained_genes = result.site_metadata_candidate.loc[:, "gene_symbol"].tolist()
    assert "CONGENE" not in retained_genes
    assert "DECOY" not in retained_genes

    report = result.quality_report
    assert report.row_count_status == IMPORTER_QUALITY_STATUS_REPORTED
    assert report.rows_read == 6
    assert report.rows_retained == 4
    assert report.rows_dropped == 2
    assert [
        (column.source_column, column.sample_id)
        for column in report.detected_intensity_columns
    ] == [
        ("Intensity Control", "Control"),
        ("Intensity Stim", "Stim"),
    ]
    assert report.missing_intensity.total_missing_values == 0
    assert report.localisation_confidence.source_column == "PTMProphet Probability"
    assert report.localisation_confidence.row_count == 4
    assert report.flagged_rows.contaminant.status == (IMPORTER_QUALITY_STATUS_REPORTED)
    assert report.flagged_rows.contaminant.count == 1
    assert report.flagged_rows.contaminant.source_column == "Contaminant"
    assert report.flagged_rows.contaminant.policy == "remove"
    assert "prefix" in str(report.flagged_rows.contaminant.reason)
    assert report.flagged_rows.reverse.status == IMPORTER_QUALITY_STATUS_NOT_APPLICABLE
    assert report.flagged_rows.decoy.count == 1
    assert report.flagged_rows.decoy.source_column == "Decoy"
    assert report.flagged_rows.decoy.policy == "remove"
    assert _fragpipe_report(result)["filtering"]["decoy_prefix_rows"] == 1
    assert _fragpipe_report(result)["adaptation"]["ambiguous_localisation_rows"] == 1
    assert report.warnings == result.warnings


def test_fragpipe_can_flag_decoys_and_contaminants_when_requested() -> None:
    result = _import_fixture(
        "ptmprophet_explicit_site_edge_cases.tsv",
        ptmprophet_position_reference="protein",
        contaminant_policy="flag",
        decoy_policy="flag",
    )

    metadata = result.site_metadata_candidate
    assert metadata.shape[0] == 6
    assert metadata.loc[:, "fragpipe_contaminant"].tolist() == [
        False,
        False,
        False,
        False,
        True,
        False,
    ]
    assert metadata.loc[:, "fragpipe_decoy"].tolist() == [
        False,
        False,
        False,
        False,
        False,
        True,
    ]
    assert _fragpipe_diagnostics(result)["filtering"]["removed_rows"] == 0
    assert result.quality_report.rows_retained == 6
    assert result.quality_report.rows_dropped == 0
    assert result.quality_report.flagged_rows.contaminant.policy == "flag"
    assert result.quality_report.flagged_rows.decoy.policy == "flag"


def test_fragpipe_missing_required_protein_start_rejects_peptide_positions() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="FragPipe Protein Start must contain non-empty values",
    ):
        _import_fixture("ptmprophet_missing_required_start.tsv")
