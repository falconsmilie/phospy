from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from phospy.api import (
    IMPORTER_QUALITY_STATUS_NOT_APPLICABLE,
    IMPORTER_QUALITY_STATUS_REPORTED,
    PhosphositeImportResult,
)
from phospy.errors import PhosPyInputError
from phospy.io.readers import (
    FragPipeColumnMapping,
    FragPipePTMProphetImporter,
    FragPipePTMProphetImportRequest,
)
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset

ROOT = Path(__file__).resolve().parents[2]
FIXTURES = ROOT / "tests" / "fixtures" / "fragpipe"


def _import_fixture(
    filename: str = "ptmprophet_sites.tsv",
    **kwargs: object,
) -> PhosphositeImportResult:
    return FragPipePTMProphetImporter().run(
        FragPipePTMProphetImportRequest(
            source=FIXTURES / filename,
            **kwargs,
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
    assert result.diagnostics["fragpipe"]["filtering"]["removed_rows"] == 1


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
        result.diagnostics["fragpipe"]["adaptation"]["ambiguous_localisation_rows"] == 1
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
    assert result.diagnostics["fragpipe"]["adaptation"]["multi_site_rows"] == 1


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
        result.diagnostics["fragpipe"]["adaptation"]["ambiguous_localisation_rows"] == 1
    )
    assert any("ambiguous localisation" in warning for warning in result.warnings)


def test_fragpipe_reports_protein_group_collapse_and_sequence_mismatch() -> None:
    result = _import_fixture(
        "ptmprophet_explicit_site_edge_cases.tsv",
        ptmprophet_position_reference="protein",
    )

    diagnostics = result.diagnostics["fragpipe"]["adaptation"]
    assert diagnostics["protein_group_rows_collapsed_to_first_accession"] == 1
    assert diagnostics["peptide_sequence_mismatch_rows"] == 1
    assert any("protein-group rows" in warning for warning in result.warnings)
    assert any("peptide sequence" in warning for warning in result.warnings)


def test_fragpipe_excludes_decoys_and_contaminants_by_default() -> None:
    result = _import_fixture(
        "ptmprophet_explicit_site_edge_cases.tsv",
        ptmprophet_position_reference="protein",
    )

    filtering = result.diagnostics["fragpipe"]["filtering"]
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
    assert (
        report.format_specific["fragpipe_ptmprophet"]["filtering"]["decoy_prefix_rows"]
        == 1
    )
    assert (
        report.format_specific["fragpipe_ptmprophet"]["adaptation"][
            "ambiguous_localisation_rows"
        ]
        == 1
    )
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
    assert result.diagnostics["fragpipe"]["filtering"]["removed_rows"] == 0
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
