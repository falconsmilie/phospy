from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from phospy.api import PhosphositeImportResult
from phospy.api.results import (
    IMPORTER_QUALITY_STATUS_NOT_APPLICABLE,
    IMPORTER_QUALITY_STATUS_REPORTED,
)
from phospy.errors import PhosPyInputError
from phospy.io.readers import (
    MaxQuantColumnMapping,
    MaxQuantPhosphositeImporter,
    MaxQuantPhosphositeImportRequest,
)

ROOT = Path(__file__).resolve().parents[2]
FIXTURES = ROOT / "tests" / "fixtures" / "maxquant"


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
