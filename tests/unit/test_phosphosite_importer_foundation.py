from __future__ import annotations

from dataclasses import replace

import pandas as pd
import pytest

from phospy.api import (
    IMPORTER_QUALITY_STATUS_NOT_APPLICABLE,
    IMPORTER_QUALITY_STATUS_REPORTED,
    PhosphositeImportRequest,
    PhosphositeImportResult,
)
from phospy.errors import PhosPyInputError
from phospy.io.readers import MappedPhosphositeTableImporter
from phospy.science.evidence import DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE


def _source_table() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "feature_id": ["f1", "f2", "f3"],
            "gene": ["MAPK1", "MAPK1", "MAPK1"],
            "site": ["S10", "S10", "S10,T12"],
            "protein": ["P28482", "P28482", "P28482"],
            "sequence_window": [
                ("A" * 15) + "S" + ("A" * 15),
                ("A" * 15) + "S" + ("A" * 15),
                ("A" * 15) + "S" + ("A" * 15),
            ],
            "loc_percent": ["95", "91", "88"],
            "sample A raw": ["10.0", "12.0", "14.0"],
            "sample B raw": ["11.0", "13.0", "15.0"],
            "peptide": ["AAAA", "BBBB", "CCCC"],
            "modified_peptide": ["AA[pS]AA", "BB[pS]BB", "CC[pS,pT]CC"],
            "site_string": ["S10", "S10", "S10;T12"],
        }
    )


def _request(source: pd.DataFrame) -> PhosphositeImportRequest:
    return PhosphositeImportRequest(
        source=source,
        sample_intensity_columns={
            "sample A raw": "sample_a",
            "sample B raw": "sample_b",
        },
        gene_symbol_column="gene",
        site_column="site",
        protein_id_column="protein",
        site_sequence_column="sequence_window",
        localisation_confidence_column="loc_percent",
        localisation_confidence_scale="percent",
        unique_feature_id_column="feature_id",
        peptide_sequence_column="peptide",
        modified_peptide_sequence_column="modified_peptide",
        peptide_site_string_column="site_string",
        source_name="synthetic_search_output",
    )


def test_importer_requires_explicit_sample_column_mapping() -> None:
    source = _source_table().drop(columns=["sample B raw"])

    with pytest.raises(
        PhosPyInputError,
        match="missing required columns: sample B raw",
    ):
        MappedPhosphositeTableImporter().run(_request(source))


def test_importer_extracts_sample_intensities_and_normalises_localisation() -> None:
    source = _source_table().copy(deep=True)
    source.loc[1, "loc_percent"] = "not_numeric"

    result = MappedPhosphositeTableImporter().run(_request(source))

    assert isinstance(result, PhosphositeImportResult)
    phospho = result.phospho_matrix_candidate
    assert list(phospho.columns) == ["sample_a", "sample_b"]
    assert phospho.iloc[0, 0] == pytest.approx(10.0)
    assert phospho.iloc[2, 1] == pytest.approx(15.0)
    metadata = result.site_metadata_candidate
    assert result.localisation_confidence_column == "localisation_confidence"
    assert metadata.iloc[0]["localisation_confidence"] == pytest.approx(0.95)
    assert pd.isna(metadata.iloc[1]["localisation_confidence"])
    localisation = result.diagnostics["localisation_confidence"]
    assert localisation["scale"] == "percent"
    assert localisation["invalid_count"] == 1
    assert any("invalid values" in warning for warning in result.warnings)


def test_importer_populates_quality_report_from_parsing_facts() -> None:
    source = _source_table().copy(deep=True)
    source.loc[1, "sample B raw"] = ""
    source["site_key"] = [
        "protein:P28482:S10",
        "protein:P28482:S10",
        "protein:P28482:S10,T12",
    ]
    source["display_id"] = ["MAPK1;S10;", "MAPK1;S10;", "MAPK1;S10,T12;"]
    request = replace(
        _request(source),
        site_key_column="site_key",
        display_id_column="display_id",
    )

    result = MappedPhosphositeTableImporter().run(request)

    report = result.quality_report
    assert report.source_name == "synthetic_search_output"
    assert report.row_count_status == IMPORTER_QUALITY_STATUS_REPORTED
    assert report.rows_read == 3
    assert report.rows_retained == 3
    assert report.rows_dropped == 0
    assert [
        (column.source_column, column.sample_id)
        for column in report.detected_intensity_columns
    ] == [
        ("sample A raw", "sample_a"),
        ("sample B raw", "sample_b"),
    ]
    assert report.missing_intensity.total_missing_values == 1
    assert report.missing_intensity.rows_with_any_missing_intensity == 1
    assert report.missing_intensity.missing_values_by_source_column == {
        "sample A raw": 0,
        "sample B raw": 1,
    }
    assert report.localisation_confidence.status == IMPORTER_QUALITY_STATUS_REPORTED
    assert report.localisation_confidence.source_column == "loc_percent"
    assert report.localisation_confidence.missing_count == 0
    assert report.flagged_rows.contaminant.status == (
        IMPORTER_QUALITY_STATUS_NOT_APPLICABLE
    )
    assert report.duplicate_keys.site_key.count == 1
    assert report.duplicate_keys.display_key.count == 1
    assert report.duplicate_keys.duplicate_site_candidate_rows == 1
    assert report.warnings == result.warnings


def test_importer_quality_report_marks_unmapped_localisation_not_applicable() -> None:
    source = _source_table().drop(columns=["loc_percent"])
    request = replace(_request(source), localisation_confidence_column=None)

    result = MappedPhosphositeTableImporter().run(request)

    report = result.quality_report
    assert result.localisation_confidence_column is None
    assert "localisation_confidence" not in result.site_metadata_candidate.columns
    assert report.localisation_confidence.status == (
        IMPORTER_QUALITY_STATUS_NOT_APPLICABLE
    )
    assert report.localisation_confidence.reason == (
        "localisation confidence column was not mapped"
    )


def test_importer_preserves_duplicate_and_multisite_peptide_evidence() -> None:
    result = MappedPhosphositeTableImporter().run(_request(_source_table()))

    evidence = result.peptide_evidence
    assert evidence is not None
    assert evidence.shape[0] == 3
    assert evidence.loc[:, "site_id"].tolist() == [
        "MAPK1;S10;",
        "MAPK1;S10;",
        "MAPK1;S10,T12;",
    ]
    assert evidence.loc[:, "multi_site"].tolist() == [False, False, True]
    assert evidence.loc[:, "peptide_sequence"].tolist() == ["AAAA", "BBBB", "CCCC"]
    assert {"sample_a", "sample_b"}.issubset(set(evidence.columns))
    assert result.diagnostics["duplicate_site_candidate_rows"] == 1
    assert result.diagnostics["multi_site_candidate_rows"] == 1
    assert any("duplicate site candidates" in warning for warning in result.warnings)
    assert any("multi-site candidates" in warning for warning in result.warnings)


def test_importer_result_can_handoff_peptide_evidence_to_builder_request() -> None:
    result = MappedPhosphositeTableImporter().run(_request(_source_table()))

    request = result.to_dataset_build_request(
        site_resolution_mode=DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE,
        multi_site_policy="split",
        input_intensity_scale="linear",
    )

    assert request.site_resolution_mode == DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE
    assert request.phospho is None
    assert request.site_metadata is None
    assert request.peptide_evidence is not None
    assert request.peptide_evidence_sample_intensity_columns == (
        "sample_a",
        "sample_b",
    )
    assert request.multi_site_policy == "split"
