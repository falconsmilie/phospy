from __future__ import annotations

import pandas as pd
import pytest

from phospy.api import (
    AnalysisReadyDatasetBuilder,
    DatasetBuildRequest,
    KinaseWorkflow,
    KinaseWorkflowRequest,
    Organism,
    PhosphositeImportRequest,
    PhosphositeImportResult,
)
from phospy.errors import PhosPyInputError, WorkflowValidationError
from phospy.io.readers import (
    FragPipePTMProphetImporter,
    FragPipePTMProphetImportRequest,
    MappedPhosphositeTableImporter,
    MaxQuantPhosphositeImporter,
    MaxQuantPhosphositeImportRequest,
)
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset


class _RecordingMappedImporter:
    def __init__(self) -> None:
        self.requests: list[PhosphositeImportRequest] = []
        self._delegate = MappedPhosphositeTableImporter()

    def run(self, request: PhosphositeImportRequest) -> PhosphositeImportResult:
        self.requests.append(request)
        return self._delegate.run(request)


def _boundary_source() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "feature_id": ["f1", "f2", "f3"],
            "gene": ["MAPK1", "MAPK1", "MAPK1"],
            "site": ["S10", "S10", "S10,T12"],
            "protein": ["P28482", "P28482", "P28482"],
            "sample A raw": ["10.0", "12.0", "14.0"],
            "sample B raw": ["11.0", "13.0", "15.0"],
            "peptide": ["AAAA", "BBBB", "CCCC"],
            "modified_peptide": ["AA[pS]AA", "BB[pS]BB", "CC[pS,pT]CC"],
            "site_string": ["S10", "S10", "S10;T12"],
        }
    )


def _boundary_request(source: pd.DataFrame) -> PhosphositeImportRequest:
    return PhosphositeImportRequest(
        source=source,
        sample_intensity_columns={
            "sample A raw": "sample_a",
            "sample B raw": "sample_b",
        },
        gene_symbol_column="gene",
        site_column="site",
        protein_id_column="protein",
        unique_feature_id_column="feature_id",
        peptide_sequence_column="peptide",
        modified_peptide_sequence_column="modified_peptide",
        peptide_site_string_column="site_string",
        source_name="boundary_fixture",
    )


def _maxquant_source() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "Proteins": ["P28482"],
            "Gene names": ["MAPK1"],
            "Modified site": ["S10"],
            "Localization prob": ["0.95"],
            "Sequence": ["AAAAASAAAA"],
            "Modified sequence": ["AAAAA(ph)SAAAA"],
            "Intensity Control": ["10.0"],
            "Intensity Stim": ["12.0"],
            "Potential contaminant": [""],
            "Reverse": [""],
        }
    )


def _fragpipe_source() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "Protein": ["sp|P28482|MK01_HUMAN"],
            "Gene": ["MAPK1"],
            "Peptide": ["AAAAASAAAA"],
            "Modified Peptide": ["AAAAA[pS]AAAA"],
            "PTMProphet Probability": ["S10(0.95)"],
            "Site": ["S10"],
            "Intensity Control": ["10.0"],
            "Intensity Stim": ["12.0"],
        }
    )


def test_mapped_importer_returns_candidate_result_not_analysis_ready_dataset() -> None:
    result = MappedPhosphositeTableImporter().run(_boundary_request(_boundary_source()))

    assert isinstance(result, PhosphositeImportResult)
    assert not isinstance(result, AnalysisReadyPhosphoDataset)

    build_request = result.to_dataset_build_request(input_intensity_scale="linear")
    assert isinstance(build_request, DatasetBuildRequest)
    assert not isinstance(build_request, AnalysisReadyPhosphoDataset)
    assert result.site_metadata_candidate.loc[:, "site"].tolist() == [
        "S10",
        "S10",
        "S10,T12",
    ]
    assert result.diagnostics["duplicate_site_candidate_rows"] == 1
    assert result.diagnostics["multi_site_candidate_rows"] == 1


def test_import_result_defers_strict_dataset_validation_to_builder() -> None:
    result = MappedPhosphositeTableImporter().run(_boundary_request(_boundary_source()))

    request = result.to_dataset_build_request(
        input_intensity_scale="linear",
        organism=Organism.HUMAN,
    )
    with pytest.raises(
        PhosPyInputError,
        match="requires strict residue/position site token",
    ):
        AnalysisReadyDatasetBuilder().run(request)


def test_maxquant_importer_delegates_to_shared_mapped_importer_path() -> None:
    mapped_importer = _RecordingMappedImporter()

    result = MaxQuantPhosphositeImporter(mapped_importer=mapped_importer).run(
        MaxQuantPhosphositeImportRequest(
            source=_maxquant_source(),
            source_name="maxquant_boundary",
        )
    )

    assert isinstance(result, PhosphositeImportResult)
    assert not isinstance(result, AnalysisReadyPhosphoDataset)
    assert len(mapped_importer.requests) == 1
    delegated_request = mapped_importer.requests[0]
    assert isinstance(delegated_request, PhosphositeImportRequest)
    assert delegated_request.source_name == "maxquant_boundary"
    assert isinstance(delegated_request.source, pd.DataFrame)


def test_fragpipe_importer_delegates_to_shared_mapped_importer_path() -> None:
    mapped_importer = _RecordingMappedImporter()

    result = FragPipePTMProphetImporter(mapped_importer=mapped_importer).run(
        FragPipePTMProphetImportRequest(
            source=_fragpipe_source(),
            ptmprophet_position_reference="protein",
            source_name="fragpipe_boundary",
        )
    )

    assert isinstance(result, PhosphositeImportResult)
    assert not isinstance(result, AnalysisReadyPhosphoDataset)
    assert len(mapped_importer.requests) == 1
    delegated_request = mapped_importer.requests[0]
    assert isinstance(delegated_request, PhosphositeImportRequest)
    assert delegated_request.source_name == "fragpipe_boundary"
    assert isinstance(delegated_request.source, pd.DataFrame)


def test_kinase_workflow_rejects_raw_vendor_tables_at_dataset_boundary() -> None:
    request = KinaseWorkflowRequest(dataset=_maxquant_source())  # type: ignore[arg-type]

    with pytest.raises(
        WorkflowValidationError,
        match="dataset must be AnalysisReadyPhosphoDataset",
    ):
        KinaseWorkflow().run(request)
