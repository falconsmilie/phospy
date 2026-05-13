from __future__ import annotations

import pandas as pd
import pytest

from phospy import AnalysisReadyDatasetBuilder
from phospy.api import (
    DatasetBuildRequest,
    DatasetIntensityTransformConfig,
    DatasetPreprocessingConfig,
)

pytestmark = pytest.mark.integration


def _site_metadata(index: pd.Index) -> pd.DataFrame:
    gene_symbols: list[str] = []
    sites: list[str] = []
    for site_id in index.astype(str):
        parts = site_id.split(";")
        gene_symbols.append(parts[0])
        sites.append(parts[1] if len(parts) > 1 else "S1")
    rows = len(gene_symbols)
    return pd.DataFrame(
        {
            "gene_symbol": gene_symbols,
            "site": sites,
            "site_sequence": ["SEQ_A"] * rows,
            "localisation_confidence": [0.95] * rows,
        },
        index=index.copy(),
    )


def _workflow_establishment_payload(dataset) -> dict[str, object]:
    assert dataset.provenance is not None
    payload = dataset.provenance.workflow_parameters.get(
        "intensity_scale_establishment"
    )
    assert isinstance(payload, dict)
    return payload


def _final_stage_establishment_payload(dataset) -> dict[str, object]:
    assert dataset.preprocessing_report is not None
    final_stage = dataset.preprocessing_report.operations.loc[
        dataset.preprocessing_report.operations.loc[:, "stage"]
        == "final_dataset_construction"
    ].iloc[0]
    parameters = final_stage["parameters"]
    assert isinstance(parameters, dict)
    payload = parameters.get("intensity_scale_establishment")
    assert isinstance(payload, dict)
    return payload


def test_builder_declared_log2_records_declared_establishment_mode_in_provenance() -> (
    None
):
    phospho = pd.DataFrame(
        {"sample_a": [1.1, 2.0], "sample_b": [2.0, 3.1]},
        index=pd.Index(["MAPK14;Y182;", "GSK3B;S9;"], name="site_id"),
    )
    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=_site_metadata(phospho.index),
            input_intensity_scale="log2",
            preprocessing_config=DatasetPreprocessingConfig(
                intensity_transform=DatasetIntensityTransformConfig(policy="identity")
            ),
        )
    )

    workflow_payload = _workflow_establishment_payload(built)
    final_stage_payload = _final_stage_establishment_payload(built)
    assert workflow_payload["establishment_mode"] == "declared"
    assert workflow_payload["input_declaration_source"] == (
        "dataset_build_request.input_intensity_scale"
    )
    assert final_stage_payload["establishment_mode"] == "declared"


def test_builder_log2_transformation_records_transformed_mode_in_provenance() -> None:
    phospho = pd.DataFrame(
        {"sample_a": [3.0, 7.0], "sample_b": [15.0, 31.0]},
        index=pd.Index(["MAPK14;Y182;", "GSK3B;S9;"], name="site_id"),
    )
    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=_site_metadata(phospho.index),
            preprocessing_config=DatasetPreprocessingConfig(
                intensity_transform=DatasetIntensityTransformConfig(
                    policy="log2",
                    pseudocount=1.0,
                )
            ),
        )
    )

    workflow_payload = _workflow_establishment_payload(built)
    assert workflow_payload["establishment_mode"] == "transformed"
    assert workflow_payload["transformer_name"] == (
        "phospy.science.datasets.preprocessing.stages.intensity_transform.log2"
    )
    assert workflow_payload["parameters"]["operation"] == "log2"
    assert workflow_payload["parameters"]["pseudocount"] == 1.0


def test_builder_identity_pass_through_records_identity_mode_in_provenance() -> None:
    phospho = pd.DataFrame(
        {"sample_a": [100.0], "sample_b": [200.0]},
        index=pd.Index(["MAPK14;Y182;"], name="site_id"),
    )
    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=_site_metadata(phospho.index),
        )
    )

    workflow_payload = _workflow_establishment_payload(built)
    assert workflow_payload["establishment_mode"] == "identity"
    assert "IdentityTransformer" in str(workflow_payload["transformer_name"])


def test_builder_records_suspicious_declared_log2_warning_in_provenance() -> None:
    phospho = pd.DataFrame(
        {"sample_a": [12000.0, 14000.0], "sample_b": [18000.0, 22000.0]},
        index=pd.Index(["MAPK14;Y182;", "GSK3B;S9;"], name="site_id"),
    )
    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=_site_metadata(phospho.index),
            input_intensity_scale="log2",
            preprocessing_config=DatasetPreprocessingConfig(
                intensity_transform=DatasetIntensityTransformConfig(policy="identity")
            ),
        )
    )

    workflow_payload = _workflow_establishment_payload(built)
    warnings = workflow_payload["diagnostic_warnings"]
    assert isinstance(warnings, list)
    assert any("declared log2 scale is suspicious" in warning for warning in warnings)


def test_builder_records_suspicious_declared_linear_warning_in_provenance() -> None:
    phospho = pd.DataFrame(
        {"sample_a": [-1.0, 3.0], "sample_b": [2.0, 4.0]},
        index=pd.Index(["MAPK14;Y182;", "GSK3B;S9;"], name="site_id"),
    )
    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=_site_metadata(phospho.index),
            input_intensity_scale="linear",
            preprocessing_config=DatasetPreprocessingConfig(
                intensity_transform=DatasetIntensityTransformConfig(policy="identity")
            ),
        )
    )

    workflow_payload = _workflow_establishment_payload(built)
    warnings = workflow_payload["diagnostic_warnings"]
    assert isinstance(warnings, list)
    assert any(
        "declared linear scale contains negative values" in warning
        for warning in warnings
    )


def test_builder_plausible_declared_scale_has_no_warning() -> None:
    phospho = pd.DataFrame(
        {"sample_a": [0.2, 1.0], "sample_b": [2.2, 4.1]},
        index=pd.Index(["MAPK14;Y182;", "GSK3B;S9;"], name="site_id"),
    )
    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=_site_metadata(phospho.index),
            input_intensity_scale="log2",
            preprocessing_config=DatasetPreprocessingConfig(
                intensity_transform=DatasetIntensityTransformConfig(policy="identity")
            ),
        )
    )

    workflow_payload = _workflow_establishment_payload(built)
    warnings = workflow_payload["diagnostic_warnings"]
    assert isinstance(warnings, list)
    assert warnings == []
