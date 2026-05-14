from __future__ import annotations

import pandas as pd
import pytest

from phospy import (
    AnalysisReadyDatasetBuilder,
    AnalysisReadyPhosphoDataset,
    KinaseWorkflow,
)
from phospy.api import (
    DatasetBuildRequest,
    KinaseWorkflowRequest,
    Organism,
    ReferencePreset,
)

pytestmark = pytest.mark.integration


def test_readme_style_kinase_workflow_builds_and_runs() -> None:
    phospho = pd.DataFrame(
        {
            "control_rep1": [8200.0, 9100.0, 6000.0],
            "control_rep2": [8000.0, 9000.0, 5900.0],
            "treatment_rep1": [16200.0, 9150.0, 13000.0],
            "treatment_rep2": [15800.0, 9050.0, 12800.0],
        },
        index=["MAPK14;Y182;", "GSK3B;S9;", "TSC2;S939;"],
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "GSK3B", "TSC2"],
            "site": ["Y182", "S9", "S939"],
            "site_sequence": [
                "AAAAAAAAAAAAAAAYAAAAAAAAAAAAAAA",
                "AAAAAAAAAAAAAAASAAAAAAAAAAAAAAA",
                "AAAAAAAAAAAAAAASAAAAAAAAAAAAAAA",
            ],
            "protein_id": ["MAPK14", "GSK3B", "TSC2"],
            "localisation_confidence": [0.95, 0.94, 0.96],
        },
        index=phospho.index.copy(),
    )
    sample_metadata = pd.DataFrame(
        {"condition": ["control", "control", "treatment", "treatment"]},
        index=phospho.columns.copy(),
    )

    dataset = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            sample_metadata=sample_metadata,
            organism=Organism.RAT,
            input_intensity_scale="linear",
        )
    )
    assert isinstance(dataset, AnalysisReadyPhosphoDataset)
    assert {"gene_symbol", "site", "site_sequence"} <= set(
        dataset.site_metadata.columns
    )
    assert set(dataset.phospho.index.astype(str)) == {
        "MAPK14;Y182;",
        "GSK3B;S9;",
        "TSC2;S939;",
    }

    result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=ReferencePreset.AUTO,
            activity_config=None,
        )
    )

    assert result.prediction_result.pred_mat.shape[0] == phospho.shape[0]
    assert result.prediction_result.pred_mat.shape[1] > 0
    assert set(result.prediction_result.pred_mat.index.astype(str)) == set(
        phospho.index.astype(str)
    )
    assert result.prediction_result.substrate_list is not None
