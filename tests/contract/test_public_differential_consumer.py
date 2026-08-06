from __future__ import annotations

import pandas as pd
import pytest

from phospy import AnalysisReadyDatasetBuilder, DifferentialAnalysisWorkflow
from phospy.api import (
    AnalysisReadyPhosphoDataset,
    Contrast,
    DatasetBuildRequest,
    DifferentialAnalysisRequest,
    DifferentialAnalysisResult,
    ExperimentalDesign,
    Organism,
    SampleDesignRecord,
)


def test_public_consumer_builds_and_runs_differential_by_site_key() -> None:
    phospho = pd.DataFrame(
        {
            "control_rep1": [10.0, 8.00],
            "control_rep2": [10.2, 8.15],
            "treated_rep1": [11.0, 8.10],
            "treated_rep2": [11.3, 8.05],
        },
        index=pd.Index(["MAPK14;Y182;", "GSK3B;S9;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "GSK3B"],
            "protein_id": ["MAPK14", "GSK3B"],
            "site": ["Y182", "S9"],
            "site_sequence": [
                ("A" * 15) + "Y" + ("A" * 15),
                ("A" * 15) + "S" + ("A" * 15),
            ],
            "localisation_confidence": [0.95, 0.93],
        },
        index=phospho.index.copy(),
    )

    dataset = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
            input_intensity_scale="log2",
        )
    )

    assert isinstance(dataset, AnalysisReadyPhosphoDataset)
    assert dataset.phospho.index.name == "site_key"
    assert dataset.site_metadata.loc[:, "site_key"].tolist() == (
        dataset.phospho.index.astype(str).tolist()
    )

    result = DifferentialAnalysisWorkflow().run(
        DifferentialAnalysisRequest(
            dataset=dataset,
            design=ExperimentalDesign(
                samples=(
                    SampleDesignRecord(
                        sample_id="control_rep1",
                        condition="control",
                        biological_replicate_id="control_r1",
                    ),
                    SampleDesignRecord(
                        sample_id="control_rep2",
                        condition="control",
                        biological_replicate_id="control_r2",
                    ),
                    SampleDesignRecord(
                        sample_id="treated_rep1",
                        condition="treated",
                        biological_replicate_id="treated_r1",
                    ),
                    SampleDesignRecord(
                        sample_id="treated_rep2",
                        condition="treated",
                        biological_replicate_id="treated_r2",
                    ),
                )
            ),
            contrasts=(
                Contrast(
                    name="treated_vs_control",
                    numerator_condition="treated",
                    denominator_condition="control",
                ),
            ),
        )
    )

    assert isinstance(result, DifferentialAnalysisResult)
    table = result.table_for("treated_vs_control")
    expected_site_keys = dataset.phospho.index.astype(str).tolist()

    assert table.index.name == "site_key"
    assert table.index.astype(str).tolist() == expected_site_keys
    assert table.index.is_unique
    assert table.loc[:, "site_key"].astype(str).tolist() == expected_site_keys
    assert {
        "site_key",
        "display_id",
        "gene_symbol",
        "site",
        "logFC",
        "P.Value",
        "adj.P.Val",
    } <= set(table.columns)

    first_site_key = expected_site_keys[0]
    matched = table.loc[[first_site_key]]
    assert matched.shape[0] == 1
    assert matched.iloc[0]["display_id"] == "MAPK14;Y182;"

    with pytest.raises(KeyError):
        table.loc["phospy:v1:missing-site-key"]
