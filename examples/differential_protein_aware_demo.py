#!/usr/bin/env python3
"""Run the experimental protein-aware differential workflow."""

from __future__ import annotations

import pandas as pd

from phospy import AnalysisReadyDatasetBuilder, DifferentialAnalysisWorkflow
from phospy.advanced import (
    DatasetProteinAwarePreparationConfig,
    DifferentialAnalysisConfig,
    DifferentialProteinAwareModelConfig,
)
from phospy.api import (
    Contrast,
    DatasetBuildRequest,
    DatasetPreprocessingConfig,
    DifferentialAnalysisRequest,
    ExperimentalDesign,
    Organism,
    SampleDesignRecord,
)

METHOD_ID = "protein_covariate_adjusted_moderated_linear_model_v1"


def _example_tables() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    phospho = pd.DataFrame(
        {
            "A_1": [1.00, 2.05, 1.48],
            "A_2": [1.15, 2.10, 1.50],
            "A_3": [0.95, 1.92, 1.46],
            "B_1": [1.75, 2.48, 1.55],
            "B_2": [1.83, 2.57, 1.58],
            "B_3": [1.69, 2.41, 1.62],
        },
        index=["MAPK14;Y182;", "AKT1;T308;", "GSK3B;S9;"],
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "AKT1", "GSK3B"],
            "site": ["Y182", "T308", "S9"],
            "protein_id": ["P53778", "P31749", "P49841"],
            "site_sequence": [
                "AAAAAAAAAAAAAAAYAAAAAAAAAAAAAAA",
                "AAAAAAAAAAAAAAATAAAAAAAAAAAAAAA",
                "AAAAAAAAAAAAAAASAAAAAAAAAAAAAAA",
            ],
            "localisation_confidence": [0.95, 0.96, 0.97],
        },
        index=phospho.index.copy(),
    )
    total = pd.DataFrame(
        {
            "A_1": [10.0, 8.0, 12.5],
            "A_2": [10.4, 8.5, 12.5],
            "A_3": [9.7, 7.8, 12.5],
            "B_1": [11.8, 8.9, 12.5],
            "B_2": [12.2, 9.7, 12.5],
            "B_3": [11.4, 8.8, 12.5],
        },
        index=pd.Index(["P53778", "P31749", "P49841"], name="protein_id"),
    )
    return phospho, site_metadata, total


def _build_prepared_dataset():
    phospho, site_metadata, total = _example_tables()
    return AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            total=total,
            organism=Organism.RAT,
            input_intensity_scale="log2",
            preprocessing_config=DatasetPreprocessingConfig(
                protein_aware_preparation=DatasetProteinAwarePreparationConfig(
                    policy="prepare_model_inputs"
                )
            ),
        )
    )


def _design() -> ExperimentalDesign:
    return ExperimentalDesign(
        samples=(
            SampleDesignRecord(
                sample_id="A_1",
                condition="A",
                biological_replicate_id="A_1_bio",
            ),
            SampleDesignRecord(
                sample_id="A_2",
                condition="A",
                biological_replicate_id="A_2_bio",
            ),
            SampleDesignRecord(
                sample_id="A_3",
                condition="A",
                biological_replicate_id="A_3_bio",
            ),
            SampleDesignRecord(
                sample_id="B_1",
                condition="B",
                biological_replicate_id="B_1_bio",
            ),
            SampleDesignRecord(
                sample_id="B_2",
                condition="B",
                biological_replicate_id="B_2_bio",
            ),
            SampleDesignRecord(
                sample_id="B_3",
                condition="B",
                biological_replicate_id="B_3_bio",
            ),
        )
    )


def main() -> None:
    dataset = _build_prepared_dataset()
    result = DifferentialAnalysisWorkflow().run(
        DifferentialAnalysisRequest(
            dataset=dataset,
            design=_design(),
            contrasts=(
                Contrast(
                    name="B_vs_A",
                    numerator_condition="B",
                    denominator_condition="A",
                ),
            ),
            config=DifferentialAnalysisConfig(
                protein_aware_model=DifferentialProteinAwareModelConfig(
                    method=METHOD_ID
                )
            ),
        )
    )

    table = result.table_for("B_vs_A")
    tested = table.loc[table["result_status"] == "tested"]
    withheld = table.loc[table["result_status"] != "tested"]
    diagnostics = result.protein_aware_diagnostics
    if diagnostics is None:
        raise AssertionError("protein-aware diagnostics were not returned")

    print("Experimental protein-aware differential workflow")
    print("Method:", diagnostics.method_id)
    print("Claim:", diagnostics.claim_status)
    print("Tested sites:", diagnostics.tested_site_count)
    print("Withheld sites:", diagnostics.withheld_site_count)
    print("Fallback policy:", diagnostics.fallback_policy)
    print("Tested rows:")
    print(tested.loc[:, ["display_id", "logFC", "P.Value", "adj.P.Val"]])
    print("Withheld rows:")
    print(withheld.loc[:, ["display_id", "result_status", "result_status_reason"]])
    print("Protein-aware diagnostics:")
    print(
        diagnostics.per_site_diagnostics_dataframe().loc[
            :, ["total_protein_row_key", "protein_covariate_coefficient"]
        ]
    )
    print("Caveat codes:", [caveat.code for caveat in result.caveats])


if __name__ == "__main__":
    main()
