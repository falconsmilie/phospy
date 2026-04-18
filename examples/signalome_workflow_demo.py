#!/usr/bin/env python3
"""Run the supported signalome workflow route over a kinase workflow result."""

from __future__ import annotations

import pandas as pd

from phospy import (
    AnalysisReadyDatasetBuilder,
    DatasetBuildRequest,
    KinaseScoringConfig,
    KinaseWorkflow,
    KinaseWorkflowRequest,
    KinaseWorkflowResult,
    Organism,
    ReferenceBundle,
    SignalomeConfig,
    SignalomeWorkflow,
    SignalomeWorkflowRequest,
    SignalomeWorkflowResult,
)


def _build_kinase_result() -> KinaseWorkflowResult:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.00, 0.70],
            "sample_b": [1.15, 0.80],
            "sample_c": [0.95, 0.75],
        },
        index=["MAPK14;Y182;", "GSK3B;S9;"],
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "GSK3B"],
            "site": ["Y182", "S9"],
            "site_sequence": [
                "LDFGLARHTDDEMTGYVATRWYRAPEIMLNW",
                "RARTSSFAEPGGGGGGGGGPGGSASPARPAR",
            ],
        },
        index=phospho.index.copy(),
    )
    dataset = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
        )
    )
    references = ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["MAP2K6", "MAP2K6", "AKT1", "AKT1"],
                "substrate_site": [
                    "MAPK14;Y182;",
                    "GSK3B;S9;",
                    "MAPK14;Y182;",
                    "GSK3B;S9;",
                ],
            }
        ),
        site_sequences=pd.DataFrame(
            {"site_sequence": dataset.site_metadata.loc[:, "site_sequence"]},
            index=pd.Index(dataset.site_metadata.index, name="site_id"),
        ),
    )
    return KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=references,
            scoring_config=KinaseScoringConfig(min_substrates=2),
            activity_config=None,
        )
    )


def run_demo() -> SignalomeWorkflowResult:
    kinase_result = _build_kinase_result()
    return SignalomeWorkflow().run(
        SignalomeWorkflowRequest(
            kinase_result=kinase_result,
            config=SignalomeConfig(
                substrate_support_cutoff=0.5,
                network_correlation_threshold=0.5,
            ),
        )
    )


def main() -> None:
    result = run_demo()
    print("Signalome workflow demo")
    print(
        "Module assignment shape:",
        result.module_assignments.table.shape,
    )
    print(
        "Signalome module shape:",
        result.signalome_modules.table.shape,
    )
    print(
        "Kinase network edge shape:",
        result.kinase_network.edges.shape,
    )


if __name__ == "__main__":
    main()
