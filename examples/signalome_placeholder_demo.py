#!/usr/bin/env python3
"""Run the signalome workflow over a real kinase result."""

from __future__ import annotations

import pandas as pd

from phospy import (
    AnalysisReadyDatasetBuilder,
    DatasetBuildRequest,
    Organism,
    ReferencePreset,
    SignalomeConfig,
    SignalomeWorkflow,
    SignalomeWorkflowRequest,
    SignalomeWorkflowResult,
    SimpleKinaseWorkflow,
    SimpleKinaseWorkflowRequest,
    SimpleKinaseWorkflowResult,
)


def _build_kinase_result() -> SimpleKinaseWorkflowResult:
    phospho = pd.DataFrame(
        {"sample_a": [1.0], "sample_b": [1.2]},
        index=["MAPK14;Y182;"],
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14"],
            "site": ["Y182"],
            "site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"],
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
    return SimpleKinaseWorkflow().run(
        SimpleKinaseWorkflowRequest(dataset=dataset, references=ReferencePreset.AUTO)
    )


def run_demo() -> SignalomeWorkflowResult:
    kinase_result = _build_kinase_result()
    return SignalomeWorkflow().run(
        SignalomeWorkflowRequest(
            kinase_result=kinase_result,
            config=SignalomeConfig(signalome_cutoff=0.5),
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
