#!/usr/bin/env python3
"""Run the signalome workflow over a real kinase result."""

from __future__ import annotations

import pandas as pd

from phospy import (
    AnalysisReadyDatasetBuilder,
    DatasetBuildRequest,
    KinaseWorkflow,
    Organism,
    ReferencePreset,
    SignalomeConfig,
    SignalomeWorkflow,
    SignalomeWorkflowRequest,
    SignalomeWorkflowResult,
    SimpleKinaseWorkflowRequest,
    SimpleKinaseWorkflowResult,
)
from phospy.transformations.models import TransformationState


def _build_kinase_result() -> SimpleKinaseWorkflowResult:
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
            transformation_state=TransformationState.raw(has_total_matrix=False),
        )
    )
    return KinaseWorkflow().run(
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
