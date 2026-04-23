#!/usr/bin/env python3
"""Run the preferred 1.5.0 signalome workflow lane.

Signalome contract note: explicit ``site_metadata.protein_id`` is required.
Gene-symbol site-ID prefixes are not treated as protein-identity fallback.
"""

from __future__ import annotations

import pandas as pd

from phospy import AnalysisReadyDatasetBuilder, KinaseWorkflow, SignalomeWorkflow
from phospy.api import (
    DatasetBuildRequest,
    KinaseWorkflowRequest,
    KinaseWorkflowResult,
    Organism,
    ReferencePreset,
    SignalomeWorkflowRequest,
    SignalomeWorkflowResult,
)


def _build_kinase_result() -> KinaseWorkflowResult:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.00, 0.70],
            "sample_b": [1.10, 0.80],
            "sample_c": [0.95, 0.75],
        },
        index=["TSC2;S939;", "GSK3B;S9;"],
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["TSC2", "GSK3B"],
            "site": ["S939", "S9"],
            "site_sequence": [
                "FDDTPEKDSFRARSTSLNERPKSLRIARAPK",
                "_______MSGRPRTTSFAESCKPVQQPSAFG",
            ],
            # Required in the supported signalome lane.
            "protein_id": ["TSC2", "GSK3B"],
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
    return KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=ReferencePreset.AUTO,
            activity_config=None,
        )
    )


def run_demo() -> SignalomeWorkflowResult:
    kinase_result = _build_kinase_result()
    return SignalomeWorkflow().run(
        SignalomeWorkflowRequest(kinase_result=kinase_result)
    )


def main() -> None:
    result = run_demo()
    protein_ids = result.kinase_result.dataset.site_metadata["protein_id"]
    print("Preferred 1.5.0 signalome workflow lane")
    print(
        "protein_id present for all sites:",
        bool(protein_ids.astype("string").str.strip().ne("").all()),
    )
    print(
        "Resolved reference organism:",
        result.kinase_result.references.organism.value,
    )
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
    if result.expanded_signalome is None:
        raise RuntimeError(
            "expanded_signalome was not materialized in the supported signalome lane"
        )
    print("Expanded signalome shape:", result.expanded_signalome.shape)


if __name__ == "__main__":
    main()
