#!/usr/bin/env python3
"""Run the supported signalome workflow route from bundled-reference kinase output.

Signalome contract note: explicit ``site_metadata.protein_id`` is required.
Gene-symbol site-ID prefixes are not treated as protein-identity fallback.
"""

from __future__ import annotations

import pandas as pd

from phospy import (
    AnalysisReadyDatasetBuilder,
    KinaseWorkflow,
    SignalomeWorkflow,
)
from phospy.api import (
    DatasetBuildRequest,
    KinaseWorkflowRequest,
    KinaseWorkflowResult,
    Organism,
    ReferencePreset,
    SignalomeConfig,
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
        SignalomeWorkflowRequest(
            kinase_result=kinase_result,
            config=SignalomeConfig(
                substrate_support_cutoff=0.5,
                network_correlation_threshold=0.5,
                score_preconditioning_policy="allow_and_report",
            ),
        )
    )


def main() -> None:
    result = run_demo()
    print("Signalome workflow demo")
    print(
        "Upstream profile score shape:",
        result.kinase_result.scoring_result.profile_scores.shape,
    )
    print(
        "Upstream prediction shape:",
        result.kinase_result.prediction_result.pred_mat.shape,
    )
    if result.kinase_result.activity_result is not None:
        print(
            "Upstream weighted activity shape:",
            result.kinase_result.activity_result.weighted_activity.shape,
        )
    else:
        print("Upstream activity output: disabled")
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
    print(
        "Module selection strategy/count:",
        result.module_selection_diagnostics.strategy,
        result.module_selection_diagnostics.selected_module_count,
    )
    print(
        "Score preconditioning policy/counts:",
        result.score_preconditioning_diagnostics.policy,
        {
            "input_rows": result.score_preconditioning_diagnostics.input_row_count,
            "dropped_all_missing_rows": (
                result.score_preconditioning_diagnostics.dropped_all_missing_row_count
            ),
            "retained_rows": result.score_preconditioning_diagnostics.retained_row_count,
        },
    )
    if result.expanded_signalome is None:
        raise RuntimeError(
            "expanded_signalome was not materialized in the supported signalome lane"
        )
    print("Expanded signalome shape:", result.expanded_signalome.shape)


if __name__ == "__main__":
    main()
