#!/usr/bin/env python3
"""Run the supported signalome workflow route over a kinase workflow result."""

from __future__ import annotations

import pandas as pd

from phospy import (
    AnalysisReadyDatasetBuilder,
    DatasetBuildRequest,
    DatasetMissingDataConfig,
    DatasetPreprocessingConfig,
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
            "protein_id": ["MAPK14", "GSK3B"],
        },
        index=phospho.index.copy(),
    )
    dataset = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
            preprocessing_config=DatasetPreprocessingConfig(
                missing_data=DatasetMissingDataConfig(policy="forbid"),
            ),
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
    if result.expanded_signalome is None:
        raise RuntimeError(
            "expanded_signalome was not materialized in the supported signalome lane"
        )
    print("Expanded signalome shape:", result.expanded_signalome.shape)


if __name__ == "__main__":
    main()
