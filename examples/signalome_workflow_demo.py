#!/usr/bin/env python3
"""Run the preferred 1.5.0 signalome workflow lane.

Signalome contract note: explicit ``site_metadata.protein_group_id`` is
required. Legacy ``protein_id`` is accepted only as a migration alias.
Gene-symbol site-ID prefixes are not treated as grouping or protein-identity
fallback.
"""

from __future__ import annotations

from dataclasses import replace

import pandas as pd

from phospy import AnalysisReadyDatasetBuilder, KinaseWorkflow, SignalomeWorkflow
from phospy.advanced import (
    KinaseReliabilityProfile,
    KinaseScoringConfig,
    ReferenceContextCompatibilityPolicy,
    SignalomeConfig,
)
from phospy.api import (
    DatasetBuildRequest,
    DatasetLocalisationConfig,
    DatasetPreprocessingConfig,
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
            "sample_a": [1.00, 0.70, 0.85, 0.92, 0.66],
            "sample_b": [1.10, 0.80, 0.88, 0.96, 0.69],
            "sample_c": [0.95, 0.75, 0.92, 0.90, 0.72],
        },
        index=[
            "TSC2;S939;",
            "GSK3A;S21;",
            "MAPK14;Y182;",
            "AKT1;T308;",
            "SRC;Y416;",
        ],
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["TSC2", "GSK3A", "MAPK14", "AKT1", "SRC"],
            "site": ["S939", "S21", "Y182", "T308", "Y416"],
            "site_sequence": [
                "FDDTPEKDSFRARSTSLNERPKSLRIARAPK",
                "PSGGGPGGSGRARTSSFAEPGGGGGGGGGGP",
                "LDFGLARHTDDEMTGYVATRWYRAPEIMLNW",
                ("A" * 15) + "T" + ("A" * 15),
                ("A" * 15) + "Y" + ("A" * 15),
            ],
            "display_id": [
                "TSC2;S939;",
                "GSK3A;S21;",
                "MAPK14;Y182;",
                "AKT1;T308;",
                "SRC;Y416;",
            ],
            "organism": ["rat"] * phospho.shape[0],
            "protein_namespace": ["protein_id"] * phospho.shape[0],
            "protein_identifier": ["TSC2", "GSK3A", "MAPK14", "AKT1", "SRC"],
            "localisation_confidence": [0.95] * phospho.shape[0],
            # Required in the supported signalome lane.
            "protein_group_id": ["TSC2", "GSK3A", "MAPK14", "AKT1", "SRC"],
        },
        index=phospho.index.copy(),
    )
    dataset = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
            input_intensity_scale="linear",
            preprocessing_config=DatasetPreprocessingConfig(
                # Fail fast on missing/below-threshold localisation confidence
                # to keep site-level module interpretation scientifically safe.
                localisation=DatasetLocalisationConfig(
                    mode="require_threshold",
                    confidence_column="localisation_confidence",
                    min_confidence=0.75,
                )
            ),
        )
    )
    protein_group_ids = (
        dataset.site_metadata["protein_group_id"].astype("string").str.strip()
    )
    if not protein_group_ids.ne("").all():
        raise RuntimeError(
            "Supported signalome lane requires explicit non-empty "
            "site_metadata.protein_group_id for every interpreted site."
        )
    return KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=ReferencePreset.AUTO,
            scoring_config=KinaseScoringConfig(
                reliability_profile=KinaseReliabilityProfile.CUSTOM,
                reference_context_compatibility_policy=(
                    ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT
                ),
            ),
            activity_config=None,
            site_sequence_conflict_policy="prefer_reference",
        )
    )


def run_demo() -> SignalomeWorkflowResult:
    kinase_result = _build_kinase_result()
    production_config = SignalomeConfig.production()
    return SignalomeWorkflow().run(
        SignalomeWorkflowRequest(
            kinase_result=kinase_result,
            config=replace(
                production_config,
                validation=replace(
                    production_config.validation,
                    reference_context_compatibility_policy=(
                        ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT
                    ),
                ),
            ),
        )
    )


def main() -> None:
    result = run_demo()
    protein_group_ids = result.kinase_result.dataset.site_metadata["protein_group_id"]
    print("Preferred 1.5.0 signalome workflow lane")
    print(
        result.kinase_result.dataset.site_metadata.loc[
            :,
            [
                "site_key",
                "display_id",
                "gene_symbol",
                "site",
                "protein_namespace",
                "protein_identifier",
                "protein_group_id",
            ],
        ]
    )
    print(
        "protein_group_id present for all sites:",
        bool(protein_group_ids.astype("string").str.strip().ne("").all()),
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
