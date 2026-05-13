from __future__ import annotations

import pandas as pd
import pytest

from phospy import (
    AnalysisReadyDatasetBuilder,
    AnalysisReadyPhosphoDataset,
    DifferentialAnalysisWorkflow,
    KinaseWorkflow,
    SignalomeWorkflow,
)
from phospy.api import (
    Contrast,
    DatasetBuildRequest,
    DatasetPreprocessingConfig,
    DatasetSiteMatrixConfig,
    DifferentialAnalysisRequest,
    ExperimentalDesign,
    KinasePredictionConfig,
    KinaseScoringConfig,
    KinaseWorkflowRequest,
    Organism,
    ReferencePreset,
    SampleDesignRecord,
    SignalomeWorkflowRequest,
)
from phospy.errors import PhosPyInputError, WorkflowValidationError
from tests.support.intensity_scale_states import (
    supported_log2_intensity_scale_state,
    supported_log2_processing_state,
)
from tests.support.rewrite_fixture_data import build_rat_l6_dataset
from tests.support.signalome_config import build_signalome_config


def test_builder_rejects_duplicate_display_ids_with_conflicting_protein_identity() -> (
    None
):
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 2.0],
            "sample_b": [1.5, 2.5],
        },
        index=pd.Index(["row_a", "row_b"], name="source_row"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "MAPK14"],
            "site": ["Y182", "Y182"],
            "site_sequence": [
                ("A" * 15) + str(site).strip().upper()[0] + ("A" * 15)
                for site in ["Y182", "Y182"]
            ],
            "localisation_confidence": [0.95, 0.95],
            "protein_id": ["P28482", "P28482"],
            "protein_accession": ["P28482-1", "P28482-2"],
        },
        index=phospho.index.copy(),
    )

    with pytest.raises(
        PhosPyInputError,
        match="conflicting scientific identities",
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=site_metadata,
                organism=Organism.RAT,
                preprocessing_config=DatasetPreprocessingConfig(
                    site_matrix=DatasetSiteMatrixConfig(
                        policy="build_from_metadata",
                        duplicate_site_policy="first",
                    )
                ),
            )
        )


def test_builder_rejects_duplicate_display_ids_with_conflicting_protein_id() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 2.0],
            "sample_b": [1.5, 2.5],
        },
        index=pd.Index(["row_a", "row_b"], name="source_row"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "MAPK14"],
            "site": ["Y182", "Y182"],
            "site_sequence": [
                ("A" * 15) + str(site).strip().upper()[0] + ("A" * 15)
                for site in ["Y182", "Y182"]
            ],
            "localisation_confidence": [0.95, 0.95],
            "protein_id": ["P28482", "Q5S007"],
            "protein_accession": ["P28482-1", "P28482-1"],
        },
        index=phospho.index.copy(),
    )

    with pytest.raises(
        PhosPyInputError,
        match="protein_id=\\['P28482', 'Q5S007'\\]",
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=site_metadata,
                organism=Organism.RAT,
                preprocessing_config=DatasetPreprocessingConfig(
                    site_matrix=DatasetSiteMatrixConfig(
                        policy="build_from_metadata",
                        duplicate_site_policy="first",
                    )
                ),
            )
        )


def test_signalome_still_requires_explicit_protein_identity() -> None:
    base_dataset = build_rat_l6_dataset(n_sites=260)
    dataset_without_protein = AnalysisReadyPhosphoDataset(
        phospho=base_dataset.phospho,
        site_metadata=base_dataset.site_metadata.drop(columns=["protein_id"]),
        sample_metadata=base_dataset.sample_metadata,
        total=base_dataset.total,
        organism=base_dataset.organism,
        intensity_scale_state=base_dataset.intensity_scale_state,
        processing_state=base_dataset.processing_state,
    )
    kinase_result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset_without_protein,
            references=ReferencePreset.AUTO,
            scoring_config=KinaseScoringConfig(min_substrates=2),
            prediction_config=KinasePredictionConfig(
                top_k=6,
                deterministic_max_selected_kinases=12,
                adaptive_ensemble_runs=12,
            ),
            activity_config=None,
        )
    )

    with pytest.raises(
        WorkflowValidationError,
        match="site_metadata is missing required columns: protein_id",
    ):
        SignalomeWorkflow().run(
            SignalomeWorkflowRequest(
                kinase_result=kinase_result,
                config=build_signalome_config(substrate_support_cutoff=0.5),
            )
        )


def test_differential_workflow_accepts_gene_site_only_dataset() -> None:
    dataset = AnalysisReadyPhosphoDataset(
        phospho=pd.DataFrame(
            {
                "A_1": [1.0, 2.0],
                "A_2": [1.1, 2.1],
                "B_1": [2.1, 2.0],
                "B_2": [2.0, 2.2],
            },
            index=pd.Index(["MAPK14;Y182;", "AKT1;T308;"], name="site_id"),
        ),
        site_metadata=pd.DataFrame(
            {
                "gene_symbol": ["MAPK14", "AKT1"],
                "site": ["Y182", "T308"],
                "site_sequence": [
                    ("A" * 15) + str(site).strip().upper()[0] + ("A" * 15)
                    for site in ["Y182", "T308"]
                ],
                "localisation_confidence": [0.95, 0.95],
            },
            index=pd.Index(["MAPK14;Y182;", "AKT1;T308;"], name="site_id"),
        ),
        organism=Organism.RAT,
        intensity_scale_state=supported_log2_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_log2_processing_state(has_total_matrix=False),
    )
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(sample_id="A_1", condition="A"),
            SampleDesignRecord(sample_id="A_2", condition="A"),
            SampleDesignRecord(sample_id="B_1", condition="B"),
            SampleDesignRecord(sample_id="B_2", condition="B"),
        )
    )
    result = DifferentialAnalysisWorkflow().run(
        DifferentialAnalysisRequest(
            dataset=dataset,
            design=design,
            contrasts=(
                Contrast(
                    name="B_vs_A", numerator_condition="B", denominator_condition="A"
                ),
            ),
        )
    )

    assert list(result.table_for("B_vs_A").index) == ["MAPK14;Y182;", "AKT1;T308;"]
