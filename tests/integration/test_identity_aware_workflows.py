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
from phospy.api.results import (
    KinasePredictionResult,
    KinaseScoringResult,
    KinaseWorkflowResult,
)
from phospy.errors import WorkflowValidationError
from phospy.science.references.models import ReferenceBundle
from phospy.workflows.differential.validator import DifferentialAnalysisValidator
from phospy.workflows.kinase.validator import KinaseWorkflowValidator
from phospy.workflows.signalome.validator import SignalomeWorkflowValidator
from tests.support.rewrite_fixture_data import build_rat_l6_dataset
from tests.support.signalome_config import build_signalome_config
from tests.support.unsafe_dataset_states import (
    unsafe_drop_dataset_site_metadata_columns,
    unsafe_set_dataset_site_metadata_columns,
    unsafe_set_dataset_site_metadata_index,
)


def _minimal_site_key_dataset(
    *, include_protein_id: bool
) -> AnalysisReadyPhosphoDataset:
    site_metadata = {
        "gene_symbol": ["MAPK14", "AKT1"],
        "site": ["Y182", "T308"],
        "site_sequence": [
            "AAAAAAAAAAAAAAAYAAAAAAAAAAAAAAA",
            "AAAAAAAAAAAAAAATAAAAAAAAAAAAAAA",
        ],
        "localisation_confidence": [0.95, 0.95],
    }
    site_metadata["protein_id"] = ["P28482", "P31749"]
    dataset = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=pd.DataFrame(
                {
                    "A_1": [1.0, 2.0],
                    "A_2": [1.1, 2.1],
                    "B_1": [2.1, 2.0],
                    "B_2": [2.0, 2.2],
                },
                index=pd.Index(["row_a", "row_b"], name="source_row"),
            ),
            site_metadata=pd.DataFrame(
                site_metadata,
                index=pd.Index(["row_a", "row_b"], name="source_row"),
            ),
            organism=Organism.RAT,
            input_intensity_scale="log2",
        )
    )
    if include_protein_id:
        return dataset
    without_protein = dataset.site_metadata.drop(columns=["protein_id"])
    return AnalysisReadyPhosphoDataset(
        phospho=dataset.phospho,
        site_metadata=without_protein,
        sample_metadata=dataset.sample_metadata,
        total=dataset.total,
        organism=dataset.organism,
        intensity_scale_state=dataset.intensity_scale_state,
        processing_state=dataset.processing_state,
        preprocessing_report=dataset.preprocessing_report,
        provenance=dataset.provenance,
        allow_opaque_site_values=dataset.allow_opaque_site_values,
    )


def _build_differential_request(
    *, dataset: AnalysisReadyPhosphoDataset
) -> DifferentialAnalysisRequest:
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(sample_id="A_1", condition="A"),
            SampleDesignRecord(sample_id="A_2", condition="A"),
            SampleDesignRecord(sample_id="B_1", condition="B"),
            SampleDesignRecord(sample_id="B_2", condition="B"),
        )
    )
    return DifferentialAnalysisRequest(
        dataset=dataset,
        design=design,
        contrasts=(
            Contrast(name="B_vs_A", numerator_condition="B", denominator_condition="A"),
        ),
    )


def _build_kinase_request(
    *, dataset: AnalysisReadyPhosphoDataset
) -> KinaseWorkflowRequest:
    return KinaseWorkflowRequest(
        dataset=dataset,
        references=ReferencePreset.AUTO,
        scoring_config=KinaseScoringConfig(min_substrates=2),
        prediction_config=KinasePredictionConfig(
            top_k=6,
            deterministic_max_selected_kinases=12,
            adaptive_ensemble_runs=12,
        ),
        activity_config=None,
    )


def _signalome_reference_bundle_for(
    dataset: AnalysisReadyPhosphoDataset,
) -> ReferenceBundle:
    metadata = dataset.site_metadata
    display_ids = metadata.loc[:, "display_id"].astype(str).tolist()
    return ReferenceBundle(
        organism=dataset.organism or Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["MAP2K6"] * len(display_ids),
                "substrate_site": display_ids,
            }
        ),
        site_sequences=pd.DataFrame(
            {"site_sequence": metadata.loc[:, "site_sequence"].astype(str).tolist()},
            index=pd.Index(display_ids, name="site_id"),
        ),
    )


def _build_signalome_request(
    *, dataset: AnalysisReadyPhosphoDataset
) -> SignalomeWorkflowRequest:
    site_index = dataset.site_metadata.index.copy()
    score_matrix = pd.DataFrame({"MAP2K6": [1.0] * len(site_index)}, index=site_index)
    prediction_matrix = pd.DataFrame(
        {"MAP2K6": [0.8] * len(site_index)},
        index=site_index,
    )
    return SignalomeWorkflowRequest(
        kinase_result=KinaseWorkflowResult(
            dataset=dataset,
            references=_signalome_reference_bundle_for(dataset),
            scoring_result=KinaseScoringResult(
                profile_scores=score_matrix,
                rank_weighted_fusion_scores=score_matrix,
            ),
            prediction_result=KinasePredictionResult(pred_mat=prediction_matrix),
            activity_result=None,
        ),
        config=build_signalome_config(substrate_support_cutoff=0.5),
    )


def _site_metadata_dataset(request: object) -> AnalysisReadyPhosphoDataset:
    if isinstance(request, SignalomeWorkflowRequest):
        return request.kinase_result.dataset
    if isinstance(request, DifferentialAnalysisRequest):
        return request.dataset
    if isinstance(request, KinaseWorkflowRequest):
        return request.dataset
    raise TypeError(f"unexpected request type: {type(request)!r}")


def test_builder_allows_duplicate_display_ids_across_distinct_protein_accessions() -> (
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

    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
            input_intensity_scale="linear",
            preprocessing_config=DatasetPreprocessingConfig(
                site_matrix=DatasetSiteMatrixConfig(
                    policy="build_from_metadata",
                    duplicate_site_policy="first",
                )
            ),
        )
    )
    assert built.phospho.shape[0] == 2


def test_builder_allows_duplicate_display_ids_across_distinct_protein_ids() -> None:
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
        },
        index=phospho.index.copy(),
    )

    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
            input_intensity_scale="linear",
            preprocessing_config=DatasetPreprocessingConfig(
                site_matrix=DatasetSiteMatrixConfig(
                    policy="build_from_metadata",
                    duplicate_site_policy="first",
                )
            ),
        )
    )
    assert built.phospho.shape[0] == 2


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
    dataset = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=pd.DataFrame(
                {
                    "A_1": [1.0, 2.0],
                    "A_2": [1.1, 2.1],
                    "B_1": [2.1, 2.0],
                    "B_2": [2.0, 2.2],
                },
                index=pd.Index(["row_a", "row_b"], name="source_row"),
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
                    "protein_id": ["P28482", "P31749"],
                },
                index=pd.Index(["row_a", "row_b"], name="source_row"),
            ),
            organism=Organism.RAT,
            input_intensity_scale="log2",
        )
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

    assert list(result.table_for("B_vs_A").index) == dataset.phospho.index.tolist()


def test_workflow_validators_require_display_id_column() -> None:
    cases: tuple[tuple[object, object], ...] = (
        (
            DifferentialAnalysisValidator(),
            _build_differential_request(
                dataset=_minimal_site_key_dataset(include_protein_id=False)
            ),
        ),
        (
            KinaseWorkflowValidator(),
            _build_kinase_request(
                dataset=_minimal_site_key_dataset(include_protein_id=True)
            ),
        ),
        (
            SignalomeWorkflowValidator(),
            _build_signalome_request(
                dataset=_minimal_site_key_dataset(include_protein_id=True)
            ),
        ),
    )
    for validator, request in cases:
        dataset = _site_metadata_dataset(request)
        unsafe_drop_dataset_site_metadata_columns(dataset, "display_id")
        with pytest.raises(
            WorkflowValidationError,
            match="missing required columns: display_id",
        ):
            validator.run(request)


def test_workflow_validators_reject_non_site_key_indexed_site_metadata() -> None:
    cases: tuple[tuple[object, object], ...] = (
        (
            DifferentialAnalysisValidator(),
            _build_differential_request(
                dataset=_minimal_site_key_dataset(include_protein_id=False)
            ),
        ),
        (
            KinaseWorkflowValidator(),
            _build_kinase_request(
                dataset=_minimal_site_key_dataset(include_protein_id=True)
            ),
        ),
        (
            SignalomeWorkflowValidator(),
            _build_signalome_request(
                dataset=_minimal_site_key_dataset(include_protein_id=True)
            ),
        ),
    )
    for validator, request in cases:
        dataset = _site_metadata_dataset(request)
        site_metadata = dataset._borrow_site_metadata_frame().copy(deep=True)
        unsafe_set_dataset_site_metadata_index(
            dataset,
            pd.Index(
                site_metadata.loc[:, "display_id"].astype(str).tolist(),
                name=site_metadata.index.name,
            ),
        )
        with pytest.raises(
            WorkflowValidationError,
            match=(
                "display-indexed direct construction|"
                "index must match .*site_key|"
                "must exactly match .*dataset\\.phospho\\.index"
            ),
        ):
            validator.run(request)


def test_workflow_validators_require_site_key_column_to_match_index() -> None:
    cases: tuple[tuple[object, object], ...] = (
        (
            DifferentialAnalysisValidator(),
            _build_differential_request(
                dataset=_minimal_site_key_dataset(include_protein_id=False)
            ),
        ),
        (
            KinaseWorkflowValidator(),
            _build_kinase_request(
                dataset=_minimal_site_key_dataset(include_protein_id=True)
            ),
        ),
        (
            SignalomeWorkflowValidator(),
            _build_signalome_request(
                dataset=_minimal_site_key_dataset(include_protein_id=True)
            ),
        ),
    )
    for validator, request in cases:
        dataset = _site_metadata_dataset(request)
        site_metadata = dataset._borrow_site_metadata_frame().copy(deep=True)
        mutated_site_keys = site_metadata.loc[:, "site_key"].tolist()
        mutated_site_keys[0] = mutated_site_keys[-1]
        unsafe_set_dataset_site_metadata_columns(
            dataset,
            {"site_key": mutated_site_keys},
        )
        with pytest.raises(
            WorkflowValidationError,
            match="site_key must exactly match .*index values",
        ):
            validator.run(request)
