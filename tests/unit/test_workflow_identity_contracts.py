from __future__ import annotations

import pandas as pd
import pytest

from phospy import AnalysisReadyDatasetBuilder, AnalysisReadyPhosphoDataset
from phospy.api import (
    Contrast,
    DatasetBuildRequest,
    DatasetPreprocessingConfig,
    DatasetSiteMatrixConfig,
    DifferentialAnalysisRequest,
    ExperimentalDesign,
    KinaseWorkflow,
    KinaseWorkflowRequest,
    Organism,
    ReferencePreset,
    SampleDesignRecord,
    SignalomeWorkflowRequest,
)
from phospy.api.configs import KinasePredictionConfig, KinaseScoringConfig
from phospy.api.results import (
    KinasePredictionResult,
    KinaseScoringResult,
    KinaseWorkflowResult,
)
from phospy.errors import (
    PhosPyInputError,
    ReferenceCompatibilityError,
    WorkflowValidationError,
)
from phospy.science.references.models import ReferenceBundle
from phospy.workflows.differential.validator import DifferentialAnalysisValidator
from phospy.workflows.kinase.validator import KinaseWorkflowValidator
from phospy.workflows.signalome.validator import SignalomeWorkflowValidator
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
    supported_log2_intensity_scale_state,
    supported_log2_processing_state,
)
from tests.support.signalome_config import build_signalome_config
from tests.support.site_keys import site_key_index_from_display_ids


def _site_keys(display_ids: list[str]) -> pd.Index:
    return site_key_index_from_display_ids(
        display_ids,
        protein_namespace="gene_symbol",
    )


def _differential_dataset(
    *, allow_opaque_site_values: bool = False
) -> AnalysisReadyPhosphoDataset:
    display_ids = ["MAPK14;Y182;", "AKT1;T308;"]
    site_ids = _site_keys(display_ids)
    metadata_sites = ["Y182", "T308"]
    sequence_values = [
        ("A" * 15) + str(site).strip().upper()[0] + ("A" * 15)
        for site in metadata_sites
    ]
    if allow_opaque_site_values:
        display_ids = ["MAPK14;FOO;", "AKT1;BAR;"]
        site_ids = _site_keys(["MAPK14;Y182;", "AKT1;T308;"])
        metadata_sites = ["FOO", "BAR"]
        sequence_values = [("A" * 15) + "Y" + ("A" * 15)] * 2
    return AnalysisReadyPhosphoDataset(
        phospho=pd.DataFrame(
            {
                "A_1": [1.0, 2.0],
                "A_2": [1.1, 2.1],
                "B_1": [2.1, 2.0],
                "B_2": [2.0, 2.2],
            },
            index=site_ids.copy(),
        ),
        site_metadata=pd.DataFrame(
            {
                "site_key": site_ids.tolist(),
                "display_id": display_ids,
                "gene_symbol": ["MAPK14", "AKT1"],
                "site": metadata_sites,
                "site_sequence": sequence_values,
            },
            index=site_ids.copy(),
        ),
        organism=Organism.RAT,
        intensity_scale_state=supported_log2_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_log2_processing_state(has_total_matrix=False),
        allow_opaque_site_values=allow_opaque_site_values,
    )


def _differential_request(
    *, dataset: AnalysisReadyPhosphoDataset | None = None
) -> DifferentialAnalysisRequest:
    resolved = _differential_dataset() if dataset is None else dataset
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(sample_id="A_1", condition="A"),
            SampleDesignRecord(sample_id="A_2", condition="A"),
            SampleDesignRecord(sample_id="B_1", condition="B"),
            SampleDesignRecord(sample_id="B_2", condition="B"),
        )
    )
    return DifferentialAnalysisRequest(
        dataset=resolved,
        design=design,
        contrasts=(
            Contrast(name="B_vs_A", numerator_condition="B", denominator_condition="A"),
        ),
    )


def _kinase_dataset() -> AnalysisReadyPhosphoDataset:
    display_ids = ["MAPK14;Y182;", "AKT1;T308;"]
    site_ids = _site_keys(display_ids)
    return AnalysisReadyPhosphoDataset(
        phospho=pd.DataFrame(
            {"sample_a": [1.0, 2.0], "sample_b": [1.1, 2.1]},
            index=site_ids.copy(),
        ),
        site_metadata=pd.DataFrame(
            {
                "site_key": site_ids.tolist(),
                "display_id": display_ids,
                "gene_symbol": ["MAPK14", "AKT1"],
                "site": ["Y182", "T308"],
                "site_sequence": ["AAAAAAAYAAAAAAA", "AAAAAAATAAAAAAA"],
                "protein_id": ["P28482", "P31749"],
                "protein_accession": ["P28482-1", "P31749-1"],
            },
            index=site_ids.copy(),
        ),
        organism=Organism.RAT,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
    )


def _references() -> ReferenceBundle:
    return ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["MAP2K6", "AKT1"],
                "substrate_site": ["MAPK14;Y182;", "AKT1;T308;"],
            }
        ),
        site_sequences=pd.DataFrame(
            {"site_sequence": ["AAAAAAAYAAAAAAA", "AAAAAAATAAAAAAA"]},
            index=pd.Index(["MAPK14;Y182;", "AKT1;T308;"], name="site_id"),
        ),
    )


def _signalome_request(
    *, dataset: AnalysisReadyPhosphoDataset
) -> SignalomeWorkflowRequest:
    site_ids = dataset.phospho.index.copy()
    score_matrix = pd.DataFrame({"MAP2K6": [1.0] * len(site_ids)}, index=site_ids)
    prediction_matrix = pd.DataFrame({"MAP2K6": [0.8] * len(site_ids)}, index=site_ids)
    return SignalomeWorkflowRequest(
        kinase_result=KinaseWorkflowResult(
            dataset=dataset,
            references=_references(),
            scoring_result=KinaseScoringResult(
                profile_scores=score_matrix,
                rank_weighted_fusion_scores=score_matrix,
            ),
            prediction_result=KinasePredictionResult(pred_mat=prediction_matrix),
            activity_result=None,
        ),
        config=build_signalome_config(substrate_support_cutoff=0.5),
    )


def test_differential_identity_contract_accepts_display_level_minimum() -> None:
    validated = DifferentialAnalysisValidator().run(
        _differential_request(dataset=_differential_dataset())
    )
    assert validated.dataset._borrow_site_metadata_frame().loc[
        :, "display_id"
    ].tolist() == [
        "MAPK14;Y182;",
        "AKT1;T308;",
    ]


def test_differential_identity_contract_allows_opaque_sites_with_explicit_opt_in() -> (
    None
):
    dataset = _differential_dataset(allow_opaque_site_values=True)
    validated = DifferentialAnalysisValidator().run(
        _differential_request(dataset=dataset)
    )
    assert validated.dataset.opaque_site_values_allowed is True


def test_kinase_identity_contract_rejects_conflicting_duplicate_display_ids() -> None:
    dataset = _kinase_dataset()
    dataset._site_metadata.loc[:, "display_id"] = ["MAPK14;Y182;", "MAPK14;Y182;"]
    dataset._site_metadata.loc[:, "gene_symbol"] = ["MAPK14", "MAPK14"]
    dataset._site_metadata.loc[:, "site"] = ["Y182", "Y182"]
    dataset._site_metadata.loc[:, "protein_id"] = ["P28482", "P28482"]
    dataset._site_metadata.loc[:, "protein_accession"] = ["P28482-1", "P28482-2"]
    request = KinaseWorkflowRequest(
        dataset=dataset,
        references=_references(),
        scoring_config=KinaseScoringConfig(min_substrates=2),
        prediction_config=KinasePredictionConfig(
            top_k=5,
            deterministic_max_selected_kinases=5,
            adaptive_ensemble_runs=5,
        ),
        activity_config=None,
    )

    with pytest.raises(
        WorkflowValidationError,
        match=(
            "identity requirement failed "
            "\\(contract=sty_site_identity_plus_sequence_context\\)"
        ),
    ) as exc_info:
        KinaseWorkflowValidator().run(request)
    assert "conflicting scientific identities for duplicate display site IDs" in str(
        exc_info.value
    )


def test_signalome_identity_contract_rejects_missing_protein_identity() -> None:
    dataset = _kinase_dataset()
    object.__setattr__(
        dataset,
        "_site_metadata",
        dataset._site_metadata.drop(columns=["protein_id"]),
    )
    request = _signalome_request(dataset=dataset)

    with pytest.raises(
        WorkflowValidationError,
        match=(
            "identity requirement failed \\(contract=protein_scoped_site_identity\\)"
        ),
    ):
        SignalomeWorkflowValidator().run(request)


def test_signalome_inputs_require_unique_display_ids_before_workflow_validation() -> (
    None
):
    phospho = pd.DataFrame(
        {"sample_a": [1.0, 2.0], "sample_b": [1.1, 2.1]},
        index=pd.Index(["row_a", "row_b"], name="source_row"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "MAPK14"],
            "site": ["Y182", "Y182"],
            "site_sequence": ["AAAAAAAYAAAAAAA", "AAAAAAAYAAAAAAA"],
            "protein_id": ["P28482", "P28482"],
            "localisation_confidence": [0.95, 0.95],
        },
        index=phospho.index.copy(),
    )

    with pytest.raises(
        PhosPyInputError,
        match="one analysis-ready row per normalised display-site identifier",
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=site_metadata,
                organism=Organism.RAT,
                preprocessing_config=DatasetPreprocessingConfig(
                    site_matrix=DatasetSiteMatrixConfig(policy="build_from_metadata")
                ),
                input_intensity_scale="linear",
            )
        )


def test_kinase_workflow_rejects_reference_organism_mismatch_where_applicable() -> None:
    request = KinaseWorkflowRequest(
        dataset=_kinase_dataset(),
        references=ReferencePreset.HUMAN,
        scoring_config=KinaseScoringConfig(min_substrates=2),
        prediction_config=KinasePredictionConfig(
            top_k=5,
            deterministic_max_selected_kinases=5,
            adaptive_ensemble_runs=5,
        ),
        activity_config=None,
    )
    with pytest.raises(ReferenceCompatibilityError):
        KinaseWorkflow().run(request)
