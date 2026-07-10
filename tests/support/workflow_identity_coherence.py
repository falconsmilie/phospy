from __future__ import annotations

import pandas as pd

from phospy import AnalysisReadyPhosphoDataset
from phospy.api import (
    Contrast,
    DifferentialAnalysisRequest,
    ExperimentalDesign,
    KinaseWorkflowRequest,
    Organism,
    ReferenceBundle,
    SampleDesignRecord,
    SignalomeWorkflowRequest,
)
from phospy.api.configs import (
    KINASE_REFERENCE_DISPLAY_AMBIGUITY_POLICY_ALLOW_WITH_DIAGNOSTICS,
    KinasePredictionConfig,
    KinaseReferenceDisplayAmbiguityPolicy,
    KinaseScoringConfig,
    ReferenceContextCompatibilityPolicy,
)
from phospy.api.results import (
    KinasePredictionResult,
    KinaseScoringResult,
    KinaseWorkflowResult,
)
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
    supported_log2_intensity_scale_state,
    supported_log2_processing_state,
)
from tests.support.signalome_config import build_signalome_config
from tests.support.site_keys import protein_site_key_index, site_key_context_columns
from tests.support.unsafe_dataset_states import (
    unsafe_corrupt_dataset_to_display_index,
    unsafe_drop_dataset_site_metadata_columns,
    unsafe_set_dataset_site_metadata_columns,
)

DUPLICATE_DISPLAY_ID = "MAPK14;Y182;"


def duplicate_display_site_index() -> pd.Index:
    return protein_site_key_index(
        protein_identifiers=["P28482", "Q99999"],
        sites=["Y182", "Y182"],
    )


def duplicate_display_site_metadata(site_index: pd.Index) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "site_key": site_index.astype(str).tolist(),
            "display_id": [DUPLICATE_DISPLAY_ID, DUPLICATE_DISPLAY_ID],
            **site_key_context_columns(site_index),
            "gene_symbol": ["MAPK14", "MAPK14"],
            "site": ["Y182", "Y182"],
            "site_sequence": [("A" * 15) + "Y" + ("A" * 15)] * 2,
            "protein_id": ["P28482", "Q99999"],
            "protein_accession": ["P28482-1", "Q99999-1"],
            "localisation_confidence": [0.95, 0.95],
        },
        index=site_index.copy(),
    )


def build_duplicate_display_differential_dataset() -> AnalysisReadyPhosphoDataset:
    site_index = duplicate_display_site_index()
    return AnalysisReadyPhosphoDataset(
        phospho=pd.DataFrame(
            {
                "A_1": [1.0, 2.0],
                "A_2": [1.1, 2.1],
                "B_1": [2.1, 2.0],
                "B_2": [2.0, 2.2],
            },
            index=site_index.copy(),
        ),
        site_metadata=duplicate_display_site_metadata(site_index),
        organism=Organism.RAT,
        intensity_scale_state=supported_log2_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_log2_processing_state(has_total_matrix=False),
    )


def build_duplicate_display_kinase_dataset() -> AnalysisReadyPhosphoDataset:
    site_index = duplicate_display_site_index()
    return AnalysisReadyPhosphoDataset(
        phospho=pd.DataFrame(
            {"sample_a": [1.0, 2.0], "sample_b": [1.1, 2.1]},
            index=site_index.copy(),
        ),
        site_metadata=duplicate_display_site_metadata(site_index),
        organism=Organism.RAT,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
    )


def build_duplicate_display_reference_bundle() -> ReferenceBundle:
    return ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["K1"],
                "substrate_site": [DUPLICATE_DISPLAY_ID],
            }
        ),
        site_sequences=pd.DataFrame(
            {"site_sequence": [("A" * 15) + "Y" + ("A" * 15)]},
            index=pd.Index([DUPLICATE_DISPLAY_ID], name="site_id"),
        ),
    )


def build_duplicate_display_differential_request(
    *,
    dataset: AnalysisReadyPhosphoDataset | None = None,
) -> DifferentialAnalysisRequest:
    resolved = (
        build_duplicate_display_differential_dataset() if dataset is None else dataset
    )
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


def build_duplicate_display_kinase_request(
    *,
    dataset: AnalysisReadyPhosphoDataset | None = None,
    reference_display_ambiguity_policy: KinaseReferenceDisplayAmbiguityPolicy = (
        KINASE_REFERENCE_DISPLAY_AMBIGUITY_POLICY_ALLOW_WITH_DIAGNOSTICS
    ),
) -> KinaseWorkflowRequest:
    return KinaseWorkflowRequest(
        dataset=(
            build_duplicate_display_kinase_dataset() if dataset is None else dataset
        ),
        references=build_duplicate_display_reference_bundle(),
        scoring_config=KinaseScoringConfig(
            min_substrates=2,
            reference_context_compatibility_policy=(
                ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT
            ),
        ),
        prediction_config=KinasePredictionConfig(
            top_k=2,
            deterministic_max_selected_kinases=2,
            adaptive_ensemble_runs=2,
        ),
        activity_config=None,
        reference_display_ambiguity_policy=reference_display_ambiguity_policy,
    )


def build_duplicate_display_signalome_request(
    *,
    dataset: AnalysisReadyPhosphoDataset | None = None,
) -> SignalomeWorkflowRequest:
    resolved = build_duplicate_display_kinase_dataset() if dataset is None else dataset
    site_index = resolved.phospho.index.copy()
    score_matrix = pd.DataFrame(
        {
            "K1": [1.0, 0.8],
            "K2": [0.2, 0.9],
        },
        index=site_index.copy(),
    )
    prediction_matrix = pd.DataFrame(
        {
            "K1": [0.8, 0.7],
            "K2": [0.7, 0.8],
        },
        index=site_index.copy(),
    )
    return SignalomeWorkflowRequest(
        kinase_result=KinaseWorkflowResult(
            dataset=resolved,
            references=build_duplicate_display_reference_bundle(),
            scoring_result=KinaseScoringResult(
                profile_scores=score_matrix,
                rank_weighted_fusion_scores=score_matrix,
            ),
            prediction_result=KinasePredictionResult(pred_mat=prediction_matrix),
            activity_result=None,
        ),
        config=build_signalome_config(
            substrate_support_cutoff=0.5,
            network_correlation_threshold=0.1,
            module_count=2,
        ),
    )


def corrupt_dataset_to_display_index(dataset: AnalysisReadyPhosphoDataset) -> None:
    unsafe_corrupt_dataset_to_display_index(dataset)


def drop_site_metadata_column(
    dataset: AnalysisReadyPhosphoDataset,
    column_name: str,
) -> None:
    unsafe_drop_dataset_site_metadata_columns(dataset, column_name)


def set_site_sequence_values(
    dataset: AnalysisReadyPhosphoDataset,
    site_sequences: list[str],
) -> None:
    unsafe_set_dataset_site_metadata_columns(
        dataset,
        {"site_sequence": site_sequences},
    )
