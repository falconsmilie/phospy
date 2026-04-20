from __future__ import annotations

import numpy as np
import pandas as pd

import phospy.api.results as result_models
from phospy import (
    AnalysisReadyPhosphoDataset,
    KinaseActivityConfig,
    KinasePredictionConfig,
    KinaseScoringConfig,
    KinaseWorkflow,
    KinaseWorkflowRequest,
    Organism,
    ReferenceBundle,
    SignalomeConfig,
    SignalomeWorkflow,
    SignalomeWorkflowRequest,
)
from phospy.api.results import (
    KinasePredictionResult,
    KinaseScoringResult,
    KinaseWorkflowResult,
    SignalomeWorkflowResult,
)
from phospy.transformations.models import TransformationState


def _dataset() -> AnalysisReadyPhosphoDataset:
    site_ids = ["MAPK14;Y182;", "GSK3B;S9;", "AKT1;T308;"]
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 2.0, 4.0],
            "sample_b": [2.0, 4.0, 1.0],
        },
        index=site_ids,
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "GSK3B", "AKT1"],
            "site": ["Y182", "S9", "T308"],
            "site_sequence": ["A" * 31, "B" * 31, "C" * 31],
            "protein_id": ["MAPK14", "GSK3B", "AKT1"],
        },
        index=site_ids,
    )
    return AnalysisReadyPhosphoDataset(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        transformation_state=TransformationState.raw(has_total_matrix=False),
    )


def _references() -> ReferenceBundle:
    site_ids = pd.Index(["MAPK14;Y182;", "GSK3B;S9;", "AKT1;T308;"], name="site_id")
    return ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["MAP2K6", "MAP2K6", "AKT1", "AKT1"],
                "substrate_site": [
                    "MAPK14;Y182;",
                    "GSK3B;S9;",
                    "GSK3B;S9;",
                    "AKT1;T308;",
                ],
            }
        ),
        site_sequences=pd.DataFrame(
            {"site_sequence": ["A" * 31, "B" * 31, "C" * 31]},
            index=site_ids,
        ),
    )


def _kinase_result() -> KinaseWorkflowResult:
    return KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=_dataset(),
            references=_references(),
            scoring_config=KinaseScoringConfig(min_substrates=2),
            prediction_config=KinasePredictionConfig(top_k=1, ensemble_size=2),
            activity_config=None,
        )
    )


def test_public_result_exports_match_contract() -> None:
    assert set(result_models.__all__) == {
        "KinaseActivityResult",
        "KinasePredictionResult",
        "KinaseScoringResult",
        "KinaseWorkflowResult",
        "SignalomeWorkflowResult",
    }


def test_kinase_result_stays_nested_and_honest_for_supported_lane() -> None:
    result = _kinase_result()
    assert isinstance(result, KinaseWorkflowResult)
    assert isinstance(result.scoring_result, KinaseScoringResult)
    assert isinstance(result.prediction_result, KinasePredictionResult)
    assert not result.scoring_result.profile_scores.empty
    assert result.scoring_result.motif_scores is None
    assert result.scoring_result.combined_scores is not None
    assert result.scoring_result.weights is None
    assert result.activity_result is None
    assert result.prediction_result.substrate_list is not None
    assert set(result.prediction_result.substrate_list.columns) == {
        "kinase",
        "substrate_site",
        "score",
        "rank",
    }
    pred_values = result.prediction_result.pred_mat.to_numpy(dtype=float)
    finite_values = pred_values[np.isfinite(pred_values)]
    assert (finite_values >= 0.0).all()
    assert (finite_values <= 1.0).all()
    assert not hasattr(result, "profile_scores")
    assert not hasattr(result, "pred_mat")
    assert not hasattr(result, "substrate_list")


def test_signalome_result_keeps_nested_kinase_result_contract() -> None:
    kinase_result = _kinase_result()
    signalome_result = SignalomeWorkflow().run(
        SignalomeWorkflowRequest(
            kinase_result=kinase_result,
            config=SignalomeConfig(substrate_support_cutoff=0.5),
        )
    )
    assert isinstance(signalome_result, SignalomeWorkflowResult)
    assert signalome_result.kinase_result is kinase_result
    assert not signalome_result.module_assignments.table.empty
    assert not signalome_result.signalome_modules.table.empty
    assert signalome_result.module_selection_diagnostics.selected_module_count >= 1
    assert signalome_result.module_selection_diagnostics.reason
    assert signalome_result.score_preconditioning_diagnostics.input_row_count >= 1
    assert (
        signalome_result.score_preconditioning_diagnostics.dropped_all_missing_row_count
        >= 0
    )
    assert signalome_result.expanded_signalome is not None
    assert not signalome_result.expanded_signalome.empty
    assert signalome_result.kinase_result.scoring_result.motif_scores is None
    assert signalome_result.kinase_result.scoring_result.combined_scores is not None
    assert signalome_result.kinase_result.scoring_result.weights is None
    assert not hasattr(signalome_result, "pred_mat")
    assert not hasattr(signalome_result, "profile_scores")


def test_kinase_result_exposes_supported_activity_stage_outputs_when_enabled() -> None:
    result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=_dataset(),
            references=_references(),
            scoring_config=KinaseScoringConfig(min_substrates=2),
            prediction_config=KinasePredictionConfig(top_k=1, ensemble_size=2),
            activity_config=KinaseActivityConfig(
                enabled=True,
                threshold=0.0,
                min_substrates=1,
                top_n_substrates=2,
            ),
        )
    )

    assert result.activity_result is not None
    assert set(result.activity_result.target_table.columns) == {
        "site_id",
        "kinase",
        "score",
    }
    assert result.activity_result.ksea_counts.name == "n_substrates"
    assert result.activity_result.target_counts.name == "n_targets"


def test_kinase_result_can_include_opt_in_diagnostic_scoring_tables() -> None:
    result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=_dataset(),
            references=_references(),
            scoring_config=KinaseScoringConfig(
                min_substrates=2,
                include_diagnostic_scoring_tables=True,
            ),
            prediction_config=KinasePredictionConfig(top_k=1, ensemble_size=2),
            activity_config=None,
        )
    )

    assert result.scoring_result.motif_scores is not None
    assert result.scoring_result.weights is not None
