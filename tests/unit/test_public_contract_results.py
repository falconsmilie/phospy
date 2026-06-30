from __future__ import annotations

import numpy as np
import pandas as pd

import phospy.api.results as result_models
from phospy import (
    AnalysisReadyPhosphoDataset,
    KinaseWorkflow,
    SignalomeWorkflow,
)
from phospy.api import (
    KinaseActivityConfig,
    KinasePredictionConfig,
    KinaseScoringConfig,
    KinaseWorkflowRequest,
    Organism,
    ReferenceBundle,
    SignalomeWorkflowRequest,
)
from phospy.api.results import (
    KinasePredictionResult,
    KinaseScoringResult,
    KinaseWorkflowResult,
    SignalomeWorkflowResult,
)
from phospy.provenance.scientific_policy_models import ScientificPolicyId
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
)
from tests.support.signalome_config import build_signalome_config
from tests.support.site_keys import (
    site_key_context_columns,
    site_key_index_from_display_ids,
)

INTENTIONAL_RESULT_COMPATIBILITY_ALIASES = {
    "KinaseEligibilityReport",
    "KinaseWorkflowPreprocessingAttritionSummary",
    "KinaseWorkflowScoringAttritionSummary",
    "KinaseWorkflowSiteAttritionSummary",
}

EXPECTED_RESULT_EXPORTS = {
    "BatchCorrectionDiagnostics",
    "BatchCorrectionPolicy",
    "BatchCorrectionReport",
    "DifferentialAnalysisResult",
    "DifferentialContrastDefinition",
    "DifferentialDesignMatrixSummary",
    "DifferentialEmpiricalBayesProvenance",
    "DifferentialFixedEffectCovariateProvenance",
    "DifferentialMissingValuePolicyProvenance",
    "DifferentialPolicyProvenance",
    "DifferentialReplicatePolicyProvenance",
    "DifferentialStatisticalTestingProvenance",
    "DifferentialTechnicalReplicateGroup",
    "DifferentialUnsupportedDesignPolicyProvenance",
    "EnrichmentResultRecord",
    "EnrichmentWorkflowResult",
    "IMPORTER_QUALITY_STATUS_NOT_APPLICABLE",
    "IMPORTER_QUALITY_STATUS_NOT_REPORTED",
    "IMPORTER_QUALITY_STATUS_REPORTED",
    "ImporterDetectedIntensityColumn",
    "ImporterDuplicateKeySummary",
    "ImporterFlaggedRowSummary",
    "ImporterLocalisationConfidenceSummary",
    "ImporterMissingIntensitySummary",
    "ImporterQualityCount",
    "ImporterQualityReport",
    "ImporterQualityStatus",
    "KinaseActivityResult",
    "KinasePredictionResult",
    "KinaseScoringResult",
    "KinaseWorkflowResult",
    "PhosphositeImportResult",
    "ProteinAwareMappingDiagnostics",
    "ProteinAwarePreparationReport",
    "ProteinAwarePreparationResult",
    "ProteinAwareSiteEligibility",
    "SignalomeWorkflowResult",
} | INTENTIONAL_RESULT_COMPATIBILITY_ALIASES


def _dataset() -> AnalysisReadyPhosphoDataset:
    site_ids = ["MAPK14;Y182;", "GSK3B;S9;", "AKT1;T308;"]
    site_index = site_key_index_from_display_ids(site_ids)
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 2.0, 4.0],
            "sample_b": [2.0, 4.0, 1.0],
        },
        index=site_index,
    )
    site_metadata = pd.DataFrame(
        {
            "site_key": site_index.astype(str).tolist(),
            "display_id": site_ids,
            **site_key_context_columns(site_index),
            "gene_symbol": ["MAPK14", "GSK3B", "AKT1"],
            "site": ["Y182", "S9", "T308"],
            "site_sequence": [
                ("A" * 15) + "Y" + ("A" * 15),
                ("A" * 15) + "S" + ("A" * 15),
                ("A" * 15) + "T" + ("A" * 15),
            ],
            "protein_id": ["MAPK14", "GSK3B", "AKT1"],
        },
        index=site_index.copy(),
    )
    return AnalysisReadyPhosphoDataset(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
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
            {
                "site_sequence": [
                    ("A" * 15) + "Y" + ("A" * 15),
                    ("A" * 15) + "S" + ("A" * 15),
                    ("A" * 15) + "T" + ("A" * 15),
                ]
            },
            index=site_ids,
        ),
    )


def _kinase_result() -> KinaseWorkflowResult:
    return KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=_dataset(),
            references=_references(),
            scoring_config=KinaseScoringConfig(min_substrates=2),
            prediction_config=KinasePredictionConfig(
                top_k=1,
                deterministic_max_selected_kinases=2,
                adaptive_ensemble_runs=2,
            ),
            activity_config=None,
        )
    )


def test_public_result_exports_match_contract() -> None:
    assert set(result_models.__all__) == EXPECTED_RESULT_EXPORTS


def test_contract_result_domain_modules_preserve_legacy_import_identity() -> None:
    import phospy.contracts.results as contracts_results
    from phospy.contracts.results.base import (
        ImporterQualityReport,
        PhosphositeImportResult,
    )
    from phospy.contracts.results.differential import DifferentialAnalysisResult
    from phospy.contracts.results.enrichment import EnrichmentWorkflowResult
    from phospy.contracts.results.kinase import (
        KinaseEligibilityReport,
        KinaseWorkflowResult,
    )
    from phospy.contracts.results.preprocessing import BatchCorrectionReport
    from phospy.contracts.results.signalome import SignalomeWorkflowResult

    assert contracts_results.ImporterQualityReport is ImporterQualityReport
    assert contracts_results.PhosphositeImportResult is PhosphositeImportResult
    assert contracts_results.DifferentialAnalysisResult is DifferentialAnalysisResult
    assert contracts_results.EnrichmentWorkflowResult is EnrichmentWorkflowResult
    assert contracts_results.KinaseEligibilityReport is KinaseEligibilityReport
    assert contracts_results.KinaseWorkflowResult is KinaseWorkflowResult
    assert contracts_results.BatchCorrectionReport is BatchCorrectionReport
    assert contracts_results.SignalomeWorkflowResult is SignalomeWorkflowResult


def test_result_compatibility_aliases_are_public_exports() -> None:
    assert INTENTIONAL_RESULT_COMPATIBILITY_ALIASES <= set(result_models.__all__)


def test_result_star_import_exposes_public_contract_without_internals() -> None:
    namespace: dict[str, object] = {}

    exec("from phospy.api.results import *", namespace)

    exported_names = {name for name in namespace if name != "__builtins__"}
    assert exported_names == EXPECTED_RESULT_EXPORTS
    for name in EXPECTED_RESULT_EXPORTS:
        assert namespace[name] is getattr(result_models, name)
    assert "_result_contracts" not in namespace


def test_kinase_result_stays_nested_and_honest_for_supported_lane() -> None:
    result = _kinase_result()
    assert isinstance(result, KinaseWorkflowResult)
    assert isinstance(result.scoring_result, KinaseScoringResult)
    assert isinstance(result.prediction_result, KinasePredictionResult)
    assert not result.scoring_result.profile_scores.empty
    assert result.scoring_result.motif_scores is None
    assert result.scoring_result.rank_weighted_fusion_scores is not None
    assert result.scoring_result.score_fusion_weights is None
    assert result.scoring_result.score_source_summary is not None
    assert result.scoring_result.score_source_matrix is None
    assert result.activity_result is None
    assert result.substrate_contributions is None
    assert result.prediction_result.substrate_list is not None
    assert {
        "kinase",
        "substrate_site",
        "score",
        "rank",
        "site_key",
        "display_id",
    }.issubset(set(result.prediction_result.substrate_list.columns))
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
            config=build_signalome_config(substrate_support_cutoff=0.5),
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
    assert signalome_result.site_membership is not None
    assert signalome_result.protein_site_context is not None
    assert signalome_result.expanded_signalome is not None
    assert not signalome_result.expanded_signalome.empty
    assert not signalome_result.site_membership.empty
    assert not signalome_result.protein_site_context.empty
    assert {"site_key", "display_id"} <= set(
        signalome_result.module_assignments.table.columns
    )
    assert {"site_key", "display_id"} <= set(signalome_result.site_membership.columns)
    assert {"site_key", "display_id"} <= set(
        signalome_result.expanded_signalome.columns
    )
    assert (
        signalome_result.module_assignments.table.loc[:, "site_key"].tolist()
        == signalome_result.module_assignments.table.index.astype(str).tolist()
    )
    assert (
        signalome_result.module_assignments.table.loc[:, "display_id"]
        .astype(str)
        .str.strip()
        .ne("")
        .all()
    )
    assert signalome_result.kinase_result.scoring_result.motif_scores is None
    assert (
        signalome_result.kinase_result.scoring_result.rank_weighted_fusion_scores
        is not None
    )
    assert signalome_result.kinase_result.scoring_result.score_fusion_weights is None
    assert not hasattr(signalome_result, "pred_mat")
    assert not hasattr(signalome_result, "profile_scores")


def test_kinase_result_exposes_supported_activity_stage_outputs_when_enabled() -> None:
    result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=_dataset(),
            references=_references(),
            scoring_config=KinaseScoringConfig(min_substrates=2),
            prediction_config=KinasePredictionConfig(
                top_k=1,
                deterministic_max_selected_kinases=2,
                adaptive_ensemble_runs=2,
            ),
            activity_config=KinaseActivityConfig(
                enabled=True,
                threshold=0.0,
                min_substrates=1,
                top_n_substrates=2,
            ),
        )
    )

    assert result.activity_result is not None
    assert {"site_id", "site_key", "display_id", "kinase", "score"} <= set(
        result.activity_result.target_table.columns
    )
    assert result.activity_result.thresholded_substrate_counts.name == "n_substrates"
    assert result.activity_result.target_counts.name == "n_targets"
    assert result.activity_result.activity_method.activity_method_id == (
        "simplified_weighted_substrate_activity_v1"
    )
    assert result.activity_result.activity_method.is_ksea is False
    assert (
        result.activity_result.activity_method.is_phosr_kinase_activity_equivalent
        is False
    )
    assert hasattr(result.activity_result, "activity_matrix")
    assert hasattr(result.activity_result, "thresholded_substrate_mean_activity")
    assert hasattr(result.activity_result, "thresholded_substrate_counts")
    assert hasattr(result.activity_result, "activity_substrate_counts")
    assert result.activity_result.activity_substrate_counts is None
    assert not hasattr(result.activity_result, "ksea_scores")
    assert not hasattr(result.activity_result, "ksea_counts")
    assert result.provenance is not None
    policy_ids = {policy.id for policy in result.provenance.scientific_policies}
    assert ScientificPolicyId.SIMPLIFIED_WEIGHTED_SUBSTRATE_ACTIVITY in policy_ids


def test_kinase_result_can_include_opt_in_diagnostic_scoring_tables() -> None:
    result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=_dataset(),
            references=_references(),
            scoring_config=KinaseScoringConfig(
                min_substrates=2,
                include_diagnostic_scoring_tables=True,
            ),
            prediction_config=KinasePredictionConfig(
                top_k=1,
                deterministic_max_selected_kinases=2,
                adaptive_ensemble_runs=2,
            ),
            activity_config=None,
        )
    )

    assert result.scoring_result.motif_scores is not None
    assert result.scoring_result.score_fusion_weights is not None
    assert result.scoring_result.score_source_summary is not None
    assert result.scoring_result.score_source_matrix is not None
    assert hasattr(result.scoring_result, "rank_weighted_fusion_scores")
    assert hasattr(result.scoring_result, "score_fusion_weights")
    assert not hasattr(result.scoring_result, "combined_scores")
    assert not hasattr(result.scoring_result, "weights")


def test_kinase_provenance_uses_renamed_scoring_and_activity_output_tables() -> None:
    result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=_dataset(),
            references=_references(),
            scoring_config=KinaseScoringConfig(
                min_substrates=2,
                include_diagnostic_scoring_tables=True,
            ),
            prediction_config=KinasePredictionConfig(
                top_k=1,
                deterministic_max_selected_kinases=2,
                adaptive_ensemble_runs=2,
            ),
            activity_config=KinaseActivityConfig(
                enabled=True,
                threshold=0.0,
                min_substrates=1,
                top_n_substrates=2,
            ),
        )
    )
    assert result.provenance is not None

    output_names = {table.name for table in result.provenance.output_tables}
    assert "outputs.scoring.rank_weighted_fusion_scores" in output_names
    assert "outputs.scoring.score_fusion_weights" in output_names
    assert "outputs.activity.thresholded_substrate_mean_activity" in output_names
    assert "outputs.activity.thresholded_substrate_counts" in output_names
    assert "outputs.scoring.combined_scores" not in output_names
    assert "outputs.scoring.weights" not in output_names
    assert "outputs.activity.ksea_scores" not in output_names
    assert "outputs.activity.ksea_counts" not in output_names
