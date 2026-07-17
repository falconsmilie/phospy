from __future__ import annotations

from collections.abc import Mapping

import pandas as pd
import pytest

from phospy.api import (
    KinaseAttritionPolicy,
    KinasePredictionConfig,
    KinaseScoringConfig,
    KinaseWorkflow,
    KinaseWorkflowRequest,
    Organism,
    ReferenceBundle,
    ReferenceContextCompatibilityPolicy,
)
from phospy.api.results import KinaseWorkflowResult
from phospy.errors.workflows import WorkflowBoundaryError
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.workflows.intensity_scale_evidence import (
    INPUT_INTENSITY_SCALE_DECLARED_CAVEAT_CODE,
)
from phospy.workflows.kinase.caveats import (
    KINASE_ATTRITION_POLICY_CAVEAT_CODE,
    KINASE_PERMISSIVE_LOCALISATION_POLICY_CAVEAT_CODE,
)
from phospy.workflows.kinase.contracts import ResolvedKinaseWorkflowRequest
from phospy.workflows.kinase.interpreter import KinaseWorkflowInterpreter
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
)
from tests.support.site_keys import (
    site_key_context_columns,
    site_key_index_from_display_ids,
)


def _window(display_id: str) -> str:
    residue = display_id.split(";")[1][0].upper()
    return ("A" * 15) + residue + ("A" * 15)


def _dataset(display_ids: list[str]) -> AnalysisReadyPhosphoDataset:
    site_index = site_key_index_from_display_ids(display_ids)
    return AnalysisReadyPhosphoDataset(
        phospho=pd.DataFrame(
            {
                "sample_a": [1.0 + index for index, _ in enumerate(display_ids)],
                "sample_b": [2.0 + index for index, _ in enumerate(display_ids)],
            },
            index=site_index.copy(),
        ),
        site_metadata=pd.DataFrame(
            {
                "site_key": site_index.astype(str).tolist(),
                "display_id": display_ids,
                **site_key_context_columns(site_index),
                "gene_symbol": [
                    display_id.split(";", 1)[0] for display_id in display_ids
                ],
                "protein_id": [
                    display_id.split(";", 1)[0] for display_id in display_ids
                ],
                "site": [display_id.split(";")[1] for display_id in display_ids],
                "site_sequence": [_window(display_id) for display_id in display_ids],
            },
            index=site_index.copy(),
        ),
        organism=Organism.RAT,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
    )


def _references(substrate_display_ids: list[str]) -> ReferenceBundle:
    return ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["KINASE_A" for _ in substrate_display_ids],
                "substrate_site": substrate_display_ids,
            }
        ),
        site_sequences=pd.DataFrame(
            {
                "site_sequence": [
                    _window(display_id) for display_id in substrate_display_ids
                ]
            },
            index=pd.Index(substrate_display_ids, name="site_id"),
        ),
    )


def _request(policy: KinaseAttritionPolicy) -> KinaseWorkflowRequest:
    display_ids = ["KIN1;S1;", "KIN2;T2;", "KIN3;Y3;", "KIN4;S4;"]
    return KinaseWorkflowRequest(
        dataset=_dataset(display_ids),
        references=_references(display_ids[:2]),
        scoring_config=KinaseScoringConfig(
            min_substrates=2,
            attrition_policy=policy,
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
    )


def _strict_scored_fraction_policy(*, on_violation: str) -> KinaseAttritionPolicy:
    return KinaseAttritionPolicy(
        minimum_scored_fraction=0.75,
        on_violation=on_violation,  # type: ignore[arg-type]
    )


def _caveat_by_code(result: KinaseWorkflowResult, code: str):
    matches = [caveat for caveat in result.caveats if caveat.code == code]
    assert len(matches) == 1
    return matches[0]


def test_kinase_attrition_policy_error_blocks_scoring() -> None:
    class _ExecutorMustNotRun:
        def run(
            self,
            request: ResolvedKinaseWorkflowRequest,
        ) -> KinaseWorkflowResult:
            _ = request
            raise AssertionError("kinase scoring must not run after policy error")

    with pytest.raises(WorkflowBoundaryError) as exc_info:
        KinaseWorkflow(executor=_ExecutorMustNotRun()).run(
            _request(_strict_scored_fraction_policy(on_violation="error"))
        )

    assert exc_info.value.seam == "kinase.interpreter.attrition_policy"


def test_kinase_result_caveats_include_attrition_warning() -> None:
    result = KinaseWorkflow().run(
        _request(_strict_scored_fraction_policy(on_violation="warn"))
    )

    assert not result.prediction_result.pred_mat.empty
    assert result.attrition_provenance is not None
    assert result.attrition_provenance.policy_outcome == "warned"
    assert result.attrition_provenance.metrics["scored_fraction"] == pytest.approx(0.5)
    assert result.attrition_provenance.policy["minimum_scored_fraction"] == (
        pytest.approx(0.75)
    )
    caveat = _caveat_by_code(result, KINASE_ATTRITION_POLICY_CAVEAT_CODE)
    assert caveat.severity == "warning"
    assert caveat.details["threshold_name"] == "minimum_scored_fraction"
    assert caveat.details["configured_threshold"] == pytest.approx(0.75)
    assert caveat.details["observed_value"] == pytest.approx(0.5)
    assert result.attrition_provenance.warning_messages == (caveat.message,)


def test_kinase_result_caveats_include_permissive_localisation_policy() -> None:
    result = KinaseWorkflow().run(
        _request(
            KinaseAttritionPolicy(
                minimum_reference_overlap_fraction=0.25,
                minimum_sequence_supported_fraction=0.25,
                minimum_scored_fraction=0.25,
                on_violation="warn",
            )
        )
    )

    caveat = _caveat_by_code(
        result,
        KINASE_PERMISSIVE_LOCALISATION_POLICY_CAVEAT_CODE,
    )

    assert caveat.severity == "warning"
    assert caveat.details["policy"] == "allow_unknown"
    assert caveat.details["workflow_scope"] == "kinase_scoring"
    assert caveat.details["minimum_probability"] is None
    assert caveat.details["site_count"] == 4


def test_kinase_result_caveats_include_declared_input_scale_evidence() -> None:
    result = KinaseWorkflow().run(
        _request(
            KinaseAttritionPolicy(
                minimum_reference_overlap_fraction=0.25,
                minimum_sequence_supported_fraction=0.25,
                minimum_scored_fraction=0.25,
                on_violation="warn",
            )
        )
    )

    caveat = _caveat_by_code(result, INPUT_INTENSITY_SCALE_DECLARED_CAVEAT_CODE)

    assert caveat.severity == "warning"
    assert caveat.details["input_intensity_scale"] == "linear"
    assert caveat.details["input_intensity_scale_evidence_level"] == (
        "declared_by_user"
    )
    assert caveat.details["input_intensity_scale_source"] == "declared_by_user"
    assert caveat.details["workflow_scope"] == "kinase_scoring"


def test_kinase_attrition_policy_error_message_contains_counts_and_threshold() -> None:
    with pytest.raises(WorkflowBoundaryError) as exc_info:
        KinaseWorkflowInterpreter().run(
            _request(_strict_scored_fraction_policy(on_violation="error"))
        )

    message = str(exc_info.value)
    assert "Kinase scoring retained 50.0% of dataset sites" in message
    assert "minimum_scored_fraction=75.0%" in message
    assert "scored_sites=2" in message
    assert "total_dataset_sites=4" in message
    assert exc_info.value.details["scored_sites"] == 2
    assert exc_info.value.details["policy_outcome"] == "failed"
    violations = exc_info.value.details["violations"]
    assert isinstance(violations, list)
    assert violations[0]["threshold_name"] == "minimum_scored_fraction"
    attrition_provenance = exc_info.value.details["attrition_provenance"]
    assert isinstance(attrition_provenance, dict)
    assert attrition_provenance["policy_outcome"] == "failed"


def test_kinase_result_exposes_attrition_metrics() -> None:
    result = KinaseWorkflow().run(
        _request(
            KinaseAttritionPolicy(
                minimum_reference_overlap_fraction=0.25,
                minimum_sequence_supported_fraction=0.25,
                minimum_scored_fraction=0.25,
                on_violation="warn",
            )
        )
    )

    attrition = result.attrition_provenance
    assert attrition is not None
    assert attrition.policy_outcome == "passed"
    assert attrition.metrics["total_dataset_sites"] == 4
    assert attrition.metrics["reference_overlap_sites"] == 2
    assert attrition.metrics["sequence_supported_sites"] == 4
    assert attrition.metrics["scored_sites"] == 2
    assert attrition.policy["minimum_scored_fraction"] == pytest.approx(0.25)
    assert attrition.policy["on_violation"] == "warn"
    assert attrition.policy_violations == ()
    assert all(
        caveat.code != KINASE_ATTRITION_POLICY_CAVEAT_CODE for caveat in result.caveats
    )


def test_kinase_provenance_records_attrition_policy_and_metrics() -> None:
    result = KinaseWorkflow().run(
        _request(_strict_scored_fraction_policy(on_violation="warn"))
    )

    summary = result.site_attrition_summary
    assert summary is not None
    caveat = _caveat_by_code(result, KINASE_ATTRITION_POLICY_CAVEAT_CODE)
    assert caveat.details["total_dataset_sites"] == (
        summary.scoring.final_quantitative_sites_entering_scoring
    )
    assert caveat.details["scored_sites"] == (
        summary.scoring.sites_with_kinase_substrate_reference_profile_evidence
    )
    assert caveat.details["observed_value"] == pytest.approx(
        summary.scoring.sites_with_kinase_substrate_reference_profile_evidence
        / summary.scoring.final_quantitative_sites_entering_scoring
    )

    assert result.provenance is not None
    assert result.attrition_provenance is not None
    workflow_parameters = result.provenance.workflow_parameters
    assert workflow_parameters["input_intensity_scale"] == "linear"
    assert workflow_parameters["input_intensity_scale_evidence_level"] == (
        "declared_by_user"
    )
    assert workflow_parameters["input_intensity_scale_source"] == "declared_by_user"
    attrition_provenance = workflow_parameters["attrition_provenance"]
    assert isinstance(attrition_provenance, Mapping)
    assert attrition_provenance["policy_outcome"] == "warned"
    assert attrition_provenance["metrics"] == result.attrition_provenance.metrics
    assert attrition_provenance["policy"] == result.attrition_provenance.policy
    row_attrition_metrics = workflow_parameters["row_attrition_metrics"]
    assert isinstance(row_attrition_metrics, Mapping)
    assert row_attrition_metrics["input_sites"] == 4
    assert row_attrition_metrics["sites_missing_valid_centered_sequence"] == 0
    assert row_attrition_metrics["sites_not_present_in_reference_resource"] == 2
    assert row_attrition_metrics["sites_with_reference_and_sequence_support"] == 2
    assert (
        row_attrition_metrics["site_kinase_pairs_considered"]
        >= (row_attrition_metrics["site_kinase_pairs_scored"])
    )
    assert "row_attrition" not in workflow_parameters
    assert (
        result.scoring_result.authoritative_scores.shape[0]
        == (row_attrition_metrics["input_sites"])
    )
    scoring_diagnostics = result.provenance.workflow_parameters["scoring_diagnostics"]
    assert isinstance(scoring_diagnostics, Mapping)
    attrition_metrics = scoring_diagnostics["attrition_metrics"]
    assert isinstance(attrition_metrics, Mapping)
    assert attrition_metrics["scored_fraction"] == pytest.approx(
        caveat.details["observed_value"]
    )
    policy_violations = scoring_diagnostics["attrition_policy_violations"]
    assert isinstance(policy_violations, tuple)
    assert policy_violations[0]["scored_sites"] == caveat.details["scored_sites"]
    scoring_config = result.provenance.workflow_parameters["scoring_config"]
    assert isinstance(scoring_config, Mapping)
    assert scoring_config["attrition_policy"] == result.attrition_provenance.policy
