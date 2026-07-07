from __future__ import annotations

import pandas as pd
import pandas.testing as pdt

from phospy.api import (
    AnalysisReadyPhosphoDataset,
    KinasePredictionConfig,
    KinaseWorkflow,
    KinaseWorkflowRequest,
    Organism,
    ReferenceBundle,
)
from phospy.api.configs import KinaseScoringConfig, ProfileSelfInclusionPolicy
from phospy.provenance.scientific_policy_models import ScientificPolicyId
from phospy.workflows.kinase.caveats import (
    KINASE_PROFILE_SELF_INCLUSION_CAVEAT_CODE,
)
from phospy.workflows.kinase.interpreter import KinaseWorkflowInterpreter
from phospy.workflows.kinase.science import (
    build_kinase_profiles,
    score_profile_correlations,
)
from phospy.workflows.kinase.scoring_runner import KinaseScoringRunner
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
)
from tests.support.site_keys import (
    site_key_context_columns,
    site_key_index_from_display_ids,
)


def _dataset() -> AnalysisReadyPhosphoDataset:
    display_ids = ["AKT1;S1;", "MAPK1;S2;", "GSK3B;S3;"]
    site_index = site_key_index_from_display_ids(display_ids)
    return AnalysisReadyPhosphoDataset(
        phospho=pd.DataFrame(
            {
                "sample_a": [1.0, 2.0, 8.0],
                "sample_b": [2.0, 3.0, 1.0],
                "sample_c": [3.0, 4.0, 5.0],
            },
            index=site_index,
        ),
        site_metadata=pd.DataFrame(
            {
                "site_key": site_index.astype(str).tolist(),
                "display_id": display_ids,
                **site_key_context_columns(site_index),
                "gene_symbol": ["AKT1", "MAPK1", "GSK3B"],
                "protein_id": ["AKT1", "MAPK1", "GSK3B"],
                "site": ["S1", "S2", "S3"],
                "site_sequence": [
                    ("A" * 15) + display_id.split(";")[1][0] + ("A" * 15)
                    for display_id in display_ids
                ],
            },
            index=site_index.copy(),
        ),
        organism=Organism.RAT,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
    )


def _references() -> ReferenceBundle:
    display_ids = ["AKT1;S1;", "MAPK1;S2;", "GSK3B;S3;"]
    return ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["KINASE_A", "KINASE_A", "KINASE_A"],
                "substrate_site": display_ids,
            }
        ),
        site_sequences=pd.DataFrame(
            {
                "site_sequence": [
                    ("A" * 15) + display_id.split(";")[1][0] + ("A" * 15)
                    for display_id in display_ids
                ]
            },
            index=pd.Index(display_ids, name="site_id"),
        ),
    )


def _request(
    scoring_config: KinaseScoringConfig | None = None,
) -> KinaseWorkflowRequest:
    return KinaseWorkflowRequest(
        dataset=_dataset(),
        references=_references(),
        scoring_config=scoring_config or KinaseScoringConfig(min_substrates=2),
        prediction_config=KinasePredictionConfig(
            top_k=3,
            deterministic_max_selected_kinases=1,
            adaptive_ensemble_runs=2,
        ),
        activity_config=None,
    )


def test_kinase_result_records_default_profile_self_inclusion_policy() -> None:
    result = KinaseWorkflow().run(_request())

    assert (
        result.scoring_result.profile_self_inclusion_policy
        is ProfileSelfInclusionPolicy.ALLOW
    )
    assert result.provenance is not None
    scoring_config = result.provenance.workflow_parameters["scoring_config"]
    assert isinstance(scoring_config, dict)
    assert scoring_config["profile_self_inclusion_policy"] == "allow"
    policy_record = next(
        policy
        for policy in result.provenance.scientific_policies
        if policy.id is ScientificPolicyId.KINASE_PROFILE_SCORING
    )
    assert policy_record.parameters["profile_self_inclusion_policy"] == "allow"
    assert policy_record.parameters["self_inclusion_allowed"] is True


def test_allow_profile_self_inclusion_policy_produces_caveat() -> None:
    result = KinaseWorkflow().run(_request())

    caveat = next(
        caveat
        for caveat in result.caveats
        if caveat.code == KINASE_PROFILE_SELF_INCLUSION_CAVEAT_CODE
    )
    assert caveat.severity == "warning"
    assert "allowed self-inclusion" in caveat.message
    assert "Scores are exploratory and may be inflated" in caveat.message
    assert caveat.details["profile_self_inclusion_policy"] == "allow"
    assert caveat.details["self_inclusion_allowed"] is True


def test_default_profile_scores_match_existing_profile_scoring_path() -> None:
    resolved = KinaseWorkflowInterpreter().run(_request())
    profile_build = build_kinase_profiles(
        phospho=resolved.activity_phospho_matrix,
        kinase_substrate_map=resolved.kinase_substrate_map,
        min_substrates=resolved.execution_config.scoring_min_substrates,
        allow_single_substrate_profiles=False,
        profile_missing_value_strategy=(
            resolved.execution_config.profile_missing_value_strategy
        ),
    )
    expected_scores = score_profile_correlations(
        phospho=resolved.activity_phospho_matrix,
        profile_matrix=profile_build.profile_matrix,
    )

    scoring_execution = KinaseScoringRunner().run(
        request=resolved,
        config=resolved.execution_config,
    )

    pdt.assert_frame_equal(
        scoring_execution.scoring_result.profile_scores,
        expected_scores,
    )
