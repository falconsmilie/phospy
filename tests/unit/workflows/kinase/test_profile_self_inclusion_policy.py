from __future__ import annotations

from collections.abc import Mapping

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
from phospy.api.configs import (
    KinaseScoringConfig,
    ProfileSelfInclusionPolicy,
    ReferenceContextCompatibilityPolicy,
)
from phospy.provenance.scientific_policy_models import ScientificPolicyId
from phospy.tables.kinase import (
    KINASE_PROFILE_SCORE_DIAGNOSTIC_REASON_INSUFFICIENT_SUBSTRATES_AFTER_LEAVE_ONE_OUT,
    KINASE_PROFILE_SCORE_DIAGNOSTIC_STATUS_SCORED,
    KINASE_PROFILE_SCORE_DIAGNOSTIC_STATUS_UNSCORED,
)
from phospy.workflows.kinase.caveats import (
    KINASE_PROFILE_LEAVE_ONE_OUT_CAVEAT_CODE,
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


def _references_with_short_and_long_kinases() -> ReferenceBundle:
    display_ids = ["AKT1;S1;", "MAPK1;S2;", "GSK3B;S3;"]
    return ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": [
                    "KINASE_SHORT",
                    "KINASE_SHORT",
                    "KINASE_LONG",
                    "KINASE_LONG",
                    "KINASE_LONG",
                ],
                "substrate_site": [
                    display_ids[0],
                    display_ids[1],
                    display_ids[0],
                    display_ids[1],
                    display_ids[2],
                ],
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
    references: ReferenceBundle | None = None,
) -> KinaseWorkflowRequest:
    return KinaseWorkflowRequest(
        dataset=_dataset(),
        references=references or _references(),
        scoring_config=scoring_config
        or KinaseScoringConfig(
            min_substrates=2,
            reference_context_compatibility_policy=(
                ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT
            ),
        ),
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
    assert isinstance(scoring_config, Mapping)
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


def test_leave_one_out_changes_known_substrate_profile_score() -> None:
    allow_result = KinaseWorkflow().run(_request())
    leave_one_out_result = KinaseWorkflow().run(
        _request(
            KinaseScoringConfig(
                min_substrates=2,
                profile_self_inclusion_policy=(
                    ProfileSelfInclusionPolicy.LEAVE_ONE_OUT
                ),
                reference_context_compatibility_policy=(
                    ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT
                ),
            )
        )
    )
    scored_site = allow_result.scoring_result.profile_scores.index[0]

    assert allow_result.scoring_result.profile_score_diagnostics is None
    assert (
        leave_one_out_result.scoring_result.profile_self_inclusion_policy
        is ProfileSelfInclusionPolicy.LEAVE_ONE_OUT
    )
    assert (
        leave_one_out_result.scoring_result.profile_scores.at[scored_site, "KINASE_A"]
        != allow_result.scoring_result.profile_scores.at[scored_site, "KINASE_A"]
    )
    diagnostics = leave_one_out_result.scoring_result.profile_score_diagnostics
    assert diagnostics is not None
    assert set(diagnostics.loc[:, "status"]) == {
        KINASE_PROFILE_SCORE_DIAGNOSTIC_STATUS_SCORED
    }


def test_leave_one_out_insufficient_substrates_are_explicitly_diagnosed() -> None:
    result = KinaseWorkflow().run(
        _request(
            KinaseScoringConfig(
                min_substrates=2,
                profile_self_inclusion_policy=(
                    ProfileSelfInclusionPolicy.LEAVE_ONE_OUT
                ),
                reference_context_compatibility_policy=(
                    ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT
                ),
            ),
            references=_references_with_short_and_long_kinases(),
        )
    )

    diagnostics = result.scoring_result.profile_score_diagnostics
    assert diagnostics is not None
    short_rows = diagnostics.loc[
        diagnostics.loc[:, "kinase"] == "KINASE_SHORT",
        :,
    ]
    assert set(short_rows.loc[:, "status"]) == {
        KINASE_PROFILE_SCORE_DIAGNOSTIC_STATUS_UNSCORED
    }
    assert set(short_rows.loc[:, "reason"]) == {
        KINASE_PROFILE_SCORE_DIAGNOSTIC_REASON_INSUFFICIENT_SUBSTRATES_AFTER_LEAVE_ONE_OUT
    }
    for substrate_site in short_rows.loc[:, "substrate_site"].astype(str):
        assert pd.isna(
            result.scoring_result.profile_scores.at[substrate_site, "KINASE_SHORT"]
        )
    non_known_short_sites = [
        site
        for site in result.scoring_result.profile_scores.index.astype(str)
        if site not in set(short_rows.loc[:, "substrate_site"].astype(str))
    ]
    assert (
        result.scoring_result.profile_scores.loc[non_known_short_sites, "KINASE_SHORT"]
        .notna()
        .any()
    )
    assert result.scoring_result.profile_scores.loc[:, "KINASE_LONG"].notna().any()


def test_leave_one_out_result_records_caveat_and_provenance() -> None:
    result = KinaseWorkflow().run(
        _request(
            KinaseScoringConfig(
                min_substrates=2,
                profile_self_inclusion_policy=(
                    ProfileSelfInclusionPolicy.LEAVE_ONE_OUT
                ),
                reference_context_compatibility_policy=(
                    ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT
                ),
            ),
            references=_references_with_short_and_long_kinases(),
        )
    )

    assert result.provenance is not None
    scoring_config = result.provenance.workflow_parameters["scoring_config"]
    assert isinstance(scoring_config, Mapping)
    assert scoring_config["profile_self_inclusion_policy"] == "leave_one_out"
    scoring_diagnostics = result.provenance.workflow_parameters["scoring_diagnostics"]
    assert isinstance(scoring_diagnostics, Mapping)
    profile_diagnostics = scoring_diagnostics["profile_score_diagnostics"]
    assert isinstance(profile_diagnostics, Mapping)
    assert profile_diagnostics["unscored_reason_counts"] == {
        KINASE_PROFILE_SCORE_DIAGNOSTIC_REASON_INSUFFICIENT_SUBSTRATES_AFTER_LEAVE_ONE_OUT: 2
    }
    policy_record = next(
        policy
        for policy in result.provenance.scientific_policies
        if policy.id is ScientificPolicyId.KINASE_PROFILE_SCORING
    )
    assert policy_record.parameters["profile_self_inclusion_policy"] == "leave_one_out"
    assert policy_record.parameters["self_inclusion_allowed"] is False
    assert policy_record.parameters["leave_one_out_enabled"] is True

    caveat = next(
        caveat
        for caveat in result.caveats
        if caveat.code == KINASE_PROFILE_LEAVE_ONE_OUT_CAVEAT_CODE
    )
    assert caveat.severity == "info"
    assert "Leave-one-out profile scoring was used" in caveat.message
    assert caveat.details["profile_self_inclusion_policy"] == "leave_one_out"
    assert caveat.details["self_inclusion_allowed"] is False
    assert caveat.details["leave_one_out_enabled"] is True
    assert caveat.details["leave_one_out_unscored_cell_count"] == 2
