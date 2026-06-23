from __future__ import annotations

import pandas as pd
import pytest

from phospy import AnalysisReadyPhosphoDataset
from phospy.api import (
    KINASE_REFERENCE_DISPLAY_AMBIGUITY_POLICY_ALLOW_WITH_DIAGNOSTICS,
    KinaseActivityConfig,
    KinasePredictionConfig,
    KinaseScoringConfig,
    KinaseWorkflow,
    KinaseWorkflowRequest,
    Organism,
    ReferenceBundle,
)
from phospy.science.prediction.scoring import (
    KINASE_SCORE_SOURCE_FUSED_MOTIF_PROFILE_EVIDENCE,
    KINASE_SCORE_SOURCE_PROFILE_ONLY_MOTIF_MISSING_OR_CONSTANT,
)
from phospy.workflows.kinase import scoring_runner as scoring_runner_module
from phospy.workflows.kinase.contributions import (
    KINASE_SUBSTRATE_CONTRIBUTION_COLUMNS,
    KINASE_SUBSTRATE_CONTRIBUTION_EXCLUDED_BELOW_MIN_SUBSTRATES,
    KINASE_SUBSTRATE_CONTRIBUTION_EXCLUDED_MISSING_SCORE_VALUE,
    KINASE_SUBSTRATE_CONTRIBUTION_EXCLUDED_NOT_QUANTIFIED,
    KINASE_SUBSTRATE_CONTRIBUTION_STATUS_EXCLUDED,
    KINASE_SUBSTRATE_CONTRIBUTION_STATUS_INCLUDED,
    KinaseSubstrateContributionReferenceSource,
    build_kinase_substrate_contribution_table,
)
from phospy.workflows.kinase.interpreter import KinaseWorkflowInterpreter
from phospy.workflows.kinase.scoring_runner import KinaseScoringRunner
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
)
from tests.support.site_keys import (
    protein_site_key_index,
    site_key_context_columns,
    site_key_index_from_display_ids,
)


def test_contribution_table_records_statuses_sources_and_ambiguity() -> None:
    contribution_table = build_kinase_substrate_contribution_table(
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["K1", "K1", "K2", "K3"],
                "substrate_site": ["S1", "S2", "S3", "S_MISSING"],
                "display_id": ["D_SHARED", "D_SHARED", "D3", "DMISS"],
            }
        ),
        scoring_values=pd.DataFrame(
            {"K1": [0.8, float("nan"), 0.2]},
            index=pd.Index(["S1", "S2", "S3"], name="site_key"),
        ),
        score_component="rank_weighted_fusion_scores",
        quantified_substrates={"K1": ["S1", "S2"]},
        substrate_counts=pd.Series(
            {"K1": 2, "K2": 1, "K3": 1},
            dtype="int64",
            name="NumSub",
        ),
        min_substrates=2,
        score_source_matrix=pd.DataFrame(
            {
                "K1": [
                    KINASE_SCORE_SOURCE_FUSED_MOTIF_PROFILE_EVIDENCE,
                    KINASE_SCORE_SOURCE_PROFILE_ONLY_MOTIF_MISSING_OR_CONSTANT,
                    "unavailable_no_score",
                ],
            },
            index=pd.Index(["S1", "S2", "S3"], name="site_key"),
        ),
        reference_source=KinaseSubstrateContributionReferenceSource(
            source_name="fixture_reference",
            source_version="v1",
            bundle_id="fixture_bundle",
            identifier_namespace="fixture_display_id",
        ),
        display_reference_matching={
            "one_to_many_display_reference_matches": [
                {"display_id": "D_SHARED", "site_keys": ["S1", "S2"]}
            ],
        },
    )

    assert tuple(contribution_table.columns) == KINASE_SUBSTRATE_CONTRIBUTION_COLUMNS
    rows = contribution_table.to_dict(orient="records")
    assert rows[0]["kinase"] == "K1"
    assert rows[0]["substrate_site"] == "S1"
    assert rows[0]["substrate_identifier"] == "D_SHARED"
    assert rows[0]["value_used_in_scoring"] == pytest.approx(0.8)
    assert rows[0]["score_component"] == "rank_weighted_fusion_scores"
    assert rows[0]["score_source"] == KINASE_SCORE_SOURCE_FUSED_MOTIF_PROFILE_EVIDENCE
    assert rows[0]["reference_source_name"] == "fixture_reference"
    assert rows[0]["reference_source_version"] == "v1"
    assert rows[0]["reference_bundle_id"] == "fixture_bundle"
    assert rows[0]["reference_identifier_namespace"] == "fixture_display_id"
    assert rows[0]["status"] == KINASE_SUBSTRATE_CONTRIBUTION_STATUS_INCLUDED
    assert pd.isna(rows[0]["exclusion_reason"])
    assert rows[0]["ambiguous"] is True
    assert rows[1]["kinase"] == "K1"
    assert rows[1]["substrate_site"] == "S2"
    assert rows[1]["substrate_identifier"] == "D_SHARED"
    assert pd.isna(rows[1]["value_used_in_scoring"])
    assert rows[1]["score_component"] == "rank_weighted_fusion_scores"
    assert (
        rows[1]["score_source"]
        == KINASE_SCORE_SOURCE_PROFILE_ONLY_MOTIF_MISSING_OR_CONSTANT
    )
    assert rows[1]["reference_source_name"] == "fixture_reference"
    assert rows[1]["reference_source_version"] == "v1"
    assert rows[1]["reference_bundle_id"] == "fixture_bundle"
    assert rows[1]["reference_identifier_namespace"] == "fixture_display_id"
    assert rows[1]["status"] == KINASE_SUBSTRATE_CONTRIBUTION_STATUS_EXCLUDED
    assert (
        rows[1]["exclusion_reason"]
        == KINASE_SUBSTRATE_CONTRIBUTION_EXCLUDED_MISSING_SCORE_VALUE
    )
    assert rows[1]["ambiguous"] is True
    assert rows[2]["kinase"] == "K2"
    assert rows[2]["substrate_site"] == "S3"
    assert rows[2]["substrate_identifier"] == "D3"
    assert pd.isna(rows[2]["value_used_in_scoring"])
    assert rows[2]["score_component"] == "rank_weighted_fusion_scores"
    assert pd.isna(rows[2]["score_source"])
    assert rows[2]["reference_source_name"] == "fixture_reference"
    assert rows[2]["reference_source_version"] == "v1"
    assert rows[2]["reference_bundle_id"] == "fixture_bundle"
    assert rows[2]["reference_identifier_namespace"] == "fixture_display_id"
    assert rows[2]["status"] == KINASE_SUBSTRATE_CONTRIBUTION_STATUS_EXCLUDED
    assert (
        rows[2]["exclusion_reason"]
        == KINASE_SUBSTRATE_CONTRIBUTION_EXCLUDED_BELOW_MIN_SUBSTRATES
    )
    assert rows[2]["ambiguous"] is False
    assert rows[3]["kinase"] == "K3"
    assert rows[3]["substrate_site"] == "S_MISSING"
    assert rows[3]["substrate_identifier"] == "DMISS"
    assert pd.isna(rows[3]["value_used_in_scoring"])
    assert rows[3]["score_component"] == "rank_weighted_fusion_scores"
    assert pd.isna(rows[3]["score_source"])
    assert rows[3]["reference_source_name"] == "fixture_reference"
    assert rows[3]["reference_source_version"] == "v1"
    assert rows[3]["reference_bundle_id"] == "fixture_bundle"
    assert rows[3]["reference_identifier_namespace"] == "fixture_display_id"
    assert rows[3]["status"] == KINASE_SUBSTRATE_CONTRIBUTION_STATUS_EXCLUDED
    assert (
        rows[3]["exclusion_reason"]
        == KINASE_SUBSTRATE_CONTRIBUTION_EXCLUDED_NOT_QUANTIFIED
    )
    assert rows[3]["ambiguous"] is False


def test_scoring_runner_skips_contribution_builder_when_collection_disabled(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def fail_if_called(**kwargs: object) -> pd.DataFrame:
        raise AssertionError("substrate contribution table should not be built")

    monkeypatch.setattr(
        scoring_runner_module,
        "build_kinase_substrate_contribution_table",
        fail_if_called,
    )
    resolved = KinaseWorkflowInterpreter().run(
        KinaseWorkflowRequest(
            dataset=_dataset(),
            references=_references(),
            scoring_config=KinaseScoringConfig(min_substrates=2),
            prediction_config=KinasePredictionConfig(
                top_k=2,
                deterministic_max_selected_kinases=2,
                adaptive_ensemble_runs=2,
            ),
            activity_config=None,
        )
    )

    scoring_execution = KinaseScoringRunner().run(
        request=resolved,
        config=resolved.execution_config,
        collect_substrate_contributions=False,
    )

    assert scoring_execution.substrate_contributions is None
    assert not scoring_execution.scoring_result.profile_scores.empty
    assert scoring_execution.scoring_result.rank_weighted_fusion_scores is not None


def test_scoring_runner_collects_internal_contributions_when_requested() -> None:
    dataset = _dataset()
    references = _references()
    request = KinaseWorkflowRequest(
        dataset=dataset,
        references=references,
        scoring_config=KinaseScoringConfig(
            min_substrates=2,
            include_diagnostic_scoring_tables=True,
            include_substrate_contributions=True,
        ),
        prediction_config=KinasePredictionConfig(
            top_k=2,
            deterministic_max_selected_kinases=2,
            adaptive_ensemble_runs=2,
        ),
        activity_config=None,
    )
    resolved = KinaseWorkflowInterpreter().run(request)

    scoring_execution = KinaseScoringRunner().run(
        request=resolved,
        config=resolved.execution_config,
        collect_substrate_contributions=True,
    )

    assert scoring_execution.substrate_contributions is not None
    assert tuple(scoring_execution.substrate_contributions.columns) == (
        KINASE_SUBSTRATE_CONTRIBUTION_COLUMNS
    )
    assert scoring_execution.substrate_contributions.loc[:, "status"].tolist() == [
        KINASE_SUBSTRATE_CONTRIBUTION_STATUS_INCLUDED,
        KINASE_SUBSTRATE_CONTRIBUTION_STATUS_INCLUDED,
    ]
    assert (
        scoring_execution.substrate_contributions.loc[:, "score_source"].notna().all()
    )
    assert not hasattr(scoring_execution.scoring_result, "substrate_contributions")
    assert not scoring_execution.scoring_result.profile_scores.empty
    assert scoring_execution.scoring_result.rank_weighted_fusion_scores is not None


def test_kinase_workflow_skips_substrate_contributions_by_default(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def fail_if_called(**kwargs: object) -> pd.DataFrame:
        raise AssertionError("substrate contribution table should not be built")

    monkeypatch.setattr(
        scoring_runner_module,
        "build_kinase_substrate_contribution_table",
        fail_if_called,
    )

    result = KinaseWorkflow().run(
        _contribution_request(include_substrate_contributions=False)
    )

    assert result.substrate_contributions is None


def test_kinase_workflow_exposes_substrate_contributions_when_enabled() -> None:
    default_result = KinaseWorkflow().run(
        _contribution_request(include_substrate_contributions=False)
    )
    result = KinaseWorkflow().run(
        _contribution_request(include_substrate_contributions=True)
    )

    table = result.substrate_contributions
    assert table is not None
    assert tuple(table.columns) == KINASE_SUBSTRATE_CONTRIBUTION_COLUMNS
    assert table.loc[:, "status"].tolist() == [
        KINASE_SUBSTRATE_CONTRIBUTION_STATUS_INCLUDED,
        KINASE_SUBSTRATE_CONTRIBUTION_STATUS_INCLUDED,
        KINASE_SUBSTRATE_CONTRIBUTION_STATUS_EXCLUDED,
    ]
    included = table.loc[
        table.loc[:, "status"] == KINASE_SUBSTRATE_CONTRIBUTION_STATUS_INCLUDED,
        :,
    ]
    assert included.loc[:, "kinase"].tolist() == ["MAP2K6", "MAP2K6"]
    assert included.loc[:, "substrate_identifier"].tolist() == [
        "MAPK14;Y182;",
        "MAPK14;Y182;",
    ]
    assert included.loc[:, "ambiguous"].tolist() == [True, True]

    excluded = table.loc[
        table.loc[:, "status"] == KINASE_SUBSTRATE_CONTRIBUTION_STATUS_EXCLUDED,
        :,
    ]
    assert excluded.loc[:, "kinase"].tolist() == ["LOW_SUPPORT_KINASE"]
    assert excluded.loc[:, "substrate_identifier"].tolist() == ["GSK3B;S9;"]
    assert excluded.loc[:, "exclusion_reason"].tolist() == [
        KINASE_SUBSTRATE_CONTRIBUTION_EXCLUDED_BELOW_MIN_SUBSTRATES
    ]
    assert excluded.loc[:, "ambiguous"].tolist() == [False]

    assert default_result.substrate_contributions is None
    pd.testing.assert_frame_equal(
        default_result.scoring_result.profile_scores,
        result.scoring_result.profile_scores,
    )
    pd.testing.assert_frame_equal(
        default_result.scoring_result.authoritative_scores,
        result.scoring_result.authoritative_scores,
    )
    pd.testing.assert_frame_equal(
        default_result.prediction_result.pred_mat,
        result.prediction_result.pred_mat,
    )
    assert default_result.prediction_result.substrate_list is not None
    assert result.prediction_result.substrate_list is not None
    pd.testing.assert_frame_equal(
        default_result.prediction_result.substrate_list,
        result.prediction_result.substrate_list,
    )


def test_substrate_contribution_flag_does_not_change_activity_outputs() -> None:
    activity_config = KinaseActivityConfig(
        enabled=True,
        threshold=0.0,
        min_substrates=2,
        top_n_substrates=3,
    )
    default_result = KinaseWorkflow().run(
        _contribution_request(
            include_substrate_contributions=False,
            activity_config=activity_config,
        )
    )
    result = KinaseWorkflow().run(
        _contribution_request(
            include_substrate_contributions=True,
            activity_config=activity_config,
        )
    )

    assert default_result.activity_result is not None
    assert result.activity_result is not None
    pd.testing.assert_frame_equal(
        default_result.activity_result.activity_matrix,
        result.activity_result.activity_matrix,
    )
    pd.testing.assert_frame_equal(
        default_result.activity_result.substrate_count_matrix,
        result.activity_result.substrate_count_matrix,
    )
    pd.testing.assert_frame_equal(
        default_result.activity_result.target_table,
        result.activity_result.target_table,
    )


def _dataset() -> AnalysisReadyPhosphoDataset:
    display_ids = ["MAPK14;Y182;", "GSK3B;S9;"]
    site_index = site_key_index_from_display_ids(display_ids)
    return AnalysisReadyPhosphoDataset(
        phospho=pd.DataFrame(
            {"sample_a": [1.0, 2.0], "sample_b": [2.0, 3.0]},
            index=site_index.copy(),
        ),
        site_metadata=pd.DataFrame(
            {
                "site_key": site_index.astype(str).tolist(),
                "display_id": display_ids,
                **site_key_context_columns(site_index),
                "gene_symbol": ["MAPK14", "GSK3B"],
                "protein_id": ["MAPK14", "GSK3B"],
                "site": ["Y182", "S9"],
                "site_sequence": [
                    "AAAAAAAAAAAAAAAYAAAAAAAAAAAAAAA",
                    "AAAAAAAAAAAAAAASAAAAAAAAAAAAAAA",
                ],
            },
            index=site_index.copy(),
        ),
        organism=Organism.RAT,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False,
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
    )


def _references() -> ReferenceBundle:
    return ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["MAP2K6", "MAP2K6"],
                "substrate_site": ["MAPK14;Y182;", "GSK3B;S9;"],
            }
        ),
        site_sequences=pd.DataFrame(
            {
                "site_sequence": [
                    "AAAAAAAAAAAAAAAYAAAAAAAAAAAAAAA",
                    "AAAAAAAAAAAAAAASAAAAAAAAAAAAAAA",
                ]
            },
            index=pd.Index(["MAPK14;Y182;", "GSK3B;S9;"], name="site_id"),
        ),
    )


def _contribution_request(
    *,
    include_substrate_contributions: bool,
    activity_config: KinaseActivityConfig | None = None,
) -> KinaseWorkflowRequest:
    return KinaseWorkflowRequest(
        dataset=_ambiguous_dataset(),
        references=_ambiguous_references(),
        scoring_config=KinaseScoringConfig(
            min_substrates=2,
            include_substrate_contributions=include_substrate_contributions,
        ),
        prediction_config=KinasePredictionConfig(
            top_k=3,
            deterministic_max_selected_kinases=3,
            adaptive_ensemble_runs=2,
        ),
        activity_config=activity_config,
        reference_display_ambiguity_policy=(
            KINASE_REFERENCE_DISPLAY_AMBIGUITY_POLICY_ALLOW_WITH_DIAGNOSTICS
        ),
    )


def _ambiguous_dataset() -> AnalysisReadyPhosphoDataset:
    display_ids = ["MAPK14;Y182;", "MAPK14;Y182;", "GSK3B;S9;"]
    site_index = protein_site_key_index(
        protein_identifiers=["MAPK14_A", "MAPK14_B", "GSK3B"],
        sites=["Y182", "Y182", "S9"],
    )
    return AnalysisReadyPhosphoDataset(
        phospho=pd.DataFrame(
            {
                "sample_a": [1.0, 1.2, 2.0],
                "sample_b": [2.0, 2.2, 3.0],
            },
            index=site_index.copy(),
        ),
        site_metadata=pd.DataFrame(
            {
                "site_key": site_index.astype(str).tolist(),
                "display_id": display_ids,
                **site_key_context_columns(site_index),
                "gene_symbol": ["MAPK14", "MAPK14", "GSK3B"],
                "protein_id": ["MAPK14_A", "MAPK14_B", "GSK3B"],
                "site": ["Y182", "Y182", "S9"],
                "site_sequence": [
                    "AAAAAAAAAAAAAAAYAAAAAAAAAAAAAAA",
                    "AAAAAAAAAAAAAAAYAAAAAAAAAAAAAAA",
                    "AAAAAAAAAAAAAAASAAAAAAAAAAAAAAA",
                ],
            },
            index=site_index.copy(),
        ),
        organism=Organism.RAT,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False,
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
    )


def _ambiguous_references() -> ReferenceBundle:
    return ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["MAP2K6", "LOW_SUPPORT_KINASE"],
                "substrate_site": ["MAPK14;Y182;", "GSK3B;S9;"],
            }
        ),
        site_sequences=pd.DataFrame(
            {
                "site_sequence": [
                    "AAAAAAAAAAAAAAAYAAAAAAAAAAAAAAA",
                    "AAAAAAAAAAAAAAASAAAAAAAAAAAAAAA",
                ]
            },
            index=pd.Index(["MAPK14;Y182;", "GSK3B;S9;"], name="site_id"),
        ),
    )
