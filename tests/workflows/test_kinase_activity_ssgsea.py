from __future__ import annotations

import pandas as pd
import pytest

from phospy import AnalysisReadyPhosphoDataset, KinaseWorkflow
from phospy.api import Organism, ReferenceBundle
from phospy.api.configs import (
    KinaseActivityConfig,
    KinasePredictionConfig,
    KinaseScoringConfig,
    ReferenceContextCompatibilityPolicy,
)
from phospy.api.requests import KinaseWorkflowRequest
from phospy.provenance.scientific_policy_models import ScientificPolicyId
from phospy.science.activities.methods import SSGSEA_SIGNIFICANCE_STATUS_AVAILABLE
from tests.support.analysis_ready_dataset_factories import (
    trusted_analysis_ready_dataset_from_tables,
)
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
)
from tests.support.site_keys import (
    site_key_context_columns,
    site_key_index_from_display_ids,
)


def _dataset() -> AnalysisReadyPhosphoDataset:
    display_ids = ["S1;S1;", "S2;S2;", "S3;S3;", "S4;S4;"]
    site_index = site_key_index_from_display_ids(display_ids)
    phospho = pd.DataFrame(
        {
            "condition_positive": [4.0, 3.0, 2.0, 1.0],
            "condition_negative": [1.0, 2.0, 3.0, 4.0],
        },
        index=site_index.copy(),
    )
    site_metadata = pd.DataFrame(
        {
            "site_key": site_index.astype(str).tolist(),
            "display_id": display_ids,
            **site_key_context_columns(site_index),
            "gene_symbol": ["S1", "S2", "S3", "S4"],
            "site": ["S1", "S2", "S3", "S4"],
            "site_sequence": ["A" * 15 + "S" + "A" * 15 for _ in display_ids],
            "protein_id": ["S1", "S2", "S3", "S4"],
            "localisation_confidence": [0.95, 0.95, 0.95, 0.95],
        },
        index=site_index.copy(),
    )
    return trusted_analysis_ready_dataset_from_tables(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
    )


def _references() -> ReferenceBundle:
    display_ids = ["S1;S1;", "S2;S2;", "S3;S3;", "S4;S4;"]
    return ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["K_TOP", "K_TOP", "K_BOTTOM", "K_BOTTOM"],
                "substrate_site": display_ids,
            }
        ),
        site_sequences=pd.DataFrame(
            {"site_sequence": ["A" * 15 + "S" + "A" * 15 for _ in display_ids]},
            index=pd.Index(display_ids, name="site_id"),
        ),
    )


def test_kinase_workflow_runs_ssgsea_substrate_enrichment_activity() -> None:
    result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=_dataset(),
            references=_references(),
            scoring_config=KinaseScoringConfig(
                reliability_profile="custom",
                min_substrates=2,
                reference_context_compatibility_policy=(
                    ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT
                ),
            ),
            prediction_config=KinasePredictionConfig(
                top_k=4,
                deterministic_max_selected_kinases=2,
                adaptive_ensemble_runs=2,
            ),
            activity_config=KinaseActivityConfig.ssgsea_with_permutation_significance(
                ssgsea_min_substrates=2,
                permutations=12,
                random_seed=19,
            ),
        )
    )

    activity = result.activity_result
    assert activity is not None
    assert activity.activity_method.activity_method_id == (
        "ssgsea_substrate_enrichment_activity_v1"
    )
    assert activity.activity_matrix.at["K_TOP", "condition_positive"] == pytest.approx(
        0.5
    )
    assert activity.activity_matrix.at["K_TOP", "condition_negative"] == pytest.approx(
        -0.5
    )
    assert activity.p_value_matrix is not None
    assert activity.q_value_matrix is not None
    assert activity.statistics_table is not None
    assert set(activity.statistics_table["significance_status"]) == {
        SSGSEA_SIGNIFICANCE_STATUS_AVAILABLE
    }
    assert activity.substrate_count_matrix.at["K_TOP", "condition_positive"] == 2
    assert {"site_key", "display_id"} <= set(activity.target_table.columns)

    assert result.provenance is not None
    activity_payload = result.provenance.workflow_parameters["activity_config"]
    assert isinstance(activity_payload, dict)
    assert activity_payload["ssgsea_permutations"] == 12
    assert activity_payload["ssgsea_random_seed"] == 19
    assert activity_payload["ssgsea_significance_status"] == (
        SSGSEA_SIGNIFICANCE_STATUS_AVAILABLE
    )
    policy_ids = {policy.id for policy in result.provenance.scientific_policies}
    assert ScientificPolicyId.SSGSEA_SUBSTRATE_ENRICHMENT_ACTIVITY in policy_ids
