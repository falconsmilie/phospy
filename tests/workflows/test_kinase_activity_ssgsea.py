from __future__ import annotations

import pandas as pd
import pytest

from phospy import AnalysisReadyPhosphoDataset, KinaseWorkflow
from phospy.api import Organism, ReferenceBundle
from phospy.api.configs import (
    KINASE_SCORING_MODE_KINASE_LIBRARY_MOTIF_ONLY,
    KinaseActivityConfig,
    KinasePredictionConfig,
    KinaseScoringConfig,
    ReferenceContextCompatibilityPolicy,
)
from phospy.api.requests import KinaseWorkflowRequest
from phospy.provenance.hashing import fingerprint_table
from phospy.provenance.models import KinaseLibraryResourceProvenance
from phospy.provenance.scientific_policy_models import ScientificPolicyId
from phospy.science.activities.methods import SSGSEA_SIGNIFICANCE_STATUS_AVAILABLE
from phospy.science.prediction.motif_scoring.models import AMINO_ACIDS
from phospy.science.references.kinase_library import (
    KinaseLibraryMatrix,
    KinaseLibraryResidueClass,
    KinaseLibraryResource,
)
from phospy.science.references.models import SequenceWindowDefinition
from phospy.science.transformations.models import QuantitativeMeaning
from tests.support.analysis_ready_dataset_factories import (
    trusted_analysis_ready_dataset_from_tables,
)
from tests.support.intensity_scale_states import (
    supported_log2_intensity_scale_state_with_meaning,
    supported_log2_processing_state_with_meaning,
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
            "site_sequence": [
                _sequence_with_left_flank(left_flank)
                for left_flank in ("A", "R", "A", "R")
            ],
            "protein_id": ["S1", "S2", "S3", "S4"],
            "localisation_confidence": [0.95, 0.95, 0.95, 0.95],
        },
        index=site_index.copy(),
    )
    return trusted_analysis_ready_dataset_from_tables(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        intensity_scale_state=supported_log2_intensity_scale_state_with_meaning(
            has_total_matrix=False,
            meaning=QuantitativeMeaning.DIFFERENTIAL_EFFECT_SIZE,
        ),
        processing_state=supported_log2_processing_state_with_meaning(
            has_total_matrix=False,
            meaning=QuantitativeMeaning.DIFFERENTIAL_EFFECT_SIZE,
        ),
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
            {
                "site_sequence": [
                    _sequence_with_left_flank(left_flank)
                    for left_flank in ("A", "R", "A", "R")
                ]
            },
            index=pd.Index(display_ids, name="site_id"),
        ),
    )


def _sequence_with_left_flank(left_flank: str) -> str:
    return ("A" * 14) + left_flank + "S" + ("A" * 15)


def _kinase_library_resource() -> KinaseLibraryResource:
    positions = tuple(range(-15, 16))
    matrices: list[KinaseLibraryMatrix] = []
    fingerprints = []
    for kinase, score in (("K_TOP", 2.0), ("K_BOTTOM", 1.0)):
        score_table = pd.DataFrame(
            0.0,
            index=pd.Index(AMINO_ACIDS, name="amino_acid"),
            columns=pd.Index(positions, name="position"),
        )
        score_table.loc["S", 0] = score
        score_table.loc["A", -1] = score
        matrices.append(
            KinaseLibraryMatrix(
                kinase=kinase,
                residue_class=KinaseLibraryResidueClass.SER_THR,
                score_table=score_table,
            )
        )
        fingerprints.append(
            fingerprint_table(
                score_table,
                name=f"references.kinase_library.score_table.{kinase.lower()}",
            )
        )
    sequence_window = SequenceWindowDefinition(
        upstream_residues=15,
        downstream_residues=15,
        central_residue_required=True,
    )
    provenance = KinaseLibraryResourceProvenance(
        source_type="local",
        source_name="synthetic_ssgsea_test_kinase_library",
        source_version="test",
        license="test-only",
        score_scale="synthetic_raw_position_sum",
        organisms=(Organism.RAT.value,),
        sequence_window=sequence_window.to_payload(),
        source_files={"kinase_library": {"path": "synthetic"}},
        table_fingerprints=tuple(fingerprints),
        manifest={
            "resource_type": "kinase_library",
            "source_name": "synthetic_ssgsea_test_kinase_library",
            "source_version": "test",
            "score_scale": "synthetic_raw_position_sum",
            "organisms": (Organism.RAT.value,),
            "sequence_window": sequence_window.to_payload(),
        },
    )
    return KinaseLibraryResource(
        matrices=tuple(matrices),
        source_name="synthetic_ssgsea_test_kinase_library",
        source_version="test",
        score_scale="synthetic_raw_position_sum",
        sequence_window=sequence_window,
        organisms=(Organism.RAT.value,),
        license="test-only",
        provenance=provenance,
    )


def test_kinase_workflow_runs_ssgsea_substrate_enrichment_activity() -> None:
    result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=_dataset(),
            references=_references(),
            scoring_config=KinaseScoringConfig(
                reliability_profile="custom",
                min_substrates=2,
                scoring_mode=KINASE_SCORING_MODE_KINASE_LIBRARY_MOTIF_ONLY,
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
            kinase_library_resource=_kinase_library_resource(),
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
