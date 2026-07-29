from __future__ import annotations

import pandas as pd
import pytest

from phospy.api import (
    KinasePredictionConfig,
    KinaseScoringConfig,
    KinaseWorkflowRequest,
    Organism,
    ReferenceBundle,
    ReferenceContextCompatibilityPolicy,
)
from phospy.api.configs import (
    KINASE_SCORING_MODE_KINASE_LIBRARY_CONTEXTUAL_MOTIF,
    KINASE_SCORING_MODE_PHOSR_RANK_WEIGHTED,
)
from phospy.errors import WorkflowBoundaryError, WorkflowValidationError
from phospy.provenance.hashing import fingerprint_table
from phospy.provenance.models import KinaseLibraryResourceProvenance
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.science.prediction.motif_scoring.models import AMINO_ACIDS
from phospy.science.references.kinase_library import (
    KinaseLibraryMatrix,
    KinaseLibraryResidueClass,
    KinaseLibraryResource,
)
from phospy.science.references.models import SequenceWindowDefinition
from phospy.workflows.kinase.kinase_library_scoring import (
    KINASE_LIBRARY_WORKFLOW_SCORE_SCALE,
)
from phospy.workflows.kinase.public import KinaseWorkflow
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

pytestmark = pytest.mark.integration


def _dataset() -> AnalysisReadyPhosphoDataset:
    display_ids = ["GENE1;S1;", "GENE2;T2;", "GENE3;Y3;"]
    site_index = site_key_index_from_display_ids(display_ids)
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 3.0, 5.0],
            "sample_b": [2.0, 4.0, 4.0],
        },
        index=site_index.copy(),
    )
    site_metadata = pd.DataFrame(
        {
            "site_key": site_index.tolist(),
            "display_id": display_ids,
            **site_key_context_columns(site_index),
            "gene_symbol": ["GENE1", "GENE2", "GENE3"],
            "protein_id": ["GENE1", "GENE2", "GENE3"],
            "site": ["S1", "T2", "Y3"],
            "site_sequence": ["ASA", "ATA", "AYA"],
            "localisation_confidence": [0.99, 0.98, 0.97],
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
    return ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["KPROFILE", "KPROFILE"],
                "substrate_site": ["GENE1;S1;", "GENE2;T2;"],
            }
        ),
        site_sequences=pd.DataFrame(
            {"site_sequence": ["ASA", "ATA", "AYA"]},
            index=pd.Index(
                ["GENE1;S1;", "GENE2;T2;", "GENE3;Y3;"],
                name="site_id",
            ),
        ),
    )


def _kinase_library_resource(
    *, organism: Organism = Organism.RAT
) -> KinaseLibraryResource:
    score_table = pd.DataFrame(
        0.0,
        index=pd.Index(AMINO_ACIDS, name="amino_acid"),
        columns=pd.Index([-1, 0, 1], name="position"),
    )
    score_table.loc["S", 0] = 2.0
    score_table.loc["T", 0] = 1.0
    matrix = KinaseLibraryMatrix(
        kinase="KLIB1",
        residue_class=KinaseLibraryResidueClass.SER_THR,
        score_table=score_table,
    )
    sequence_window = SequenceWindowDefinition(
        upstream_residues=1,
        downstream_residues=1,
        central_residue_required=True,
    )
    provenance = KinaseLibraryResourceProvenance(
        source_type="local",
        source_name="synthetic_kinase_library",
        source_version="test",
        license="test-only",
        score_scale="synthetic_raw_position_sum",
        organisms=(organism.value,),
        sequence_window=sequence_window.to_payload(),
        source_files={"kinase_library": {"path": "synthetic"}},
        table_fingerprints=(
            fingerprint_table(
                score_table,
                name="references.kinase_library.score_table.klib1.ser_thr",
            ),
        ),
        manifest={
            "resource_type": "kinase_library",
            "source_name": "synthetic_kinase_library",
            "source_version": "test",
            "score_scale": "synthetic_raw_position_sum",
            "organisms": (organism.value,),
            "sequence_window": sequence_window.to_payload(),
        },
    )
    return KinaseLibraryResource(
        matrices=(matrix,),
        source_name="synthetic_kinase_library",
        source_version="test",
        score_scale="synthetic_raw_position_sum",
        sequence_window=sequence_window,
        organisms=(organism.value,),
        license="test-only",
        provenance=provenance,
    )


def _request(
    *,
    scoring_mode: str = KINASE_SCORING_MODE_PHOSR_RANK_WEIGHTED,
    kinase_library_resource: KinaseLibraryResource | None = None,
) -> KinaseWorkflowRequest:
    return KinaseWorkflowRequest(
        dataset=_dataset(),
        references=_references(),
        scoring_config=KinaseScoringConfig(
            reliability_profile="custom",
            min_substrates=2,
            scoring_mode=scoring_mode,
            include_diagnostic_scoring_tables=True,
            reference_context_compatibility_policy=(
                ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT
            ),
        ),
        prediction_config=KinasePredictionConfig(
            top_k=2,
            deterministic_max_selected_kinases=2,
            adaptive_ensemble_runs=2,
            mode="deterministic_ranking",
        ),
        activity_config=None,
        kinase_library_resource=kinase_library_resource,
    )


def test_existing_default_scoring_mode_is_unchanged() -> None:
    result = KinaseWorkflow().run(_request())

    assert (
        KinaseScoringConfig.exploratory().scoring_mode
        == KINASE_SCORING_MODE_PHOSR_RANK_WEIGHTED
    )
    assert result.scoring_result.scoring_mode == KINASE_SCORING_MODE_PHOSR_RANK_WEIGHTED
    assert result.scoring_result.score_source == "rank_weighted_fusion_scores"
    assert result.scoring_result.kinase_library_motif_scores is None
    assert result.scoring_result.rank_weighted_fusion_scores is not None


def test_kinase_library_mode_requires_resource_in_validator() -> None:
    with pytest.raises(WorkflowValidationError, match="kinase_library_resource"):
        KinaseWorkflow().run(
            _request(scoring_mode=KINASE_SCORING_MODE_KINASE_LIBRARY_CONTEXTUAL_MOTIF)
        )


def test_kinase_library_resource_organism_mismatch_fails_in_interpreter() -> None:
    with pytest.raises(
        WorkflowBoundaryError,
        match="kinase_library_resource_organism",
    ):
        KinaseWorkflow().run(
            _request(
                scoring_mode=KINASE_SCORING_MODE_KINASE_LIBRARY_CONTEXTUAL_MOTIF,
                kinase_library_resource=_kinase_library_resource(
                    organism=Organism.MOUSE
                ),
            )
        )


def test_kinase_library_mode_exposes_scores_source_provenance_and_attrition() -> None:
    result = KinaseWorkflow().run(
        _request(
            scoring_mode=KINASE_SCORING_MODE_KINASE_LIBRARY_CONTEXTUAL_MOTIF,
            kinase_library_resource=_kinase_library_resource(),
        )
    )

    scores = result.scoring_result.kinase_library_motif_scores
    assert scores is not None
    assert scores.shape == (3, 1)
    assert scores.columns.astype(str).tolist() == ["KLIB1"]
    assert scores.index.astype(str).tolist() == result.dataset.phospho.index.tolist()
    assert result.scoring_result.score_source == "kinase_library_motif_scores"
    assert result.scoring_result.score_scale == KINASE_LIBRARY_WORKFLOW_SCORE_SCALE
    assert result.prediction_result.pred_mat.columns.astype(str).tolist() == ["KLIB1"]

    site_diagnostics = result.scoring_result.kinase_library_site_diagnostics
    assert site_diagnostics is not None
    assert site_diagnostics.loc[:, "status"].value_counts().to_dict() == {
        "valid_scored_site": 2,
        "wrong_residue_class": 1,
    }

    assert result.site_attrition_summary is not None
    scoring_attrition = result.site_attrition_summary.scoring
    assert scoring_attrition.final_quantitative_sites_entering_scoring == 3
    assert scoring_attrition.sites_with_valid_site_sequence == 2
    assert scoring_attrition.sites_without_usable_site_sequence == 1
    assert scoring_attrition.sites_eligible_for_motif_scoring == 2

    assert result.provenance is not None
    scoring_config = result.provenance.workflow_parameters["scoring_config"]
    assert (
        scoring_config["scoring_mode"]
        == KINASE_SCORING_MODE_KINASE_LIBRARY_CONTEXTUAL_MOTIF
    )
    assert scoring_config["score_source"] == "kinase_library_motif_scores"
    assert scoring_config["score_scale"] == KINASE_LIBRARY_WORKFLOW_SCORE_SCALE
    policy_ids = {policy.id.value for policy in result.provenance.scientific_policies}
    assert "kinase_library_motif_scoring_v1" in policy_ids
