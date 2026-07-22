from __future__ import annotations

from typing import NoReturn

import pandas as pd
import pytest

from phospy.api import (
    KinasePredictionConfig,
    KinaseScoringConfig,
    KinaseWorkflow,
    KinaseWorkflowRequest,
    Organism,
    ReferenceBundle,
    ReferenceContextCompatibilityPolicy,
)
from phospy.api.configs import (
    KINASE_SCORING_MODE_KINASE_LIBRARY_CONTEXTUAL_MOTIF,
)
from phospy.errors import WorkflowValidationError
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
from phospy.workflows.kinase.contracts import ResolvedKinaseWorkflowRequest
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


class _ExecutorMustNotRun:
    def run(self, request: ResolvedKinaseWorkflowRequest) -> NoReturn:
        _ = request
        raise AssertionError("invalid sequence context reached kinase executor")


def _window(residue: str, *, flank: int = 7) -> str:
    return ("A" * flank) + residue + ("A" * flank)


def _dataset() -> AnalysisReadyPhosphoDataset:
    display_ids = ["GENE1;S10;"]
    site_index = site_key_index_from_display_ids(display_ids)
    return trusted_analysis_ready_dataset_from_tables(
        phospho=pd.DataFrame({"sample_a": [1.0]}, index=site_index.copy()),
        site_metadata=pd.DataFrame(
            {
                "site_key": site_index.astype(str).tolist(),
                "display_id": display_ids,
                **site_key_context_columns(site_index),
                "gene_symbol": ["GENE1"],
                "site": ["S10"],
                "protein_id": ["GENE1"],
                "site_sequence": [_window("S")],
            },
            index=site_index.copy(),
        ),
        organism=Organism.RAT,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
    )


def _references_with_invalid_window() -> ReferenceBundle:
    return ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {"kinase": ["KLIB_ST"], "substrate_site": ["GENE1;S10;"]}
        ),
        site_sequences=pd.DataFrame(
            {"site_sequence": [_window("S", flank=15)]},
            index=pd.Index(["GENE1;S10;"], name="site_id"),
        ),
    )


def _kinase_library_resource() -> KinaseLibraryResource:
    positions = tuple(range(-7, 8))
    score_table = pd.DataFrame(
        0.0,
        index=pd.Index(AMINO_ACIDS, name="amino_acid"),
        columns=pd.Index(positions, name="position"),
    )
    score_table.loc["S", 0] = 1.0
    matrix = KinaseLibraryMatrix(
        kinase="KLIB_ST",
        residue_class=KinaseLibraryResidueClass.SER_THR,
        score_table=score_table,
    )
    sequence_window = SequenceWindowDefinition(
        upstream_residues=7,
        downstream_residues=7,
        central_residue_required=True,
    )
    provenance = KinaseLibraryResourceProvenance(
        source_type="local",
        source_name="synthetic_kinase_library",
        source_version="test",
        license="test-only",
        score_scale="synthetic_raw_position_sum",
        organisms=(Organism.RAT.value,),
        sequence_window=sequence_window.to_payload(),
        source_files={"kinase_library": {"path": "synthetic"}},
        table_fingerprints=(
            fingerprint_table(
                score_table,
                name="references.kinase_library.score_table.klib_st.ser_thr",
            ),
        ),
        manifest={
            "resource_type": "kinase_library",
            "source_name": "synthetic_kinase_library",
            "source_version": "test",
            "score_scale": "synthetic_raw_position_sum",
            "organisms": (Organism.RAT.value,),
            "sequence_window": sequence_window.to_payload(),
        },
    )
    return KinaseLibraryResource(
        matrices=(matrix,),
        source_name="synthetic_kinase_library",
        source_version="test",
        score_scale="synthetic_raw_position_sum",
        sequence_window=sequence_window,
        organisms=(Organism.RAT.value,),
        license="test-only",
        provenance=provenance,
    )


def test_invalid_fixed_window_sequence_is_rejected_before_scoring() -> None:
    request = KinaseWorkflowRequest(
        dataset=_dataset(),
        references=_references_with_invalid_window(),
        scoring_config=KinaseScoringConfig(
            min_substrates=2,
            scoring_mode=KINASE_SCORING_MODE_KINASE_LIBRARY_CONTEXTUAL_MOTIF,
            reference_context_compatibility_policy=(
                ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT
            ),
        ),
        prediction_config=KinasePredictionConfig(
            top_k=1,
            deterministic_max_selected_kinases=1,
            adaptive_ensemble_runs=1,
            mode="deterministic_ranking",
        ),
        activity_config=None,
        kinase_library_resource=_kinase_library_resource(),
    )

    with pytest.raises(WorkflowValidationError) as exc_info:
        KinaseWorkflow._with_components(executor=_ExecutorMustNotRun()).run(request)

    message = str(exc_info.value)
    assert "workflow-specific sequence context contract failed" in message
    assert "expected_length=15" in message
    assert "observed_length=31" in message
