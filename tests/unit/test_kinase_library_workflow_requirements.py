from __future__ import annotations

import inspect

import pandas as pd
import pytest

from phospy.api import (
    KINASE_SCORING_MODE_KINASE_LIBRARY_MOTIF,
    KinasePredictionConfig,
    KinaseScoringConfig,
    KinaseWorkflowRequest,
    Organism,
    ReferenceBundle,
)
from phospy.errors import WorkflowBoundaryError, WorkflowValidationError
from phospy.provenance.hashing import fingerprint_table
from phospy.provenance.models import KinaseLibraryResourceProvenance
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.science.prediction.motif_scoring import (
    KinaseLibraryMotifScorer,
    score_kinase_library_motifs,
)
from phospy.science.prediction.motif_scoring.models import AMINO_ACIDS
from phospy.science.references.kinase_library import (
    KinaseLibraryMatrix,
    KinaseLibraryResidueClass,
    KinaseLibraryResource,
)
from phospy.science.references.models import SequenceWindowDefinition
from phospy.workflows.kinase.executor import KinaseWorkflowExecutor
from phospy.workflows.kinase.interpreter import KinaseWorkflowInterpreter
from phospy.workflows.kinase.kinase_library_scoring import (
    KINASE_LIBRARY_WORKFLOW_SCORE_SCALE,
    KinaseLibraryWorkflowScorer,
)
from phospy.workflows.kinase.public import KinaseWorkflow
from phospy.workflows.kinase.scoring_runner import KinaseScoringRunner
from phospy.workflows.kinase.site_sequence_support import (
    KinaseSiteSequenceSupportResult,
)
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
)
from tests.support.site_keys import (
    site_key_context_columns,
    site_key_index_from_display_ids,
)

_SEQUENCES_BY_DISPLAY_ID = {
    "GENE1;S1;": "ASA",
    "GENE2;T2;": "ATA",
    "GENE3;Y3;": "AYA",
    "OTHER;S1;": "ASA",
}


def test_kinase_library_workflow_scoring_succeeds_with_all_required_inputs() -> None:
    result = KinaseWorkflow().run(
        _request(
            scoring_mode=KINASE_SCORING_MODE_KINASE_LIBRARY_MOTIF,
            kinase_library_resource=_kinase_library_resource(),
        )
    )

    scores = result.scoring_result.kinase_library_motif_scores
    assert scores is not None
    assert result.scoring_result.profile_scores is not None
    assert result.scoring_result.rank_weighted_fusion_scores is None
    assert result.scoring_result.motif_scores is None
    assert result.scoring_result.score_source == "kinase_library_motif_scores"
    assert result.scoring_result.score_scale == KINASE_LIBRARY_WORKFLOW_SCORE_SCALE
    assert result.prediction_result.pred_mat.columns.astype(str).tolist() == ["KLIB_ST"]
    assert scores.columns.astype(str).tolist() == ["KLIB_ST"]


def test_kinase_library_workflow_mode_requires_local_resource() -> None:
    with pytest.raises(WorkflowValidationError, match="kinase_library_resource"):
        KinaseWorkflow().run(
            _request(scoring_mode=KINASE_SCORING_MODE_KINASE_LIBRARY_MOTIF)
        )


def test_kinase_library_workflow_mode_still_requires_eligible_reference_map() -> None:
    with pytest.raises(WorkflowBoundaryError) as exc_info:
        KinaseWorkflow().run(
            _request(
                references=_references(mapped_display_ids=("GENE1;S1;",)),
                scoring_mode=KINASE_SCORING_MODE_KINASE_LIBRARY_MOTIF,
                kinase_library_resource=_kinase_library_resource(),
            )
        )

    error = exc_info.value
    assert error.seam == "kinase.interpreter.eligible_kinases"
    assert error.details["scoring_config_min_substrates"] == 2
    assert error.details["max_quantified_sites_per_kinase"] == 1


def test_kinase_library_workflow_mode_fails_when_no_site_sequences_resolve() -> None:
    class _NoSequenceSupportBuilder:
        def run(self, **kwargs: object) -> KinaseSiteSequenceSupportResult:
            dataset = kwargs["dataset"]
            assert isinstance(dataset, pd.DataFrame)
            empty = pd.DataFrame(columns=["site_sequence", "display_id"], dtype=object)
            empty.index = pd.Index([], name=dataset.index.name)
            return KinaseSiteSequenceSupportResult(
                site_sequences=empty,
                dataset_sequences_added=0,
                dataset_sequences_missing=int(dataset.shape[0]),
                dataset_sequences_available=0,
                conflict_policy="prefer_reference",
                conflicts=(),
                display_reference_multi_matches=(),
            )

    interpreter = KinaseWorkflowInterpreter(
        site_sequence_support_builder=_NoSequenceSupportBuilder()
    )

    with pytest.raises(WorkflowBoundaryError) as exc_info:
        interpreter.run(
            _request(
                scoring_mode=KINASE_SCORING_MODE_KINASE_LIBRARY_MOTIF,
                kinase_library_resource=_kinase_library_resource(),
            )
        )

    error = exc_info.value
    assert error.seam == "kinase.interpreter.sequence_support"
    assert error.details["sequence_supported_sites"] == 0


def test_kinase_library_workflow_mode_fails_when_residue_class_resources_are_missing() -> (
    None
):
    dataset = _dataset(
        display_ids=("GENE1;S1;", "GENE2;T2;"),
        sequences=("ASA", "ATA"),
    )

    with pytest.raises(WorkflowBoundaryError) as exc_info:
        KinaseWorkflow().run(
            _request(
                dataset=dataset,
                references=_references(
                    mapped_display_ids=("GENE1;S1;", "GENE2;T2;"),
                    sequence_display_ids=("GENE1;S1;", "GENE2;T2;"),
                ),
                scoring_mode=KINASE_SCORING_MODE_KINASE_LIBRARY_MOTIF,
                kinase_library_resource=_kinase_library_resource(
                    kinase="KLIB_TYR",
                    residue_class=KinaseLibraryResidueClass.TYR,
                    center_scores={"Y": 5.0},
                ),
            )
        )

    error = exc_info.value
    assert error.seam == "kinase.interpreter.kinase_library_resource_usability"
    assert error.details["resource_residue_classes"] == ("tyr",)
    assert error.details["scoring_site_residue_classes"] == ("ser_thr",)


def test_kinase_library_workflow_mode_does_not_call_phosr_inspired_fallback() -> None:
    interpreted = KinaseWorkflowInterpreter().run(
        _request(
            scoring_mode=KINASE_SCORING_MODE_KINASE_LIBRARY_MOTIF,
            kinase_library_resource=_kinase_library_resource(),
        )
    )

    def _fail_if_phosr_inspired_motif_scorer_runs(**_: object) -> object:
        raise AssertionError("PhosR-inspired motif scorer must not run")

    scoring = KinaseScoringRunner(
        score_motifs=_fail_if_phosr_inspired_motif_scorer_runs
    ).run(
        request=interpreted,
        config=interpreted.execution_config,
    )

    assert scoring.downstream_score_source == "kinase_library_motif_scores"
    assert scoring.scoring_result.kinase_library_motif_scores is not None
    assert scoring.scoring_result.motif_scores is None
    assert scoring.scoring_result.rank_weighted_fusion_scores is None


def test_raw_science_scores_and_workflow_support_scores_remain_distinct() -> None:
    resource = _kinase_library_resource(center_scores={"S": 10.0, "T": 20.0})
    raw_result = score_kinase_library_motifs(
        site_sequences={"GENE1;S1;": "ASA", "GENE2;T2;": "ATA"},
        matrices=resource,
    )

    assert raw_result.raw_scores.at["GENE1;S1;", "KLIB_ST"] == pytest.approx(10.0)
    assert raw_result.raw_scores.at["GENE2;T2;", "KLIB_ST"] == pytest.approx(20.0)
    assert raw_result.score_scale_metadata.score_scale == ("synthetic_raw_position_sum")

    display_ids = ["GENE1;S1;", "GENE2;T2;"]
    site_index = site_key_index_from_display_ids(display_ids)
    workflow_result = KinaseLibraryWorkflowScorer().run(
        resource=resource,
        site_sequences=pd.Series(["ASA", "ATA"], index=site_index),
        site_identities=pd.Series(display_ids, index=site_index),
        site_index=site_index,
    )

    first_site, second_site = site_index.astype(str).tolist()
    assert workflow_result.raw_scores.at[first_site, "KLIB_ST"] == pytest.approx(10.0)
    assert workflow_result.raw_scores.at[second_site, "KLIB_ST"] == pytest.approx(20.0)
    assert workflow_result.scores.at[first_site, "KLIB_ST"] == pytest.approx(0.0)
    assert workflow_result.scores.at[second_site, "KLIB_ST"] == pytest.approx(1.0)
    assert workflow_result.score_scale_metadata["resource_score_scale"] == (
        "synthetic_raw_position_sum"
    )
    assert workflow_result.score_scale_metadata["workflow_score_scale"] == (
        KINASE_LIBRARY_WORKFLOW_SCORE_SCALE
    )


def test_kinase_library_responsibility_boundaries_remain_layered() -> None:
    interpreter_source = inspect.getsource(KinaseWorkflowInterpreter)
    executor_source = inspect.getsource(KinaseWorkflowExecutor)
    raw_scorer_source = inspect.getsource(KinaseLibraryMotifScorer)
    workflow_scorer_source = inspect.getsource(KinaseLibraryWorkflowScorer)

    assert "ReferenceResolver" in interpreter_source
    assert "BundledReferenceProvider" in interpreter_source
    assert "ReferenceResolver" not in executor_source
    assert "BundledReferenceProvider" not in executor_source
    assert "minmax_scale_columns" not in raw_scorer_source
    assert "minmax_scale_columns" in workflow_scorer_source
    assert "score_kinase_library_motifs" in workflow_scorer_source


def _dataset(
    *,
    display_ids: tuple[str, ...] = ("GENE1;S1;", "GENE2;T2;", "GENE3;Y3;"),
    sequences: tuple[str, ...] = ("ASA", "ATA", "AYA"),
) -> AnalysisReadyPhosphoDataset:
    site_index = site_key_index_from_display_ids(list(display_ids))
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0 + float(index) for index, _ in enumerate(display_ids)],
            "sample_b": [2.0 + float(index) for index, _ in enumerate(display_ids)],
        },
        index=site_index.copy(),
    )
    site_metadata = pd.DataFrame(
        {
            "site_key": site_index.astype(str).tolist(),
            "display_id": list(display_ids),
            **site_key_context_columns(site_index),
            "gene_symbol": [site.split(";", 1)[0] for site in display_ids],
            "protein_id": [site.split(";", 1)[0] for site in display_ids],
            "site": [site.split(";")[1] for site in display_ids],
            "site_sequence": list(sequences),
            "localisation_confidence": [0.99 for _ in display_ids],
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


def _references(
    *,
    mapped_display_ids: tuple[str, ...] = ("GENE1;S1;", "GENE2;T2;"),
    sequence_display_ids: tuple[str, ...] = (
        "GENE1;S1;",
        "GENE2;T2;",
        "GENE3;Y3;",
    ),
) -> ReferenceBundle:
    return ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["KPROFILE" for _ in mapped_display_ids],
                "substrate_site": list(mapped_display_ids),
            }
        ),
        site_sequences=pd.DataFrame(
            {
                "site_sequence": [
                    _SEQUENCES_BY_DISPLAY_ID[display_id]
                    for display_id in sequence_display_ids
                ]
            },
            index=pd.Index(list(sequence_display_ids), name="site_id"),
        ),
    )


def _kinase_library_resource(
    *,
    kinase: str = "KLIB_ST",
    residue_class: KinaseLibraryResidueClass = KinaseLibraryResidueClass.SER_THR,
    center_scores: dict[str, float] | None = None,
) -> KinaseLibraryResource:
    center_scores = center_scores or {"S": 2.0, "T": 1.0}
    score_table = pd.DataFrame(
        0.0,
        index=pd.Index(AMINO_ACIDS, name="amino_acid"),
        columns=pd.Index([-1, 0, 1], name="position"),
    )
    for amino_acid, score in center_scores.items():
        score_table.loc[amino_acid, 0] = float(score)
    matrix = KinaseLibraryMatrix(
        kinase=kinase,
        residue_class=residue_class,
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
        organisms=(Organism.RAT.value,),
        sequence_window=sequence_window.to_payload(),
        source_files={"kinase_library": {"path": "synthetic"}},
        table_fingerprints=(
            fingerprint_table(
                score_table,
                name=(
                    "references.kinase_library.score_table."
                    f"{kinase}.{residue_class.value}"
                ),
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


def _request(
    *,
    dataset: AnalysisReadyPhosphoDataset | None = None,
    references: ReferenceBundle | None = None,
    scoring_mode: str = KINASE_SCORING_MODE_KINASE_LIBRARY_MOTIF,
    kinase_library_resource: KinaseLibraryResource | None = None,
) -> KinaseWorkflowRequest:
    return KinaseWorkflowRequest(
        dataset=dataset or _dataset(),
        references=references or _references(),
        scoring_config=KinaseScoringConfig(
            min_substrates=2,
            scoring_mode=scoring_mode,
            include_diagnostic_scoring_tables=True,
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
