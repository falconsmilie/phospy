from __future__ import annotations

import pandas as pd

from phospy import AnalysisReadyPhosphoDataset
from phospy.api import (
    KinasePredictionConfig,
    KinaseScoringConfig,
    KinaseWorkflowRequest,
    Organism,
    ReferenceBundle,
)
from phospy.prediction.motif_scoring import (
    DEFAULT_MOTIF_FLANK_SIZE,
    build_motif_library_from_sequences,
    score_phosphosite_motifs,
)
from phospy.prediction.sequence_validation import (
    SEQUENCE_VALIDATION_STATUS_MISSING_SEQUENCE,
    SEQUENCE_VALIDATION_STATUS_OFF_CENTRE_SEQUENCE,
    SEQUENCE_VALIDATION_STATUS_SHORT_SEQUENCE,
    SEQUENCE_VALIDATION_STATUS_SITE_RESIDUE_MISMATCH,
    SEQUENCE_VALIDATION_STATUS_UNSUPPORTED_RESIDUE_CHARACTER,
    MotifSequenceValidator,
    SequenceValidationInput,
)
from phospy.workflows.kinase.executor import KinaseWorkflowExecutor
from phospy.workflows.kinase.interpreter import KinaseWorkflowInterpreter
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
)


def _validator() -> MotifSequenceValidator:
    return MotifSequenceValidator(
        expected_window_size=(2 * DEFAULT_MOTIF_FLANK_SIZE) + 1
    )


def test_valid_centred_site_sequence_passes_validation() -> None:
    result = _validator().run(
        rows=[
            SequenceValidationInput(
                site_id="MAPK1;S202;",
                site_sequence="AAAAAAASAAAAAAA",
            )
        ]
    )

    assert result.total_sequences == 1
    assert result.valid_sequences == 1
    assert result.invalid_sequences == 0


def test_missing_sequence_is_reported_and_excluded() -> None:
    result = _validator().run(
        rows=[
            SequenceValidationInput(
                site_id="MAPK1;S202;",
                site_sequence=None,
            )
        ]
    )

    assert result.invalid_sequences == 1
    assert result.rows[0].status == SEQUENCE_VALIDATION_STATUS_MISSING_SEQUENCE


def test_short_sequence_is_reported_and_excluded() -> None:
    result = _validator().run(
        rows=[
            SequenceValidationInput(
                site_id="MAPK1;S202;",
                site_sequence="ASAA",
            )
        ]
    )

    assert result.short_sequences == 1
    assert result.rows[0].status == SEQUENCE_VALIDATION_STATUS_SHORT_SEQUENCE


def test_off_centre_sequence_is_reported_and_excluded() -> None:
    result = _validator().run(
        rows=[
            SequenceValidationInput(
                site_id="MAPK1;S202;",
                site_sequence="SAAAAAAAAAAAAAA",
            )
        ]
    )

    assert result.off_centre_sequences == 1
    assert result.rows[0].status == SEQUENCE_VALIDATION_STATUS_OFF_CENTRE_SEQUENCE


def test_site_residue_mismatch_is_reported_and_excluded() -> None:
    result = _validator().run(
        rows=[
            SequenceValidationInput(
                site_id="MAPK1;S202;",
                site_sequence="AAAAAAATAAAAAAA",
            )
        ]
    )

    assert result.site_residue_mismatches == 1
    assert result.rows[0].status == SEQUENCE_VALIDATION_STATUS_SITE_RESIDUE_MISMATCH


def test_unsupported_residue_characters_are_reported_and_excluded() -> None:
    result = _validator().run(
        rows=[
            SequenceValidationInput(
                site_id="MAPK1;S202;",
                site_sequence="AAAAAAAXAAAAAAA",
            )
        ]
    )

    assert result.unsupported_residue_characters == 1
    assert (
        result.rows[0].status
        == SEQUENCE_VALIDATION_STATUS_UNSUPPORTED_RESIDUE_CHARACTER
    )


def test_motif_scoring_excludes_invalid_sequences_and_reports_diagnostics() -> None:
    motif_frequency_matrices, motif_sizes = build_motif_library_from_sequences(
        motif_sequences={
            "K1": [
                "AAAAAAASAAAAAAA",
                "VVVVVVVTVVVVVVV",
                "SSSSSSSSSSSSSSS",
            ]
        },
        flank_size=DEFAULT_MOTIF_FLANK_SIZE,
    )
    site_index = pd.Index(
        ["MAPK1;S202;", "MAPK1;T205;", "MAPK1;S210;"],
        name="site_id",
    )
    result = score_phosphosite_motifs(
        site_sequences={
            "MAPK1;S202;": "AAAAAAASAAAAAAA",
            "MAPK1;T205;": "VVVVVVVTVVVVVVV",
            "MAPK1;S210;": "AAAAAAAXAAAAAAA",
        },
        motif_frequency_matrices=motif_frequency_matrices,
        motif_sizes=motif_sizes,
        site_index=site_index,
        min_motif_size=1,
        flank_size=DEFAULT_MOTIF_FLANK_SIZE,
    )

    assert "MAPK1;S210;" in result.sequence_validation.excluded_site_ids
    invalid_row = next(
        row for row in result.sequence_validation.rows if row.site_id == "MAPK1;S210;"
    )
    assert (
        invalid_row.status == SEQUENCE_VALIDATION_STATUS_UNSUPPORTED_RESIDUE_CHARACTER
    )
    assert result.sequence_validation.sequences_excluded_from_motif_scoring == 1
    assert pd.isna(result.motif_scores.loc["MAPK1;S210;", "K1"])
    assert result.motif_scores.loc[["MAPK1;S202;", "MAPK1;T205;"], "K1"].notna().all()


def test_kinase_workflow_exposes_sequence_validation_diagnostics() -> None:
    site_ids = pd.Index(["MAPK1;S202;", "MAPK1;T205;"], name="site_id")
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 2.0],
            "sample_b": [2.0, 1.0],
        },
        index=site_ids,
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK1", "MAPK1"],
            "site": ["S202", "T205"],
            "site_sequence": ["A" * 31, "A" * 31],
        },
        index=site_ids,
    )
    dataset = AnalysisReadyPhosphoDataset(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
    )
    references = ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["KINASE_A", "KINASE_A"],
                "substrate_site": site_ids.tolist(),
            }
        ),
        site_sequences=pd.DataFrame(
            {
                "site_sequence": ["AAAAAAASAAAAAAA", "AAAAAAAXAAAAAAA"],
            },
            index=site_ids,
        ),
    )

    request = KinaseWorkflowRequest(
        dataset=dataset,
        references=references,
        scoring_config=KinaseScoringConfig(
            min_substrates=2,
            include_diagnostic_scoring_tables=True,
        ),
        prediction_config=KinasePredictionConfig(
            top_k=1,
            deterministic_max_selected_kinases=1,
            adaptive_ensemble_runs=1,
        ),
        activity_config=None,
    )
    interpreted = KinaseWorkflowInterpreter().run(request)
    scoring_execution = KinaseWorkflowExecutor()._run_scoring_stage(
        request=interpreted,
        config=interpreted.execution_config,
    )
    scoring_result = scoring_execution.scoring_result

    assert scoring_result.motif_sequence_validation is not None
    summary = scoring_result.motif_sequence_validation.summary()
    assert summary["total_sequences"] == 2
    assert summary["valid_sequences"] == 1
    assert summary["invalid_sequences"] == 1
    assert summary["short_sequences"] == 0
    assert summary["off_centre_sequences"] == 0
    assert summary["site_residue_mismatches"] == 0
    assert summary["unsupported_residue_characters"] == 1
    assert summary["sequences_excluded_from_motif_scoring"] == 1

    assert scoring_result.motif_scores is not None
    assert pd.isna(scoring_result.motif_scores.loc["MAPK1;T205;", "KINASE_A"])
