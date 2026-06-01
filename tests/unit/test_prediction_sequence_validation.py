from __future__ import annotations

import pandas as pd

from phospy import AnalysisReadyPhosphoDataset, KinaseWorkflow
from phospy.api import (
    KinasePredictionConfig,
    KinaseScoringConfig,
    KinaseWorkflowRequest,
    Organism,
    ReferenceBundle,
)
from phospy.science.prediction.motif_scoring import (
    DEFAULT_MOTIF_FLANK_SIZE,
    SEQUENCE_SEMANTICS_CENTRED_SEQUENCE,
    SEQUENCE_SEMANTICS_CENTRED_WINDOW,
    ExplicitMotifSequence,
    build_motif_library,
    build_motif_library_from_sequences,
    get_motif_library_validation,
    score_phosphosite_motifs,
)
from phospy.science.prediction.sequence_validation import (
    SEQUENCE_VALIDATION_STATUS_INVALID_SITE_ID,
    SEQUENCE_VALIDATION_STATUS_MISSING_SEQUENCE,
    SEQUENCE_VALIDATION_STATUS_NON_PHOSPHO_CENTRE_RESIDUE,
    SEQUENCE_VALIDATION_STATUS_OFF_CENTRE_SEQUENCE,
    SEQUENCE_VALIDATION_STATUS_SHORT_SEQUENCE,
    SEQUENCE_VALIDATION_STATUS_SITE_RESIDUE_MISMATCH,
    SEQUENCE_VALIDATION_STATUS_UNSUPPORTED_RESIDUE_CHARACTER,
    MotifSequenceValidator,
    SequenceValidationInput,
)
from phospy.workflows.kinase.interpreter import KinaseWorkflowInterpreter
from phospy.workflows.kinase.scoring_runner import KinaseScoringRunner
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
)
from tests.support.site_keys import (
    site_key_from_display_id,
    site_key_index_from_display_ids,
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


def test_valid_centred_t_and_y_site_sequences_pass_validation() -> None:
    result = _validator().run(
        rows=[
            SequenceValidationInput(
                site_id="MAPK1;T205;",
                site_sequence="AAAAAAATAAAAAAA",
            ),
            SequenceValidationInput(
                site_id="MAPK1;Y207;",
                site_sequence="AAAAAAAYAAAAAAA",
            ),
        ]
    )

    assert result.total_sequences == 2
    assert result.valid_sequences == 2
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


def test_even_length_sequence_is_reported_as_off_centre_and_excluded() -> None:
    result = _validator().run(
        rows=[
            SequenceValidationInput(
                site_id="MAPK1;S202;",
                site_sequence="AAAAAAASAAAAAAAA",
            )
        ]
    )

    assert result.off_centre_sequences == 1
    assert result.rows[0].status == SEQUENCE_VALIDATION_STATUS_OFF_CENTRE_SEQUENCE


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
    coverage = result.sequence_validation.site_sequence_coverage_summary()
    assert coverage["total_sites_considered"] == 3
    assert coverage["sites_with_valid_site_sequence"] == 2
    assert coverage["sites_without_valid_site_sequence"] == 1
    assert coverage["site_sequence_coverage_fraction"] == 2 / 3
    assert coverage["sites_used_for_motif_scoring"] == 2
    assert coverage["sites_excluded_from_motif_scoring_due_to_sequence"] == 1
    assert pd.isna(result.motif_scores.loc["MAPK1;S210;", "K1"])
    assert result.motif_scores.loc[["MAPK1;S202;", "MAPK1;T205;"], "K1"].notna().all()


def test_motif_scoring_requires_exact_centred_windows_by_default() -> None:
    motif_frequency_matrices, motif_sizes = build_motif_library_from_sequences(
        motif_sequences={"K1": ["AAAAAAASAAAAAAA"]},
        flank_size=DEFAULT_MOTIF_FLANK_SIZE,
    )
    result = score_phosphosite_motifs(
        site_sequences={"MAPK1;S202;": "AAAAAAASAAAAAAASAAAAAAASAAAAAAA"},
        motif_frequency_matrices=motif_frequency_matrices,
        motif_sizes=motif_sizes,
        site_index=pd.Index(["MAPK1;S202;"], name="site_id"),
        min_motif_size=1,
        flank_size=DEFAULT_MOTIF_FLANK_SIZE,
        sequence_semantics=SEQUENCE_SEMANTICS_CENTRED_WINDOW,
    )

    row = result.sequence_validation.rows[0]
    assert row.status == SEQUENCE_VALIDATION_STATUS_OFF_CENTRE_SEQUENCE
    assert "do not supply full-protein-like sequences" in str(row.reason)
    assert pd.isna(result.motif_scores.loc["MAPK1;S202;", "K1"])


def test_motif_scoring_can_explicitly_extract_window_from_centred_long_sequence() -> (
    None
):
    motif_frequency_matrices, motif_sizes = build_motif_library_from_sequences(
        motif_sequences={"K1": ["AAAAAAASAAAAAAA"]},
        flank_size=DEFAULT_MOTIF_FLANK_SIZE,
    )
    result = score_phosphosite_motifs(
        site_sequences={"MAPK1;S202;": "GGGGGGGAAAAAAASAAAAAAAGGGGGGG"},
        motif_frequency_matrices=motif_frequency_matrices,
        motif_sizes=motif_sizes,
        site_index=pd.Index(["MAPK1;S202;"], name="site_id"),
        min_motif_size=1,
        flank_size=DEFAULT_MOTIF_FLANK_SIZE,
        sequence_semantics=SEQUENCE_SEMANTICS_CENTRED_SEQUENCE,
    )

    assert result.sequence_validation.valid_sequences == 1
    assert result.sequence_windows.loc["MAPK1;S202;"] == "AAAAAAASAAAAAAA"
    assert "MAPK1;S202;" not in result.sequence_validation.excluded_site_ids


def test_centre_extraction_mode_ignores_unsupported_tail_characters_outside_window() -> (
    None
):
    motif_frequency_matrices, motif_sizes = build_motif_library_from_sequences(
        motif_sequences={"K1": ["AAAAAAASAAAAAAA"]},
        flank_size=DEFAULT_MOTIF_FLANK_SIZE,
    )
    result = score_phosphosite_motifs(
        site_sequences={"MAPK1;S202;": "__AAAAAAASAAAAAAA__"},
        motif_frequency_matrices=motif_frequency_matrices,
        motif_sizes=motif_sizes,
        site_index=pd.Index(["MAPK1;S202;"], name="site_id"),
        min_motif_size=1,
        flank_size=DEFAULT_MOTIF_FLANK_SIZE,
        sequence_semantics=SEQUENCE_SEMANTICS_CENTRED_SEQUENCE,
    )

    assert result.sequence_validation.valid_sequences == 1
    assert result.sequence_windows.loc["MAPK1;S202;"] == "AAAAAAASAAAAAAA"
    assert "MAPK1;S202;" not in result.sequence_validation.excluded_site_ids


def test_motif_scoring_rejects_long_off_centre_sequence_even_with_centre_extraction_mode() -> (
    None
):
    motif_frequency_matrices, motif_sizes = build_motif_library_from_sequences(
        motif_sequences={"K1": ["AAAAAAASAAAAAAA"]},
        flank_size=DEFAULT_MOTIF_FLANK_SIZE,
    )
    result = score_phosphosite_motifs(
        site_sequences={"MAPK1;S202;": "SAAAAAAAAAAAAAATGGGGGGGGGGGGGGG"},
        motif_frequency_matrices=motif_frequency_matrices,
        motif_sizes=motif_sizes,
        site_index=pd.Index(["MAPK1;S202;"], name="site_id"),
        min_motif_size=1,
        flank_size=DEFAULT_MOTIF_FLANK_SIZE,
        sequence_semantics=SEQUENCE_SEMANTICS_CENTRED_SEQUENCE,
    )

    row = result.sequence_validation.rows[0]
    assert row.status == SEQUENCE_VALIDATION_STATUS_OFF_CENTRE_SEQUENCE
    assert pd.isna(result.motif_scores.loc["MAPK1;S202;", "K1"])


def test_motif_scoring_rejects_non_site_identifier_without_explicit_position_support() -> (
    None
):
    motif_frequency_matrices, motif_sizes = build_motif_library_from_sequences(
        motif_sequences={"K1": ["AAAAAAASAAAAAAA"]},
        flank_size=DEFAULT_MOTIF_FLANK_SIZE,
    )
    result = score_phosphosite_motifs(
        site_sequences={"SITE_1": "GGGGGGGAAAAAAASAAAAAAAGGGGGGG"},
        motif_frequency_matrices=motif_frequency_matrices,
        motif_sizes=motif_sizes,
        site_index=pd.Index(["SITE_1"], name="site_id"),
        min_motif_size=1,
        flank_size=DEFAULT_MOTIF_FLANK_SIZE,
        sequence_semantics=SEQUENCE_SEMANTICS_CENTRED_WINDOW,
    )

    row = result.sequence_validation.rows[0]
    assert row.status == SEQUENCE_VALIDATION_STATUS_INVALID_SITE_ID
    assert pd.isna(result.motif_scores.loc["SITE_1", "K1"])


def test_motif_scoring_reports_missing_sequence_when_required_for_scoring() -> None:
    motif_frequency_matrices, motif_sizes = build_motif_library_from_sequences(
        motif_sequences={"K1": ["AAAAAAASAAAAAAA"]},
        flank_size=DEFAULT_MOTIF_FLANK_SIZE,
    )
    result = score_phosphosite_motifs(
        site_sequences={"MAPK1;S202;": None},
        motif_frequency_matrices=motif_frequency_matrices,
        motif_sizes=motif_sizes,
        site_index=pd.Index(["MAPK1;S202;"], name="site_id"),
        min_motif_size=1,
        flank_size=DEFAULT_MOTIF_FLANK_SIZE,
        sequence_semantics=SEQUENCE_SEMANTICS_CENTRED_WINDOW,
    )

    row = result.sequence_validation.rows[0]
    assert row.status == SEQUENCE_VALIDATION_STATUS_MISSING_SEQUENCE
    assert pd.isna(result.motif_scores.loc["MAPK1;S202;", "K1"])


def test_motif_scoring_rejects_non_phospho_centre_residue() -> None:
    motif_frequency_matrices, motif_sizes = build_motif_library_from_sequences(
        motif_sequences={"K1": ["AAAAAAASAAAAAAA"]},
        flank_size=DEFAULT_MOTIF_FLANK_SIZE,
    )
    result = score_phosphosite_motifs(
        site_sequences={"MAPK1;S202;": "AAAAAAAQAAAAAAA"},
        motif_frequency_matrices=motif_frequency_matrices,
        motif_sizes=motif_sizes,
        site_index=pd.Index(["MAPK1;S202;"], name="site_id"),
        min_motif_size=1,
        flank_size=DEFAULT_MOTIF_FLANK_SIZE,
        sequence_semantics=SEQUENCE_SEMANTICS_CENTRED_WINDOW,
    )

    row = result.sequence_validation.rows[0]
    assert row.status == SEQUENCE_VALIDATION_STATUS_NON_PHOSPHO_CENTRE_RESIDUE
    assert pd.isna(result.motif_scores.loc["MAPK1;S202;", "K1"])


def test_kinase_workflow_exposes_sequence_validation_diagnostics() -> None:
    display_ids = ["MAPK1;S202;", "MAPK1;T205;", "MAPK1;S210;"]
    site_ids = site_key_index_from_display_ids(display_ids)
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 2.0, 1.5],
            "sample_b": [2.0, 1.0, 1.2],
        },
        index=site_ids,
    )
    site_metadata = pd.DataFrame(
        {
            "site_key": site_ids.astype(str).tolist(),
            "display_id": display_ids,
            "gene_symbol": ["MAPK1", "MAPK1", "MAPK1"],
            "site": ["S202", "T205", "S210"],
            "site_sequence": [
                ("A" * 15) + str(site).strip().upper()[0] + ("A" * 15)
                for site in ["S202", "T205", "S210"]
            ],
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
                "substrate_site": ["MAPK1;S202;", "MAPK1;T205;"],
            }
        ),
        site_sequences=pd.DataFrame(
            {
                "site_sequence": [
                    "AAAAAAASAAAAAAA",
                    "VVVVVVVTVVVVVVV",
                    "AAAAAAAXAAAAAAA",
                ],
            },
            index=pd.Index(display_ids, name="site_id"),
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
    scoring_result = (
        KinaseScoringRunner()
        .run(
            request=interpreted,
            config=interpreted.execution_config,
        )
        .scoring_result
    )

    assert scoring_result.motif_sequence_validation is not None
    summary = scoring_result.motif_sequence_validation.summary()
    assert summary["total_sequences"] == 3
    assert summary["valid_sequences"] == 2
    assert summary["invalid_sequences"] == 1
    assert summary["short_sequences"] == 0
    assert summary["off_centre_sequences"] == 0
    assert summary["site_residue_mismatches"] == 0
    assert summary["unsupported_residue_characters"] == 1
    assert summary["sequences_excluded_from_motif_scoring"] == 1
    assert summary["total_sites_considered"] == 3
    assert summary["sites_with_valid_site_sequence"] == 2
    assert summary["sites_without_valid_site_sequence"] == 1
    assert summary["site_sequence_coverage_fraction"] == 2 / 3
    assert summary["sites_used_for_motif_scoring"] == 2
    assert summary["sites_excluded_from_motif_scoring_due_to_sequence"] == 1

    assert scoring_result.motif_library_validation is not None
    library_summary = scoring_result.motif_library_validation.summary()
    assert library_summary["reference_sequences_provided"] == 2
    assert library_summary["reference_sequences_accepted"] == 2
    assert library_summary["reference_sequences_excluded"] == 0
    assert library_summary["excluded_missing_sequence"] == 0
    assert library_summary["excluded_short_window"] == 0
    assert library_summary["excluded_unsupported_residue"] == 0
    assert library_summary["excluded_off_centre_residue"] == 0
    assert library_summary["excluded_site_residue_mismatch"] == 0
    assert library_summary["expected_window_size"] == 15

    assert scoring_result.motif_scores is not None
    s210_site_key = site_key_from_display_id("MAPK1;S210;")
    assert pd.isna(scoring_result.motif_scores.loc[s210_site_key, "KINASE_A"])


def test_sequence_coverage_summary_reports_full_coverage_for_all_valid_sites() -> None:
    result = _validator().run(
        rows=[
            SequenceValidationInput(
                site_id="MAPK1;S202;",
                site_sequence="AAAAAAASAAAAAAA",
            ),
            SequenceValidationInput(
                site_id="MAPK1;T205;",
                site_sequence="VVVVVVVTVVVVVVV",
            ),
        ]
    )

    coverage = result.site_sequence_coverage_summary()
    assert coverage["total_sites_considered"] == 2
    assert coverage["sites_with_valid_site_sequence"] == 2
    assert coverage["sites_without_valid_site_sequence"] == 0
    assert coverage["site_sequence_coverage_fraction"] == 1.0
    assert coverage["sites_used_for_motif_scoring"] == 2
    assert coverage["sites_excluded_from_motif_scoring_due_to_sequence"] == 0


def test_sequence_coverage_summary_reports_zero_coverage_for_all_invalid_sites() -> (
    None
):
    result = _validator().run(
        rows=[
            SequenceValidationInput(
                site_id="MAPK1;S202;",
                site_sequence=None,
            ),
            SequenceValidationInput(
                site_id="MAPK1;S203;",
                site_sequence="AAAAAAAXAAAAAAA",
            ),
        ]
    )

    coverage = result.site_sequence_coverage_summary()
    assert coverage["total_sites_considered"] == 2
    assert coverage["sites_with_valid_site_sequence"] == 0
    assert coverage["sites_without_valid_site_sequence"] == 2
    assert coverage["site_sequence_coverage_fraction"] == 0.0
    assert coverage["sites_used_for_motif_scoring"] == 0
    assert coverage["sites_excluded_from_motif_scoring_due_to_sequence"] == 2


def test_kinase_workflow_reports_partial_sequence_coverage_in_provenance() -> None:
    display_ids = ["MAPK1;S202;", "MAPK1;T205;", "MAPK1;S210;"]
    site_ids = site_key_index_from_display_ids(display_ids)
    dataset = AnalysisReadyPhosphoDataset(
        phospho=pd.DataFrame(
            {
                "sample_a": [1.0, 2.0, 3.0],
                "sample_b": [2.0, 3.0, 4.0],
            },
            index=site_ids,
        ),
        site_metadata=pd.DataFrame(
            {
                "site_key": site_ids.astype(str).tolist(),
                "display_id": display_ids,
                "gene_symbol": ["MAPK1", "MAPK1", "MAPK1"],
                "site": ["S202", "T205", "S210"],
                "site_sequence": [
                    ("A" * 15) + str(site).strip().upper()[0] + ("A" * 15)
                    for site in ["S202", "T205", "S210"]
                ],
            },
            index=site_ids,
        ),
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
                "substrate_site": ["MAPK1;S202;", "MAPK1;T205;"],
            }
        ),
        site_sequences=pd.DataFrame(
            {
                "site_sequence": [
                    "AAAAAAASAAAAAAA",
                    "VVVVVVVTVVVVVVV",
                    "AAAAAAAXAAAAAAA",
                ],
            },
            index=pd.Index(display_ids, name="site_id"),
        ),
    )

    result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=references,
            scoring_config=KinaseScoringConfig(min_substrates=2),
            prediction_config=KinasePredictionConfig(
                top_k=1,
                deterministic_max_selected_kinases=1,
                adaptive_ensemble_runs=1,
            ),
            activity_config=None,
        )
    )

    assert not result.prediction_result.pred_mat.empty
    assert result.provenance is not None
    diagnostics = result.provenance.workflow_parameters["scoring_diagnostics"]
    assert diagnostics["total_sites_considered"] == 3
    assert diagnostics["sites_with_valid_site_sequence"] == 2
    assert diagnostics["sites_without_valid_site_sequence"] == 1
    assert diagnostics["site_sequence_coverage_fraction"] == 2 / 3
    assert diagnostics["sites_used_for_motif_scoring"] == 2
    assert diagnostics["sites_excluded_from_motif_scoring_due_to_sequence"] == 1
    nested_coverage = diagnostics["motif_site_sequence_coverage"]
    assert nested_coverage["total_sites_considered"] == 3
    assert nested_coverage["sites_with_valid_site_sequence"] == 2
    assert nested_coverage["sites_without_valid_site_sequence"] == 1


def test_kinase_workflow_continues_when_no_sites_have_valid_sequence() -> None:
    display_ids = ["MAPK1;S202;", "MAPK1;T205;"]
    site_ids = site_key_index_from_display_ids(display_ids)
    dataset = AnalysisReadyPhosphoDataset(
        phospho=pd.DataFrame(
            {"sample_a": [1.0, 2.0], "sample_b": [2.0, 3.0]},
            index=site_ids,
        ),
        site_metadata=pd.DataFrame(
            {
                "site_key": site_ids.astype(str).tolist(),
                "display_id": display_ids,
                "gene_symbol": ["MAPK1", "MAPK1"],
                "site": ["S202", "T205"],
                "site_sequence": [
                    ("A" * 15) + str(site).strip().upper()[0] + ("A" * 15)
                    for site in ["S202", "T205"]
                ],
            },
            index=site_ids,
        ),
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
                "substrate_site": ["MAPK1;S202;", "MAPK1;T205;"],
            }
        ),
        site_sequences=pd.DataFrame(
            {"site_sequence": ["AAAAAAAXAAAAAAA", "ASAA"]},
            index=pd.Index(display_ids, name="site_id"),
        ),
    )

    result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
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
    )

    assert not result.prediction_result.pred_mat.empty
    assert result.scoring_result.motif_scores is not None
    assert result.scoring_result.motif_scores.isna().all().all()
    assert result.scoring_result.rank_weighted_fusion_scores is not None
    assert result.provenance is not None
    diagnostics = result.provenance.workflow_parameters["scoring_diagnostics"]
    assert diagnostics["total_sites_considered"] == 2
    assert diagnostics["sites_with_valid_site_sequence"] == 0
    assert diagnostics["sites_without_valid_site_sequence"] == 2
    assert diagnostics["site_sequence_coverage_fraction"] == 0.0


def test_motif_library_accepts_valid_reference_sequence_and_builds_profile() -> None:
    motif_frequency_matrices, motif_sizes = build_motif_library(
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["K1"],
                "substrate_site": ["MAPK1;S202;"],
            }
        ),
        site_sequences=pd.Series(
            {"MAPK1;S202;": "AAAAAAASAAAAAAA"},
            dtype=object,
        ),
        flank_size=DEFAULT_MOTIF_FLANK_SIZE,
    )

    assert motif_sizes.loc["K1"] == 1.0
    matrix = motif_frequency_matrices["K1"]
    assert matrix.loc["S", "p8"] == 1.0

    validation = get_motif_library_validation(motif_sizes)
    assert validation is not None
    assert validation.accepted_reference_sequences == 1
    assert validation.excluded_reference_sequences == 0


def test_motif_library_excludes_invalid_reference_windows_and_reports_diagnostics() -> (
    None
):
    motif_frequency_matrices, motif_sizes = build_motif_library(
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["K1"] * 6,
                "substrate_site": [
                    "MAPK1;S202;",
                    "MAPK1;S203;",
                    "MAPK1;S204;",
                    "MAPK1;S205;",
                    "MAPK1;S206;",
                    "MAPK1;S207;",
                ],
            }
        ),
        site_sequences=pd.Series(
            {
                "MAPK1;S202;": "AAAAAAASAAAAAAA",
                "MAPK1;S203;": "ASAA",
                "MAPK1;S204;": "AAAAAAAXAAAAAAA",
                "MAPK1;S205;": "SAAAAAAAAAAAAAA",
                "MAPK1;S206;": "AAAAAAATAAAAAAA",
            },
            dtype=object,
        ),
        flank_size=DEFAULT_MOTIF_FLANK_SIZE,
    )

    validation = get_motif_library_validation(motif_sizes)
    assert validation is not None
    summary = validation.summary()
    assert summary["reference_sequences_provided"] == 6
    assert summary["reference_sequences_accepted"] == 1
    assert summary["reference_sequences_excluded"] == 5
    assert summary["excluded_missing_sequence"] == 1
    assert summary["excluded_short_window"] == 1
    assert summary["excluded_unsupported_residue"] == 1
    assert summary["excluded_off_centre_residue"] == 0
    assert summary["excluded_non_phospho_centre_residue"] == 1
    assert summary["excluded_site_residue_mismatch"] == 1
    assert summary["expected_window_size"] == 15
    assert "minimum length 15" in summary["accepted_window_length_policy"]
    assert "supported residues" in summary["unsupported_residue_policy"]
    assert motif_sizes.loc["K1"] == 1.0
    assert motif_frequency_matrices["K1"].loc["S", "p8"] == 1.0


def test_invalid_reference_sequences_do_not_change_motif_profile() -> None:
    valid_map = pd.DataFrame(
        {
            "kinase": ["K1"],
            "substrate_site": ["MAPK1;S202;"],
        }
    )
    valid_sequences = pd.Series(
        {"MAPK1;S202;": "AAAAAAASAAAAAAA"},
        dtype=object,
    )
    baseline_matrices, baseline_sizes = build_motif_library(
        kinase_substrate_map=valid_map,
        site_sequences=valid_sequences,
        flank_size=DEFAULT_MOTIF_FLANK_SIZE,
    )

    mixed_map = pd.DataFrame(
        {
            "kinase": ["K1", "K1", "K1"],
            "substrate_site": ["MAPK1;S202;", "MAPK1;S203;", "MAPK1;S204;"],
        }
    )
    mixed_sequences = pd.Series(
        {
            "MAPK1;S202;": "AAAAAAASAAAAAAA",
            "MAPK1;S203;": "AAAAAAAXAAAAAAA",
            "MAPK1;S204;": "ASAA",
        },
        dtype=object,
    )
    mixed_matrices, mixed_sizes = build_motif_library(
        kinase_substrate_map=mixed_map,
        site_sequences=mixed_sequences,
        flank_size=DEFAULT_MOTIF_FLANK_SIZE,
    )

    pd.testing.assert_series_equal(baseline_sizes, mixed_sizes)
    pd.testing.assert_frame_equal(baseline_matrices["K1"], mixed_matrices["K1"])


def test_structured_explicit_sequences_support_site_aware_validation_and_exclusion() -> (
    None
):
    motif_frequency_matrices, motif_sizes = build_motif_library_from_sequences(
        motif_sequences={
            "K1": [
                ExplicitMotifSequence(
                    reference_id="REF_VALID",
                    site_id="MAPK1;S202;",
                    kinase="K1",
                    sequence="AAAAAAASAAAAAAA",
                ),
                {
                    "reference_id": "REF_SHORT",
                    "site_id": "MAPK1;S203;",
                    "kinase": "K1",
                    "sequence": "ASAA",
                },
                {
                    "reference_id": "REF_UNSUPPORTED",
                    "site_id": "MAPK1;S204;",
                    "kinase": "K1",
                    "sequence": "AAAAAAAXAAAAAAA",
                },
                {
                    "reference_id": "REF_NON_PHOSPHO",
                    "site_id": None,
                    "kinase": "K1",
                    "sequence": "AAAAAAAAAAAAAAA",
                },
                {
                    "reference_id": "REF_MISMATCH",
                    "site_id": "MAPK1;S205;",
                    "kinase": "K1",
                    "sequence": "AAAAAAATAAAAAAA",
                },
            ]
        },
        flank_size=DEFAULT_MOTIF_FLANK_SIZE,
    )

    validation = get_motif_library_validation(motif_sizes)
    assert validation is not None
    summary = validation.summary()
    assert summary["reference_sequences_provided"] == 5
    assert summary["reference_sequences_accepted"] == 1
    assert summary["reference_sequences_excluded"] == 4
    assert summary["excluded_short_window"] == 1
    assert summary["excluded_unsupported_residue"] == 1
    assert summary["excluded_non_phospho_centre_residue"] == 1
    assert summary["excluded_site_residue_mismatch"] == 1
    assert summary["sequences_excluded_from_motif_profile_construction"] == 4

    rows_by_reference = {row.reference_id: row for row in validation.rows}
    assert rows_by_reference["REF_VALID"].status == "valid"
    assert rows_by_reference["REF_VALID"].site_id == "MAPK1;S202;"
    assert (
        rows_by_reference["REF_NON_PHOSPHO"].status
        == SEQUENCE_VALIDATION_STATUS_NON_PHOSPHO_CENTRE_RESIDUE
    )
    assert (
        rows_by_reference["REF_MISMATCH"].status
        == SEQUENCE_VALIDATION_STATUS_SITE_RESIDUE_MISMATCH
    )
    assert motif_sizes.loc["K1"] == 1.0
    assert motif_frequency_matrices["K1"].loc["S", "p8"] == 1.0


def test_structured_explicit_sequence_rejects_invalid_site_id_format() -> None:
    _, motif_sizes = build_motif_library_from_sequences(
        motif_sequences={
            "K1": [
                {
                    "reference_id": "REF_BAD_SITE",
                    "site_id": "S202",
                    "kinase": "K1",
                    "sequence": "AAAAAAASAAAAAAA",
                }
            ]
        },
        flank_size=DEFAULT_MOTIF_FLANK_SIZE,
    )

    validation = get_motif_library_validation(motif_sizes)
    assert validation is not None
    summary = validation.summary()
    assert summary["reference_sequences_accepted"] == 0
    assert summary["excluded_invalid_site_id"] == 1
    assert "K1" not in motif_sizes.index
    assert validation.rows[0].status == SEQUENCE_VALIDATION_STATUS_INVALID_SITE_ID


def test_bare_explicit_sequences_remain_supported_without_site_mismatch_claims() -> (
    None
):
    _, motif_sizes = build_motif_library_from_sequences(
        motif_sequences={"K1": ["AAAAAAASAAAAAAA", "AAAAAAATAAAAAAA"]},
        flank_size=DEFAULT_MOTIF_FLANK_SIZE,
    )

    validation = get_motif_library_validation(motif_sizes)
    assert validation is not None
    summary = validation.summary()
    assert summary["reference_sequences_provided"] == 2
    assert summary["reference_sequences_accepted"] == 2
    assert summary["excluded_site_residue_mismatch"] == 0
    assert summary["excluded_invalid_site_id"] == 0
    for row in validation.rows:
        assert row.status == "valid"
        assert row.site_id is None
        assert row.expected_centre_residue is None


def test_bare_explicit_invalid_sequence_is_rejected_and_reported() -> None:
    _, motif_sizes = build_motif_library_from_sequences(
        motif_sequences={"K1": ["AAAAAAAXAAAAAAA"]},
        flank_size=DEFAULT_MOTIF_FLANK_SIZE,
    )

    validation = get_motif_library_validation(motif_sizes)
    assert validation is not None
    summary = validation.summary()
    assert summary["reference_sequences_accepted"] == 0
    assert summary["reference_sequences_excluded"] == 1
    assert summary["excluded_unsupported_residue"] == 1
    assert (
        validation.rows[0].status
        == SEQUENCE_VALIDATION_STATUS_UNSUPPORTED_RESIDUE_CHARACTER
    )


def test_motif_library_validation_rows_preserve_reference_provenance() -> None:
    _, motif_sizes = build_motif_library(
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["K1", "K1"],
                "substrate_site": ["MAPK1;S202;", "MAPK1;S203;"],
            }
        ),
        site_sequences=pd.Series(
            {
                "MAPK1;S202;": "AAAAAAASAAAAAAA",
                "MAPK1;S203;": "AAAAAAATAAAAAAA",
            },
            dtype=object,
        ),
        flank_size=DEFAULT_MOTIF_FLANK_SIZE,
    )

    validation = get_motif_library_validation(motif_sizes)
    assert validation is not None
    assert len(validation.rows) == 2

    rows = {(row.kinase, row.reference_id): row for row in validation.rows}
    assert rows[("K1", "MAPK1;S202;")].status == "valid"
    mismatch = rows[("K1", "MAPK1;S203;")]
    assert mismatch.status == SEQUENCE_VALIDATION_STATUS_SITE_RESIDUE_MISMATCH
    assert mismatch.expected_centre_residue == "S"
    assert mismatch.observed_centre_residue == "T"
    assert mismatch.reason is not None


def test_query_sequence_validation_behavior_remains_unchanged_with_library_validation() -> (
    None
):
    motif_frequency_matrices, motif_sizes = build_motif_library(
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["K1", "K1"],
                "substrate_site": ["MAPK1;S202;", "MAPK1;S203;"],
            }
        ),
        site_sequences=pd.Series(
            {
                "MAPK1;S202;": "AAAAAAASAAAAAAA",
                "MAPK1;S203;": "AAAAAAAXAAAAAAA",
            },
            dtype=object,
        ),
        flank_size=DEFAULT_MOTIF_FLANK_SIZE,
    )
    motif_result = score_phosphosite_motifs(
        site_sequences={
            "MAPK1;S202;": "AAAAAAASAAAAAAA",
            "MAPK1;S210;": "AAAAAAAXAAAAAAA",
        },
        motif_frequency_matrices=motif_frequency_matrices,
        motif_sizes=motif_sizes,
        site_index=pd.Index(["MAPK1;S202;", "MAPK1;S210;"], name="site_id"),
        min_motif_size=1,
        flank_size=DEFAULT_MOTIF_FLANK_SIZE,
        library_validation=get_motif_library_validation(motif_sizes),
    )

    summary = motif_result.sequence_validation.summary()
    assert summary["total_sequences"] == 2
    assert summary["valid_sequences"] == 1
    assert summary["invalid_sequences"] == 1
    assert summary["short_sequences"] == 0
    assert summary["off_centre_sequences"] == 0
    assert summary["site_residue_mismatches"] == 0
    assert summary["unsupported_residue_characters"] == 1
    assert summary["sequences_excluded_from_motif_scoring"] == 1
