from __future__ import annotations

import pandas as pd
import pytest

from phospy.provenance.scientific_policy_models import ScientificPolicyId
from phospy.science.prediction.motif_scoring import (
    KINASE_LIBRARY_MATRIX_STATUS_FILTERED_RESIDUE_CLASS,
    KINASE_LIBRARY_SITE_STATUS_MISSING_SEQUENCE,
    KINASE_LIBRARY_SITE_STATUS_UNSUPPORTED_SEQUENCE_LENGTH,
    KINASE_LIBRARY_SITE_STATUS_VALID_SCORED_SITE,
    KINASE_LIBRARY_SITE_STATUS_WRONG_CENTRAL_RESIDUE,
    KINASE_LIBRARY_SITE_STATUS_WRONG_RESIDUE_CLASS,
    KinaseLibraryMotifMatrix,
    score_kinase_library_motifs,
)
from phospy.science.prediction.policies import KinaseLibraryMotifScoringPolicy


def test_kinase_library_motif_scoring_calculates_exact_raw_scores() -> None:
    result = score_kinase_library_motifs(
        site_sequences={
            "P1;S10;": "AST",
            "P2;Y20;": "AYP",
        },
        matrices=[
            _matrix(
                "K_ST",
                "ser_thr",
                {
                    ("A", -1): 1.0,
                    ("S", 0): 2.0,
                    ("T", 1): 3.0,
                },
            ),
            _matrix(
                "K_TYR",
                "tyr",
                {
                    ("A", -1): 4.0,
                    ("Y", 0): 5.0,
                    ("P", 1): 6.0,
                },
            ),
        ],
        flank_size=1,
        score_scale="raw_log2_enrichment",
    )

    assert list(result.raw_scores.index) == ["P1;S10;", "P2;Y20;"]
    assert list(result.raw_scores.columns) == ["K_ST", "K_TYR"]
    assert result.raw_scores.at["P1;S10;", "K_ST"] == pytest.approx(6.0)
    assert pd.isna(result.raw_scores.at["P1;S10;", "K_TYR"])
    assert pd.isna(result.raw_scores.at["P2;Y20;", "K_ST"])
    assert result.raw_scores.at["P2;Y20;", "K_TYR"] == pytest.approx(15.0)
    assert (
        result.site_diagnostics.at["P1;S10;", "status"]
        == KINASE_LIBRARY_SITE_STATUS_VALID_SCORED_SITE
    )
    assert result.score_scale_metadata.score_scale == "raw_log2_enrichment"
    assert result.score_scale_metadata.sequence_window["window_size"] == 3


def test_kinase_library_motif_scoring_builds_percentiles_and_reference_ranks() -> None:
    result = score_kinase_library_motifs(
        site_sequences={"P1;S10;": "AST", "P2;Y20;": "AYP"},
        matrices=[
            _matrix(
                "K_ST",
                "ser_thr",
                {("A", -1): 1.0, ("S", 0): 2.0, ("T", 1): 3.0},
            ),
            _matrix(
                "K_TYR",
                "tyr",
                {("A", -1): 4.0, ("Y", 0): 5.0, ("P", 1): 6.0},
            ),
        ],
        flank_size=1,
        reference_distributions={
            ("K_ST", "ser_thr"): [0.0, 4.0, 6.0, 8.0],
            ("K_TYR", "tyr"): [10.0, 15.0, 20.0],
        },
    )

    assert result.percentile_ranks is not None
    assert result.reference_ranks is not None
    assert result.percentile_ranks.at["P1;S10;", "K_ST"] == pytest.approx(75.0)
    assert result.reference_ranks.at["P1;S10;", "K_ST"] == pytest.approx(2.0)
    assert result.percentile_ranks.at["P2;Y20;", "K_TYR"] == pytest.approx(
        100.0 * 2.0 / 3.0
    )
    assert result.reference_ranks.at["P2;Y20;", "K_TYR"] == pytest.approx(2.0)
    assert result.score_scale_metadata.percentile_method is not None
    assert result.score_scale_metadata.rank_method is not None


def test_kinase_library_motif_scoring_filters_wrong_residue_class() -> None:
    result = score_kinase_library_motifs(
        site_sequences={"P1;S10;": "AST", "P2;Y20;": "AYP"},
        matrices=[
            _matrix("K_ST", "ser_thr", {("A", -1): 1.0, ("S", 0): 2.0}),
            _matrix("K_TYR", "tyr", {("A", -1): 4.0, ("Y", 0): 5.0}),
        ],
        residue_classes=("ser_thr",),
        flank_size=1,
    )

    assert list(result.raw_scores.columns) == ["K_ST"]
    assert (
        result.site_diagnostics.at["P2;Y20;", "status"]
        == KINASE_LIBRARY_SITE_STATUS_WRONG_RESIDUE_CLASS
    )
    assert pd.isna(result.raw_scores.at["P2;Y20;", "K_ST"])
    assert (
        result.kinase_diagnostics.loc[("K_TYR", "tyr"), "status"]
        == KINASE_LIBRARY_MATRIX_STATUS_FILTERED_RESIDUE_CLASS
    )


def test_kinase_library_motif_scoring_reports_required_site_attrition_statuses() -> (
    None
):
    result = score_kinase_library_motifs(
        site_sequences={
            "P0;S1;": None,
            "P1;S10;": "AQT",
            "P2;S20;": "ASTA",
            "P3;S30;": "ATT",
            "P4;S40;": "AST",
        },
        matrices=[_matrix("K_ST", "ser_thr", {("A", -1): 1.0, ("S", 0): 2.0})],
        residue_classes=("ser_thr",),
        flank_size=1,
    )

    statuses = result.site_diagnostics.loc[:, "status"].to_dict()
    assert statuses["P0;S1;"] == KINASE_LIBRARY_SITE_STATUS_MISSING_SEQUENCE
    assert statuses["P1;S10;"] == KINASE_LIBRARY_SITE_STATUS_WRONG_CENTRAL_RESIDUE
    assert statuses["P2;S20;"] == KINASE_LIBRARY_SITE_STATUS_UNSUPPORTED_SEQUENCE_LENGTH
    assert statuses["P3;S30;"] == KINASE_LIBRARY_SITE_STATUS_WRONG_CENTRAL_RESIDUE
    assert statuses["P4;S40;"] == KINASE_LIBRARY_SITE_STATUS_VALID_SCORED_SITE
    assert bool(result.site_diagnostics.at["P0;S1;", "excluded_from_scoring"]) is True
    assert result.site_diagnostics.at["P4;S40;", "scored_kinase_count"] == 1


def test_kinase_library_motif_scoring_output_order_is_deterministic() -> None:
    site_sequences = pd.Series(
        ["ATT", "AST"],
        index=pd.Index(["P2;T20;", "P1;S10;"], name="site_id"),
    )
    matrices = [
        _matrix("K2", "ser_thr", {("A", -1): 10.0, ("T", 0): 1.0}),
        _matrix("K1", "ser_thr", {("A", -1): 3.0, ("S", 0): 2.0}),
    ]

    first = score_kinase_library_motifs(
        site_sequences=site_sequences,
        matrices=matrices,
        residue_classes=("ser_thr",),
        flank_size=1,
    )
    second = score_kinase_library_motifs(
        site_sequences=site_sequences,
        matrices=matrices,
        residue_classes=("ser_thr",),
        flank_size=1,
    )

    assert list(first.raw_scores.index) == ["P2;T20;", "P1;S10;"]
    assert list(first.raw_scores.columns) == ["K2", "K1"]
    pd.testing.assert_frame_equal(first.raw_scores, second.raw_scores)
    pd.testing.assert_frame_equal(first.site_diagnostics, second.site_diagnostics)


def test_kinase_library_motif_scoring_accepts_reference_matrix_model() -> None:
    from phospy.science.references import KinaseLibraryMatrix

    reference_matrix = KinaseLibraryMatrix(
        kinase="K_REF",
        residue_class="ser_thr",
        score_table=_matrix(
            "K_REF",
            "ser_thr",
            {("A", -1): 2.0, ("S", 0): 3.0, ("T", 1): 4.0},
        ).score_table,
    )

    result = score_kinase_library_motifs(
        site_sequences={"P1;S10;": "AST"},
        matrices=[reference_matrix],
        flank_size=1,
    )

    assert result.raw_scores.at["P1;S10;", "K_REF"] == pytest.approx(9.0)


def test_kinase_library_motif_scoring_policy_exposes_score_scale_metadata() -> None:
    policy = KinaseLibraryMotifScoringPolicy(
        score_scale="raw_log2_enrichment",
        residue_classes=("ser_thr", "tyr"),
        upstream_residues=1,
        downstream_residues=1,
        sequence_semantics="centred_window",
        reference_distributions_supplied=True,
    )

    record = policy.record

    assert record.id == ScientificPolicyId.KINASE_LIBRARY_MOTIF_SCORING
    assert record.parameters["score_scale"] == "raw_log2_enrichment"
    assert record.parameters["reference_distributions_supplied"] is True
    assert record.parameters["percentile_method"] is not None


def _matrix(
    kinase: str,
    residue_class: str,
    scores: dict[tuple[str, int], float],
) -> KinaseLibraryMotifMatrix:
    score_table = pd.DataFrame(
        0.0,
        index=pd.Index(["A", "S", "T", "Y", "P"], name="amino_acid"),
        columns=pd.Index([-1, 0, 1], name="position"),
    )
    for (amino_acid, position), score in scores.items():
        score_table.at[amino_acid, position] = score
    return KinaseLibraryMotifMatrix(
        kinase=kinase,
        residue_class=residue_class,
        score_table=score_table,
    )
