from __future__ import annotations

from collections.abc import Mapping, Sequence

import pandas as pd
import pandas.testing as pdt
import pytest

from phospy.api.configs import (
    SIGNALOME_ASSIGNMENT_POLICY_CUTOFF_BINARY,
    SIGNALOME_ASSIGNMENT_POLICY_WEIGHTED_TOP,
)
from phospy.errors import WorkflowStageError
from phospy.science.signalomes.constants import (
    EXPANDED_SIGNALOME_ROW_KIND_COLUMN,
    EXPANDED_SIGNALOME_ROW_KIND_SITE,
    EXPANDED_SIGNALOME_ROW_KIND_SUMMARY,
)
from phospy.science.signalomes.science import (
    build_expanded_signalome_table,
    build_kinase_network,
    build_kinase_network_with_diagnostics,
    build_signalome_module_table,
)


def _with_site_identity(module_assignments: pd.DataFrame) -> pd.DataFrame:
    resolved = module_assignments.copy(deep=True)
    site_keys = resolved.index.astype(str).tolist()
    resolved.index = pd.Index(site_keys, name="site_key")
    resolved.loc[:, "site_key"] = site_keys
    resolved.loc[:, "display_id"] = site_keys
    resolved.loc[:, "gene_symbol"] = (
        resolved.loc[:, "protein_group_id"].astype(str).tolist()
    )
    resolved.loc[:, "site"] = site_keys
    resolved.loc[:, "protein_accession"] = ""
    resolved.loc[:, "isoform_id"] = ""
    return resolved


def _historical_baseline_build_signalome_module_table(
    *,
    module_assignments: pd.DataFrame,
    kinase_substrates: Mapping[str, Sequence[str]],
    kinase_order: Sequence[str],
) -> pd.DataFrame:
    module_index = pd.Index(
        sorted(
            {
                int(value)
                for value in module_assignments.loc[:, "module_id"]
                if int(value) > 0
            }
        ),
        name="module_id",
    )
    kinase_index = pd.Index([str(kinase) for kinase in kinase_order], name="kinase")
    module_table = pd.DataFrame(
        0.0,
        index=module_index.copy(),
        columns=kinase_index.copy(),
    )

    protein_to_module = (
        module_assignments.loc[:, ["protein_group_id", "module_id"]]
        .drop_duplicates(subset=["protein_group_id"])
        .set_index("protein_group_id")
        .loc[:, "module_id"]
        .astype("int64")
    )
    protein_to_module = protein_to_module.loc[protein_to_module > 0]
    site_to_protein = module_assignments.loc[:, "protein_group_id"].astype(str)
    site_to_protein.index = pd.Index(
        site_to_protein.index.astype(str),
        name="site_id",
    )

    for kinase in kinase_index:
        substrate_sites = pd.Index(
            [str(site_id) for site_id in kinase_substrates.get(str(kinase), ())],
            name="site_id",
        )
        if substrate_sites.empty:
            continue
        substrate_proteins = (
            site_to_protein.reindex(substrate_sites).dropna().astype(str)
        )
        if substrate_proteins.empty:
            continue
        unique_proteins = pd.Index(sorted(set(substrate_proteins.tolist())))
        module_hits = (
            protein_to_module.reindex(unique_proteins).dropna().astype("int64")
        )
        if module_hits.empty:
            continue
        counts = module_hits.value_counts().astype(float)
        module_table.loc[counts.index.astype(int), kinase] = counts.to_numpy(
            dtype=float, copy=False
        )

    row_totals = module_table.sum(axis=1)
    non_zero_rows = row_totals > 0.0
    if non_zero_rows.any():
        module_table.loc[non_zero_rows] = (
            module_table.loc[non_zero_rows].div(row_totals.loc[non_zero_rows], axis=0)
            * 100.0
        )
    return module_table.astype(float).round(3)


def test_build_signalome_module_table_matches_historical_baseline_semantics() -> None:
    module_assignments = pd.DataFrame(
        {
            "protein_group_id": ["P1", "P1", "P2", "P3", "P4", "P5", "P6", "P3"],
            "module_id": [1, 1, 2, 2, 0, 3, 3, 3],
        },
        index=pd.Index(
            ["S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8"],
            name="site_id",
        ),
    )
    kinase_substrates = {
        "K1": ("S1", "S2", "S4", "S8"),
        "K2": ("S3", "S5", "S6", "MISSING"),
        "K3": ("S7",),
    }
    kinase_order = ["K2", "K1", "K3", "K1"]

    expected = _historical_baseline_build_signalome_module_table(
        module_assignments=module_assignments,
        kinase_substrates=kinase_substrates,
        kinase_order=kinase_order,
    )
    observed = build_signalome_module_table(
        module_assignments=module_assignments,
        kinase_substrates=kinase_substrates,
        kinase_order=kinase_order,
    )

    pdt.assert_frame_equal(observed, expected, check_dtype=False)


def test_build_signalome_module_table_weighted_top_propagates_fractional_support() -> (
    None
):
    module_assignments = pd.DataFrame(
        {
            "protein_group_id": ["P1", "P1", "P2", "P2"],
            "module_id": [1, 1, 2, 2],
            "top_kinase_weights": [
                (("K1", 0.5), ("K2", 0.5)),
                (("K1", 1.0),),
                (("K2", 1.0),),
                (("K2", 1.0),),
            ],
        },
        index=pd.Index(["S1", "S2", "S3", "S4"], name="site_id"),
    )

    weighted = build_signalome_module_table(
        module_assignments=module_assignments,
        kinase_substrates={"K1": (), "K2": ()},
        kinase_order=["K1", "K2"],
        assignment_policy=SIGNALOME_ASSIGNMENT_POLICY_WEIGHTED_TOP,
    )
    binary = build_signalome_module_table(
        module_assignments=module_assignments,
        kinase_substrates={"K1": (), "K2": ()},
        kinase_order=["K1", "K2"],
        assignment_policy=SIGNALOME_ASSIGNMENT_POLICY_CUTOFF_BINARY,
    )

    assert weighted.loc[1, "K1"] == 66.667
    assert weighted.loc[1, "K2"] == 33.333
    assert weighted.loc[2, "K2"] == 100.0
    assert binary.loc[1, "K1"] == 0.0
    assert binary.loc[1, "K2"] == 0.0


def test_build_kinase_network_positive_only_excludes_negative_edges() -> None:
    downstream_scores = pd.DataFrame(
        {
            "K1": [1.0, 2.0, 3.0, 4.0],
            "K2": [4.0, 3.0, 2.0, 1.0],
            "K3": [1.0, 2.0, 2.0, 3.0],
        },
        index=pd.Index(["S1", "S2", "S3", "S4"], name="site_id"),
    )

    edges, nodes = build_kinase_network(
        downstream_score_matrix=downstream_scores,
        kinase_order=["K1", "K2", "K3"],
        kinase_substrates={"K1": ("S1", "S2"), "K2": ("S3",), "K3": ()},
        threshold=0.9,
        network_policy="positive_only",
        min_paired_observations=3,
    )

    pdt.assert_frame_equal(
        edges,
        pd.DataFrame(
            {
                "source_kinase": ["K1"],
                "target_kinase": ["K3"],
                "correlation": [0.9486832980505138],
                "valid_observations": [4],
            }
        ).astype(
            {
                "source_kinase": str,
                "target_kinase": str,
                "correlation": float,
                "valid_observations": "int64",
            }
        ),
    )
    pdt.assert_frame_equal(
        nodes,
        pd.DataFrame(
            {
                "degree": [1, 0, 1],
                "n_substrates": [2, 1, 0],
            },
            index=pd.Index(["K1", "K2", "K3"], name="kinase"),
        ).astype({"degree": "int64", "n_substrates": "int64"}),
    )


def test_build_kinase_network_absolute_threshold_uses_unsigned_correlation_values() -> (
    None
):
    downstream_scores = pd.DataFrame(
        {
            "K1": [1.0, 2.0, 3.0, 4.0],
            "K2": [4.0, 3.0, 2.0, 1.0],
            "K3": [1.0, 2.0, 2.0, 3.0],
        },
        index=pd.Index(["S1", "S2", "S3", "S4"], name="site_id"),
    )

    edges, _ = build_kinase_network(
        downstream_score_matrix=downstream_scores,
        kinase_order=["K1", "K2", "K3"],
        kinase_substrates={"K1": (), "K2": (), "K3": ()},
        threshold=0.9,
        network_policy="absolute_threshold",
        min_paired_observations=3,
    )

    expected = pd.DataFrame(
        {
            "source_kinase": ["K1", "K1", "K2"],
            "target_kinase": ["K2", "K3", "K3"],
            "correlation": [
                1.0,
                0.9486832980505138,
                0.9486832980505138,
            ],
            "valid_observations": [4, 4, 4],
        }
    ).astype(
        {
            "source_kinase": str,
            "target_kinase": str,
            "correlation": float,
            "valid_observations": "int64",
        }
    )
    pdt.assert_frame_equal(edges, expected)


def test_build_kinase_network_signed_policy_retains_negative_edges_with_sign() -> None:
    downstream_scores = pd.DataFrame(
        {
            "K1": [1.0, 2.0, 3.0, 4.0],
            "K2": [4.0, 3.0, 2.0, 1.0],
            "K3": [1.0, 2.0, 2.0, 3.0],
        },
        index=pd.Index(["S1", "S2", "S3", "S4"], name="site_id"),
    )

    edges, _ = build_kinase_network(
        downstream_score_matrix=downstream_scores,
        kinase_order=["K1", "K2", "K3"],
        kinase_substrates={"K1": (), "K2": (), "K3": ()},
        threshold=0.9,
        network_policy="signed",
        min_paired_observations=3,
    )

    expected = pd.DataFrame(
        {
            "source_kinase": ["K1", "K1", "K2"],
            "target_kinase": ["K2", "K3", "K3"],
            "correlation": [
                -1.0,
                0.9486832980505138,
                -0.9486832980505138,
            ],
            "valid_observations": [4, 4, 4],
        }
    ).astype(
        {
            "source_kinase": str,
            "target_kinase": str,
            "correlation": float,
            "valid_observations": "int64",
        }
    )
    pdt.assert_frame_equal(edges, expected)


def test_build_kinase_network_rejects_unsupported_network_policy() -> None:
    downstream_scores = pd.DataFrame(
        {
            "K1": [1.0, 2.0],
            "K2": [2.0, 1.0],
        },
        index=pd.Index(["S1", "S2"], name="site_id"),
    )

    with pytest.raises(
        WorkflowStageError,
        match="unsupported network_policy 'unsupported'",
    ):
        build_kinase_network(
            downstream_score_matrix=downstream_scores,
            kinase_order=["K1", "K2"],
            kinase_substrates={"K1": (), "K2": ()},
            threshold=0.5,
            network_policy="unsupported",  # type: ignore[arg-type]
        )


def test_candidate_correlation_true_zero_is_finite_zero() -> None:
    downstream_scores = pd.DataFrame(
        {
            "K1": [1.0, -1.0, 1.0, -1.0],
            "K2": [1.0, 1.0, -1.0, -1.0],
        },
        index=pd.Index(["S1", "S2", "S3", "S4"], name="site_id"),
    )

    _, _, candidates, _ = build_kinase_network_with_diagnostics(
        downstream_score_matrix=downstream_scores,
        kinase_order=["K1", "K2"],
        kinase_substrates={"K1": (), "K2": ()},
        threshold=0.0,
        network_policy="signed",
        min_paired_observations=3,
    )

    assert candidates.shape[0] == 1
    assert candidates.at[0, "correlation"] == 0.0
    assert candidates.at[0, "correlation_status"] == "finite"


def test_candidate_correlation_constant_profile_remains_nan_and_creates_no_edge() -> (
    None
):
    downstream_scores = pd.DataFrame(
        {
            "K1": [1.0, 1.0, 1.0, 1.0],
            "K2": [1.0, 2.0, 3.0, 4.0],
        },
        index=pd.Index(["S1", "S2", "S3", "S4"], name="site_id"),
    )

    edges, _, candidates, diagnostics = build_kinase_network_with_diagnostics(
        downstream_score_matrix=downstream_scores,
        kinase_order=["K1", "K2"],
        kinase_substrates={"K1": (), "K2": ()},
        threshold=0.0,
        network_policy="signed",
        min_paired_observations=3,
    )

    assert edges.empty
    assert candidates.shape[0] == 1
    assert pd.isna(candidates.at[0, "correlation"])
    assert candidates.at[0, "correlation_status"] == "constant_profile"
    assert diagnostics.edges_skipped_constant_profile == 1


def test_candidate_correlation_insufficient_observations_remains_nan() -> None:
    downstream_scores = pd.DataFrame(
        {
            "K1": [1.0],
            "K2": [2.0],
        },
        index=pd.Index(["S1"], name="site_id"),
    )

    edges, _, candidates, _ = build_kinase_network_with_diagnostics(
        downstream_score_matrix=downstream_scores,
        kinase_order=["K1", "K2"],
        kinase_substrates={"K1": (), "K2": ()},
        threshold=0.0,
        network_policy="signed",
    )

    assert edges.empty
    assert candidates.shape[0] == 1
    assert pd.isna(candidates.at[0, "correlation"])
    assert candidates.at[0, "correlation_status"] == "insufficient_observations"


def test_candidate_correlation_two_paired_values_is_insufficient_for_new_edges() -> (
    None
):
    downstream_scores = pd.DataFrame(
        {
            "K1": [1.0, 2.0],
            "K2": [3.0, 4.0],
        },
        index=pd.Index(["S1", "S2"], name="site_id"),
    )

    edges, _, candidates, diagnostics = build_kinase_network_with_diagnostics(
        downstream_score_matrix=downstream_scores,
        kinase_order=["K1", "K2"],
        kinase_substrates={"K1": (), "K2": ()},
        threshold=0.0,
        network_policy="signed",
        min_paired_observations=3,
    )

    assert edges.empty
    assert candidates.at[0, "valid_observations"] == 2
    assert pd.isna(candidates.at[0, "correlation"])
    assert candidates.at[0, "correlation_status"] == "insufficient_observations"
    assert diagnostics.edges_skipped_insufficient_paired_observations == 1


def test_candidate_correlation_rejects_legacy_two_observation_threshold() -> None:
    downstream_scores = pd.DataFrame(
        {
            "K1": [1.0, 2.0],
            "K2": [3.0, 4.0],
        },
        index=pd.Index(["S1", "S2"], name="site_id"),
    )

    with pytest.raises(
        WorkflowStageError,
        match="Legacy threshold 2 cannot be used for new signalome network execution",
    ):
        build_kinase_network_with_diagnostics(
            downstream_score_matrix=downstream_scores,
            kinase_order=["K1", "K2"],
            kinase_substrates={"K1": (), "K2": ()},
            threshold=0.0,
            network_policy="signed",
            min_paired_observations=2,
        )


def test_candidate_correlation_default_minimum_requires_five_observations() -> None:
    downstream_scores = pd.DataFrame(
        {
            "K1": [1.0, 2.0, 3.0, 4.0],
            "K2": [1.0, 2.0, 3.0, 4.0],
        },
        index=pd.Index(["S1", "S2", "S3", "S4"], name="site_id"),
    )

    edges, _, candidates, diagnostics = build_kinase_network_with_diagnostics(
        downstream_score_matrix=downstream_scores,
        kinase_order=["K1", "K2"],
        kinase_substrates={"K1": (), "K2": ()},
        threshold=0.0,
        network_policy="signed",
    )

    assert edges.empty
    assert candidates.at[0, "valid_observations"] == 4
    assert candidates.at[0, "correlation_status"] == "insufficient_observations"
    assert diagnostics.edges_skipped_insufficient_paired_observations == 1


def test_candidate_correlation_respects_minimum_paired_observations() -> None:
    downstream_scores = pd.DataFrame(
        {
            "K1": [1.0, 2.0, 3.0],
            "K2": [1.0, 2.0, 3.0],
        },
        index=pd.Index(["S1", "S2", "S3"], name="site_id"),
    )

    edges, _, candidates, diagnostics = build_kinase_network_with_diagnostics(
        downstream_score_matrix=downstream_scores,
        kinase_order=["K1", "K2"],
        kinase_substrates={"K1": (), "K2": ()},
        threshold=0.0,
        network_policy="signed",
        min_paired_observations=4,
    )

    assert edges.empty
    assert candidates.at[0, "valid_observations"] == 3
    assert candidates.at[0, "correlation_status"] == "insufficient_observations"
    assert diagnostics.insufficient_observation_correlations == 1
    assert diagnostics.edges_skipped_insufficient_paired_observations == 1


def test_candidate_correlation_threshold_three_permits_three_observation_edge() -> None:
    downstream_scores = pd.DataFrame(
        {
            "K1": [1.0, 2.0, 3.0],
            "K2": [1.0, 2.0, 3.0],
        },
        index=pd.Index(["S1", "S2", "S3"], name="site_id"),
    )

    edges, _, candidates, diagnostics = build_kinase_network_with_diagnostics(
        downstream_score_matrix=downstream_scores,
        kinase_order=["K1", "K2"],
        kinase_substrates={"K1": (), "K2": ()},
        threshold=0.9,
        network_policy="signed",
        min_paired_observations=3,
    )

    assert candidates.at[0, "correlation_status"] == "finite"
    assert candidates.at[0, "valid_observations"] == 3
    assert edges.to_dict("records") == [
        {
            "source_kinase": "K1",
            "target_kinase": "K2",
            "correlation": pytest.approx(1.0),
            "valid_observations": 3,
        }
    ]
    assert diagnostics.edges_created == 1


def test_candidate_correlation_classifies_missing_and_non_finite_inputs() -> None:
    missing_scores = pd.DataFrame(
        {
            "K1": [1.0, float("nan"), float("nan")],
            "K2": [2.0, 3.0, float("nan")],
        },
        index=pd.Index(["S1", "S2", "S3"], name="site_id"),
    )
    non_finite_scores = pd.DataFrame(
        {
            "K1": [1.0, float("inf"), 3.0],
            "K2": [1.0, 2.0, 3.0],
        },
        index=pd.Index(["S1", "S2", "S3"], name="site_id"),
    )

    _, _, missing_candidates, missing_diagnostics = (
        build_kinase_network_with_diagnostics(
            downstream_score_matrix=missing_scores,
            kinase_order=["K1", "K2"],
            kinase_substrates={"K1": (), "K2": ()},
            threshold=0.0,
            network_policy="signed",
        )
    )
    _, _, non_finite_candidates, non_finite_diagnostics = (
        build_kinase_network_with_diagnostics(
            downstream_score_matrix=non_finite_scores,
            kinase_order=["K1", "K2"],
            kinase_substrates={"K1": (), "K2": ()},
            threshold=0.0,
            network_policy="signed",
        )
    )

    assert missing_candidates.at[0, "correlation_status"] == "missing_values"
    assert pd.isna(missing_candidates.at[0, "correlation"])
    assert missing_diagnostics.edges_skipped_missing_score == 1
    assert non_finite_candidates.at[0, "correlation_status"] == "non_finite_values"
    assert pd.isna(non_finite_candidates.at[0, "correlation"])
    assert non_finite_diagnostics.edges_skipped_non_finite_score == 1


def test_network_edge_creation_uses_only_finite_candidate_correlations() -> None:
    downstream_scores = pd.DataFrame(
        {
            "K1": [1.0, -1.0, 1.0, -1.0],
            "K2": [1.0, 1.0, -1.0, -1.0],
            "K3": [1.0, -1.0, 1.0, -1.0],
            "K4": [5.0, 5.0, 5.0, 5.0],
            "K5": [1.0, float("nan"), float("nan"), float("nan")],
        },
        index=pd.Index(["S1", "S2", "S3", "S4"], name="site_id"),
    )

    edges, _, candidates, _ = build_kinase_network_with_diagnostics(
        downstream_score_matrix=downstream_scores,
        kinase_order=["K1", "K2", "K3", "K4", "K5"],
        kinase_substrates={"K1": (), "K2": (), "K3": (), "K4": (), "K5": ()},
        threshold=0.8,
        network_policy="signed",
        min_paired_observations=3,
    )

    assert edges.to_dict("records") == [
        {
            "source_kinase": "K1",
            "target_kinase": "K3",
            "correlation": pytest.approx(1.0),
            "valid_observations": 4,
        }
    ]
    assert "finite" in set(candidates.loc[:, "correlation_status"])
    assert "constant_profile" in set(candidates.loc[:, "correlation_status"])
    assert "missing_values" in set(candidates.loc[:, "correlation_status"])


def test_network_correlation_diagnostics_count_statuses_and_skips() -> None:
    downstream_scores = pd.DataFrame(
        {
            "K1": [1.0, -1.0, 1.0, -1.0],
            "K2": [1.0, 1.0, -1.0, -1.0],
            "K3": [1.0, -1.0, 1.0, -1.0],
            "K4": [5.0, 5.0, 5.0, 5.0],
            "K5": [1.0, float("nan"), float("nan"), float("nan")],
        },
        index=pd.Index(["S1", "S2", "S3", "S4"], name="site_id"),
    )

    edges, _, _, diagnostics = build_kinase_network_with_diagnostics(
        downstream_score_matrix=downstream_scores,
        kinase_order=["K1", "K2", "K3", "K4", "K5"],
        kinase_substrates={"K1": (), "K2": (), "K3": (), "K4": (), "K5": ()},
        threshold=0.8,
        network_policy="signed",
        min_paired_observations=3,
    )

    assert diagnostics.total_candidate_correlations == 10
    assert diagnostics.finite_correlations == 3
    assert diagnostics.undefined_correlations == 7
    assert diagnostics.constant_profile_correlations == 3
    assert diagnostics.insufficient_observation_correlations == 0
    assert diagnostics.missing_value_correlations == 4
    assert diagnostics.non_finite_value_correlations == 0
    assert diagnostics.edges_skipped_non_finite_correlation == 7
    assert diagnostics.edges_skipped_below_threshold == 2
    assert diagnostics.edges_skipped_constant_profile == 3
    assert diagnostics.edges_skipped_missing_score == 4
    assert diagnostics.edges_skipped_non_finite_score == 0
    assert diagnostics.edges_skipped_undefined_correlation == 0
    assert diagnostics.edges_created == int(edges.shape[0]) == 1


@pytest.mark.parametrize(
    ("network_policy", "expected_edges", "expected_threshold_skips"),
    [
        ("positive_only", 1, 2),
        ("absolute_threshold", 3, 0),
        ("signed", 3, 0),
    ],
)
def test_network_edge_diagnostics_follow_network_policy(
    network_policy: str,
    expected_edges: int,
    expected_threshold_skips: int,
) -> None:
    downstream_scores = pd.DataFrame(
        {
            "K1": [1.0, 2.0, 3.0, 4.0],
            "K2": [4.0, 3.0, 2.0, 1.0],
            "K3": [1.0, 2.0, 2.0, 3.0],
        },
        index=pd.Index(["S1", "S2", "S3", "S4"], name="site_id"),
    )

    edges, _, candidates, diagnostics = build_kinase_network_with_diagnostics(
        downstream_score_matrix=downstream_scores,
        kinase_order=["K1", "K2", "K3"],
        kinase_substrates={"K1": (), "K2": (), "K3": ()},
        threshold=0.9,
        network_policy=network_policy,
        min_paired_observations=3,
    )

    assert int(edges.shape[0]) == expected_edges
    assert diagnostics.edges_created == expected_edges
    assert diagnostics.finite_correlations == int(candidates.shape[0]) == 3
    assert diagnostics.edges_skipped_below_threshold == expected_threshold_skips
    assert diagnostics.edges_skipped_non_finite_correlation == 0


def test_network_regression_undefined_correlations_are_not_zero_imputed() -> None:
    downstream_scores = pd.DataFrame(
        {
            "K1": [7.0, 7.0, 7.0, 7.0],
            "K2": [1.0, 2.0, 3.0, 4.0],
        },
        index=pd.Index(["S1", "S2", "S3", "S4"], name="site_id"),
    )

    edges, _, candidates, _ = build_kinase_network_with_diagnostics(
        downstream_score_matrix=downstream_scores,
        kinase_order=["K1", "K2"],
        kinase_substrates={"K1": (), "K2": ()},
        threshold=0.0,
        network_policy="signed",
        min_paired_observations=3,
    )

    assert edges.empty
    assert pd.isna(candidates.at[0, "correlation"])
    assert candidates.at[0, "correlation"] != 0.0
    assert candidates.at[0, "correlation_status"] != "finite"


def test_build_expanded_signalome_table_tracks_membership_and_site_order() -> None:
    module_assignments = _with_site_identity(
        pd.DataFrame(
            {
                "protein_group_id": ["P1", "P2", "P3", "P4"],
                "module_id": [2, 1, 2, 3],
                "top_kinase": ["K2", "K1", "K2", "K3"],
                "top_score": [0.91, 0.93, 0.92, 0.88],
                "top_kinase_weights": [
                    (("K2", 1.0),),
                    (("K1", 1.0),),
                    (("K2", 1.0),),
                    (("K3", 1.0),),
                ],
            },
            index=pd.Index(["S3", "S1", "S4", "S2"], name="site_key"),
        )
    )
    signalome_modules = pd.DataFrame(
        {
            "K1": [5.0, 85.0, 0.0],
            "K2": [90.0, 10.0, 0.0],
            "K3": [0.0, 0.0, 100.0],
        },
        index=pd.Index([1, 2, 3], name="module_id"),
    )
    kinase_network_edges = pd.DataFrame(
        {
            "source_kinase": ["K1"],
            "target_kinase": ["K2"],
            "correlation": [0.95],
            "valid_observations": [4],
        }
    )
    kinase_substrates = {
        "K1": ("S1", "S2"),
        "K2": ("S3", "S4"),
    }

    expanded = build_expanded_signalome_table(
        module_assignments=module_assignments,
        signalome_modules=signalome_modules,
        kinase_network_edges=kinase_network_edges,
        kinase_substrates=kinase_substrates,
        assignment_policy=SIGNALOME_ASSIGNMENT_POLICY_CUTOFF_BINARY,
    )

    k1_sites = expanded.loc[
        (expanded.loc[:, "kinase"] == "K1")
        & (
            expanded.loc[:, EXPANDED_SIGNALOME_ROW_KIND_COLUMN]
            == EXPANDED_SIGNALOME_ROW_KIND_SITE
        )
        & (expanded.loc[:, "site_id"] != ""),
        :,
    ]
    assert k1_sites.loc[:, "site_id"].tolist() == ["S3", "S1", "S4"]
    assert k1_sites.loc[:, "site_order"].tolist() == [0, 1, 2]
    assert k1_sites.loc[:, "linked_kinases"].iloc[0] == '["K1","K2"]'
    assert k1_sites.loc[:, "regulated_module_ids"].iloc[0] == "[1,2]"


def test_build_expanded_signalome_table_weighted_top_uses_fractional_site_support() -> (
    None
):
    module_assignments = _with_site_identity(
        pd.DataFrame(
            {
                "protein_group_id": ["P1", "P1", "P2"],
                "module_id": [1, 1, 2],
                "top_kinase": ["K1", "K1", "K2"],
                "top_score": [0.95, 0.96, 0.92],
                "top_kinase_weights": [
                    (("K1", 0.5), ("K2", 0.5)),
                    (("K1", 1.0),),
                    (("K2", 1.0),),
                ],
            },
            index=pd.Index(["S1", "S2", "S3"], name="site_key"),
        )
    )
    signalome_modules = pd.DataFrame(
        {"K1": [90.0, 0.0], "K2": [10.0, 100.0]},
        index=pd.Index([1, 2], name="module_id"),
    )
    kinase_network_edges = pd.DataFrame(
        {
            "source_kinase": ["K1"],
            "target_kinase": ["K2"],
            "correlation": [0.91],
            "valid_observations": [3],
        }
    )

    expanded = build_expanded_signalome_table(
        module_assignments=module_assignments,
        signalome_modules=signalome_modules,
        kinase_network_edges=kinase_network_edges,
        kinase_substrates={"K1": (), "K2": ()},
        assignment_policy=SIGNALOME_ASSIGNMENT_POLICY_WEIGHTED_TOP,
    )

    site_s1 = expanded.loc[
        (expanded.loc[:, "kinase"] == "K1") & (expanded.loc[:, "site_id"] == "S1"),
        :,
    ].iloc[0]
    assert site_s1["support_kinases"] == '["K1","K2"]'
    assert float(site_s1["support_weight"]) == 1.0


def test_build_expanded_signalome_table_emits_expected_shape_for_site_and_summary_rows() -> (
    None
):
    module_assignments = _with_site_identity(
        pd.DataFrame(
            {
                "protein_group_id": ["P1", "P2"],
                "module_id": [1, 2],
                "top_kinase": ["K1", "K2"],
                "top_score": [0.9, 0.8],
                "top_kinase_weights": [
                    (("K1", 1.0),),
                    (("K2", 1.0),),
                ],
            },
            index=pd.Index(["S1", "S2"], name="site_key"),
        )
    )
    signalome_modules = pd.DataFrame(
        {
            "K1": [90.0, 0.0, 0.0],
            "K2": [0.0, 90.0, 0.0],
            "K3": [0.0, 0.0, 0.0],
        },
        index=pd.Index([1, 2, 3], name="module_id"),
    )
    kinase_network_edges = pd.DataFrame(
        columns=[
            "source_kinase",
            "target_kinase",
            "correlation",
            "valid_observations",
        ]
    )

    expanded = build_expanded_signalome_table(
        module_assignments=module_assignments,
        signalome_modules=signalome_modules,
        kinase_network_edges=kinase_network_edges,
        kinase_substrates={"K1": ("S1",), "K2": ("S2",), "K3": ()},
        assignment_policy=SIGNALOME_ASSIGNMENT_POLICY_CUTOFF_BINARY,
    )

    assert expanded.loc[:, "kinase"].tolist() == ["K1", "K2", "K3"]
    assert expanded.loc[:, EXPANDED_SIGNALOME_ROW_KIND_COLUMN].tolist() == [
        EXPANDED_SIGNALOME_ROW_KIND_SITE,
        EXPANDED_SIGNALOME_ROW_KIND_SITE,
        EXPANDED_SIGNALOME_ROW_KIND_SUMMARY,
    ]
    assert expanded.loc[:, "site_id"].tolist() == ["S1", "S2", ""]
    assert expanded.loc[:, "support_kinases"].tolist() == ['["K1"]', '["K2"]', "[]"]
    assert expanded.loc[:, "regulated_module_ids"].tolist() == ["[1]", "[2]", "[]"]
