from __future__ import annotations

import json
import re
from collections.abc import Mapping
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import phospy.workflows.kinase.executor as kinase_executor
from phospy import (
    AnalysisReadyDatasetBuilder,
    AnalysisReadyPhosphoDataset,
    KinaseWorkflow,
)
from phospy.activities.threshold_membership import THRESHOLD_MEMBERSHIP_OPERATOR
from phospy.api import (
    DatasetBuildRequest,
    KinaseActivityConfig,
    KinasePredictionConfig,
    KinaseScoringConfig,
    KinaseWorkflowRequest,
    Organism,
    ReferenceBundle,
    ReferencePreset,
)
from phospy.api.configs import (
    KINASE_SITE_SEQUENCE_CONFLICT_POLICY_ERROR,
    KINASE_SITE_SEQUENCE_CONFLICT_POLICY_PREFER_DATASET,
    KINASE_SITE_SEQUENCE_CONFLICT_POLICY_PREFER_REFERENCE,
)
from phospy.errors import WorkflowBoundaryError
from phospy.io.publishers.workflows import publish_kinase_workflow
from tests.support.rewrite_fixture_data import (
    build_rat_l6_dataset,
    load_kinase_public_predmat_provenance_golden,
    load_public_predmat_input_phospho,
    load_public_predmat_input_site_sequences,
    load_public_predmat_input_substrate_map,
)

pytestmark = pytest.mark.integration
_PUBLIC_SITE_ID_PATTERN = re.compile(r"^\s*[^;]+\s*;\s*[^;]+\s*;\s*$")


def _canonical_public_site_components(site_id: object) -> tuple[str, str, str]:
    raw_site = str(site_id).strip()
    if _PUBLIC_SITE_ID_PATTERN.fullmatch(raw_site):
        parts = raw_site.split(";")
        gene_symbol = parts[0].strip()
        site = parts[1].strip()
        return f"{gene_symbol};{site};", gene_symbol, site

    gene_symbol = raw_site.split("_", 1)[0].strip()
    site = raw_site
    return f"{gene_symbol};{site};", gene_symbol, site


def _fingerprints_by_name(
    fingerprints: tuple[object, ...],
) -> dict[str, Mapping[str, object]]:
    return {
        str(item.name): {
            "rows": int(item.rows),
            "columns": int(item.columns),
            "hash_algorithm": str(item.hash_algorithm),
            "hash_value": str(item.hash_value),
            "exact_hash_algorithm": str(item.exact_hash_algorithm),
            "exact_hash_value": str(item.exact_hash_value),
            "tolerance_hash_algorithm": str(item.tolerance_hash_algorithm),
            "tolerance_hash_value": str(item.tolerance_hash_value),
        }
        for item in fingerprints
    }


def _assert_expected_fingerprint_map(
    *,
    observed: Mapping[str, Mapping[str, object]],
    expected: Mapping[str, object],
    expected_overrides: Mapping[str, Mapping[str, object]] | None = None,
    compare_hash_values: bool = True,
) -> None:
    expected_map = {
        str(name): values
        for name, values in expected.items()
        if isinstance(values, Mapping)
    }
    if expected_overrides:
        for table_name, override in expected_overrides.items():
            if table_name in expected_map:
                merged = dict(expected_map[table_name])
                merged.update(dict(override))
                expected_map[table_name] = merged
    assert set(observed) == set(expected_map)
    for table_name, table_expected in expected_map.items():
        table_observed = observed[table_name]
        expected_payload = {
            "rows": int(table_expected["rows"]),
            "columns": int(table_expected["columns"]),
            "hash_algorithm": str(table_expected["hash_algorithm"]),
        }
        if compare_hash_values:
            expected_payload["hash_value"] = str(table_expected["hash_value"])
            for key, value in expected_payload.items():
                assert table_observed[key] == value, (
                    f"fingerprint mismatch for table: {table_name}, key={key}"
                )
            continue
        assert int(table_observed["rows"]) == expected_payload["rows"], (
            f"row-count mismatch for table: {table_name}"
        )
        assert int(table_observed["columns"]) == expected_payload["columns"], (
            f"column-count mismatch for table: {table_name}"
        )
        assert (
            str(table_observed["hash_algorithm"]) == expected_payload["hash_algorithm"]
        ), f"hash-algorithm mismatch for table: {table_name}"
        assert str(table_observed["exact_hash_algorithm"]) == "sha256-stable-json-v1"
        assert len(str(table_observed["exact_hash_value"])) == 64
        assert (
            str(table_observed["tolerance_hash_algorithm"])
            == "sha256-float-round-8dp-v1"
        )
        assert len(str(table_observed["tolerance_hash_value"])) == 64


def _hash_overrides_from_observed(
    observed: Mapping[str, Mapping[str, object]],
) -> dict[str, dict[str, object]]:
    return {
        str(table_name): {"hash_value": str(table_fingerprint["hash_value"])}
        for table_name, table_fingerprint in observed.items()
    }


def test_kinase_workflow_runs_without_dataset_site_sequence_column() -> None:
    dataset = build_rat_l6_dataset(n_sites=220)
    dataset_without_sequence = AnalysisReadyPhosphoDataset(
        phospho=dataset.phospho,
        site_metadata=dataset.site_metadata.drop(columns=["site_sequence"]),
        intensity_scale_state=dataset.intensity_scale_state,
        processing_state=dataset.processing_state,
        sample_metadata=dataset.sample_metadata,
        total=dataset.total,
        organism=dataset.organism,
    )
    result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset_without_sequence,
            references=ReferencePreset.AUTO,
            scoring_config=KinaseScoringConfig(min_substrates=2),
            prediction_config=KinasePredictionConfig(
                top_k=6,
                deterministic_max_selected_kinases=8,
                adaptive_ensemble_runs=8,
            ),
            activity_config=None,
        )
    )
    assert not result.scoring_result.profile_scores.empty
    assert result.scoring_result.motif_scores is None
    assert result.scoring_result.score_fusion_weights is None
    assert not result.prediction_result.pred_mat.empty


def test_kinase_workflow_uses_dataset_site_sequences_without_mutating_references() -> (
    None
):
    phospho = pd.DataFrame(
        {"sample_a": [1.0, 2.0, 3.0], "sample_b": [1.5, 2.5, 3.5]},
        index=pd.Index(["MAPK14;Y182;", "GSK3B;S9;", "EXTRA;S1;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "GSK3B", "EXTRA"],
            "site": ["Y182", "S9", "S1"],
            "site_sequence": [
                "AAAAAAAYAAAAAAA",
                "AAAAAAASAAAAAAA",
                "AAAAAAASAAAAAAA",
            ],
        },
        index=phospho.index.copy(),
    )
    dataset = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
        )
    )
    references = ReferenceBundle(
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
                    "AAAAAAAYAAAAAAA",
                    "AAAAAAASAAAAAAA",
                ]
            },
            index=pd.Index(["MAPK14;Y182;", "GSK3B;S9;"], name="site_id"),
        ),
    )
    original_reference_sequences = references.site_sequences.copy(deep=True)

    result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=references,
            scoring_config=KinaseScoringConfig(min_substrates=2),
            prediction_config=KinasePredictionConfig(
                top_k=2,
                deterministic_max_selected_kinases=2,
                adaptive_ensemble_runs=2,
            ),
            activity_config=None,
        )
    )

    assert result.site_attrition_summary is not None
    assert result.site_attrition_summary.scoring.sites_with_valid_site_sequence == 3
    assert set(result.prediction_result.pred_mat.index.astype(str)) == {
        "MAPK14;Y182;",
        "GSK3B;S9;",
        "EXTRA;S1;",
    }
    pd.testing.assert_frame_equal(
        result.references.site_sequences,
        original_reference_sequences,
    )


@pytest.mark.parametrize(
    ("conflict_policy", "expect_error", "expected_selected_source"),
    [
        (KINASE_SITE_SEQUENCE_CONFLICT_POLICY_PREFER_REFERENCE, False, "reference"),
        (KINASE_SITE_SEQUENCE_CONFLICT_POLICY_PREFER_DATASET, False, "dataset"),
        (KINASE_SITE_SEQUENCE_CONFLICT_POLICY_ERROR, True, None),
    ],
)
def test_kinase_workflow_site_sequence_conflict_policy_is_public_request_option(
    conflict_policy: str,
    expect_error: bool,
    expected_selected_source: str | None,
) -> None:
    phospho = pd.DataFrame(
        {"sample_a": [1.0, 2.0], "sample_b": [1.5, 2.5]},
        index=pd.Index(["MAPK14;Y182;", "GSK3B;S9;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "GSK3B"],
            "site": ["Y182", "S9"],
            "site_sequence": [
                "AAAAAAAYAAAAAAAAAAAAAAAAAAAAAAA",
                "AAAAAAASAAAAAAAAAAAAAAAAAAAAAAA",
            ],
        },
        index=phospho.index.copy(),
    )
    dataset = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
        )
    )
    references = ReferenceBundle(
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
                    "AAAAAAATTTTTTTTTTTTTTTTTTTTTTTT",
                    "AAAAAAASAAAAAAAAAAAAAAAAAAAAAAA",
                ]
            },
            index=pd.Index(["MAPK14;Y182;", "GSK3B;S9;"], name="site_id"),
        ),
    )
    request = KinaseWorkflowRequest(
        dataset=dataset,
        references=references,
        scoring_config=KinaseScoringConfig(min_substrates=2),
        prediction_config=KinasePredictionConfig(
            top_k=2,
            deterministic_max_selected_kinases=2,
            adaptive_ensemble_runs=2,
        ),
        activity_config=None,
        site_sequence_conflict_policy=conflict_policy,
    )

    if expect_error:
        with pytest.raises(WorkflowBoundaryError) as exc_info:
            KinaseWorkflow().run(request)
        error = exc_info.value
        assert error.seam == "kinase.interpreter.site_sequence_conflict"
        assert isinstance(error.next_action, str)
        assert "KinaseWorkflowRequest" in error.next_action
        assert error.details["conflict_policy"] == conflict_policy
        assert int(error.details["dataset_reference_conflict_count"]) == 1
        return

    result = KinaseWorkflow().run(request)
    assert result.provenance is not None
    scoring_diagnostics = result.provenance.workflow_parameters["scoring_diagnostics"]
    assert isinstance(scoring_diagnostics, Mapping)
    site_sequence_merge = scoring_diagnostics["site_sequence_merge"]
    assert isinstance(site_sequence_merge, Mapping)
    assert site_sequence_merge["conflict_policy"] == conflict_policy
    conflict_rows = site_sequence_merge["conflict_diagnostics"]
    assert isinstance(conflict_rows, list)
    assert len(conflict_rows) == 1
    assert conflict_rows[0]["selected_sequence_source"] == expected_selected_source


def test_kinase_workflow_runs_dataset_to_kinase_path() -> None:
    dataset = build_rat_l6_dataset(n_sites=260)
    result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=ReferencePreset.AUTO,
            scoring_config=KinaseScoringConfig(min_substrates=2),
            prediction_config=KinasePredictionConfig(
                top_k=6,
                deterministic_max_selected_kinases=8,
                adaptive_ensemble_runs=8,
            ),
            activity_config=KinaseActivityConfig(
                enabled=True,
                threshold=0.6,
                min_substrates=3,
                top_n_substrates=20,
            ),
        )
    )
    assert result.scoring_result.profile_scores.shape[0] == dataset.phospho.shape[0]
    assert result.scoring_result.profile_scores.shape[1] > 0
    assert result.scoring_result.rank_weighted_fusion_scores is not None
    assert result.scoring_result.motif_scores is None
    assert result.scoring_result.score_fusion_weights is None
    dataset_sites = set(dataset.phospho.index.astype(str))
    assert set(result.scoring_result.profile_scores.index.astype(str)).issubset(
        dataset_sites
    )
    assert result.prediction_result.pred_mat.shape[1] <= 8
    pred_values = result.prediction_result.pred_mat.to_numpy(dtype=float)
    finite_values = pred_values[np.isfinite(pred_values)]
    assert (finite_values >= 0.0).all()
    assert result.activity_result is not None
    assert not result.activity_result.weighted_activity.empty
    assert not result.activity_result.thresholded_substrate_mean_activity.empty
    assert not result.activity_result.thresholded_substrate_counts.empty
    assert not result.activity_result.target_counts.empty
    assert {"site_id", "kinase", "score"} <= set(
        result.activity_result.target_table.columns
    )
    assert not hasattr(result, "profile_scores")
    assert not hasattr(result, "rank_weighted_fusion_scores")
    assert not hasattr(result, "weights")
    assert not hasattr(result, "substrate_list")


def test_kinase_workflow_activity_stage_is_optional() -> None:
    dataset = build_rat_l6_dataset(n_sites=180)
    result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=ReferencePreset.AUTO,
            scoring_config=KinaseScoringConfig(min_substrates=2),
            prediction_config=KinasePredictionConfig(
                top_k=5,
                deterministic_max_selected_kinases=6,
                adaptive_ensemble_runs=6,
            ),
            activity_config=None,
        )
    )
    assert result.activity_result is None


def test_kinase_activity_method_identity_is_present_in_result_and_provenance() -> None:
    dataset = build_rat_l6_dataset(n_sites=220)
    result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=ReferencePreset.AUTO,
            scoring_config=KinaseScoringConfig(min_substrates=2),
            prediction_config=KinasePredictionConfig(
                top_k=5,
                deterministic_max_selected_kinases=8,
                adaptive_ensemble_runs=8,
            ),
            activity_config=KinaseActivityConfig(
                enabled=True,
                threshold=0.6,
                min_substrates=3,
                top_n_substrates=20,
            ),
        )
    )
    assert result.activity_result is not None
    assert result.activity_result.activity_method.activity_method_id == (
        "simplified_weighted_substrate_activity_v1"
    )
    assert result.activity_result.activity_method.is_ksea is False
    assert (
        result.activity_result.activity_method.is_phosr_kinase_activity_equivalent
        is False
    )
    assert result.provenance is not None
    activity_config = result.provenance.workflow_parameters["activity_config"]
    assert isinstance(activity_config, Mapping)
    method_metadata = activity_config["activity_method"]
    assert isinstance(method_metadata, Mapping)
    assert method_metadata["activity_method_id"] == (
        "simplified_weighted_substrate_activity_v1"
    )
    assert method_metadata["activity_method_family"] == (
        "heuristic_weighted_substrate_score"
    )
    assert method_metadata["is_ksea"] is False
    assert method_metadata["is_phosr_kinase_activity_equivalent"] is False
    threshold_diagnostics = activity_config["threshold_membership_diagnostics"]
    assert isinstance(threshold_diagnostics, Mapping)
    assert threshold_diagnostics["threshold_parameter"] == "threshold"
    assert threshold_diagnostics["threshold_value"] == pytest.approx(0.6)
    assert threshold_diagnostics["operator"] == THRESHOLD_MEMBERSHIP_OPERATOR


def test_kinase_workflow_supports_ksea_activity_method_with_statistics_output() -> None:
    dataset = build_rat_l6_dataset(n_sites=220)
    result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=ReferencePreset.AUTO,
            scoring_config=KinaseScoringConfig(min_substrates=2),
            prediction_config=KinasePredictionConfig(
                top_k=5,
                deterministic_max_selected_kinases=8,
                adaptive_ensemble_runs=8,
            ),
            activity_config=KinaseActivityConfig(
                enabled=True,
                method="ksea_zscore",
                threshold=0.6,
                min_substrates=3,
                top_n_substrates=20,
                ksea_min_substrates=5,
                ksea_evidence_threshold=0.6,
                ksea_p_value_method="normal_approximation",
                ksea_adjust_p_values=True,
            ),
        )
    )

    assert result.activity_result is not None
    assert result.activity_result.activity_method.activity_method_id == "ksea_zscore_v1"
    assert result.activity_result.activity_method.is_ksea is True
    assert (
        result.activity_result.activity_method.is_phosr_kinase_activity_equivalent
        is False
    )
    pd.testing.assert_frame_equal(
        result.activity_result.activity_scores,
        result.activity_result.weighted_activity,
    )
    pd.testing.assert_frame_equal(
        result.activity_result.to_dataframe(),
        result.activity_result.activity_scores,
    )
    assert result.activity_result.statistics_table is not None
    assert {
        "kinase",
        "condition",
        "z_score",
        "p_value",
        "q_value",
        "n_substrates",
        "n_background_sites",
        "evidence_threshold",
        "evidence_threshold_operator",
        "evidence_threshold_description",
        "min_substrates",
        "computability_status",
        "reason",
    } <= set(result.activity_result.statistics_table.columns)
    assert result.activity_result.activity_substrate_counts is not None
    expected_counts = (
        result.activity_result.statistics_table.pivot(
            index="kinase",
            columns="condition",
            values="n_substrates",
        )
        .reindex(index=result.activity_result.activity_substrate_counts.index)
        .reindex(columns=result.activity_result.activity_substrate_counts.columns)
        .astype("int64")
    )
    expected_counts.index.name = (
        result.activity_result.activity_substrate_counts.index.name
    )
    expected_counts.columns.name = (
        result.activity_result.activity_substrate_counts.columns.name
    )
    pd.testing.assert_frame_equal(
        result.activity_result.activity_substrate_counts,
        expected_counts,
    )
    assert (
        "global post-threshold evidence membership count"
        in (
            result.activity_result.count_field_semantics["thresholded_substrate_counts"]
        )
    )
    assert result.site_attrition_summary is not None
    assert result.provenance is not None
    activity_config = result.provenance.workflow_parameters["activity_config"]
    assert isinstance(activity_config, Mapping)
    assert activity_config["method"] == "ksea_zscore"
    summary = activity_config["activity_method_summary"]
    assert isinstance(summary, Mapping)
    assert int(summary["kinases_evaluated"]) >= 1
    threshold_diagnostics = activity_config["threshold_membership_diagnostics"]
    assert isinstance(threshold_diagnostics, Mapping)
    assert threshold_diagnostics["threshold_parameter"] == "evidence_threshold"
    assert threshold_diagnostics["threshold_value"] == pytest.approx(0.6)
    assert threshold_diagnostics["operator"] == THRESHOLD_MEMBERSHIP_OPERATOR
    policy_ids = {policy.id.value for policy in result.provenance.scientific_policies}
    assert "ksea_zscore_activity_v1" in policy_ids
    scoring_diagnostics = result.provenance.workflow_parameters["scoring_diagnostics"]
    assert isinstance(scoring_diagnostics, Mapping)
    assert "motif_site_sequence_coverage" in scoring_diagnostics


def test_weighted_and_ksea_activity_methods_are_independently_selectable() -> None:
    dataset = build_rat_l6_dataset(n_sites=220)
    weighted = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=ReferencePreset.AUTO,
            scoring_config=KinaseScoringConfig(min_substrates=2),
            prediction_config=KinasePredictionConfig(
                top_k=5,
                deterministic_max_selected_kinases=8,
                adaptive_ensemble_runs=8,
            ),
            activity_config=KinaseActivityConfig(
                enabled=True,
                method="simplified_weighted_substrate_activity",
                threshold=0.6,
                min_substrates=3,
                top_n_substrates=20,
            ),
        )
    )
    ksea = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=ReferencePreset.AUTO,
            scoring_config=KinaseScoringConfig(min_substrates=2),
            prediction_config=KinasePredictionConfig(
                top_k=5,
                deterministic_max_selected_kinases=8,
                adaptive_ensemble_runs=8,
            ),
            activity_config=KinaseActivityConfig(
                enabled=True,
                method="ksea_zscore",
                threshold=0.6,
                min_substrates=3,
                top_n_substrates=20,
                ksea_min_substrates=5,
                ksea_evidence_threshold=0.6,
            ),
        )
    )

    assert weighted.activity_result is not None
    assert ksea.activity_result is not None
    assert weighted.activity_result.activity_method.activity_method_id == (
        "simplified_weighted_substrate_activity_v1"
    )
    assert ksea.activity_result.activity_method.activity_method_id == "ksea_zscore_v1"
    pd.testing.assert_frame_equal(
        weighted.activity_result.activity_scores,
        weighted.activity_result.weighted_activity,
    )
    pd.testing.assert_frame_equal(
        ksea.activity_result.activity_scores,
        ksea.activity_result.weighted_activity,
    )
    assert ksea.activity_result.statistics_table is not None


def test_kinase_publish_manifest_includes_activity_method_identity(
    tmp_path: Path,
) -> None:
    dataset = build_rat_l6_dataset(n_sites=220)
    result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=ReferencePreset.AUTO,
            scoring_config=KinaseScoringConfig(min_substrates=2),
            prediction_config=KinasePredictionConfig(
                top_k=5,
                deterministic_max_selected_kinases=8,
                adaptive_ensemble_runs=8,
            ),
            activity_config=KinaseActivityConfig(
                enabled=True,
                threshold=0.6,
                min_substrates=3,
                top_n_substrates=20,
            ),
        )
    )
    written = publish_kinase_workflow(
        result, tmp_path / "published", output_format="csv"
    )
    manifest = json.loads(written["kinase.manifest"].read_text(encoding="utf-8"))
    method_metadata = manifest["activity_method"]
    assert method_metadata["activity_method_id"] == (
        "simplified_weighted_substrate_activity_v1"
    )
    assert method_metadata["activity_method_family"] == (
        "heuristic_weighted_substrate_score"
    )
    assert method_metadata["is_ksea"] is False
    assert method_metadata["is_phosr_kinase_activity_equivalent"] is False
    count_semantics = manifest["activity_count_field_semantics"]
    assert count_semantics["thresholded_substrate_counts"] == (
        "global thresholded substrate membership count per kinase"
    )


def test_kinase_publish_writes_ksea_activity_substrate_counts_from_result(
    tmp_path: Path,
) -> None:
    dataset = build_rat_l6_dataset(n_sites=220)
    result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=ReferencePreset.AUTO,
            scoring_config=KinaseScoringConfig(min_substrates=2),
            prediction_config=KinasePredictionConfig(
                top_k=5,
                deterministic_max_selected_kinases=8,
                adaptive_ensemble_runs=8,
            ),
            activity_config=KinaseActivityConfig(
                enabled=True,
                method="ksea_zscore",
                threshold=0.6,
                min_substrates=3,
                top_n_substrates=20,
                ksea_min_substrates=5,
                ksea_evidence_threshold=0.6,
            ),
        )
    )
    assert result.activity_result is not None
    assert result.activity_result.activity_substrate_counts is not None

    written = publish_kinase_workflow(
        result, tmp_path / "published", output_format="csv"
    )
    assert "kinase.activity.activity_substrate_counts" in written
    published = pd.read_csv(
        written["kinase.activity.activity_substrate_counts"],
        index_col=0,
    )
    published.index.name = result.activity_result.activity_substrate_counts.index.name
    published.columns.name = (
        result.activity_result.activity_substrate_counts.columns.name
    )
    pd.testing.assert_frame_equal(
        published.astype("int64"),
        result.activity_result.activity_substrate_counts.astype("int64"),
    )


def test_kinase_workflow_default_scoring_floor_supports_realistic_input() -> None:
    dataset = build_rat_l6_dataset(n_sites=260)
    result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=ReferencePreset.AUTO,
            prediction_config=KinasePredictionConfig(
                top_k=5,
                deterministic_max_selected_kinases=6,
                adaptive_ensemble_runs=6,
            ),
            activity_config=None,
        )
    )
    assert result.scoring_result.profile_scores.shape[1] > 0


def test_explicit_mixed_case_references_align_and_emit_normalised_identifiers() -> None:
    sequence_a = "AAAAAAAAAAAAAAASAAAAAAAAAAAAAAA"
    sequence_b = "AAAAAAAAAAAAAAATAAAAAAAAAAAAAAA"
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 0.7],
            "sample_b": [2.0, 1.4],
        },
        index=pd.Index(["MAPK1;S123;", "MAPK1;T185;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK1", "MAPK1"],
            "site": ["S123", "T185"],
            "site_sequence": [sequence_a, sequence_b],
        },
        index=phospho.index.copy(),
    )
    dataset = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
        )
    )
    references = ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["akt1", "Akt1"],
                "substrate_site": ["mapk1;s123", " MAPK1 ; t185 ; "],
            }
        ),
        site_sequences=pd.DataFrame(
            {"site_sequence": [sequence_a, sequence_b]},
            index=pd.Index(["mapk1;s123", "MAPK1;T185"], name="site_id"),
        ),
    )

    result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=references,
            scoring_config=KinaseScoringConfig(min_substrates=2),
            prediction_config=KinasePredictionConfig(
                top_k=2,
                deterministic_max_selected_kinases=1,
                adaptive_ensemble_runs=1,
            ),
            activity_config=None,
        )
    )

    assert set(result.references.kinase_substrate_map.loc[:, "kinase"]) == {"AKT1"}
    assert set(result.references.kinase_substrate_map.loc[:, "substrate_site"]) == {
        "MAPK1;S123;",
        "MAPK1;T185;",
    }
    assert "MAPK1;S123;" in result.scoring_result.profile_scores.index
    assert "MAPK1;S123;" in result.prediction_result.pred_mat.index
    assert list(result.prediction_result.pred_mat.columns) == ["AKT1"]
    substrate_list = result.prediction_result.substrate_list
    assert substrate_list is not None
    assert set(substrate_list.loc[:, "kinase"]) == {"AKT1"}


def test_prediction_changes_when_downstream_matrix_switches_profile_vs_combined(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    dataset = build_rat_l6_dataset(n_sites=260)
    request = KinaseWorkflowRequest(
        dataset=dataset,
        references=ReferencePreset.AUTO,
        scoring_config=KinaseScoringConfig(min_substrates=2),
        prediction_config=KinasePredictionConfig(
            top_k=6,
            deterministic_max_selected_kinases=12,
            adaptive_ensemble_runs=12,
        ),
        activity_config=None,
    )

    combined_lane = KinaseWorkflow().run(request)

    def _force_profile_lane(*, profile_scores, rank_weighted_fusion_scores):
        _ = rank_weighted_fusion_scores
        return profile_scores, "profile_scores"

    monkeypatch.setattr(
        kinase_executor,
        "select_downstream_score_matrix",
        _force_profile_lane,
    )
    profile_lane = KinaseWorkflow().run(request)

    combined_pred = combined_lane.prediction_result.pred_mat
    profile_pred = profile_lane.prediction_result.pred_mat
    assert not combined_pred.equals(profile_pred)

    rank_weighted_fusion_scores = (
        combined_lane.scoring_result.rank_weighted_fusion_scores
    )
    assert rank_weighted_fusion_scores is not None
    profile_scores = combined_lane.scoring_result.profile_scores
    shared_kinases = pd.Index(combined_pred.columns).intersection(profile_pred.columns)
    assert not shared_kinases.empty

    matched_a_differing_kinase = False
    for kinase in shared_kinases.astype(str):
        if combined_pred.loc[:, kinase].dropna().empty:
            continue
        if profile_pred.loc[:, kinase].dropna().empty:
            continue
        combined_top = rank_weighted_fusion_scores.loc[:, kinase].astype(float).idxmax()
        profile_top = profile_scores.loc[:, kinase].astype(float).idxmax()
        if combined_top == profile_top:
            continue
        assert (
            combined_pred.loc[:, kinase].astype(float).dropna().idxmax() == combined_top
        )
        assert (
            profile_pred.loc[:, kinase].astype(float).dropna().idxmax() == profile_top
        )
        matched_a_differing_kinase = True
        break

    assert matched_a_differing_kinase


def test_diagnostic_scoring_tables_are_opt_in_without_changing_supported_lane() -> None:
    dataset = build_rat_l6_dataset(n_sites=260)
    request_default = KinaseWorkflowRequest(
        dataset=dataset,
        references=ReferencePreset.AUTO,
        scoring_config=KinaseScoringConfig(min_substrates=2),
        prediction_config=KinasePredictionConfig(
            top_k=6,
            deterministic_max_selected_kinases=12,
            adaptive_ensemble_runs=12,
        ),
        activity_config=None,
    )
    request_with_diagnostics = KinaseWorkflowRequest(
        dataset=dataset,
        references=ReferencePreset.AUTO,
        scoring_config=KinaseScoringConfig(
            min_substrates=2,
            include_diagnostic_scoring_tables=True,
        ),
        prediction_config=KinasePredictionConfig(
            top_k=6,
            deterministic_max_selected_kinases=12,
            adaptive_ensemble_runs=12,
        ),
        activity_config=None,
    )

    default_result = KinaseWorkflow().run(request_default)
    diagnostics_result = KinaseWorkflow().run(request_with_diagnostics)

    assert default_result.scoring_result.motif_scores is None
    assert default_result.scoring_result.score_fusion_weights is None
    assert default_result.scoring_result.score_source_summary is not None
    assert default_result.scoring_result.score_source_matrix is None
    assert diagnostics_result.scoring_result.motif_scores is not None
    assert diagnostics_result.scoring_result.score_fusion_weights is not None
    assert diagnostics_result.scoring_result.score_source_summary is not None
    assert diagnostics_result.scoring_result.score_source_matrix is not None

    pd.testing.assert_frame_equal(
        default_result.scoring_result.profile_scores,
        diagnostics_result.scoring_result.profile_scores,
        check_dtype=False,
    )
    assert default_result.scoring_result.rank_weighted_fusion_scores is not None
    assert diagnostics_result.scoring_result.rank_weighted_fusion_scores is not None
    pd.testing.assert_frame_equal(
        default_result.scoring_result.rank_weighted_fusion_scores,
        diagnostics_result.scoring_result.rank_weighted_fusion_scores,
        check_dtype=False,
    )
    pd.testing.assert_frame_equal(
        default_result.prediction_result.pred_mat,
        diagnostics_result.prediction_result.pred_mat,
        check_dtype=False,
    )
    assert default_result.provenance is not None
    scoring_diagnostics = default_result.provenance.workflow_parameters[
        "scoring_diagnostics"
    ]
    assert isinstance(scoring_diagnostics, Mapping)
    by_kinase = scoring_diagnostics["kinase_score_source_counts_by_kinase"]
    totals = scoring_diagnostics["kinase_score_source_counts_total"]
    assert isinstance(by_kinase, Mapping)
    assert isinstance(totals, Mapping)
    assert int(totals["total_sites_count"]) > 0


def test_profile_missing_value_strategy_flows_from_request_to_science(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    dataset = build_rat_l6_dataset(n_sites=220)
    captured_strategies: list[str] = []
    original_build_kinase_profiles = kinase_executor.build_kinase_profiles

    def _capture_strategy(**kwargs):
        captured_strategies.append(kwargs["profile_missing_value_strategy"])
        return original_build_kinase_profiles(**kwargs)

    monkeypatch.setattr(
        kinase_executor,
        "build_kinase_profiles",
        _capture_strategy,
    )
    result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=ReferencePreset.AUTO,
            scoring_config=KinaseScoringConfig(
                min_substrates=2,
                profile_missing_value_strategy="median_skipna",
            ),
            prediction_config=KinasePredictionConfig(
                top_k=5,
                deterministic_max_selected_kinases=8,
                adaptive_ensemble_runs=8,
            ),
            activity_config=None,
        )
    )

    assert captured_strategies == ["median_skipna"]
    assert not result.scoring_result.profile_scores.empty


def test_motif_library_build_is_limited_to_profile_eligible_kinases(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    dataset = build_rat_l6_dataset(n_sites=40)
    overlap_sites = dataset.phospho.index.astype(str).tolist()[:2]
    offlane_sites = [f"OFFLANE{i};S{i};" for i in range(1, 9)]
    references = ReferenceBundle(
        organism=dataset.organism,
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["K_ELIGIBLE", "K_ELIGIBLE"]
                + ["K_OFFLANE" for _ in offlane_sites],
                "substrate_site": overlap_sites + offlane_sites,
            }
        ),
        site_sequences=pd.DataFrame(
            {"site_sequence": ["A" * 31 for _ in overlap_sites + offlane_sites]},
            index=pd.Index(overlap_sites + offlane_sites, name="site_id"),
        ),
    )
    captured_kinases: list[str] = []
    original_build_motif_library = kinase_executor.build_motif_library

    def _capture_eligible_kinases(*, kinase_substrate_map, site_sequences, flank_size):
        nonlocal captured_kinases
        captured_kinases = sorted(
            kinase_substrate_map.loc[:, "kinase"].astype(str).unique().tolist()
        )
        return original_build_motif_library(
            kinase_substrate_map=kinase_substrate_map,
            site_sequences=site_sequences,
            flank_size=flank_size,
        )

    monkeypatch.setattr(
        kinase_executor,
        "build_motif_library",
        _capture_eligible_kinases,
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
                top_k=2,
                deterministic_max_selected_kinases=1,
                adaptive_ensemble_runs=1,
            ),
            activity_config=None,
        )
    )

    assert captured_kinases == ["K_ELIGIBLE"]
    assert list(result.scoring_result.profile_scores.columns) == ["K_ELIGIBLE"]
    assert result.scoring_result.motif_scores is not None
    assert list(result.scoring_result.motif_scores.columns) == ["K_ELIGIBLE"]


def test_kinase_public_predmat_provenance_matches_golden_contract() -> None:
    input_phospho = load_public_predmat_input_phospho()
    site_sequences = load_public_predmat_input_site_sequences()
    canonical_components = [
        _canonical_public_site_components(site_id) for site_id in input_phospho.index
    ]
    phospho = input_phospho.copy(deep=True)
    phospho.index = pd.Index(
        [canonical_site_id for canonical_site_id, _, _ in canonical_components],
        name=input_phospho.index.name,
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": [gene_symbol for _, gene_symbol, _ in canonical_components],
            "site": [site for _, _, site in canonical_components],
            "site_sequence": [
                str(site_sequences[str(site_id).strip()])
                for site_id in input_phospho.index.astype(str)
            ],
        },
        index=phospho.index.copy(),
    )
    substrate_map = load_public_predmat_input_substrate_map()
    dataset = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
        )
    )
    references = ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            [
                {
                    "kinase": str(kinase),
                    "substrate_site": _canonical_public_site_components(site_id)[0],
                }
                for kinase, site_ids in substrate_map.items()
                for site_id in site_ids
            ]
        ),
        site_sequences=pd.DataFrame(
            {
                "site_sequence": [
                    str(sequence) for _, sequence in site_sequences.items()
                ]
            },
            index=pd.Index(
                [
                    _canonical_public_site_components(site_id)[0]
                    for site_id, _ in site_sequences.items()
                ],
                name="site_id",
            ),
        ),
    )

    result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=references,
            scoring_config=KinaseScoringConfig(min_substrates=2),
            prediction_config=KinasePredictionConfig(
                top_k=4,
                deterministic_max_selected_kinases=3,
                adaptive_ensemble_runs=3,
                mode="adaptive_ensemble",
                adaptive_policy="stable",
                n_iterations=2,
                random_state=17,
            ),
            activity_config=None,
        )
    )
    provenance = result.provenance
    assert provenance is not None
    golden = load_kinase_public_predmat_provenance_golden()

    assert provenance.workflow_name == golden["workflow_name"]
    assert provenance.random_state == int(golden["random_state"])
    assert provenance.random_seed_policy == golden["random_seed_policy"]
    assert (
        provenance.workflow_parameters["scoring_config"]
        == (golden["workflow_parameters"]["scoring_config"])
    )
    assert (
        provenance.workflow_parameters["prediction_config"]
        == (golden["workflow_parameters"]["prediction_config"])
    )
    assert provenance.workflow_parameters["prediction_config"]["random_state"] == 17
    assert provenance.reference is not None
    assert provenance.reference.source_type == golden["reference"]["source_type"]
    assert provenance.reference.organism == golden["reference"]["organism"]
    assert provenance.reference.bundle_id is None

    observed_input_tables = _fingerprints_by_name(provenance.input_tables)
    _assert_expected_fingerprint_map(
        observed=observed_input_tables,
        expected=golden["input_tables"],
        expected_overrides={
            **_hash_overrides_from_observed(observed_input_tables),
        },
        compare_hash_values=False,
    )
    _assert_expected_fingerprint_map(
        observed=_fingerprints_by_name(provenance.output_tables),
        expected=golden["output_tables"],
        compare_hash_values=False,
    )
    _assert_expected_fingerprint_map(
        observed=_fingerprints_by_name(provenance.reference.table_fingerprints),
        expected=golden["reference"]["table_fingerprints"],
        compare_hash_values=False,
    )
    assert [
        {"id": item.id.value, "name": item.name, "version": item.version}
        for item in provenance.scientific_policies
    ] == golden["scientific_policies"]
