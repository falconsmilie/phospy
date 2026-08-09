from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pandas.testing as pdt
import pytest

from phospy.contracts.configs.preprocessing import (
    CorrectionMaskPolicy,
    CorrectionMissingnessPolicy,
    InternalBatchCorrectionControlSiteMode,
    InternalBatchCorrectionControlSiteSource,
    InternalBatchCorrectionImputationPolicy,
    InternalBatchCorrectionMethod,
    InternalBatchCorrectionMissingValuePolicy,
    InternalBatchCorrectionRequest,
    InternalBatchCorrectionStageOrder,
    ObservationMask,
    OriginallyMissingCellTracking,
    TemporaryImputationMethod,
    TemporaryImputationPolicy,
)
from phospy.provenance.hashing import fingerprint_table_normalized_axes
from phospy.science.activities.membership import (
    ActivityMembershipSelection,
    fingerprint_ksea_tested_quantitative_matrix,
)
from phospy.science.activities.methods.ksea_zscore import (
    KSEA_STATUS_COMPUTED,
    KSEA_STATUS_INSUFFICIENT_SUBSTRATES,
    KseaZScoreActivityMethod,
)
from phospy.science.activities.models import KinaseActivityInputs, PredMatOverlapSummary
from phospy.science.activities.semantics import ActivityInputMatrix
from phospy.science.batch_correction import DeterministicSpsRuvStyleExecutor
from phospy.science.datasets.preprocessing.control_sites import (
    ControlSiteSet,
    ControlSiteSourceMetadata,
)
from phospy.science.evidence import PeptideEvidenceDatasetResolver, PeptideEvidenceTable
from phospy.science.signalomes.clustering import cluster_sites_with_diagnostics
from phospy.validation.datasets.batch_correction import ResolvedBatchDesignMetadata
from phospy.validation.workflows.batch_correction import ControlSiteEligibilityValidator
from phospy.workflows.batch_correction import BatchCorrectionPlanInterpreter
from tests.support.fixture_byte_policy import assert_text_fixture_matches_sha256
from tests.support.site_keys import site_key_index_from_display_ids

pytestmark = pytest.mark.release_gate

ROOT = Path(__file__).resolve().parents[2]
FIXTURE_ROOT = ROOT / "tests" / "fixtures" / "release_validation_regression"


def _fixture_dir(name: str) -> Path:
    return FIXTURE_ROOT / name


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _read_matrix(path: Path) -> pd.DataFrame:
    return pd.read_csv(path).set_index("site_id")


def _read_bool_text(value: object) -> bool:
    return str(value).strip().lower() == "true"


def _validate_manifest_file_hashes(fixture_dir: Path) -> dict[str, Any]:
    manifest = _read_json(fixture_dir / "MANIFEST.json")
    for file_entry in manifest["files"]:
        path = fixture_dir / str(file_entry["relative_path"])
        assert path.is_file()
        assert hashlib.sha256(path.read_bytes()).hexdigest() == file_entry["sha256"]
    return manifest


@pytest.mark.parametrize(
    ("fixture_family", "classification"),
    (
        ("sps_ruv_planted_unwanted_factor", "synthetic_validation"),
        ("peptide_site_bias_regimes", "synthetic_validation"),
        ("kinase_activity_known_membership", "synthetic_validation"),
        ("signalome_planted_modules", "synthetic_validation"),
        ("importer_edge_cases", "regression"),
    ),
)
def test_scientific_validation_fixture_manifests_classify_evidence(
    fixture_family: str,
    classification: str,
) -> None:
    manifest = _validate_manifest_file_hashes(_fixture_dir(fixture_family))

    assert manifest["fixture_family"] == fixture_family
    assert manifest["classification"] == classification
    assert manifest.get("evidence_category", classification) == classification
    assert "not external parity" in str(manifest["source_policy"])
    assert "empirical validation" not in str(manifest["notes"]).casefold()


def test_sps_ruv_synthetic_fixture_recovers_planted_factor_and_protected_signal() -> (
    None
):
    fixture_dir = _fixture_dir("sps_ruv_planted_unwanted_factor")
    truth = _read_json(fixture_dir / "known_truth.json")
    phospho = _read_matrix(fixture_dir / "phospho.csv")
    sample_metadata = pd.read_csv(fixture_dir / "sample_metadata.csv")

    result = DeterministicSpsRuvStyleExecutor().run(
        phospho=phospho,
        plan=_sps_ruv_plan(
            phospho=phospho,
            sample_metadata=sample_metadata,
            control_site_ids=tuple(str(site) for site in truth["control_site_ids"]),
        ),
    )

    planted_factor = np.asarray(
        [
            truth["planted_unwanted_factor_by_sample"][sample_id]
            for sample_id in result.estimated_unwanted_factors.index.astype(str)
        ],
        dtype=float,
    )
    estimated = result.estimated_unwanted_factors.iloc[:, 0].to_numpy(dtype=float)
    factor_correlation = float(np.corrcoef(planted_factor, estimated)[0, 1])
    recovery = truth["unwanted_factor_recovery"]
    assert abs(factor_correlation) >= recovery["minimum_abs_correlation"]

    corrected_controls = result.corrected_matrix.loc[list(truth["control_site_ids"])]
    control_span = corrected_controls.max(axis=1) - corrected_controls.min(axis=1)
    assert (
        float(control_span.max()) <= recovery["control_row_max_span_after_correction"]
    )

    contrast = truth["protected_condition_effect"]
    before = _condition_effect(phospho, sample_metadata)
    after = _condition_effect(result.corrected_matrix, sample_metadata)
    signal_site = str(truth["protected_signal_site_id"])
    assert before.loc[signal_site] == pytest.approx(
        contrast["expected_difference"],
        abs=contrast["acceptance_absolute_tolerance"],
    )
    assert after.loc[signal_site] == pytest.approx(
        contrast["expected_difference"],
        abs=contrast["acceptance_absolute_tolerance"],
    )
    assert (
        "not used to estimate unwanted factors"
        in (result.diagnostics.term_roles["batch_terms"])
    )


def _read_peptide_bias_evidence(fixture_dir: Path) -> pd.DataFrame:
    evidence = pd.read_csv(fixture_dir / "peptide_evidence.csv")
    evidence.loc[:, "multi_site"] = evidence.loc[:, "multi_site"].map(_read_bool_text)
    return evidence


def _read_peptide_bias_mapping(fixture_dir: Path) -> pd.DataFrame:
    mapping = pd.read_csv(fixture_dir / "site_mapping.csv")
    mapping.loc[:, "mapping_uncertainty"] = mapping.loc[
        :,
        "mapping_uncertainty",
    ].map(_read_bool_text)
    mapping.loc[:, "is_multi_site"] = mapping.loc[:, "is_multi_site"].map(
        _read_bool_text
    )
    return mapping


def _run_peptide_bias_resolution(
    *,
    evidence: pd.DataFrame,
    mapping: pd.DataFrame,
):
    return PeptideEvidenceDatasetResolver().run(
        evidence=PeptideEvidenceTable(
            frame=evidence,
            sample_intensity_columns=("sample_1", "sample_2"),
            site_mapping=mapping,
        ),
        multi_site_policy="split",
        input_intensity_scale="linear",
    )


def _expected_site_estimate_matrix(fixture_dir: Path) -> pd.DataFrame:
    expected = pd.read_csv(fixture_dir / "expected_site_estimates.csv")
    matrix = expected.pivot(
        index="site_id",
        columns="sample_id",
        values="expected_site_abundance",
    ).sort_index()
    matrix.columns.name = None
    matrix.index.name = "site_id"
    return matrix.loc[:, ["sample_1", "sample_2"]].astype(float)


def test_peptide_site_bias_expected_output_generator_is_independent() -> None:
    fixture_dir = _fixture_dir("peptide_site_bias_regimes")
    manifest = _validate_manifest_file_hashes(fixture_dir)
    generator_path = ROOT / str(manifest["expected_output_generator"])
    source = generator_path.read_text(encoding="utf-8")
    parsed = ast.parse(source)

    imported_modules: set[str] = set()
    for node in ast.walk(parsed):
        if isinstance(node, ast.Import):
            imported_modules.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module is not None:
            imported_modules.add(node.module)

    forbidden_import_roots = (
        "phospy",
        "src.phospy",
        "tests",
    )
    assert not any(
        imported_module == root or imported_module.startswith(f"{root}.")
        for imported_module in imported_modules
        for root in forbidden_import_roots
    )
    assert "PeptideEvidenceDatasetResolver" not in source
    assert "allocate_peptide_signals_to_resolved_sites" not in source
    assert "tests.support" not in source
    assert (
        hashlib.sha256(generator_path.read_bytes()).hexdigest()
        == (manifest["expected_output_generator_sha256"])
    )


def test_peptide_site_bias_raw_fixture_runs_through_production_resolver() -> None:
    fixture_dir = _fixture_dir("peptide_site_bias_regimes")
    assumptions = _read_json(fixture_dir / "policy_assumptions.json")
    expected_matrix = _expected_site_estimate_matrix(fixture_dir)

    resolved = _run_peptide_bias_resolution(
        evidence=_read_peptide_bias_evidence(fixture_dir),
        mapping=_read_peptide_bias_mapping(fixture_dir),
    )
    observed = resolved.phospho.sort_index().loc[:, ["sample_1", "sample_2"]]

    pdt.assert_frame_equal(
        observed,
        expected_matrix,
        check_exact=False,
        rtol=float(assumptions.get("absolute_tolerance", 1.0e-12)),
        atol=float(assumptions.get("absolute_tolerance", 1.0e-12)),
    )

    payload = resolved.summary.to_payload()
    assert payload["input_intensity_scale"] == assumptions["input_intensity_scale"]
    assert (
        payload["input_quantitative_meaning"]
        == (assumptions["input_quantitative_meaning"])
    )
    assert payload["allocation_domain"] == assumptions["allocation_domain"]
    assert payload["multi_site_policy"] == assumptions["multi_site_policy"]
    assert int(payload["fractional_mapping_rows"]) == 4
    assert int(payload["duplicate_peptide_rows"]) == 4


def test_peptide_site_bias_raw_fixture_is_row_order_invariant() -> None:
    fixture_dir = _fixture_dir("peptide_site_bias_regimes")
    evidence = _read_peptide_bias_evidence(fixture_dir)
    mapping = _read_peptide_bias_mapping(fixture_dir)
    expected = _run_peptide_bias_resolution(
        evidence=evidence,
        mapping=mapping,
    ).phospho.sort_index()

    permuted_evidence = evidence.sort_values(
        ["site_id", "peptide_row_id"],
        ascending=[False, False],
    ).reset_index(drop=True)
    permuted_mapping = mapping.iloc[::-1, :].reset_index(drop=True)
    observed = _run_peptide_bias_resolution(
        evidence=permuted_evidence,
        mapping=permuted_mapping,
    ).phospho.sort_index()

    pdt.assert_frame_equal(observed, expected)


def test_peptide_site_bias_fixture_quantifies_adverse_regimes() -> None:
    fixture_dir = _fixture_dir("peptide_site_bias_regimes")
    rows = pd.read_csv(fixture_dir / "bias_regimes.csv")
    expected = _read_json(fixture_dir / "expected_bias_summary.json")["regimes"]

    observed = rows.groupby("regime")["absolute_bias"].mean().to_dict()
    for regime, expected_payload in expected.items():
        assert observed[regime] == pytest.approx(
            float(expected_payload["mean_absolute_bias"])
        )
        assert (
            int(
                rows.loc[
                    (rows["regime"] == regime)
                    & (rows["bias_status"] == "nonestimable_missing_output"),
                    :,
                ].shape[0]
            )
            == expected_payload["missing_estimate_count"]
        )

    assert set(rows["evidence_classification"]) == {"synthetic_validation"}
    concordant = rows.loc[rows["regime"] == "duplicate_concordant_control"]
    assert set(concordant["signed_bias"].astype(float)) == {0.0}
    discordant = rows.loc[rows["regime"] == "duplicate_discordant"]
    assert set(discordant["signed_bias"].astype(float)) == {15.0}
    ambiguous = rows.loc[rows["regime"] == "ambiguous_equal_split"]
    assert set(ambiguous["signed_bias"].astype(float)) == {-50.0, 50.0}
    unequal = rows.loc[rows["regime"] == "ambiguous_unequal_fraction"]
    assert set(unequal["signed_bias"].astype(float)) == {-60.0, 60.0}
    missing_partial = rows.loc[
        rows["case_id"] == "missing_value_in_one_sample",
    ].set_index("sample_id")
    assert float(missing_partial.loc["sample_1", "signed_bias"]) == pytest.approx(0.0)
    assert float(missing_partial.loc["sample_2", "signed_bias"]) == pytest.approx(-20.0)
    missing_entire = rows.loc[
        rows["case_id"] == "entirely_missing_peptide_observation",
    ]
    assert missing_entire["expected_site_abundance"].isna().all()
    assert set(missing_entire["bias_status"]) == {"nonestimable_missing_output"}
    localisation = rows.loc[rows["regime"] == "localisation_error"]
    assert set(localisation["signed_bias"].astype(float)) == {-100.0, 100.0}
    assert set(rows["regime"]) == {
        "duplicate_concordant_control",
        "duplicate_discordant",
        "ambiguous_equal_split",
        "ambiguous_unequal_fraction",
        "missing_observation",
        "localisation_error",
    }


def test_peptide_site_bias_fixture_bias_matches_production_output_and_truth() -> None:
    fixture_dir = _fixture_dir("peptide_site_bias_regimes")
    resolved = _run_peptide_bias_resolution(
        evidence=_read_peptide_bias_evidence(fixture_dir),
        mapping=_read_peptide_bias_mapping(fixture_dir),
    )
    production_long = resolved.phospho.reset_index().melt(
        id_vars="site_id",
        var_name="sample_id",
        value_name="production_site_abundance",
    )
    bias = pd.read_csv(fixture_dir / "bias_regimes.csv")
    merged = bias.merge(production_long, on=("site_id", "sample_id"), how="left")

    finite = merged.loc[merged["bias_status"] == "finite", :]
    assert not finite.empty
    np.testing.assert_allclose(
        finite["production_site_abundance"].to_numpy(dtype=float),
        finite["expected_site_abundance"].to_numpy(dtype=float),
        rtol=1.0e-12,
        atol=1.0e-12,
    )
    np.testing.assert_allclose(
        finite["production_site_abundance"].to_numpy(dtype=float)
        - finite["true_site_abundance"].to_numpy(dtype=float),
        finite["signed_bias"].to_numpy(dtype=float),
        rtol=1.0e-12,
        atol=1.0e-12,
    )

    missing = merged.loc[merged["bias_status"] == "nonestimable_missing_output", :]
    assert set(missing["case_id"]) == {"entirely_missing_peptide_observation"}
    assert missing["production_site_abundance"].isna().all()


def test_kinase_activity_known_membership_fixture_direction_and_coverage() -> None:
    fixture_dir = _fixture_dir("kinase_activity_known_membership")
    truth = _read_json(fixture_dir / "known_truth.json")
    effects = _read_matrix(fixture_dir / "phospho_effects.csv")
    membership = _read_matrix(fixture_dir / "membership_scores.csv")
    threshold = float(truth["membership_threshold"])

    min_two = _ksea_result(
        pred_mat=membership,
        effects=effects,
        threshold=threshold,
        min_substrates=2,
    )
    min_three = _ksea_result(
        pred_mat=membership,
        effects=effects,
        threshold=threshold,
        min_substrates=3,
    )

    assert float(min_two.activity_matrix.at["K_UP", "stim_effect"]) > 0.0
    assert float(min_two.activity_matrix.at["K_DOWN", "stim_effect"]) < 0.0

    coverage = truth["substrate_coverage_sensitivity"]
    min_two_stats = _activity_statistics_for_profile(min_two.statistics_table)
    min_three_stats = _activity_statistics_for_profile(min_three.statistics_table)
    kinase = str(coverage["kinase"])
    assert min_two_stats.at[kinase, "computability_status"] == KSEA_STATUS_COMPUTED
    assert min_three_stats.at[kinase, "computability_status"] == (
        KSEA_STATUS_INSUFFICIENT_SUBSTRATES
    )
    sparse = str(truth["sparse_membership"]["kinase"])
    assert min_two_stats.at[sparse, "computability_status"] == (
        KSEA_STATUS_INSUFFICIENT_SUBSTRATES
    )


def test_signalome_planted_module_fixture_recovers_modules_and_perturbation_stability() -> (
    None
):
    fixture_dir = _fixture_dir("signalome_planted_modules")
    truth = _read_json(fixture_dir / "known_truth.json")
    baseline = _read_matrix(fixture_dir / "score_matrix.csv")
    perturbed = _read_matrix(fixture_dir / "score_matrix_perturbed.csv")
    planted = pd.read_csv(fixture_dir / "planted_modules.csv").set_index("site_id")
    planted_labels = planted.loc[baseline.index, "planted_module"].astype(str)

    baseline_result = cluster_sites_with_diagnostics(
        scoring_matrix=baseline,
        requested_module_count=None,
        max_clusters=4,
    )
    perturbed_result = cluster_sites_with_diagnostics(
        scoring_matrix=perturbed,
        requested_module_count=None,
        max_clusters=4,
    )

    assert (
        baseline_result.module_selection_diagnostics.selected_module_count
        == (truth["expected_module_count"])
    )
    assert baseline_result.module_selection_diagnostics.stability_report.status == (
        "stable"
    )
    assert _pairwise_coassignment_accuracy(
        planted_labels,
        baseline_result.site_clusters,
    ) == pytest.approx(truth["required_pairwise_accuracy"])
    assert _pairwise_coassignment_accuracy(
        baseline_result.site_clusters,
        perturbed_result.site_clusters,
    ) == pytest.approx(truth["required_perturbation_stability"])


def test_importer_edge_case_manifest_references_existing_fixture_bytes() -> None:
    fixture_dir = _fixture_dir("importer_edge_cases")
    index = _read_json(fixture_dir / "fixture_index.json")

    assert index["classification"] == "regression"
    assert "not broad vendor parity" in index["supported_interpretation"]
    assert "Spectronaut/DIA-NN" in index["supported_interpretation"]
    for file_entry in index["referenced_fixture_files"]:
        path = ROOT / str(file_entry["relative_path"])
        assert path.is_file()
        assert_text_fixture_matches_sha256(
            path,
            expected_sha256=str(file_entry["sha256"]),
            repo_root=ROOT,
        )

    coverage = index["edge_case_coverage"]
    assert "raw/LFQ intensity ambiguity" in coverage["maxquant"]
    assert "ambiguous localisation diagnostics" in coverage["fragpipe_ptmprophet"]


def _sps_ruv_plan(
    *,
    phospho: pd.DataFrame,
    sample_metadata: pd.DataFrame,
    control_site_ids: tuple[str, ...],
):
    sample_metadata = sample_metadata.set_index("sample_id", drop=False)
    sample_ids = tuple(phospho.columns.astype(str))
    metadata = ResolvedBatchDesignMetadata(
        batch_by_sample={
            sample_id: str(sample_metadata.at[sample_id, "batch"])
            for sample_id in sample_ids
        },
        condition_by_sample={
            sample_id: str(sample_metadata.at[sample_id, "condition"])
            for sample_id in sample_ids
        },
        replicate_by_sample={
            sample_id: str(sample_metadata.at[sample_id, "replicate"])
            for sample_id in sample_ids
        },
        sample_order=sample_ids,
    )
    control_mapping = ControlSiteEligibilityValidator().run(
        control_set=ControlSiteSet.from_site_keys(
            control_site_ids,
            source_metadata=ControlSiteSourceMetadata(
                organism="rat",
                identifier_namespace="site_key",
                source_name="synthetic-known-truth",
                source_version="v1",
                license="project synthetic data",
                redistribution="project test fixture",
            ),
        ),
        site_keys=tuple(phospho.index.astype(str)),
        method="sps_ruv_style",
        min_eligible_controls=2,
        n_unwanted_factors=1,
    )
    missingness_policy = CorrectionMissingnessPolicy(
        temporary_imputation=TemporaryImputationPolicy(
            allowed=True,
            method=TemporaryImputationMethod.ROW_MEDIAN_TEMPORARY,
            method_parameters={"min_observed_values": 2},  # type: ignore[arg-type]
        ),
        originally_missing_cells_tracked_by=(
            OriginallyMissingCellTracking.OBSERVATION_MASK
        ),
        correction_mask_policy=CorrectionMaskPolicy(),
        observation_mask=ObservationMask(
            feature_ids=tuple(phospho.index.astype(str)),
            sample_ids=sample_ids,
            originally_missing_cells=(),
        ),
    )
    return BatchCorrectionPlanInterpreter().run(
        config=InternalBatchCorrectionRequest(
            method=InternalBatchCorrectionMethod.SPS_RUV_STYLE,
            batch_column="batch",
            condition_columns=("condition",),
            replicate_column="replicate",
            control_site_source=(
                InternalBatchCorrectionControlSiteSource.CALLER_SUPPLIED
            ),
            control_site_mode=InternalBatchCorrectionControlSiteMode.SITE_KEY_LIST,
            missing_value_policy=(
                InternalBatchCorrectionMissingValuePolicy.ALLOW_TEMPORARY_IMPUTATION
            ),
            imputation_policy=(
                InternalBatchCorrectionImputationPolicy.ROW_MEDIAN_TEMPORARY
            ),
            n_unwanted_factors=1,
            stage_order=(
                InternalBatchCorrectionStageOrder.AFTER_MISSING_DATA_BEFORE_DOWNSTREAM
            ),
            diagnostics_enabled=True,
        ),
        dataset_metadata=metadata,
        control_site_mapping=control_mapping,
        missingness_policy=missingness_policy,
    )


def _condition_effect(
    matrix: pd.DataFrame,
    sample_metadata: pd.DataFrame,
) -> pd.Series:
    sample_metadata = sample_metadata.set_index("sample_id", drop=False)
    condition_a_samples = sample_metadata.loc[
        sample_metadata["condition"].astype(str) == "A",
        "sample_id",
    ].astype(str)
    condition_b_samples = sample_metadata.loc[
        sample_metadata["condition"].astype(str) == "B",
        "sample_id",
    ].astype(str)
    return matrix.loc[:, condition_b_samples.tolist()].mean(axis=1) - matrix.loc[
        :, condition_a_samples.tolist()
    ].mean(axis=1)


def _ksea_result(
    *,
    pred_mat: pd.DataFrame,
    effects: pd.DataFrame,
    threshold: float,
    min_substrates: int,
):
    pred_mat = _with_site_key_index(pred_mat)
    effects = _with_site_key_index(effects)
    selected_substrates = pred_mat.loc[
        pred_mat.ge(float(threshold)).any(axis=1)
    ].index.astype(str)
    membership = ActivityMembershipSelection.fixed_external_reference(
        provider_method_identifier="release_validation_known_membership",
        provider_method_version="synthetic-v1",
        provider_score_source_identifier="known_membership_scores",
        threshold_top_k_policy={
            "evidence_threshold": float(threshold),
            "evidence_threshold_operator": ">=",
            "top_k": None,
        },
        source_reference_fingerprints=(
            fingerprint_table_normalized_axes(
                pred_mat,
                name="release_validation.kinase_activity_known_membership",
            ),
        ),
        tested_quantitative_matrix_fingerprint=(
            fingerprint_ksea_tested_quantitative_matrix(effects)
        ),
        selected_kinase_universe=pred_mat.columns.astype(str).tolist(),
        selected_substrate_universe=selected_substrates.tolist(),
    )
    return KseaZScoreActivityMethod(
        evidence_threshold=threshold,
        min_substrates=min_substrates,
    ).run(
        KinaseActivityInputs(
            pred_mat=pred_mat,
            phospho_matrix=effects,
            threshold=threshold,
            min_substrates=min_substrates,
            top_n_substrates=1,
            overlap_summary=PredMatOverlapSummary(
                overlap_count=int(pred_mat.index.intersection(effects.index).size),
                pred_mat_rows=int(pred_mat.index.size),
                phospho_rows=int(effects.index.size),
            ),
            activity_input=ActivityInputMatrix.standardised_effect(
                effects,
                _assume_owned=True,
            ),
            membership_selection=membership,
        )
    )


def _with_site_key_index(frame: pd.DataFrame) -> pd.DataFrame:
    labels = frame.index.astype(str).tolist()
    if all(label.startswith("phospy:v1|") for label in labels):
        return frame
    keyed = frame.copy(deep=True)
    keyed.index = site_key_index_from_display_ids(
        labels,
        protein_namespace="gene_symbol",
    )
    return keyed


def _activity_statistics_for_profile(statistics_table: pd.DataFrame) -> pd.DataFrame:
    return statistics_table.loc[
        statistics_table["profile_id"].astype(str) == "stim_effect",
        :,
    ].set_index("kinase")


def _pairwise_coassignment_accuracy(
    expected: pd.Series,
    observed: pd.Series,
) -> float:
    expected = expected.reindex(observed.index).astype(str)
    observed = observed.astype(str)
    comparisons = 0
    matches = 0
    labels = observed.index.astype(str).tolist()
    for left_position, left_label in enumerate(labels):
        for right_label in labels[left_position + 1 :]:
            expected_same = expected.at[left_label] == expected.at[right_label]
            observed_same = observed.at[left_label] == observed.at[right_label]
            comparisons += 1
            matches += int(expected_same == observed_same)
    if comparisons == 0:
        raise AssertionError("pairwise coassignment requires at least two sites")
    return float(matches) / float(comparisons)
