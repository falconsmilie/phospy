from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pandas.testing as pdt
import pytest

from phospy.advanced import (
    SpsDiscoveryConfig,
    SpsDiscoveryRequest,
    SpsDiscoveryWorkflow,
    SpsReferenceDataset,
)
from phospy.science.batch_correction import (
    RuvIIIReplicateStructure,
    run_ruv_iii,
)

ROOT = Path(__file__).resolve().parents[2]
FIXTURE_ROOT = ROOT / "tests" / "fixtures" / "rewrite_parity" / "phosr_ruv"

pytestmark = [pytest.mark.parity, pytest.mark.release_gate]


def _read_json(relative_path: str) -> dict[str, Any]:
    return json.loads((FIXTURE_ROOT / relative_path).read_text(encoding="utf-8"))


def _read_matrix(relative_path: str) -> pd.DataFrame:
    frame = pd.read_csv(FIXTURE_ROOT / relative_path)
    return frame.set_index("site_key").rename_axis(index="site_key")


def _sps_result():
    metadata = _read_json("METADATA.json")
    sample_metadata = pd.read_csv(FIXTURE_ROOT / "sps" / "sample_metadata.csv")
    references: list[SpsReferenceDataset] = []
    for dataset_id in ("reference_a", "reference_b", "reference_c"):
        matrix = _read_matrix(f"sps/{dataset_id}_condition_relative_log2.csv")
        rows = sample_metadata.loc[sample_metadata["dataset_id"] == dataset_id]
        condition_by_sample = dict(
            zip(rows["sample_id"], rows["condition"], strict=True)
        )
        references.append(
            SpsReferenceDataset.from_condition_relative_log2(
                dataset_id=dataset_id,
                intensities=matrix,
                condition_by_sample=condition_by_sample,
                log2_scale_established_by=(
                    "pinned PhosR parity fixture source matrices are log2 abundance"
                ),
                baseline_centering_established_by=(
                    "fixture generator subtracts each site's baseline-condition mean"
                ),
                organism="rat",
                baseline_context="fixture baseline condition",
                reference_context="synthetic rat PhosR parity fixture",
                source_name="synthetic PhosR parity fixture",
                source_version="phosr-ruv-parity-v1",
                source_uri=("https://example.test/phosr-parity/" + str(dataset_id)),
            )
        )
    case = metadata["sps_case"]
    return SpsDiscoveryWorkflow().run(
        SpsDiscoveryRequest(
            reference_datasets=tuple(references),
            config=SpsDiscoveryConfig(
                top_n=int(case["requested_top_n"]),
                minimum_reference_datasets=int(case["reference_dataset_count"]),
                minimum_datasets_per_site=int(case["minimum_datasets_per_site"]),
                minimum_shared_sites=200,
            ),
        )
    )


def test_sps_preparation_is_explicit_and_both_implementations_receive_same_state() -> (
    None
):
    metadata = _read_json("METADATA.json")
    quantitative = metadata["quantitative_contract"]
    assert quantitative == {
        "sps_source_quantitative_meaning": (
            "absolute processed phosphosite log2 abundance"
        ),
        "sps_source_scale": "log2",
        "selected_baseline_condition": "baseline",
        "centering_operation": (
            "for each site and reference dataset, subtract the arithmetic mean of "
            "baseline_1 and baseline_2 from every sample on the log2 scale"
        ),
        "centered_quantitative_meaning": "condition-relative log2 fold change",
        "zero_semantics": "the per-site baseline-condition mean is exactly zero",
        "matrix_passed_to_getSPS": "sps/*_condition_relative_log2.csv",
        "matrix_passed_to_phospy": (
            "the identical sps/*_condition_relative_log2.csv bytes"
        ),
        "raw_absolute_input_is_not_valid_phospy_sps_input": True,
    }

    for dataset_id in ("reference_a", "reference_b", "reference_c"):
        absolute = _read_matrix(f"sps/{dataset_id}_source_absolute_log2.csv")
        relative = _read_matrix(f"sps/{dataset_id}_condition_relative_log2.csv")
        baseline = absolute.loc[:, ["baseline_1", "baseline_2"]].mean(axis=1)
        expected = absolute.subtract(baseline, axis=0)
        pdt.assert_frame_equal(relative, expected, atol=1e-14, rtol=1e-14)
        np.testing.assert_allclose(
            relative.loc[:, ["baseline_1", "baseline_2"]].mean(axis=1),
            0.0,
            atol=1e-14,
            rtol=0.0,
        )


def test_sps_ranking_intermediates_and_top_n_match_pinned_phosr_getsps() -> None:
    metadata = _read_json("METADATA.json")
    tolerance = float(
        metadata["comparison_policy"]["sps_stability_and_consensus_absolute_tolerance"]
    )
    trace = pd.read_csv(FIXTURE_ROOT / "sps" / "phosr_getsps_trace.csv").set_index(
        "site_key"
    )
    expected_selection = pd.read_csv(FIXTURE_ROOT / "sps" / "phosr_getsps_selected.csv")
    result = _sps_result()

    assert set(record.site_key for record in result.site_ranking) == set(trace.index)
    for record in result.site_ranking:
        expected = trace.loc[record.site_key]
        assert record.consensus_stability_score == pytest.approx(
            expected["consensus_score"], abs=tolerance, rel=0.0
        )
        assert record.contributing_dataset_count == 3
        for statistic in record.dataset_statistics:
            assert statistic.stability_score == pytest.approx(
                expected[f"{statistic.dataset_id}_stability"],
                abs=tolerance,
                rel=0.0,
            )

    # This fixture's selected boundary has no cross-boundary tie, so the exact
    # top-N sequence is defined. Tie-group evidence is checked separately below.
    assert result.selected_site_keys == tuple(expected_selection["site_key"])
    assert len(result.selected_site_keys) == int(
        metadata["sps_case"]["requested_top_n"]
    )
    selected_numbers = (
        expected_selection["phosr_site_id"]
        .str.extract(r"SPS(\d{4})", expand=False)
        .astype(int)
    )
    assert int(selected_numbers.max()) <= 25
    strongly_changing = trace.loc[
        trace["phosr_site_id"].str.extract(r"SPS(\d{4})", expand=False).astype(int)
        >= 181
    ]
    assert (
        strongly_changing[
            [
                "reference_a_stability",
                "reference_b_stability",
                "reference_c_stability",
            ]
        ]
        .min()
        .min()
        > 2.0
    )


def test_sps_tie_group_matches_phosr_evidence_and_phospy_order_is_deterministic() -> (
    None
):
    result = _sps_result()
    records = {record.site_key: record for record in result.site_ranking}
    tie_keys = tuple(
        sorted(
            key
            for key in records
            if "protein_identifier=SPS0020|" in key
            or "protein_identifier=SPS0021|" in key
        )
    )
    assert len(tie_keys) == 2
    assert records[tie_keys[0]].consensus_stability_score == pytest.approx(
        records[tie_keys[1]].consensus_stability_score,
        abs=1e-15,
        rel=0.0,
    )
    positions = tuple(record.site_key for record in result.site_ranking)
    observed_tie_order = tuple(key for key in positions if key in set(tie_keys))
    assert observed_tie_order == tie_keys


def _ruv_inputs() -> tuple[
    pd.DataFrame,
    tuple[str, ...],
    pd.DataFrame,
]:
    matrix = _read_matrix("ruv_iii/input_matrix.csv")
    controls = tuple(
        pd.read_csv(FIXTURE_ROOT / "ruv_iii" / "negative_controls.csv")["site_key"]
    )
    sample_metadata = pd.read_csv(FIXTURE_ROOT / "ruv_iii" / "sample_metadata.csv")
    return matrix, controls, sample_metadata


@pytest.mark.parametrize("k", (0, 1, 2))
def test_ruv_iii_corrected_matrices_match_pinned_ruv(k: int) -> None:
    metadata = _read_json("METADATA.json")
    diagnostics = _read_json("ruv_iii/reference_diagnostics.json")["cases"][f"k{k}"]
    matrix, controls, sample_metadata = _ruv_inputs()
    sample_order = tuple(sample_metadata["sample_id"])
    replicates = RuvIIIReplicateStructure.from_assignments(
        sample_order=sample_order,
        replicate_by_sample=dict(
            zip(
                sample_metadata["sample_id"],
                sample_metadata["replicate_set"],
                strict=True,
            )
        ),
    )

    observed = run_ruv_iii(
        matrix,
        control_site_keys=controls,
        replicate_structure=replicates,
        k=k,
    )
    expected = _read_matrix(f"ruv_iii/ruv_expected_k{k}.csv")
    policy = metadata["comparison_policy"]
    pdt.assert_frame_equal(
        observed.corrected_matrix,
        expected,
        atol=float(policy["ruv_corrected_matrix_absolute_tolerance"]),
        rtol=float(policy["ruv_corrected_matrix_relative_tolerance"]),
    )
    assert observed.diagnostics.input_matrix_rank == diagnostics["input_matrix_rank"]
    assert (
        observed.diagnostics.replicate_mapping_rank
        == diagnostics["replicate_mapping_rank"]
    )
    assert (
        observed.diagnostics.replicate_residual_rank
        == diagnostics["replicate_residual_rank"]
    )
    np.testing.assert_allclose(
        observed.diagnostics.residual_singular_values,
        np.atleast_1d(diagnostics["residual_singular_values"])
        if k > 0
        else np.asarray([], dtype=float),
        atol=float(policy["diagnostics_absolute_tolerance"]),
        rtol=1e-12,
    )

    if k == 2:
        # Replicate-set mean contrasts carry the planted biological signal and
        # must remain after removing the two within-set technical factors.
        replicate_sets = list(sample_metadata["replicate_set"])
        original_means = matrix.T.groupby(replicate_sets).mean().T
        corrected_means = observed.corrected_matrix.T.groupby(replicate_sets).mean().T
        pdt.assert_frame_equal(
            corrected_means.subtract(corrected_means.iloc[:, 0], axis=0),
            original_means.subtract(original_means.iloc[:, 0], axis=0),
            atol=1e-10,
            rtol=1e-12,
        )


def test_ruv_iii_sample_permutation_matches_external_reference() -> None:
    diagnostics = _read_json("ruv_iii/reference_diagnostics.json")
    matrix, controls, sample_metadata = _ruv_inputs()
    sample_order = tuple(diagnostics["sample_permutation"])
    assignments = dict(
        zip(
            sample_metadata["sample_id"],
            sample_metadata["replicate_set"],
            strict=True,
        )
    )
    observed = run_ruv_iii(
        matrix.loc[:, list(sample_order)],
        control_site_keys=controls,
        replicate_structure=RuvIIIReplicateStructure.from_assignments(
            sample_order=sample_order,
            replicate_by_sample=assignments,
        ),
        k=2,
    )
    expected = _read_matrix("ruv_iii/ruv_expected_k2_permuted.csv")
    pdt.assert_frame_equal(
        observed.corrected_matrix,
        expected,
        atol=1e-10,
        rtol=1e-12,
    )


def test_fixture_documents_supported_scope_and_contract_differences() -> None:
    metadata = _read_json("METADATA.json")
    differences = metadata["contract_differences"]
    assert (
        "exactly minimum_datasets_per_site"
        in differences["partial_reference_contributions"]
    )
    assert (
        "requires established condition-relative log2"
        in differences["ambiguous_sps_input"]
    )
    assert "No missing-value parity is claimed" in differences["missing_values"]
    assert (
        "rejects singleton replicate groups" in differences["singleton_replicate_sets"]
    )
    assert "does not silently cap or reduce it" in differences["non_estimable_k"]
    assert "supported overlapping domain" in differences["non_estimable_k"]
    assert metadata["redistribution"]["status"] == (
        "approved_for_repository_test_fixture_redistribution"
    )
