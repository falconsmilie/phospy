from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import cast

import numpy as np
import pandas as pd
import pytest

import phospy.science.batch_correction.sps_discovery as sps_science
from phospy.advanced import (
    SpsDiscoveryConfig,
    SpsDiscoveryRequest,
    SpsDiscoveryValidationError,
    SpsDiscoveryWorkflow,
    SpsReferenceDataset,
)
from phospy.contracts.configs.preprocessing import (
    InternalBatchCorrectionControlSiteMode,
    InternalBatchCorrectionControlSiteSource,
    InternalBatchCorrectionImputationPolicy,
    InternalBatchCorrectionMethod,
    InternalBatchCorrectionMissingValuePolicy,
    InternalBatchCorrectionRequest,
    InternalBatchCorrectionStageOrder,
)
from phospy.validation.workflows.batch_correction import (
    BatchCorrectionWorkflowControlSiteValidator,
)
from phospy.workflows.batch_correction import (
    BatchCorrectionWorkflowRequest,
)
from tests.support.site_keys import protein_site_key_index


def _site_keys(count: int) -> pd.Index:
    return protein_site_key_index(
        protein_identifiers=[f"P{position:04d}" for position in range(count)],
        sites=[f"S{position + 1}" for position in range(count)],
    )


def _reference(
    dataset_id: str,
    *,
    site_keys: Sequence[str] | pd.Index,
    columns: Mapping[str, Sequence[float]],
    conditions: Mapping[str, str],
) -> SpsReferenceDataset:
    return SpsReferenceDataset.from_condition_relative_log2(
        dataset_id=dataset_id,
        intensities=pd.DataFrame(columns, index=pd.Index(site_keys, name="site_key")),
        condition_by_sample=conditions,
        log2_scale_established_by="test fixture log2 preparation",
        baseline_centering_established_by="test fixture control subtraction",
        organism="rat",
        baseline_context="control condition",
        reference_context="synthetic rat SPS science fixture",
        source_name=f"synthetic-{dataset_id}",
        source_version="1",
        source_uri=f"https://example.test/sps/{dataset_id}",
    )


def _hand_references() -> tuple[SpsReferenceDataset, SpsReferenceDataset]:
    keys = _site_keys(4)
    conditions_a = {
        "a_control_1": "control",
        "a_control_2": "control",
        "a_treated_1": "treated",
        "a_treated_2": "treated",
    }
    conditions_b = {
        "b_base_1": "base",
        "b_base_2": "base",
        "b_stim_1": "stim",
        "b_stim_2": "stim",
    }
    return (
        _reference(
            "study_a",
            site_keys=keys,
            columns={
                "a_control_1": [0.10, 0.30, 0.00, 0.40],
                "a_control_2": [0.10, 0.30, 0.00, 0.40],
                "a_treated_1": [0.20, 0.50, 4.00, 0.40],
                "a_treated_2": [0.20, 0.50, 4.00, 0.40],
            },
            conditions=conditions_a,
        ),
        _reference(
            "study_b",
            site_keys=keys,
            columns={
                "b_base_1": [0.15, 0.60, -0.20, 0.40],
                "b_base_2": [0.15, 0.60, -0.20, 0.40],
                "b_stim_1": [0.10, 1.00, 3.00, 0.40],
                "b_stim_2": [0.10, 1.00, 3.00, 0.40],
            },
            conditions=conditions_b,
        ),
    )


def _run(
    references: Sequence[SpsReferenceDataset],
    *,
    config: SpsDiscoveryConfig | None = None,
):
    return SpsDiscoveryWorkflow().run(
        SpsDiscoveryRequest(
            reference_datasets=references,
            config=config or SpsDiscoveryConfig(top_n=2),
        )
    )


def _scalar_reference_stability_evidence(
    reference: SpsReferenceDataset,
) -> sps_science._ReferenceStabilityEvidence:
    """Retain the reviewed scalar-access semantics as an equivalence oracle."""

    intensities = reference.intensities
    conditions = tuple(
        sorted({str(value) for value in reference.condition_by_sample.values()})
    )
    samples_by_condition = {
        condition: tuple(
            sorted(
                str(column)
                for column in intensities.columns
                if str(reference.condition_by_sample[str(column)]) == condition
            )
        )
        for condition in conditions
    }
    scores: dict[str, float] = {}
    for site_key in sorted(str(value) for value in intensities.index):
        condition_means: list[float] = []
        for condition in conditions:
            values = tuple(
                float(cast(float | int, intensities.at[site_key, sample_id]))
                for sample_id in samples_by_condition[condition]
                if not pd.isna(intensities.at[site_key, sample_id])
            )
            if not values or not all(math.isfinite(value) for value in values):
                condition_means = []
                break
            condition_means.append(math.fsum(values) / float(len(values)))
        if condition_means:
            scores[site_key] = max(abs(value) for value in condition_means)
    return sps_science._ReferenceStabilityEvidence(
        dataset_id=reference.dataset_id,
        scores_by_site=scores,
    )


def test_obvious_stable_sites_rank_ahead_of_changing_sites() -> None:
    references = _hand_references()
    result = _run(references, config=SpsDiscoveryConfig(top_n=4))
    keys = tuple(str(value) for value in references[0].intensities.index)

    assert result.site_ranking[0].site_key == keys[0]
    assert result.site_ranking[-1].site_key == keys[2]
    assert result.site_ranking[0].consensus_stability_score > (
        result.site_ranking[-1].consensus_stability_score
    )
    assert tuple(
        record.consensus_stability_score for record in result.site_ranking
    ) == (
        pytest.approx(0.765625),
        pytest.approx(0.390625),
        pytest.approx(0.140625),
        pytest.approx(0.015625),
    )
    first_statistics = result.site_ranking[0].dataset_statistics
    assert tuple(item.stability_score for item in first_statistics) == (0.2, 0.15)


def test_row_reference_and_replicate_order_do_not_change_result() -> None:
    reference_a, reference_b = _hand_references()

    reordered_a = SpsReferenceDataset(
        dataset_id=reference_a.dataset_id,
        intensities=reference_a.intensities.iloc[::-1, ::-1],
        condition_by_sample={
            sample_id: str(condition)
            for sample_id, condition in reversed(reference_a.sample_conditions.items())
        },
        intensity_scale_state=reference_a.intensity_scale_state,
        organism=reference_a.organism,
        baseline_context=reference_a.baseline_context,
        reference_context=reference_a.reference_context,
        source_name=reference_a.source_name,
        source_version=reference_a.source_version,
        source_uri=reference_a.source_uri,
    )
    reordered_b = SpsReferenceDataset(
        dataset_id=reference_b.dataset_id,
        intensities=reference_b.intensities.iloc[[2, 0, 3, 1], [2, 3, 0, 1]],
        condition_by_sample={
            sample_id: str(reference_b.sample_conditions[sample_id])
            for sample_id in reversed(reference_b.intensities.columns)
        },
        intensity_scale_state=reference_b.intensity_scale_state,
        organism=reference_b.organism,
        baseline_context=reference_b.baseline_context,
        reference_context=reference_b.reference_context,
        source_name=reference_b.source_name,
        source_version=reference_b.source_version,
        source_uri=reference_b.source_uri,
    )

    expected = _run((reference_a, reference_b))
    reordered = _run((reordered_b, reordered_a))

    assert reordered == expected


def test_ties_have_equal_evidence_and_use_site_key_ascending_order() -> None:
    keys = _site_keys(3)
    columns_a = {
        "a_1": [0.1, 0.1, 2.0],
        "a_2": [0.2, 0.2, 3.0],
    }
    columns_b = {
        "b_1": [0.1, 0.1, 3.0],
        "b_2": [0.2, 0.2, 4.0],
    }
    references = (
        _reference(
            "a",
            site_keys=keys,
            columns=columns_a,
            conditions={"a_1": "one", "a_2": "two"},
        ),
        _reference(
            "b",
            site_keys=keys,
            columns=columns_b,
            conditions={"b_1": "one", "b_2": "two"},
        ),
    )

    result = _run(references, config=SpsDiscoveryConfig(top_n=3))
    tied = result.site_ranking[:2]

    assert tied[0].site_key < tied[1].site_key
    assert tied[0].consensus_stability_score == tied[1].consensus_stability_score
    assert all(item.rank == 1 for record in tied for item in record.dataset_statistics)


def test_missing_data_requires_each_condition_but_not_every_replicate() -> None:
    keys = _site_keys(3)
    references = []
    for dataset_number in range(3):
        prefix = f"d{dataset_number}"
        changing_missing = np.nan if dataset_number < 2 else 0.2
        references.append(
            _reference(
                prefix,
                site_keys=keys,
                columns={
                    f"{prefix}_c1": [np.nan, 0.1, 0.1],
                    f"{prefix}_c2": [0.1, 0.1, 0.1],
                    f"{prefix}_t1": [0.2, changing_missing, 2.0],
                    f"{prefix}_t2": [0.2, changing_missing, 2.0],
                },
                conditions={
                    f"{prefix}_c1": "control",
                    f"{prefix}_c2": "control",
                    f"{prefix}_t1": "treated",
                    f"{prefix}_t2": "treated",
                },
            )
        )

    result = _run(
        references,
        config=SpsDiscoveryConfig(
            top_n=10,
            minimum_reference_datasets=3,
            minimum_datasets_per_site=2,
        ),
    )

    assert tuple(record.site_key for record in result.site_ranking) == (
        str(keys[0]),
        str(keys[2]),
    )
    assert result.site_ranking[0].contributing_dataset_count == 3
    assert result.site_ranking[1].contributing_dataset_count == 3
    assert result.provenance.actual_control_count == 2
    assert result.provenance.requested_control_count == 10


def test_insufficient_rankable_overlap_after_missingness_fails_clearly() -> None:
    keys = _site_keys(2)
    references = tuple(
        _reference(
            dataset_id,
            site_keys=keys,
            columns={
                f"{dataset_id}_c": [0.1, 0.1],
                f"{dataset_id}_t": [np.nan, np.nan],
            },
            conditions={f"{dataset_id}_c": "control", f"{dataset_id}_t": "treated"},
        )
        for dataset_id in ("a", "b")
    )

    with pytest.raises(SpsDiscoveryValidationError) as caught:
        _run(references)

    issue = caught.value.validation_result.issues[0]
    assert issue.code == "insufficient_rankable_overlap"
    assert caught.value.validation_result.potential_overlap_site_count == 2


def test_infinite_reference_measurements_fail_validation() -> None:
    keys = _site_keys(2)
    references = (
        _reference(
            "finite",
            site_keys=keys,
            columns={"f_c": [0.1, 0.2], "f_t": [0.2, 0.3]},
            conditions={"f_c": "control", "f_t": "treated"},
        ),
        _reference(
            "infinite",
            site_keys=keys,
            columns={"i_c": [0.1, np.inf], "i_t": [0.2, 0.3]},
            conditions={"i_c": "control", "i_t": "treated"},
        ),
    )

    with pytest.raises(SpsDiscoveryValidationError) as caught:
        _run(references)

    assert "non_finite_intensity_data" in {
        issue.code for issue in caught.value.validation_result.issues
    }


def test_site_key_alignment_and_attrition_are_explicit() -> None:
    reference_a, reference_b = _hand_references()
    reordered_b = SpsReferenceDataset(
        dataset_id=reference_b.dataset_id,
        intensities=reference_b.intensities.iloc[[3, 1, 0, 2]],
        condition_by_sample=reference_b.sample_conditions,  # type: ignore[arg-type]
        intensity_scale_state=reference_b.intensity_scale_state,
        organism=reference_b.organism,
        baseline_context=reference_b.baseline_context,
        reference_context=reference_b.reference_context,
        source_name=reference_b.source_name,
        source_version=reference_b.source_version,
        source_uri=reference_b.source_uri,
    )

    result = _run((reference_a, reordered_b), config=SpsDiscoveryConfig(top_n=2))

    assert result.provenance.selection_boundaries.to_payload() == {
        "total_unique_sites": 4,
        "sites_meeting_dataset_overlap": 4,
        "sites_with_valid_stability": 4,
        "sites_ranked": 4,
        "sites_selected": 2,
    }
    assert tuple(
        (
            source.dataset_id,
            source.site_count,
            source.sites_passing_required_data,
            source.sites_entering_consensus,
        )
        for source in result.provenance.source_datasets
    ) == (("study_a", 4, 4, 4), ("study_b", 4, 4, 4))
    assert tuple(
        statistic.dataset_id for statistic in result.site_ranking[0].dataset_statistics
    ) == ("study_a", "study_b")
    assert len(result.control_site_set.annotations) == 2


def test_partial_reference_contributions_attrition_and_native_handoff() -> None:
    keys = tuple(str(value) for value in _site_keys(4))
    all_references, exact_minimum, below_minimum, other_eligible = keys
    references = (
        _reference(
            "study_a",
            site_keys=keys,
            columns={
                "a_control": [0.0, 0.0, 0.0, 0.0],
                "a_treated": [0.1, 0.2, 0.3, 0.4],
            },
            conditions={"a_control": "control", "a_treated": "treated"},
        ),
        _reference(
            "study_b",
            site_keys=(other_eligible, exact_minimum, all_references, below_minimum),
            columns={
                "b_control": [0.0, 0.0, 0.0, 0.0],
                "b_treated": [0.4, 0.2, 0.1, np.nan],
            },
            conditions={"b_control": "control", "b_treated": "treated"},
        ),
        _reference(
            "study_c",
            site_keys=(below_minimum, all_references, other_eligible, exact_minimum),
            columns={
                "c_control": [0.0, 0.0, 0.0, 0.0],
                "c_treated": [np.nan, 0.15, 0.45, np.nan],
            },
            conditions={"c_control": "control", "c_treated": "treated"},
        ),
    )
    config = SpsDiscoveryConfig(
        top_n=3,
        minimum_reference_datasets=3,
        minimum_datasets_per_site=2,
        minimum_shared_sites=3,
    )

    result = _run(references, config=config)

    records = {record.site_key: record for record in result.site_ranking}
    exact_record = records[exact_minimum]
    assert exact_record.contributing_dataset_count == config.minimum_datasets_per_site
    assert tuple(
        statistic.dataset_id for statistic in exact_record.dataset_statistics
    ) == ("study_a", "study_b")
    assert "study_c" not in {
        statistic.dataset_id for statistic in exact_record.dataset_statistics
    }
    assert below_minimum not in records
    assert below_minimum not in result.selected_site_keys
    assert exact_minimum in result.selected_site_keys
    assert records[all_references].contributing_dataset_count == 3
    assert _run(tuple(reversed(references)), config=config) == result

    boundaries = result.provenance.selection_boundaries
    assert boundaries.to_payload() == {
        "total_unique_sites": 4,
        "sites_meeting_dataset_overlap": 4,
        "sites_with_valid_stability": 3,
        "sites_ranked": 3,
        "sites_selected": 3,
    }
    sources = result.provenance.source_datasets
    assert tuple(
        (
            source.dataset_id,
            source.site_count,
            source.sites_passing_required_data,
            source.sites_entering_consensus,
        )
        for source in sources
    ) == (
        ("study_a", 4, 4, 3),
        ("study_b", 4, 3, 3),
        ("study_c", 4, 2, 2),
    )
    removed_for_insufficient_contribution = (
        boundaries.sites_meeting_dataset_overlap - boundaries.sites_with_valid_stability
    )
    assert removed_for_insufficient_contribution == 1
    assert sum(source.sites_entering_consensus or 0 for source in sources) == sum(
        record.contributing_dataset_count for record in result.site_ranking
    )
    assert (
        sum(source.sites_passing_required_data or 0 for source in sources)
        - sum(source.sites_entering_consensus or 0 for source in sources)
        == removed_for_insufficient_contribution
    )

    discovered_controls = result.control_site_set
    native_request = BatchCorrectionWorkflowRequest(
        phospho=pd.DataFrame(
            {
                "target_1": [1.0, 2.0, 3.0, 4.0],
                "target_2": [1.1, 2.1, 3.1, 4.1],
            },
            index=pd.Index(keys, name="site_key"),
        ),
        config=InternalBatchCorrectionRequest(
            method=InternalBatchCorrectionMethod.SPS_RUV_STYLE,
            batch_column="batch",
            condition_columns=("condition",),
            replicate_column=None,
            control_site_source=(
                InternalBatchCorrectionControlSiteSource.CALLER_SUPPLIED
            ),
            control_site_mode=InternalBatchCorrectionControlSiteMode.SITE_KEY_LIST,
            missing_value_policy=(
                InternalBatchCorrectionMissingValuePolicy.REJECT_MISSING
            ),
            imputation_policy=InternalBatchCorrectionImputationPolicy.NONE,
            n_unwanted_factors=1,
            stage_order=(
                InternalBatchCorrectionStageOrder.AFTER_MISSING_DATA_BEFORE_DOWNSTREAM
            ),
            diagnostics_enabled=True,
        ),
        sample_metadata=None,
        control_site_set=discovered_controls,
    )
    mapping = BatchCorrectionWorkflowControlSiteValidator().run(request=native_request)

    assert (
        tuple(row.site_key for row in mapping.row_eligibility if row.is_control)
        == result.selected_site_keys
    )
    assert mapping.control_status_by_site_key[below_minimum].value == "non_control"


def test_array_backed_stability_path_is_exactly_scalar_semantics_equivalent(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    keys = tuple(str(value) for value in _site_keys(6))
    references = (
        _reference(
            "layout_three_conditions",
            site_keys=keys,
            columns={
                "a_c1_r1": [0.1, 0.1, 0.2, 0.5, 0.3, 0.1],
                "a_c1_r2": [0.2, 0.2, 0.2, 0.6, 0.4, 0.2],
                "a_c2_r1": [0.1, 0.1, 0.3, 1.5, 0.2, 0.1],
                "a_c2_r2": [0.2, 0.2, 0.3, 1.6, 0.3, 0.2],
                "a_c3_r1": [0.1, 0.1, 0.4, 2.5, 0.1, 0.1],
                "a_c3_r2": [0.2, 0.2, 0.4, 2.6, 0.2, 0.2],
            },
            conditions={
                "a_c1_r1": "one",
                "a_c1_r2": "one",
                "a_c2_r1": "two",
                "a_c2_r2": "two",
                "a_c3_r1": "three",
                "a_c3_r2": "three",
            },
        ),
        _reference(
            "layout_two_conditions_reordered",
            site_keys=(keys[3], keys[1], keys[5], keys[0], keys[4], keys[2]),
            columns={
                "b_t_r3": [2.4, 0.1, np.nan, 0.1, 0.3, 0.4],
                "b_c_r2": [0.6, 0.2, 0.2, 0.2, 0.4, 0.2],
                "b_t_r1": [2.2, 0.1, np.nan, 0.1, 0.3, 0.4],
                "b_c_r1": [0.5, 0.1, 0.1, 0.1, 0.3, 0.2],
                "b_t_r2": [2.3, 0.2, np.nan, 0.2, 0.2, 0.3],
                "b_c_r3": [0.7, 0.3, 0.3, 0.3, 0.5, 0.2],
            },
            conditions={
                "b_t_r3": "treated",
                "b_c_r2": "control",
                "b_t_r1": "treated",
                "b_c_r1": "control",
                "b_t_r2": "treated",
                "b_c_r3": "control",
            },
        ),
        _reference(
            "layout_four_conditions_partial",
            site_keys=(keys[4], keys[2], keys[0], keys[5], keys[3], keys[1]),
            columns={
                "c_four": [np.nan, 0.4, 0.1, np.nan, 2.4, 0.1],
                "c_two": [np.nan, 0.2, 0.1, np.nan, 1.4, 0.1],
                "c_one": [np.nan, 0.1, 0.1, 0.1, 0.4, 0.1],
                "c_three": [np.nan, 0.3, 0.1, np.nan, 1.9, 0.1],
            },
            conditions={
                "c_four": "four",
                "c_two": "two",
                "c_one": "one",
                "c_three": "three",
            },
        ),
    )
    config = SpsDiscoveryConfig(
        top_n=4,
        minimum_reference_datasets=3,
        minimum_datasets_per_site=2,
        minimum_shared_sites=4,
    )

    optimized_evidence = tuple(
        sps_science._reference_stability_evidence(reference) for reference in references
    )
    scalar_evidence = tuple(
        _scalar_reference_stability_evidence(reference) for reference in references
    )
    optimized_result = _run(references, config=config)
    with monkeypatch.context() as context:
        context.setattr(
            sps_science,
            "_reference_stability_evidence",
            _scalar_reference_stability_evidence,
        )
        scalar_result = _run(references, config=config)

    assert optimized_evidence == scalar_evidence
    assert optimized_result == scalar_result
    assert optimized_result.selected_site_keys == scalar_result.selected_site_keys
    assert tuple(
        (
            record.site_key,
            record.consensus_rank,
            record.consensus_stability_score,
            record.contributing_dataset_count,
            record.dataset_statistics,
        )
        for record in optimized_result.site_ranking
    ) == tuple(
        (
            record.site_key,
            record.consensus_rank,
            record.consensus_stability_score,
            record.contributing_dataset_count,
            record.dataset_statistics,
        )
        for record in scalar_result.site_ranking
    )
    assert keys[5] not in optimized_result.selected_site_keys
    assert optimized_result.provenance.selection_boundaries.sites_ranked == 5


def test_larger_deterministic_population_selects_only_known_stable_sites() -> None:
    stable_count = 120
    changing_count = 80
    site_count = stable_count + changing_count
    keys = _site_keys(site_count)
    references: list[SpsReferenceDataset] = []
    positions = np.arange(site_count, dtype=float)
    for dataset_number in range(3):
        prefix = f"ref_{dataset_number}"
        stable = 0.01 + (positions[:stable_count] % 17.0) / 1000.0
        changing = 2.0 + (positions[stable_count:] % 23.0) / 10.0
        change = np.concatenate((stable, changing)) + dataset_number * 0.001
        references.append(
            _reference(
                prefix,
                site_keys=keys,
                columns={
                    f"{prefix}_c1": np.zeros(site_count).tolist(),
                    f"{prefix}_c2": np.zeros(site_count).tolist(),
                    f"{prefix}_t1": change.tolist(),
                    f"{prefix}_t2": change.tolist(),
                },
                conditions={
                    f"{prefix}_c1": "control",
                    f"{prefix}_c2": "control",
                    f"{prefix}_t1": "treated",
                    f"{prefix}_t2": "treated",
                },
            )
        )

    result = _run(
        references,
        config=SpsDiscoveryConfig(
            top_n=100,
            minimum_reference_datasets=3,
            minimum_datasets_per_site=3,
            minimum_shared_sites=100,
        ),
    )

    stable_keys = {str(value) for value in keys[:stable_count]}
    assert len(result.selected_site_keys) == 100
    assert set(result.selected_site_keys).issubset(stable_keys)
    assert result.provenance.selection_boundaries.sites_ranked == site_count
