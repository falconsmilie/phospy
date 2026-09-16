from __future__ import annotations

from dataclasses import replace

import pandas as pd
import pytest

from phospy.advanced import (
    SpsDatasetSiteStatistic,
    SpsDiscoveryConfig,
    SpsDiscoveryProvenance,
    SpsDiscoveryRequest,
    SpsDiscoveryResult,
    SpsDiscoveryValidationError,
    SpsDiscoveryWorkflow,
    SpsReferenceDataset,
    SpsSelectionBoundaryCounts,
    SpsSiteStabilityRecord,
)
from phospy.errors.input import PhosPyInputError
from phospy.science.datasets.preprocessing.control_sites import (
    ControlSiteSet,
    ControlSiteStatus,
)
from tests.support.site_keys import protein_site_key_index


def _site_keys(*proteins: str) -> pd.Index:
    return protein_site_key_index(
        protein_identifiers=list(proteins),
        sites=[f"S{position + 1}" for position in range(len(proteins))],
    )


def _reference(
    dataset_id: str,
    site_keys: pd.Index,
    *,
    conditions: dict[str, str] | None = None,
) -> SpsReferenceDataset:
    return SpsReferenceDataset(
        dataset_id=dataset_id,
        intensities=pd.DataFrame(
            {
                f"{dataset_id}_control": [1.0 + i for i in range(len(site_keys))],
                f"{dataset_id}_treated": [2.0 + i for i in range(len(site_keys))],
            },
            index=site_keys,
        ),
        condition_by_sample=(
            conditions
            if conditions is not None
            else {
                f"{dataset_id}_control": "control",
                f"{dataset_id}_treated": "treated",
            }
        ),
        source_name=f"reference-{dataset_id}",
        source_version="2026-01",
    )


def test_valid_multi_dataset_discovery_configuration() -> None:
    keys = _site_keys("P1", "P2", "P3")
    config = SpsDiscoveryConfig(
        top_n=2,
        minimum_reference_datasets=2,
        minimum_datasets_per_site=2,
        minimum_shared_sites=2,
    )

    request = SpsDiscoveryRequest(
        reference_datasets=(
            _reference("study_a", keys[:2]),
            _reference(
                "study_b",
                keys,
                conditions={
                    "study_b_control": "baseline",
                    "study_b_treated": "stimulated",
                },
            ),
        ),
        config=config,
    )
    validation = SpsDiscoveryWorkflow().require_valid(request)

    assert validation.valid is True
    assert request.config == config
    assert len(request.reference_datasets) == 2
    assert request.reference_datasets[0].intensities.index.name == "site_key"
    assert request.reference_datasets[0].sample_conditions != (
        request.reference_datasets[1].sample_conditions
    )


def test_too_few_reference_datasets_has_structured_validation() -> None:
    request = SpsDiscoveryRequest(
        reference_datasets=(_reference("study_a", _site_keys("P1", "P2")),),
        config=SpsDiscoveryConfig(),
    )

    with pytest.raises(SpsDiscoveryValidationError) as caught:
        SpsDiscoveryWorkflow().require_valid(request)

    result = caught.value.validation_result
    assert result.valid is False
    assert result.reference_dataset_count == 1
    assert "insufficient_reference_datasets" in {issue.code for issue in result.issues}


def test_duplicate_and_invalid_site_keys_have_structured_validation() -> None:
    valid_key = _site_keys("P1")[0]
    duplicate_index = pd.Index([valid_key, valid_key], name="site_key")
    invalid_index = pd.Index(["P1;S1;"], name="site_key")

    request = SpsDiscoveryRequest(
        reference_datasets=(
            _reference("duplicates", duplicate_index),
            _reference("invalid", invalid_index),
        ),
        config=SpsDiscoveryConfig(),
    )

    with pytest.raises(SpsDiscoveryValidationError) as caught:
        SpsDiscoveryWorkflow().require_valid(request)

    codes = {issue.code for issue in caught.value.validation_result.issues}
    assert "duplicate_site_key" in codes
    assert "invalid_site_key" in codes


def test_invalid_reference_entry_has_structured_validation() -> None:
    request = SpsDiscoveryRequest(
        reference_datasets=(
            _reference("study_a", _site_keys("P1", "P2")),
            object(),  # type: ignore[arg-type]
        ),
        config=SpsDiscoveryConfig(),
    )

    with pytest.raises(SpsDiscoveryValidationError) as caught:
        SpsDiscoveryWorkflow().require_valid(request)

    assert "invalid_reference_dataset" in {
        issue.code for issue in caught.value.validation_result.issues
    }


def test_duplicate_reference_dataset_ids_have_structured_validation() -> None:
    keys = _site_keys("P1", "P2")

    request = SpsDiscoveryRequest(
        reference_datasets=(
            _reference("duplicated", keys),
            _reference("duplicated", keys),
        ),
        config=SpsDiscoveryConfig(),
    )

    with pytest.raises(SpsDiscoveryValidationError) as caught:
        SpsDiscoveryWorkflow().require_valid(request)

    assert "duplicate_reference_dataset_id" in {
        issue.code for issue in caught.value.validation_result.issues
    }


def test_missing_condition_metadata_has_structured_validation() -> None:
    keys = _site_keys("P1", "P2")

    request = SpsDiscoveryRequest(
        reference_datasets=(
            _reference("study_a", keys),
            _reference(
                "study_b",
                keys,
                conditions={"study_b_control": "control"},
            ),
        ),
        config=SpsDiscoveryConfig(),
    )

    with pytest.raises(SpsDiscoveryValidationError) as caught:
        SpsDiscoveryWorkflow().require_valid(request)

    assert "unresolvable_condition_metadata" in {
        issue.code for issue in caught.value.validation_result.issues
    }


def test_absent_condition_metadata_has_structured_validation() -> None:
    keys = _site_keys("P1", "P2")
    missing_conditions = SpsReferenceDataset(
        dataset_id="missing_conditions",
        intensities=pd.DataFrame(
            {
                "sample_a": [1.0, 2.0],
                "sample_b": [2.0, 3.0],
            },
            index=keys,
        ),
        condition_by_sample=None,  # type: ignore[arg-type]
    )

    request = SpsDiscoveryRequest(
        reference_datasets=(
            _reference("study_a", keys),
            missing_conditions,
        ),
        config=SpsDiscoveryConfig(),
    )

    with pytest.raises(SpsDiscoveryValidationError) as caught:
        SpsDiscoveryWorkflow().require_valid(request)

    assert "missing_condition_metadata" in {
        issue.code for issue in caught.value.validation_result.issues
    }


def test_invalid_condition_mapping_has_structured_validation() -> None:
    with pytest.raises(SpsDiscoveryValidationError) as caught:
        SpsReferenceDataset(
            dataset_id="invalid_conditions",
            intensities=pd.DataFrame(
                {"sample_a": [1.0], "sample_b": [2.0]},
                index=_site_keys("P1"),
            ),
            condition_by_sample={1: "control"},  # type: ignore[dict-item]
        )

    assert caught.value.validation_result.issues[0].code == (
        "invalid_condition_metadata"
    )


@pytest.mark.parametrize(
    "intensities",
    [
        object(),
        pd.DataFrame(index=_site_keys("P1", "P2")),
    ],
    ids=["malformed", "empty"],
)
def test_missing_intensity_data_has_structured_workflow_validation(
    intensities: object,
) -> None:
    invalid_reference = SpsReferenceDataset(
        dataset_id="invalid_intensities",
        intensities=intensities,  # type: ignore[arg-type]
        condition_by_sample={},
    )
    request = SpsDiscoveryRequest(
        reference_datasets=(
            _reference("valid", _site_keys("P1", "P2")),
            invalid_reference,
        ),
        config=SpsDiscoveryConfig(),
    )

    validation = SpsDiscoveryWorkflow().validate(request)

    assert validation.valid is False
    assert "missing_intensity_data" in {issue.code for issue in validation.issues}


@pytest.mark.parametrize(
    "intensity_values",
    [
        (["low", "high"], ["medium", "higher"]),
        ([True, False], [False, True]),
    ],
    ids=["non-numeric", "boolean"],
)
def test_non_numeric_intensity_data_has_structured_workflow_validation(
    intensity_values: tuple[list[object], list[object]],
) -> None:
    control_values, treated_values = intensity_values
    keys = _site_keys("P1", "P2")
    invalid_reference = SpsReferenceDataset(
        dataset_id="invalid_intensities",
        intensities=pd.DataFrame(
            {"control": control_values, "treated": treated_values},
            index=keys,
        ),
        condition_by_sample={"control": "control", "treated": "treated"},
    )
    request = SpsDiscoveryRequest(
        reference_datasets=(
            _reference("valid", keys),
            invalid_reference,
        ),
        config=SpsDiscoveryConfig(),
    )

    validation = SpsDiscoveryWorkflow().validate(request)

    assert validation.valid is False
    assert "non_numeric_intensity_data" in {issue.code for issue in validation.issues}


@pytest.mark.parametrize("top_n", [0, -1, True])
def test_invalid_top_n_is_rejected(top_n: int) -> None:
    with pytest.raises(PhosPyInputError, match="top_n"):
        SpsDiscoveryConfig(top_n=top_n)


@pytest.mark.parametrize("value", [1, 0, -1, True])
def test_invalid_minimum_reference_datasets_is_rejected(value: int) -> None:
    with pytest.raises(PhosPyInputError, match="minimum_reference_datasets"):
        SpsDiscoveryConfig(minimum_reference_datasets=value)


@pytest.mark.parametrize("value", [1, 0, -1, True])
def test_invalid_minimum_datasets_per_site_is_rejected(value: int) -> None:
    with pytest.raises(PhosPyInputError, match="minimum_datasets_per_site"):
        SpsDiscoveryConfig(minimum_datasets_per_site=value)


@pytest.mark.parametrize("value", [0, -1, True])
def test_invalid_minimum_shared_sites_is_rejected(value: int) -> None:
    with pytest.raises(PhosPyInputError, match="minimum_shared_sites"):
        SpsDiscoveryConfig(minimum_shared_sites=value)


def test_unsupported_tie_handling_is_rejected() -> None:
    with pytest.raises(PhosPyInputError, match="tie_handling"):
        SpsDiscoveryConfig(tie_handling="input_order")  # type: ignore[arg-type]


def test_unsupported_selection_method_is_rejected() -> None:
    with pytest.raises(PhosPyInputError, match="selection_method"):
        SpsDiscoveryConfig(selection_method="lowest_variance")  # type: ignore[arg-type]


def test_insufficient_potential_overlap_has_structured_validation() -> None:
    request = SpsDiscoveryRequest(
        reference_datasets=(
            _reference("study_a", _site_keys("P1", "P2")),
            _reference("study_b", _site_keys("P3", "P4")),
        ),
        config=SpsDiscoveryConfig(minimum_shared_sites=1),
    )

    with pytest.raises(SpsDiscoveryValidationError) as caught:
        SpsDiscoveryWorkflow().require_valid(request)

    validation = caught.value.validation_result
    assert validation.potential_overlap_site_count == 0
    assert "insufficient_potential_overlap" in {
        issue.code for issue in validation.issues
    }


def _provenance_and_ranking() -> tuple[
    tuple[str, str],
    SpsDiscoveryProvenance,
    tuple[SpsSiteStabilityRecord, ...],
]:
    keys = _site_keys("P1", "P2", "P3")
    references = (
        _reference("study_a", keys),
        _reference("study_b", keys),
    )
    config = SpsDiscoveryConfig(top_n=2, minimum_shared_sites=2)
    request = SpsDiscoveryRequest(reference_datasets=references, config=config)
    provenance = SpsDiscoveryWorkflow().assemble_provenance(
        request,
        actual_control_count=2,
        selection_boundaries=SpsSelectionBoundaryCounts(
            total_unique_sites=3,
            sites_meeting_dataset_overlap=3,
            sites_with_valid_stability=3,
            sites_ranked=3,
            sites_selected=2,
        ),
    )
    ranking = tuple(
        SpsSiteStabilityRecord(
            site_key=str(site_key),
            consensus_rank=rank,
            consensus_stability_score=score,
            contributing_dataset_count=2,
            dataset_statistics=(
                SpsDatasetSiteStatistic("study_a", rank, score),
                SpsDatasetSiteStatistic("study_b", rank, score + 0.01),
            ),
        )
        for rank, (site_key, score) in enumerate(
            zip(keys, (0.1, 0.2, 0.3), strict=True),
            start=1,
        )
    )
    return (str(keys[0]), str(keys[1])), provenance, ranking


def test_request_construction_is_passive_for_contextually_invalid_payload() -> None:
    request = SpsDiscoveryRequest(
        reference_datasets=None,  # type: ignore[arg-type]
        config=SpsDiscoveryConfig(),
    )

    assert request.reference_datasets is None


def test_malformed_reference_collection_has_structured_validation() -> None:
    request = SpsDiscoveryRequest(
        reference_datasets=None,  # type: ignore[arg-type]
        config=SpsDiscoveryConfig(),
    )

    with pytest.raises(SpsDiscoveryValidationError) as caught:
        SpsDiscoveryWorkflow().require_valid(request)

    issue = caught.value.validation_result.issues[0]
    assert issue.code == "invalid_reference_collection"
    assert issue.field_name == "reference_datasets"


def test_reference_matrix_isolated_from_input_and_accessor_mutation() -> None:
    keys = _site_keys("P1", "P2")
    original = pd.DataFrame(
        {"control": [1.0, 2.0], "treated": [2.0, 3.0]},
        index=keys,
    )
    reference = SpsReferenceDataset(
        dataset_id="study_a",
        intensities=original,
        condition_by_sample={"control": "control", "treated": "treated"},
    )

    original.index = pd.Index(["display-a", "display-b"], name="display_id")
    original.iloc[0, 0] = 999.0
    public_snapshot = reference.intensities
    public_snapshot.index = pd.Index(["other-a", "other-b"], name="display_id")
    public_snapshot.iloc[0, 0] = 888.0

    assert reference.intensities.index.equals(keys)
    assert reference.intensities.index.name == "site_key"
    assert reference.intensities.iloc[0, 0] == 1.0


def test_accessor_mutation_cannot_change_validated_provenance() -> None:
    keys = _site_keys("P1", "P2")
    references = (_reference("study_a", keys), _reference("study_b", keys))
    request = SpsDiscoveryRequest(
        reference_datasets=references,
        config=SpsDiscoveryConfig(top_n=1),
    )
    workflow = SpsDiscoveryWorkflow()
    workflow.require_valid(request)

    exposed = request.reference_datasets[0].intensities
    exposed.index = pd.Index(["display-a", "display-b"], name="display_id")
    provenance = workflow.assemble_provenance(
        request,
        actual_control_count=1,
        selection_boundaries=SpsSelectionBoundaryCounts(2, 2, 2, 2, 1),
    )

    assert provenance.source_datasets[0].intensity_fingerprint.index_name == "site_key"


def test_invalid_reference_cannot_produce_source_provenance() -> None:
    invalid = SpsReferenceDataset(
        dataset_id="invalid",
        intensities=pd.DataFrame(
            {"control": [1.0], "treated": [2.0]},
            index=pd.Index(["not-a-site-key"], name="display_id"),
        ),
        condition_by_sample={"control": "control", "treated": "treated"},
    )
    request = SpsDiscoveryRequest(
        reference_datasets=(invalid, _reference("valid", _site_keys("P1"))),
        config=SpsDiscoveryConfig(top_n=1),
    )

    with pytest.raises(SpsDiscoveryValidationError) as caught:
        SpsDiscoveryWorkflow().assemble_provenance(
            request,
            actual_control_count=1,
            selection_boundaries=SpsSelectionBoundaryCounts(2, 1, 1, 1, 1),
        )

    issue = next(
        issue
        for issue in caught.value.validation_result.issues
        if issue.dataset_id == "invalid"
    )
    assert issue.code == "missing_site_key_identity"
    assert issue.field_name == "intensities.index"


def test_config_and_provenance_comparison_and_serialization_are_deterministic() -> None:
    selected, provenance, _ = _provenance_and_ranking()
    selected_again, provenance_again, _ = _provenance_and_ranking()

    assert selected_again == selected
    assert provenance_again.config == provenance.config
    assert provenance_again == provenance
    assert provenance_again.to_payload() == provenance.to_payload()
    assert provenance.to_payload()["algorithm_parameters"] == {
        "top_n": 2,
        "minimum_reference_datasets": 2,
        "minimum_datasets_per_site": 2,
        "minimum_shared_sites": 2,
        "tie_handling": "site_key_ascending",
        "selection_method": "consensus_stability",
    }
    assert tuple(source.dataset_id for source in provenance.source_datasets) == (
        "study_a",
        "study_b",
    )
    assert all(
        source.intensity_fingerprint.index_name == "site_key"
        for source in provenance.source_datasets
    )
    assert provenance.requested_control_count == 2
    assert provenance.actual_control_count == 2
    assert provenance.selection_boundaries.to_payload() == {
        "total_unique_sites": 3,
        "sites_meeting_dataset_overlap": 3,
        "sites_with_valid_stability": 3,
        "sites_ranked": 3,
        "sites_selected": 2,
    }


def test_discovery_result_converts_to_existing_control_site_set() -> None:
    selected, provenance, ranking = _provenance_and_ranking()
    result = SpsDiscoveryResult(
        selected_site_keys=selected,
        site_ranking=ranking,
        provenance=provenance,
    )

    controls = result.to_control_site_set()

    assert isinstance(controls, ControlSiteSet)
    assert result.control_site_set == controls
    assert tuple(annotation.site_key for annotation in controls.annotations) == selected
    assert all(
        annotation.control_status is ControlSiteStatus.CONTROL
        for annotation in controls.annotations
    )
    assert controls.source_metadata.identifier_namespace == "site_key"
    assert controls.source_metadata.source_type == "sps_discovery"
    assert controls.source_metadata.selection_method == "consensus_stability"


def test_site_stability_record_rejects_noncanonical_site_key() -> None:
    _, _, ranking = _provenance_and_ranking()
    record = ranking[0]

    with pytest.raises(PhosPyInputError, match="exact canonical site_key"):
        SpsSiteStabilityRecord(
            site_key=f" {record.site_key}",
            consensus_rank=record.consensus_rank,
            consensus_stability_score=record.consensus_stability_score,
            contributing_dataset_count=record.contributing_dataset_count,
            dataset_statistics=record.dataset_statistics,
        )


def test_site_stability_record_rejects_invalid_nested_statistic_type() -> None:
    site_key = str(_site_keys("P1")[0])

    with pytest.raises(
        PhosPyInputError,
        match="dataset_statistics must contain only SpsDatasetSiteStatistic",
    ):
        SpsSiteStabilityRecord(
            site_key=site_key,
            consensus_rank=1,
            consensus_stability_score=0.1,
            contributing_dataset_count=1,
            dataset_statistics=(object(),),  # type: ignore[arg-type]
        )


def test_provenance_rejects_requested_count_different_from_top_n() -> None:
    _, provenance, _ = _provenance_and_ranking()

    with pytest.raises(PhosPyInputError, match="requested_control_count"):
        replace(provenance, requested_control_count=1)


def test_provenance_rejects_actual_count_different_from_boundary_count() -> None:
    _, provenance, _ = _provenance_and_ranking()
    boundaries = replace(provenance.selection_boundaries, sites_selected=1)

    with pytest.raises(PhosPyInputError, match="sites_selected"):
        replace(provenance, selection_boundaries=boundaries)


def test_result_rejects_selected_count_different_from_actual_count() -> None:
    selected, provenance, ranking = _provenance_and_ranking()
    one_selected_provenance = replace(
        provenance,
        actual_control_count=1,
        selection_boundaries=replace(
            provenance.selection_boundaries,
            sites_selected=1,
        ),
    )

    with pytest.raises(PhosPyInputError, match="actual_control_count"):
        SpsDiscoveryResult(
            selected_site_keys=selected,
            site_ranking=ranking,
            provenance=one_selected_provenance,
        )


def test_result_rejects_ranked_count_different_from_ranking_length() -> None:
    selected, provenance, ranking = _provenance_and_ranking()
    shortened_ranking = ranking[:2]

    with pytest.raises(PhosPyInputError, match="sites_ranked"):
        SpsDiscoveryResult(
            selected_site_keys=selected,
            site_ranking=shortened_ranking,
            provenance=provenance,
        )


def test_result_rejects_site_below_minimum_contributing_datasets() -> None:
    selected, provenance, ranking = _provenance_and_ranking()
    first = ranking[0]
    insufficient = replace(
        first,
        contributing_dataset_count=1,
        dataset_statistics=first.dataset_statistics[:1],
    )

    with pytest.raises(PhosPyInputError, match="minimum_datasets_per_site"):
        SpsDiscoveryResult(
            selected_site_keys=selected,
            site_ranking=(insufficient, *ranking[1:]),
            provenance=provenance,
        )


def test_result_rejects_unknown_source_dataset_id() -> None:
    selected, provenance, ranking = _provenance_and_ranking()
    first = ranking[0]
    unknown_source = replace(
        first,
        dataset_statistics=(
            first.dataset_statistics[0],
            replace(first.dataset_statistics[1], dataset_id="unknown-study"),
        ),
    )

    with pytest.raises(PhosPyInputError, match="unknown source dataset ids"):
        SpsDiscoveryResult(
            selected_site_keys=selected,
            site_ranking=(unknown_source, *ranking[1:]),
            provenance=provenance,
        )


def test_discovery_result_rejects_noncanonical_selected_site_key() -> None:
    selected, provenance, ranking = _provenance_and_ranking()

    with pytest.raises(PhosPyInputError, match="exact canonical site_key"):
        SpsDiscoveryResult(
            selected_site_keys=(f"{selected[0]} ", selected[1]),
            site_ranking=ranking,
            provenance=provenance,
        )


def test_discovery_result_rejects_reversed_equal_score_site_key_tie() -> None:
    _, provenance, ranking = _provenance_and_ranking()
    reversed_tie_ranking = (
        SpsSiteStabilityRecord(
            site_key=ranking[1].site_key,
            consensus_rank=1,
            consensus_stability_score=0.1,
            contributing_dataset_count=ranking[1].contributing_dataset_count,
            dataset_statistics=ranking[1].dataset_statistics,
        ),
        SpsSiteStabilityRecord(
            site_key=ranking[0].site_key,
            consensus_rank=2,
            consensus_stability_score=0.1,
            contributing_dataset_count=ranking[0].contributing_dataset_count,
            dataset_statistics=ranking[0].dataset_statistics,
        ),
        ranking[2],
    )

    with pytest.raises(PhosPyInputError, match="ascending site_key tie handling"):
        SpsDiscoveryResult(
            selected_site_keys=tuple(
                record.site_key for record in reversed_tie_ranking[:2]
            ),
            site_ranking=reversed_tie_ranking,
            provenance=provenance,
        )
