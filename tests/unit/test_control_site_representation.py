from __future__ import annotations

from phospy.science.datasets.preprocessing.control_sites import (
    CONTROL_SITE_REASON_CONTROL_ANNOTATION_NOT_IN_SITE_ROWS,
    CONTROL_SITE_REASON_DUPLICATE_CONTROL_ANNOTATION,
    CONTROL_SITE_REASON_INVALID_WEIGHT,
    CONTROL_SITE_REASON_MISSING_CONTROL_ANNOTATION,
    ControlSiteAnnotation,
    ControlSiteSet,
    ControlSiteSourceMetadata,
    ControlSiteStatus,
)


def test_binary_control_annotations_map_to_site_key_rows() -> None:
    control_set = ControlSiteSet.from_binary_controls(
        {
            "AKT1_S473": True,
            "MAPK1_T202": False,
        },
        source_metadata=ControlSiteSourceMetadata(
            organism="human",
            identifier_namespace="site_key",
            source_version="manual-v1",
            license="caller local use",
            redistribution="not redistributed",
        ),
    )

    mapping = control_set.map_to_site_keys(("AKT1_S473", "MAPK1_T202"))

    assert mapping.control_status_by_site_key == {
        "AKT1_S473": ControlSiteStatus.CONTROL,
        "MAPK1_T202": ControlSiteStatus.NON_CONTROL,
    }
    assert [row.is_control for row in mapping.row_eligibility] == [True, False]
    assert mapping.row_eligibility[0].organism == "human"
    assert mapping.row_eligibility[0].identifier_namespace == "site_key"
    assert mapping.row_eligibility[0].source_version == "manual-v1"
    assert mapping.row_eligibility[0].license == "caller local use"


def test_control_metadata_missing_rationale_maps_to_eligibility() -> None:
    control_set = ControlSiteSet.from_site_keys(
        ("AKT1_S473",),
        source_metadata=ControlSiteSourceMetadata(
            organism="human",
            identifier_namespace="site_key",
            metadata_missing_reason={
                "source_version": "caller-supplied local controls have no version",
            },
        ),
    )

    mapping = control_set.map_to_site_keys(("AKT1_S473",))

    assert mapping.row_eligibility[0].metadata_missing_reason == {
        "source_version": "caller-supplied local controls have no version",
    }


def test_weighted_control_annotations_preserve_weight_without_correction() -> None:
    control_set = ControlSiteSet.from_weighted_controls(
        {
            "AKT1_S473": 0.25,
            "GSK3B_S9": 1.0,
        }
    )

    mapping = control_set.map_to_site_keys(("AKT1_S473", "GSK3B_S9", "MAPK1_T202"))

    assert mapping.control_weight_by_site_key == {
        "AKT1_S473": 0.25,
        "GSK3B_S9": 1.0,
    }
    assert mapping.row_eligibility[0].is_weighted_control is True
    assert mapping.row_eligibility[2].control_status is ControlSiteStatus.NON_CONTROL
    assert mapping.row_eligibility[2].valid is True


def test_grouped_control_annotations_preserve_group_labels() -> None:
    control_set = ControlSiteSet(
        annotations=(
            ControlSiteAnnotation(
                site_key="AKT1_S473",
                control_status="control",
                group="growth_factor",
            ),
            ControlSiteAnnotation(
                site_key="GSK3B_S9",
                control_status="control",
                group="growth_factor",
            ),
            ControlSiteAnnotation(
                site_key="MAPK1_T202",
                control_status="control",
                group="stress",
            ),
        )
    )

    mapping = control_set.map_to_site_keys(("AKT1_S473", "GSK3B_S9", "MAPK1_T202"))

    assert mapping.grouped_control_site_keys == {
        "growth_factor": ("AKT1_S473", "GSK3B_S9"),
        "stress": ("MAPK1_T202",),
    }


def test_missing_controls_are_represented_for_later_validation() -> None:
    control_set = ControlSiteSet(
        annotations=(ControlSiteAnnotation("AKT1_S473", control_status=True),)
    )

    mapping = control_set.map_to_site_keys(("AKT1_S473", "MAPK1_T202"))

    missing_row = mapping.row_eligibility[1]
    assert missing_row.site_key == "MAPK1_T202"
    assert missing_row.control_status is ControlSiteStatus.UNKNOWN
    assert missing_row.valid is False
    assert missing_row.reasons == (CONTROL_SITE_REASON_MISSING_CONTROL_ANNOTATION,)


def test_caller_supplied_control_not_in_rows_is_not_silently_dropped() -> None:
    control_set = ControlSiteSet.from_site_keys(("AKT1_S473", "MISSING_SITE"))

    mapping = control_set.map_to_site_keys(("AKT1_S473", "MAPK1_T202"))

    assert mapping.row_eligibility[1].control_status is ControlSiteStatus.NON_CONTROL
    assert mapping.row_eligibility[1].valid is True
    assert len(mapping.unmapped_annotations) == 1
    assert mapping.unmapped_annotations[0].site_key == "MISSING_SITE"
    assert mapping.unmapped_annotations[0].valid is False
    assert mapping.unmapped_annotations[0].reasons == (
        CONTROL_SITE_REASON_CONTROL_ANNOTATION_NOT_IN_SITE_ROWS,
    )


def test_duplicate_control_annotations_are_invalid_mapped_states() -> None:
    control_set = ControlSiteSet(
        annotations=(
            ControlSiteAnnotation("AKT1_S473", control_status=True),
            ControlSiteAnnotation("AKT1_S473", control_status=False),
        )
    )

    mapping = control_set.map_to_site_keys(("AKT1_S473",))
    row = mapping.row_eligibility[0]

    assert row.control_status is ControlSiteStatus.INVALID
    assert row.valid is False
    assert row.annotation_indices == (0, 1)
    assert row.reasons == (CONTROL_SITE_REASON_DUPLICATE_CONTROL_ANNOTATION,)


def test_structurally_invalid_annotations_are_represented_not_rejected() -> None:
    control_set = ControlSiteSet(
        annotations=(
            ControlSiteAnnotation(
                site_key="AKT1_S473",
                control_status=True,
                weight=-1.0,
            ),
        )
    )

    mapping = control_set.map_to_site_keys(("AKT1_S473",))
    row = mapping.row_eligibility[0]

    assert row.site_key == "AKT1_S473"
    assert row.control_status is ControlSiteStatus.CONTROL
    assert row.valid is False
    assert row.reasons == (CONTROL_SITE_REASON_INVALID_WEIGHT,)
