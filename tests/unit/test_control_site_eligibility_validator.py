from __future__ import annotations

import pandas as pd
import pytest

from phospy.errors import PhosPyInputError
from phospy.science.datasets.preprocessing.control_sites import (
    ControlSiteAnnotation,
    ControlSiteSet,
    ControlSiteSourceMetadata,
)
from phospy.science.references.models import Organism
from phospy.validation.workflows.batch_correction import (
    ControlSiteEligibilityValidator,
)


def _metadata(site_keys: tuple[str, ...], *, organism: str = "rat") -> pd.DataFrame:
    return pd.DataFrame(
        {"organism": [organism] * len(site_keys)},
        index=pd.Index(site_keys, name="site_key"),
    )


def _complete_source_metadata() -> ControlSiteSourceMetadata:
    return ControlSiteSourceMetadata(
        organism="rat",
        identifier_namespace="site_key",
        source_name="manual-curated-controls",
        source_version="manual-v1",
        license="caller local use",
        redistribution="not redistributed",
    )


def _reasoned_source_metadata() -> ControlSiteSourceMetadata:
    return ControlSiteSourceMetadata(
        organism="rat",
        identifier_namespace="site_key",
        metadata_missing_reason={
            "source_version": "caller-supplied local controls have no formal version",
            "license": "caller-supplied local controls are not licensed data",
            "redistribution": "caller-supplied local controls are not redistributed",
        },
    )


def test_control_site_validator_accepts_valid_caller_supplied_controls() -> None:
    site_keys = ("AKT1_S473", "GSK3B_S9", "MAPK1_T202")
    control_set = ControlSiteSet.from_site_keys(
        ("AKT1_S473", "GSK3B_S9"),
        source_metadata=_complete_source_metadata(),
    )

    mapping = ControlSiteEligibilityValidator().run(
        control_set=control_set,
        site_keys=site_keys,
        site_metadata=_metadata(site_keys),
        dataset_organism=Organism.RAT,
        method="sps_ruv_style",
        min_eligible_controls=2,
    )

    assert [row.site_key for row in mapping.row_eligibility if row.is_control] == [
        "AKT1_S473",
        "GSK3B_S9",
    ]
    assert mapping.row_eligibility[0].source_version == "manual-v1"


def test_control_site_validator_accepts_caller_metadata_with_explicit_rationale() -> (
    None
):
    site_keys = ("AKT1_S473", "GSK3B_S9")
    control_set = ControlSiteSet.from_site_keys(
        site_keys,
        source_metadata=_reasoned_source_metadata(),
    )

    mapping = ControlSiteEligibilityValidator().run(
        control_set=control_set,
        site_keys=site_keys,
        site_metadata=_metadata(site_keys),
        dataset_organism=Organism.RAT,
        method="sps_ruv_style",
        min_eligible_controls=2,
    )

    assert [row.site_key for row in mapping.row_eligibility if row.is_control] == [
        "AKT1_S473",
        "GSK3B_S9",
    ]
    assert mapping.row_eligibility[0].metadata_missing_reason["source_version"]


def test_control_site_validator_rejects_missing_control_mappings() -> None:
    site_keys = ("AKT1_S473", "MAPK1_T202")
    control_set = ControlSiteSet(
        annotations=(ControlSiteAnnotation("AKT1_S473", control_status=True),)
    )

    with pytest.raises(
        PhosPyInputError,
        match="missing control mappings for dataset site_key values 'MAPK1_T202'",
    ):
        ControlSiteEligibilityValidator().run(
            control_set=control_set,
            site_keys=site_keys,
            method="sps_ruv_style",
            min_eligible_controls=1,
        )


def test_control_site_validator_rejects_duplicate_control_mappings() -> None:
    control_set = ControlSiteSet(
        annotations=(
            ControlSiteAnnotation("AKT1_S473", control_status=True),
            ControlSiteAnnotation("AKT1_S473", control_status=False),
        )
    )

    with pytest.raises(
        PhosPyInputError,
        match="duplicate control mappings for site_key values 'AKT1_S473'",
    ):
        ControlSiteEligibilityValidator().run(
            control_set=control_set,
            site_keys=("AKT1_S473",),
            method="sps_ruv_style",
            min_eligible_controls=1,
        )


def test_control_site_validator_rejects_control_site_keys_absent_from_dataset() -> None:
    control_set = ControlSiteSet.from_site_keys(
        ("AKT1_S473", "MISSING_SITE"),
        source_metadata=_complete_source_metadata(),
    )

    with pytest.raises(
        PhosPyInputError,
        match="control site_key values are absent from the dataset: 'MISSING_SITE'",
    ):
        ControlSiteEligibilityValidator().run(
            control_set=control_set,
            site_keys=("AKT1_S473",),
            method="sps_ruv_style",
            min_eligible_controls=1,
        )


def test_control_site_validator_rejects_incompatible_organism_metadata() -> None:
    site_keys = ("AKT1_S473",)
    control_set = ControlSiteSet.from_site_keys(
        site_keys,
        source_metadata=ControlSiteSourceMetadata(
            organism="human",
            identifier_namespace="site_key",
            source_version="manual-v1",
            license="caller local use",
            redistribution="not redistributed",
        ),
    )

    with pytest.raises(
        PhosPyInputError,
        match=(
            "incompatible organism metadata for site_key 'AKT1_S473'; control "
            "source declares 'human', dataset expects 'rat'"
        ),
    ):
        ControlSiteEligibilityValidator().run(
            control_set=control_set,
            site_keys=site_keys,
            site_metadata=_metadata(site_keys, organism="rat"),
            dataset_organism=Organism.RAT,
            method="sps_ruv_style",
            min_eligible_controls=1,
        )


def test_control_site_validator_rejects_duplicate_site_metadata_index_values() -> None:
    site_metadata = pd.DataFrame(
        {"organism": ["rat", "human"]},
        index=pd.Index(["AKT1_S473", "AKT1_S473"], name="site_key"),
    )
    control_set = ControlSiteSet.from_site_keys(
        ("AKT1_S473",),
        source_metadata=_complete_source_metadata(),
    )

    with pytest.raises(
        PhosPyInputError,
        match="site_metadata\\.index values must be unique.*'AKT1_S473'",
    ):
        ControlSiteEligibilityValidator().run(
            control_set=control_set,
            site_keys=("AKT1_S473",),
            site_metadata=site_metadata,
            dataset_organism=Organism.RAT,
            method="sps_ruv_style",
            min_eligible_controls=1,
        )


def test_control_site_validator_rejects_blank_site_metadata_index_values() -> None:
    site_metadata = pd.DataFrame(
        {"organism": ["rat", "rat"]},
        index=pd.Index(["AKT1_S473", "  "], name="site_key"),
    )
    control_set = ControlSiteSet.from_site_keys(
        ("AKT1_S473",),
        source_metadata=_complete_source_metadata(),
    )

    with pytest.raises(
        PhosPyInputError,
        match="site_metadata\\.index values must not be blank.*blank positions \\[1\\]",
    ):
        ControlSiteEligibilityValidator().run(
            control_set=control_set,
            site_keys=("AKT1_S473",),
            site_metadata=site_metadata,
            dataset_organism=Organism.RAT,
            method="sps_ruv_style",
            min_eligible_controls=1,
        )


def test_control_site_validator_rejects_missing_site_metadata_for_accepted_control() -> (
    None
):
    site_keys = ("AKT1_S473", "GSK3B_S9")
    control_set = ControlSiteSet.from_site_keys(
        site_keys,
        source_metadata=_complete_source_metadata(),
    )

    with pytest.raises(
        PhosPyInputError,
        match=(
            "missing site_metadata rows for accepted control site_key values 'GSK3B_S9'"
        ),
    ):
        ControlSiteEligibilityValidator().run(
            control_set=control_set,
            site_keys=site_keys,
            site_metadata=_metadata(("AKT1_S473",)),
            dataset_organism=Organism.RAT,
            method="sps_ruv_style",
            min_eligible_controls=2,
        )


def test_control_site_validator_allows_extra_site_metadata_rows_after_index_checks() -> (
    None
):
    site_keys = ("AKT1_S473", "GSK3B_S9")
    site_metadata = pd.DataFrame(
        {"organism": ["rat", "rat", "human"]},
        index=pd.Index(["AKT1_S473", "GSK3B_S9", "EXTRA_SITE"], name="site_key"),
    )
    control_set = ControlSiteSet.from_site_keys(
        site_keys,
        source_metadata=_complete_source_metadata(),
    )

    mapping = ControlSiteEligibilityValidator().run(
        control_set=control_set,
        site_keys=site_keys,
        site_metadata=site_metadata,
        dataset_organism=Organism.RAT,
        method="sps_ruv_style",
        min_eligible_controls=2,
    )

    assert [row.site_key for row in mapping.row_eligibility if row.is_control] == [
        "AKT1_S473",
        "GSK3B_S9",
    ]


def test_control_site_validator_rejects_incompatible_identifier_namespace() -> None:
    control_set = ControlSiteSet.from_site_keys(
        ("AKT1_S473",),
        source_metadata=ControlSiteSourceMetadata(
            organism="rat",
            identifier_namespace="display_id",
            source_version="manual-v1",
            license="caller local use",
            redistribution="not redistributed",
        ),
    )

    with pytest.raises(
        PhosPyInputError,
        match="incompatible identifier namespace metadata for site_key 'AKT1_S473'",
    ):
        ControlSiteEligibilityValidator().run(
            control_set=control_set,
            site_keys=("AKT1_S473",),
            dataset_organism=Organism.RAT,
            method="sps_ruv_style",
            min_eligible_controls=1,
        )


def test_control_site_validator_rejects_ambiguous_control_labels() -> None:
    control_set = ControlSiteSet(
        annotations=(
            ControlSiteAnnotation(
                "AKT1_S473",
                control_status=True,
                label="housekeeping",
            ),
            ControlSiteAnnotation(
                "GSK3B_S9",
                control_status=False,
                label="housekeeping",
            ),
        )
    )

    with pytest.raises(
        PhosPyInputError,
        match="ambiguous control labels 'housekeeping'",
    ):
        ControlSiteEligibilityValidator().run(
            control_set=control_set,
            site_keys=("AKT1_S473", "GSK3B_S9"),
            method="sps_ruv_style",
            min_eligible_controls=1,
        )


def test_control_site_validator_rejects_ambiguous_control_metadata() -> None:
    control_set = ControlSiteSet(
        annotations=(
            ControlSiteAnnotation(
                "AKT1_S473",
                control_status=True,
                source_name="manual-curated-controls",
            ),
            ControlSiteAnnotation(
                "GSK3B_S9",
                control_status=True,
                source_name="second-control-source",
            ),
        ),
        source_metadata=_complete_source_metadata(),
    )

    with pytest.raises(
        PhosPyInputError,
        match="ambiguous control metadata for 'source_name'",
    ):
        ControlSiteEligibilityValidator().run(
            control_set=control_set,
            site_keys=("AKT1_S473", "GSK3B_S9"),
            site_metadata=_metadata(("AKT1_S473", "GSK3B_S9")),
            dataset_organism=Organism.RAT,
            method="sps_ruv_style",
            min_eligible_controls=2,
        )


def test_control_site_validator_rejects_too_few_eligible_controls_for_method() -> None:
    control_set = ControlSiteSet.from_site_keys(
        ("AKT1_S473",),
        source_metadata=_complete_source_metadata(),
    )

    with pytest.raises(
        PhosPyInputError,
        match=(
            "too few eligible controls for method 'sps_ruv_style' and "
            "n_unwanted_factors=1; required at least 2, observed 1"
        ),
    ):
        ControlSiteEligibilityValidator().run(
            control_set=control_set,
            site_keys=("AKT1_S473", "MAPK1_T202"),
            method="sps_ruv_style",
            min_eligible_controls=1,
            n_unwanted_factors=1,
        )


def test_control_site_validator_rejects_one_zero_weight_control() -> None:
    control_set = ControlSiteSet.from_weighted_controls(
        {"AKT1_S473": 0.0},
        source_metadata=_complete_source_metadata(),
    )

    with pytest.raises(
        PhosPyInputError,
        match="control weights must be positive finite values.*'AKT1_S473'",
    ):
        ControlSiteEligibilityValidator().run(
            control_set=control_set,
            site_keys=("AKT1_S473",),
            site_metadata=_metadata(("AKT1_S473",)),
            dataset_organism=Organism.RAT,
            method="sps_ruv_style",
            min_eligible_controls=1,
            supports_weights=True,
        )


def test_control_site_validator_rejects_mixed_positive_and_zero_weights() -> None:
    site_keys = ("AKT1_S473", "GSK3B_S9")
    control_set = ControlSiteSet.from_weighted_controls(
        {"AKT1_S473": 1.0, "GSK3B_S9": 0.0},
        source_metadata=_complete_source_metadata(),
    )

    with pytest.raises(
        PhosPyInputError,
        match="control weights must be positive finite values.*'GSK3B_S9'",
    ):
        ControlSiteEligibilityValidator().run(
            control_set=control_set,
            site_keys=site_keys,
            site_metadata=_metadata(site_keys),
            dataset_organism=Organism.RAT,
            method="sps_ruv_style",
            min_eligible_controls=2,
            supports_weights=True,
        )


def test_control_site_validator_accepts_all_positive_weights() -> None:
    site_keys = ("AKT1_S473", "GSK3B_S9")
    control_set = ControlSiteSet.from_weighted_controls(
        {"AKT1_S473": 1.0, "GSK3B_S9": 0.5},
        source_metadata=_complete_source_metadata(),
    )

    mapping = ControlSiteEligibilityValidator().run(
        control_set=control_set,
        site_keys=site_keys,
        site_metadata=_metadata(site_keys),
        dataset_organism=Organism.RAT,
        method="sps_ruv_style",
        min_eligible_controls=2,
        supports_weights=True,
    )

    assert mapping.control_weight_by_site_key == {
        "AKT1_S473": 1.0,
        "GSK3B_S9": 0.5,
    }


def test_control_site_validator_rejects_unsupported_weighted_grouped_controls() -> None:
    control_set = ControlSiteSet(
        annotations=(
            ControlSiteAnnotation(
                "AKT1_S473",
                control_status=True,
                weight=1.0,
                group="reference",
            ),
            ControlSiteAnnotation(
                "GSK3B_S9",
                control_status=True,
                weight=0.5,
                group="reference",
            ),
        ),
        source_metadata=_complete_source_metadata(),
    )

    with pytest.raises(
        PhosPyInputError,
        match="unsupported control weighting/grouping combination",
    ):
        ControlSiteEligibilityValidator().run(
            control_set=control_set,
            site_keys=("AKT1_S473", "GSK3B_S9"),
            method="sps_ruv_style",
            min_eligible_controls=2,
            supports_weights=True,
            supports_groups=True,
            supports_weighted_groups=False,
        )


def test_control_site_validator_rejects_missing_metadata_without_rationale() -> None:
    control_set = ControlSiteSet.from_site_keys(("AKT1_S473", "GSK3B_S9"))

    with pytest.raises(
        PhosPyInputError,
        match="missing control metadata.*metadata_missing_reason",
    ):
        ControlSiteEligibilityValidator().run(
            control_set=control_set,
            site_keys=("AKT1_S473", "GSK3B_S9"),
            site_metadata=_metadata(("AKT1_S473", "GSK3B_S9")),
            dataset_organism=Organism.RAT,
            method="sps_ruv_style",
            min_eligible_controls=2,
        )


def test_control_site_validator_rejects_incomplete_packaged_metadata() -> None:
    control_set = ControlSiteSet.from_site_keys(
        ("AKT1_S473", "GSK3B_S9"),
        source_metadata=ControlSiteSourceMetadata(
            source_type="packaged_reference",
            organism="rat",
            identifier_namespace="site_key",
            source_name="packaged-controls",
        ),
    )

    with pytest.raises(
        PhosPyInputError,
        match=(
            "packaged/reference/external control-source metadata is incomplete"
            ".*source_version.*license.*redistribution.*selection_method"
        ),
    ):
        ControlSiteEligibilityValidator().run(
            control_set=control_set,
            site_keys=("AKT1_S473", "GSK3B_S9"),
            site_metadata=_metadata(("AKT1_S473", "GSK3B_S9")),
            dataset_organism=Organism.RAT,
            method="sps_ruv_style",
            min_eligible_controls=2,
            control_site_source_type="packaged_reference",
        )


def test_control_site_validator_rejects_row_level_packaged_source_missing_license() -> (
    None
):
    control_set = ControlSiteSet(
        annotations=(
            ControlSiteAnnotation(
                "AKT1_S473",
                control_status=True,
                source_type="packaged_reference",
                organism="rat",
                identifier_namespace="site_key",
                source_name="packaged-controls",
                source_version="2026.1",
                redistribution="redistributable with package",
                selection_method="curated stable-control list",
            ),
        ),
        unannotated_site_status=False,
    )

    with pytest.raises(
        PhosPyInputError,
        match=(
            "packaged/reference/external control-source metadata is incomplete"
            ".*missing 'license'"
        ),
    ):
        ControlSiteEligibilityValidator().run(
            control_set=control_set,
            site_keys=("AKT1_S473",),
            site_metadata=_metadata(("AKT1_S473",)),
            dataset_organism=Organism.RAT,
            method="sps_ruv_style",
            min_eligible_controls=1,
            control_site_source_type="caller_supplied",
        )


def test_control_site_validator_rejects_row_level_packaged_source_missing_version() -> (
    None
):
    control_set = ControlSiteSet(
        annotations=(
            ControlSiteAnnotation(
                "AKT1_S473",
                control_status=True,
                source_type="packaged_reference",
                organism="rat",
                identifier_namespace="site_key",
                source_name="packaged-controls",
                license="CC-BY-4.0",
                redistribution="redistributable with package",
                selection_method="curated stable-control list",
            ),
        ),
        unannotated_site_status=False,
    )

    with pytest.raises(
        PhosPyInputError,
        match=(
            "packaged/reference/external control-source metadata is incomplete"
            ".*missing 'source_version'"
        ),
    ):
        ControlSiteEligibilityValidator().run(
            control_set=control_set,
            site_keys=("AKT1_S473",),
            site_metadata=_metadata(("AKT1_S473",)),
            dataset_organism=Organism.RAT,
            method="sps_ruv_style",
            min_eligible_controls=1,
            control_site_source_type="caller_supplied",
        )


def test_control_site_validator_top_level_caller_supplied_does_not_mask_row_packaged_metadata() -> (
    None
):
    control_set = ControlSiteSet(
        annotations=(
            ControlSiteAnnotation(
                "AKT1_S473",
                control_status=True,
                source_type="packaged_reference",
                organism="rat",
                identifier_namespace="site_key",
                source_name="packaged-controls",
                license="CC-BY-4.0",
                redistribution="redistributable with package",
                selection_method="curated stable-control list",
            ),
        ),
        unannotated_site_status=False,
    )

    with pytest.raises(
        PhosPyInputError,
        match=(
            "packaged/reference/external control-source metadata is incomplete"
            ".*missing 'source_version'"
        ),
    ):
        ControlSiteEligibilityValidator().run(
            control_set=control_set,
            site_keys=("AKT1_S473",),
            site_metadata=_metadata(("AKT1_S473",)),
            dataset_organism=Organism.RAT,
            method="sps_ruv_style",
            min_eligible_controls=1,
            control_site_source_type="caller_supplied",
        )


def test_control_site_validator_accepts_fully_described_external_controls() -> None:
    control_set = ControlSiteSet(
        annotations=(
            ControlSiteAnnotation(
                "AKT1_S473",
                control_status=True,
                source_type="external",
                organism="rat",
                identifier_namespace="site_key",
                source_name="published-control-panel",
                source_version="2026-06",
                license="CC-BY-4.0",
                redistribution="redistributable under CC-BY-4.0",
                selection_method="published stable-control panel",
            ),
        ),
        unannotated_site_status=False,
    )

    mapping = ControlSiteEligibilityValidator().run(
        control_set=control_set,
        site_keys=("AKT1_S473",),
        site_metadata=_metadata(("AKT1_S473",)),
        dataset_organism=Organism.RAT,
        method="sps_ruv_style",
        min_eligible_controls=1,
        control_site_source_type="caller_supplied",
    )

    assert mapping.row_eligibility[0].is_control
    assert mapping.row_eligibility[0].source_type == "external"
