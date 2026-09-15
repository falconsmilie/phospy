from __future__ import annotations

import pandas as pd
import pytest

from phospy.errors.input import PhosPyInputError
from phospy.science.datasets.preprocessing.group_coverage_metadata import (
    GroupCoverageFilterMetadataValidator,
    ResolvedGroupCoverageFilterMetadata,
)
from phospy.science.datasets.preprocessing.sample_group_metadata import (
    ResolvedSampleGroups,
    SampleGroupMetadataResolver,
)
from phospy.validation.datasets import group_coverage_filter as validation_facade


def _phospho(*, columns: pd.Index | None = None) -> pd.DataFrame:
    sample_columns = (
        pd.Index(["sample_a", "sample_b", "sample_c", "sample_d"], name="sample")
        if columns is None
        else columns
    )
    return pd.DataFrame([[1.0, 2.0, 3.0, 4.0]], columns=sample_columns)


def _metadata(*, index: pd.Index | None = None) -> pd.DataFrame:
    sample_index = _phospho().columns if index is None else index
    groups = {
        "sample_a": "control",
        "sample_b": "treated",
        "sample_c": "control",
        "sample_d": "treated",
    }
    return pd.DataFrame(
        {"condition": [groups[str(sample).strip()] for sample in sample_index]},
        index=sample_index,
    )


def _resolve(
    *,
    phospho: pd.DataFrame | None = None,
    sample_metadata: pd.DataFrame | None = None,
    group_column: str | None = "condition",
) -> ResolvedSampleGroups:
    return SampleGroupMetadataResolver().run(
        phospho=_phospho() if phospho is None else phospho,
        sample_metadata=_metadata() if sample_metadata is None else sample_metadata,
        group_column=group_column,
    )


def test_resolver_aligns_groups_to_phospho_sample_order() -> None:
    resolved = _resolve()

    assert resolved.group_column == "condition"
    assert resolved.sample_order == (
        "sample_a",
        "sample_b",
        "sample_c",
        "sample_d",
    )
    assert resolved.group_by_sample == {
        "sample_a": "control",
        "sample_b": "treated",
        "sample_c": "control",
        "sample_d": "treated",
    }


def test_resolver_requires_sample_metadata() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="sample-group metadata resolution requires sample_metadata input data",
    ):
        SampleGroupMetadataResolver().run(
            phospho=_phospho(),
            sample_metadata=None,
            group_column="condition",
        )


@pytest.mark.parametrize("group_column", [None, "", "  "])
def test_resolver_requires_non_blank_group_column(group_column: str | None) -> None:
    with pytest.raises(
        PhosPyInputError,
        match="group_column must be a non-empty string",
    ):
        _resolve(group_column=group_column)


def test_resolver_reorders_metadata_rows_and_preserves_stable_group_order() -> None:
    metadata_order = pd.Index(
        ["sample_d", "sample_c", "sample_b", "sample_a"],
        name="sample",
    )

    resolved = _resolve(sample_metadata=_metadata(index=metadata_order))

    assert resolved.group_labels == ("control", "treated")
    assert resolved.sample_order_by_group == {
        "control": ("sample_a", "sample_c"),
        "treated": ("sample_b", "sample_d"),
    }


@pytest.mark.parametrize(
    ("metadata_index", "message"),
    [
        (
            pd.Index(["sample_a", "sample_b", "sample_c"], name="sample"),
            "missing rows for grouped samples: 'sample_d'",
        ),
        (
            pd.Index(
                ["sample_a", "sample_b", "sample_c", "sample_d", "sample_extra"],
                name="sample",
            ),
            "rows not present in phospho columns.*'sample_extra'",
        ),
    ],
)
def test_resolver_rejects_missing_or_extra_metadata_rows(
    metadata_index: pd.Index,
    message: str,
) -> None:
    groups = ["control", "treated", "control", "treated", "control"]
    metadata = pd.DataFrame(
        {"condition": groups[: len(metadata_index)]},
        index=metadata_index,
    )

    with pytest.raises(PhosPyInputError, match=message):
        _resolve(sample_metadata=metadata)


@pytest.mark.parametrize(
    ("phospho_columns", "metadata_index", "field_name"),
    [
        (
            pd.Index(["sample_a", " sample_a ", "sample_c", "sample_d"]),
            None,
            "phospho.columns",
        ),
        (
            None,
            pd.Index(["sample_a", " sample_a ", "sample_c", "sample_d"]),
            "sample_metadata.index",
        ),
    ],
)
def test_resolver_rejects_duplicate_normalized_sample_labels(
    phospho_columns: pd.Index | None,
    metadata_index: pd.Index | None,
    field_name: str,
) -> None:
    phospho = _phospho(columns=phospho_columns) if phospho_columns is not None else None
    metadata = (
        pd.DataFrame(
            {"condition": ["control", "treated", "control", "treated"]},
            index=metadata_index,
        )
        if metadata_index is not None
        else None
    )

    with pytest.raises(
        PhosPyInputError,
        match=rf"{field_name} contains duplicate sample labels.*'sample_a'",
    ):
        _resolve(phospho=phospho, sample_metadata=metadata)


def test_resolver_rejects_missing_group_column() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="group_column references missing sample_metadata column 'condition'",
    ):
        _resolve(sample_metadata=pd.DataFrame(index=_phospho().columns))


def test_resolver_rejects_duplicate_group_columns() -> None:
    metadata = pd.DataFrame(
        [
            ["control", "control"],
            ["treated", "treated"],
            ["control", "control"],
            ["treated", "treated"],
        ],
        columns=["condition", "condition"],
        index=_phospho().columns,
    )

    with pytest.raises(
        PhosPyInputError,
        match="group_column resolves to duplicate sample_metadata column 'condition'",
    ):
        _resolve(sample_metadata=metadata)


@pytest.mark.parametrize(
    ("invalid_label", "label_kind"),
    [(None, "missing"), (float("nan"), "missing"), ("  ", "blank")],
)
def test_resolver_rejects_missing_and_blank_group_labels(
    invalid_label: object,
    label_kind: str,
) -> None:
    metadata = _metadata()
    metadata.loc["sample_b", "condition"] = invalid_label

    with pytest.raises(
        PhosPyInputError,
        match=rf"{label_kind} group labels for samples: 'sample_b'",
    ):
        _resolve(sample_metadata=metadata)


def test_group_coverage_validator_composes_resolver_then_checks_threshold() -> None:
    validator = GroupCoverageFilterMetadataValidator(
        sample_group_resolver=SampleGroupMetadataResolver()
    )

    with pytest.raises(
        PhosPyInputError,
        match=r"min_groups_passing_threshold.*threshold=3, observed_groups=2",
    ):
        validator.run(
            phospho=_phospho(),
            sample_metadata=_metadata(),
            group_column="condition",
            min_groups_passing_threshold=3,
        )


def test_validation_facade_reuses_preprocessing_owned_coverage_objects() -> None:
    assert (
        validation_facade.GroupCoverageFilterMetadataValidator
        is GroupCoverageFilterMetadataValidator
    )
    assert (
        validation_facade.ResolvedGroupCoverageFilterMetadata
        is ResolvedGroupCoverageFilterMetadata
    )
