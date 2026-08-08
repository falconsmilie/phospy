from __future__ import annotations

import pandas as pd
import pandas.testing as pdt
import pytest

from phospy.contracts.configs import (
    DATASET_SITE_MATRIX_DUPLICATE_POLICY_AGGREGATE_MEAN,
    DATASET_SITE_MATRIX_DUPLICATE_POLICY_AGGREGATE_MEDIAN,
    DATASET_SITE_MATRIX_DUPLICATE_POLICY_ERROR,
    DATASET_SITE_MATRIX_DUPLICATE_POLICY_FIRST,
    DATASET_SITE_MATRIX_DUPLICATE_POLICY_MAX_MEAN_SIGNAL,
)
from phospy.errors.input import PhosPyInputError
from phospy.science.datasets.preprocessing.models import DuplicateSiteResolutionResult
from phospy.science.datasets.preprocessing.stages.site_matrix_components import (
    DuplicateSiteResolver,
)
from phospy.science.sites.site_keys import ProteinScopedPhosphositeKey, encode_site_key
from tests.support.site_keys import protein_site_key

_MAPK14_KEY = protein_site_key(protein_identifier="PROT_A", site="Y182")
_AKT1_KEY = protein_site_key(protein_identifier="PROT_C", site="T308")


def _duplicate_policy_inputs() -> tuple[
    pd.DataFrame, pd.DataFrame, pd.Series, pd.Series
]:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 3.0, 5.0],
            "sample_b": [2.0, 4.0, 6.0],
        },
        index=pd.Index(["row_a", "row_b", "row_c"], name="source_row"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "MAPK14", "AKT1"],
            "site": ["Y182", "Y182", "T308"],
            "display_id": ["MAPK14;Y182;", "MAPK14;Y182;", "AKT1;T308;"],
            "site_sequence": ["SEQ_A", "SEQ_R", "SEQ_C"],
            "protein_id": ["PROT_A", "PROT_A", "PROT_C"],
            "uid": ["A", "B", "C"],
        },
        index=phospho.index.copy(),
    )
    scientific_row_key = pd.Series(
        [_MAPK14_KEY, _MAPK14_KEY, _AKT1_KEY],
        index=phospho.index.copy(),
        name="site_key",
    )
    constructed_display_id = pd.Series(
        ["MAPK14;Y182;", "MAPK14;Y182;", "AKT1;T308;"],
        index=phospho.index.copy(),
        name="display_id",
    )
    return phospho, site_metadata, scientific_row_key, constructed_display_id


def test_duplicate_site_resolver_returns_structured_result() -> None:
    phospho, site_metadata, scientific_row_key, constructed_display_id = (
        _duplicate_policy_inputs()
    )
    result = DuplicateSiteResolver().resolve(
        phospho=phospho,
        site_metadata=site_metadata,
        scientific_row_key=scientific_row_key,
        constructed_display_id=constructed_display_id,
        duplicate_site_policy=DATASET_SITE_MATRIX_DUPLICATE_POLICY_FIRST,
    )
    assert isinstance(result, DuplicateSiteResolutionResult)


def test_duplicate_site_policy_first_reports_retained_and_dropped_rows() -> None:
    phospho, site_metadata, scientific_row_key, constructed_display_id = (
        _duplicate_policy_inputs()
    )
    result = DuplicateSiteResolver().resolve(
        phospho=phospho,
        site_metadata=site_metadata,
        scientific_row_key=scientific_row_key,
        constructed_display_id=constructed_display_id,
        duplicate_site_policy=DATASET_SITE_MATRIX_DUPLICATE_POLICY_FIRST,
    )

    table = result.duplicate_site_resolution.sort_values("source_row_id")
    assert table.shape[0] == 2
    assert {
        "site_key",
        "display_id",
        "site_id",
        "source_row_id",
        "retained",
        "resolution_policy",
        "retained_reason",
        "dropped_reason",
        "observed_values",
        "mean_signal",
        "affected_sample_columns",
        "source_protein_id",
    }.issubset(set(table.columns))
    assert table.loc[table["source_row_id"] == "row_a", "retained"].item()
    assert not table.loc[table["source_row_id"] == "row_b", "retained"].item()
    assert "input order" in str(
        table.loc[table["source_row_id"] == "row_a", "retained_reason"].item()
    )
    assert table.loc[
        table["source_row_id"] == "row_a", "affected_sample_columns"
    ].item() == ("sample_a", "sample_b")
    assert table.loc[table["source_row_id"] == "row_a", "display_id"].item() == (
        "MAPK14;Y182;"
    )
    assert table.loc[table["source_row_id"] == "row_a", "site_key"].item() == (
        _MAPK14_KEY
    )
    assert table.loc[table["source_row_id"] == "row_a", "source_protein_id"].item() == (
        "PROT_A"
    )


def test_duplicate_site_policy_max_mean_signal_reports_mean_and_selection() -> None:
    phospho, site_metadata, scientific_row_key, constructed_display_id = (
        _duplicate_policy_inputs()
    )
    result = DuplicateSiteResolver().resolve(
        phospho=phospho,
        site_metadata=site_metadata,
        scientific_row_key=scientific_row_key,
        constructed_display_id=constructed_display_id,
        duplicate_site_policy=DATASET_SITE_MATRIX_DUPLICATE_POLICY_MAX_MEAN_SIGNAL,
    )

    table = result.duplicate_site_resolution.sort_values("source_row_id")
    retained_row_id = table.loc[table["retained"], "source_row_id"].astype(str).tolist()
    assert retained_row_id == ["row_b"]
    assert table.loc[
        table["source_row_id"] == "row_a", "mean_signal"
    ].item() == pytest.approx(1.5)
    assert table.loc[
        table["source_row_id"] == "row_b", "mean_signal"
    ].item() == pytest.approx(3.5)


def test_duplicate_site_policy_aggregate_mean_reports_all_contributing_rows() -> None:
    phospho, site_metadata, scientific_row_key, constructed_display_id = (
        _duplicate_policy_inputs()
    )
    result = DuplicateSiteResolver().resolve(
        phospho=phospho,
        site_metadata=site_metadata,
        scientific_row_key=scientific_row_key,
        constructed_display_id=constructed_display_id,
        duplicate_site_policy=DATASET_SITE_MATRIX_DUPLICATE_POLICY_AGGREGATE_MEAN,
    )

    table = result.duplicate_site_resolution.sort_values("source_row_id")
    assert table["source_row_id"].tolist() == ["row_a", "row_b"]
    assert table["retained"].tolist() == [True, True]
    assert table["n_aggregated_rows"].tolist() == [2, 2]
    assert (table["resolution_policy"] == "aggregate_mean").all()


def test_duplicate_site_policy_aggregate_median_reports_all_contributing_rows() -> None:
    phospho, site_metadata, scientific_row_key, constructed_display_id = (
        _duplicate_policy_inputs()
    )
    result = DuplicateSiteResolver().resolve(
        phospho=phospho,
        site_metadata=site_metadata,
        scientific_row_key=scientific_row_key,
        constructed_display_id=constructed_display_id,
        duplicate_site_policy=DATASET_SITE_MATRIX_DUPLICATE_POLICY_AGGREGATE_MEDIAN,
    )

    table = result.duplicate_site_resolution.sort_values("source_row_id")
    assert table["source_row_id"].tolist() == ["row_a", "row_b"]
    assert table["retained"].tolist() == [True, True]
    assert table["n_aggregated_rows"].tolist() == [2, 2]
    assert (table["resolution_policy"] == "aggregate_median").all()


def test_duplicate_site_policy_reports_metadata_conflicts() -> None:
    phospho, site_metadata, scientific_row_key, constructed_display_id = (
        _duplicate_policy_inputs()
    )
    result = DuplicateSiteResolver().resolve(
        phospho=phospho,
        site_metadata=site_metadata,
        scientific_row_key=scientific_row_key,
        constructed_display_id=constructed_display_id,
        duplicate_site_policy=DATASET_SITE_MATRIX_DUPLICATE_POLICY_AGGREGATE_MEAN,
    )

    assert not result.metadata_conflicts.empty
    uid_conflict = result.metadata_conflicts.loc[
        result.metadata_conflicts["field"] == "uid"
    ]
    assert uid_conflict.shape[0] == 1
    assert uid_conflict.iloc[0]["site_key"] == _MAPK14_KEY
    assert uid_conflict.iloc[0]["display_id"] == "MAPK14;Y182;"
    assert uid_conflict.iloc[0]["site_id"] == "MAPK14;Y182;"
    assert uid_conflict.iloc[0]["n_distinct_values"] == 2
    assert uid_conflict.iloc[0]["source_row_ids"] == ("row_a", "row_b")
    assert "site_sequence" in set(result.metadata_conflicts.loc[:, "field"])
    assert "uid" in set(result.metadata_conflicts.loc[:, "field"])
    assert pd.isna(result.site_metadata.loc[_MAPK14_KEY, "site_sequence"])
    assert pd.isna(result.site_metadata.loc[_MAPK14_KEY, "uid"])
    assert (
        result.duplicate_site_resolution.loc[:, "metadata_conflict_detected"]
        .astype(bool)
        .all()
    )


@pytest.mark.parametrize(
    ("policy", "expected_phospho", "expected_site_metadata"),
    [
        (
            DATASET_SITE_MATRIX_DUPLICATE_POLICY_FIRST,
            pd.DataFrame(
                {
                    "sample_a": [1.0, 5.0],
                    "sample_b": [2.0, 6.0],
                },
                index=pd.Index([_MAPK14_KEY, _AKT1_KEY], name="site_key"),
            ),
            pd.DataFrame(
                {
                    "gene_symbol": ["MAPK14", "AKT1"],
                    "site": ["Y182", "T308"],
                    "display_id": ["MAPK14;Y182;", "AKT1;T308;"],
                    "site_sequence": ["SEQ_A", "SEQ_C"],
                    "protein_id": ["PROT_A", "PROT_C"],
                    "uid": ["A", "C"],
                },
                index=pd.Index([_MAPK14_KEY, _AKT1_KEY], name="site_key"),
            ),
        ),
        (
            DATASET_SITE_MATRIX_DUPLICATE_POLICY_MAX_MEAN_SIGNAL,
            pd.DataFrame(
                {
                    "sample_a": [3.0, 5.0],
                    "sample_b": [4.0, 6.0],
                },
                index=pd.Index([_MAPK14_KEY, _AKT1_KEY], name="site_key"),
            ),
            pd.DataFrame(
                {
                    "gene_symbol": ["MAPK14", "AKT1"],
                    "site": ["Y182", "T308"],
                    "display_id": ["MAPK14;Y182;", "AKT1;T308;"],
                    "site_sequence": ["SEQ_R", "SEQ_C"],
                    "protein_id": ["PROT_A", "PROT_C"],
                    "uid": ["B", "C"],
                },
                index=pd.Index([_MAPK14_KEY, _AKT1_KEY], name="site_key"),
            ),
        ),
        (
            DATASET_SITE_MATRIX_DUPLICATE_POLICY_AGGREGATE_MEAN,
            pd.DataFrame(
                {
                    "sample_a": [2.0, 5.0],
                    "sample_b": [3.0, 6.0],
                },
                index=pd.Index([_MAPK14_KEY, _AKT1_KEY], name="site_key"),
            ),
            pd.DataFrame(
                {
                    "gene_symbol": ["MAPK14", "AKT1"],
                    "site": ["Y182", "T308"],
                    "display_id": ["MAPK14;Y182;", "AKT1;T308;"],
                    "site_sequence": [pd.NA, "SEQ_C"],
                    "protein_id": ["PROT_A", "PROT_C"],
                    "uid": [pd.NA, "C"],
                },
                index=pd.Index([_MAPK14_KEY, _AKT1_KEY], name="site_key"),
            ),
        ),
        (
            DATASET_SITE_MATRIX_DUPLICATE_POLICY_AGGREGATE_MEDIAN,
            pd.DataFrame(
                {
                    "sample_a": [2.0, 5.0],
                    "sample_b": [3.0, 6.0],
                },
                index=pd.Index([_MAPK14_KEY, _AKT1_KEY], name="site_key"),
            ),
            pd.DataFrame(
                {
                    "gene_symbol": ["MAPK14", "AKT1"],
                    "site": ["Y182", "T308"],
                    "display_id": ["MAPK14;Y182;", "AKT1;T308;"],
                    "site_sequence": [pd.NA, "SEQ_C"],
                    "protein_id": ["PROT_A", "PROT_C"],
                    "uid": [pd.NA, "C"],
                },
                index=pd.Index([_MAPK14_KEY, _AKT1_KEY], name="site_key"),
            ),
        ),
    ],
)
def test_duplicate_site_policy_preserves_existing_outputs_for_current_policies(
    policy: str,
    expected_phospho: pd.DataFrame,
    expected_site_metadata: pd.DataFrame,
) -> None:
    phospho, site_metadata, scientific_row_key, constructed_display_id = (
        _duplicate_policy_inputs()
    )
    result = DuplicateSiteResolver().resolve(
        phospho=phospho,
        site_metadata=site_metadata,
        scientific_row_key=scientific_row_key,
        constructed_display_id=constructed_display_id,
        duplicate_site_policy=policy,
    )

    pdt.assert_frame_equal(result.phospho, expected_phospho)
    pdt.assert_frame_equal(result.site_metadata, expected_site_metadata)


@pytest.mark.parametrize(
    "policy",
    [
        DATASET_SITE_MATRIX_DUPLICATE_POLICY_AGGREGATE_MEAN,
        DATASET_SITE_MATRIX_DUPLICATE_POLICY_AGGREGATE_MEDIAN,
    ],
)
def test_duplicate_site_aggregate_complete_observations_preserve_values(
    policy: str,
) -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [2.0, 4.0],
            "sample_b": [6.0, 8.0],
        },
        index=pd.Index(["row_a", "row_b"], name="source_row"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "MAPK14"],
            "site": ["Y182", "Y182"],
            "site_sequence": ["SEQ_A", "SEQ_A"],
        },
        index=phospho.index.copy(),
    )
    scientific_row_key = pd.Series(
        [_MAPK14_KEY, _MAPK14_KEY],
        index=phospho.index.copy(),
        name="site_key",
    )
    constructed_display_id = pd.Series(
        ["MAPK14;Y182;", "MAPK14;Y182;"],
        index=phospho.index.copy(),
        name="display_id",
    )

    result = DuplicateSiteResolver().resolve(
        phospho=phospho,
        site_metadata=site_metadata,
        scientific_row_key=scientific_row_key,
        constructed_display_id=constructed_display_id,
        duplicate_site_policy=policy,
    )

    assert result.phospho.shape == (1, 2)
    assert result.phospho.loc[_MAPK14_KEY, "sample_a"] == pytest.approx(3.0)
    assert result.phospho.loc[_MAPK14_KEY, "sample_b"] == pytest.approx(7.0)
    assert result.duplicate_aggregation_diagnostics["missing_value_policy"] == (
        "skip_missing_values"
    )
    assert (
        result.duplicate_aggregation_diagnostics["missing_cells_before_aggregation"]
        == 0
    )
    assert (
        result.duplicate_aggregation_diagnostics["missing_cells_after_aggregation"] == 0
    )


@pytest.mark.parametrize(
    ("policy", "expected_sample_a", "expected_sample_b"),
    [
        (DATASET_SITE_MATRIX_DUPLICATE_POLICY_AGGREGATE_MEAN, 1.0, 2.0),
        (DATASET_SITE_MATRIX_DUPLICATE_POLICY_AGGREGATE_MEDIAN, 1.0, 2.0),
    ],
)
def test_duplicate_site_aggregate_partial_missing_uses_explicit_skip_missing_policy(
    policy: str,
    expected_sample_a: float,
    expected_sample_b: float,
) -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, float("nan"), 5.0],
            "sample_b": [float("nan"), 2.0, 7.0],
        },
        index=pd.Index(["row_a", "row_b", "row_c"], name="source_row"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "MAPK14", "AKT1"],
            "site": ["Y182", "Y182", "T308"],
            "site_sequence": ["SEQ_A", "SEQ_A", "SEQ_C"],
        },
        index=phospho.index.copy(),
    )
    scientific_row_key = pd.Series(
        [_MAPK14_KEY, _MAPK14_KEY, _AKT1_KEY],
        index=phospho.index.copy(),
        name="site_key",
    )
    constructed_display_id = pd.Series(
        ["MAPK14;Y182;", "MAPK14;Y182;", "AKT1;T308;"],
        index=phospho.index.copy(),
        name="display_id",
    )

    result = DuplicateSiteResolver().resolve(
        phospho=phospho,
        site_metadata=site_metadata,
        scientific_row_key=scientific_row_key,
        constructed_display_id=constructed_display_id,
        duplicate_site_policy=policy,
    )

    assert result.phospho.loc[_MAPK14_KEY, "sample_a"] == pytest.approx(
        expected_sample_a
    )
    assert result.phospho.loc[_MAPK14_KEY, "sample_b"] == pytest.approx(
        expected_sample_b
    )
    assert result.duplicate_aggregation_diagnostics["missing_value_policy"] == (
        "skip_missing_values"
    )
    assert (
        result.duplicate_aggregation_diagnostics["aggregation_reduced_missingness"]
        is True
    )


@pytest.mark.parametrize(
    "policy",
    [
        DATASET_SITE_MATRIX_DUPLICATE_POLICY_AGGREGATE_MEAN,
        DATASET_SITE_MATRIX_DUPLICATE_POLICY_AGGREGATE_MEDIAN,
    ],
)
def test_duplicate_site_aggregate_one_source_row_all_missing_is_not_dropped(
    policy: str,
) -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [float("nan"), 4.0],
            "sample_b": [float("nan"), 6.0],
        },
        index=pd.Index(["row_all_missing", "row_observed"], name="source_row"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "MAPK14"],
            "site": ["Y182", "Y182"],
            "site_sequence": ["SEQ_A", "SEQ_A"],
        },
        index=phospho.index.copy(),
    )
    scientific_row_key = pd.Series(
        [_MAPK14_KEY, _MAPK14_KEY],
        index=phospho.index.copy(),
        name="site_key",
    )
    constructed_display_id = pd.Series(
        ["MAPK14;Y182;", "MAPK14;Y182;"],
        index=phospho.index.copy(),
        name="display_id",
    )

    result = DuplicateSiteResolver().resolve(
        phospho=phospho,
        site_metadata=site_metadata,
        scientific_row_key=scientific_row_key,
        constructed_display_id=constructed_display_id,
        duplicate_site_policy=policy,
    )

    assert result.phospho.shape[0] == 1
    assert result.phospho.loc[_MAPK14_KEY, "sample_a"] == pytest.approx(4.0)
    assert result.phospho.loc[_MAPK14_KEY, "sample_b"] == pytest.approx(6.0)
    assert result.duplicate_aggregation_diagnostics["rows_collapsed_count"] == 1


def test_duplicate_site_aggregate_all_missing_duplicates_keep_group_row() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [float("nan"), float("nan")],
            "sample_b": [float("nan"), float("nan")],
        },
        index=pd.Index(["row_a", "row_b"], name="source_row"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "MAPK14"],
            "site": ["Y182", "Y182"],
            "site_sequence": ["SEQ_A", "SEQ_A"],
        },
        index=phospho.index.copy(),
    )
    scientific_row_key = pd.Series(
        [_MAPK14_KEY, _MAPK14_KEY],
        index=phospho.index.copy(),
        name="site_key",
    )
    constructed_display_id = pd.Series(
        ["MAPK14;Y182;", "MAPK14;Y182;"],
        index=phospho.index.copy(),
        name="display_id",
    )

    result = DuplicateSiteResolver().resolve(
        phospho=phospho,
        site_metadata=site_metadata,
        scientific_row_key=scientific_row_key,
        constructed_display_id=constructed_display_id,
        duplicate_site_policy=DATASET_SITE_MATRIX_DUPLICATE_POLICY_AGGREGATE_MEAN,
    )

    assert result.phospho.shape == (1, 2)
    assert result.phospho.index.tolist() == [_MAPK14_KEY]
    assert result.duplicate_aggregation_diagnostics["duplicate_group_count"] == 1


def test_duplicate_display_ids_with_distinct_site_keys_are_not_treated_as_duplicates() -> (
    None
):
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 2.0],
            "sample_b": [3.0, 4.0],
        },
        index=pd.Index(["row_a", "row_b"], name="source_row"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "MAPK14"],
            "site": ["Y182", "Y182"],
            "display_id": ["MAPK14;Y182;", "MAPK14;Y182;"],
            "site_sequence": ["SEQ_A", "SEQ_B"],
            "protein_id": ["P28482", "Q9Y2R5"],
        },
        index=phospho.index.copy(),
    )
    scientific_row_key = pd.Series(
        [
            protein_site_key(protein_identifier="P28482", site="Y182"),
            protein_site_key(protein_identifier="Q9Y2R5", site="Y182"),
        ],
        index=phospho.index.copy(),
        name="site_key",
    )
    constructed_display_id = pd.Series(
        ["MAPK14;Y182;", "MAPK14;Y182;"],
        index=phospho.index.copy(),
        name="display_id",
    )

    result = DuplicateSiteResolver().resolve(
        phospho=phospho,
        site_metadata=site_metadata,
        scientific_row_key=scientific_row_key,
        constructed_display_id=constructed_display_id,
        duplicate_site_policy="error",
    )

    assert result.phospho.index.tolist() == scientific_row_key.tolist()
    assert result.duplicate_site_resolution.empty


def test_duplicate_site_keys_fail_under_strict_policy() -> None:
    phospho, site_metadata, scientific_row_key, constructed_display_id = (
        _duplicate_policy_inputs()
    )

    with pytest.raises(PhosPyInputError, match="duplicate site_key values"):
        DuplicateSiteResolver().resolve(
            phospho=phospho,
            site_metadata=site_metadata,
            scientific_row_key=scientific_row_key,
            constructed_display_id=constructed_display_id,
            duplicate_site_policy=DATASET_SITE_MATRIX_DUPLICATE_POLICY_ERROR,
        )


def test_duplicate_resolution_rejects_missing_site_key() -> None:
    phospho, site_metadata, _, constructed_display_id = _duplicate_policy_inputs()

    with pytest.raises(PhosPyInputError, match="requires site_key row identity"):
        DuplicateSiteResolver().resolve(
            phospho=phospho,
            site_metadata=site_metadata,
            constructed_display_id=constructed_display_id,
            duplicate_site_policy=DATASET_SITE_MATRIX_DUPLICATE_POLICY_FIRST,
        )


def test_duplicate_resolution_rejects_raw_invalid_site_key_strings() -> None:
    phospho, site_metadata, _, constructed_display_id = _duplicate_policy_inputs()
    invalid_site_key = pd.Series(
        ["site_key_1", "site_key_1", "site_key_2"],
        index=phospho.index.copy(),
        name="site_key",
    )

    with pytest.raises(PhosPyInputError, match="must start with 'phospy:v1'"):
        DuplicateSiteResolver().resolve(
            phospho=phospho,
            site_metadata=site_metadata,
            scientific_row_key=invalid_site_key,
            constructed_display_id=constructed_display_id,
            duplicate_site_policy=DATASET_SITE_MATRIX_DUPLICATE_POLICY_FIRST,
        )


def test_missing_and_specified_isoform_scope_collision_fails() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 2.0],
            "sample_b": [3.0, 4.0],
        },
        index=pd.Index(["row_missing_isoform", "row_explicit_isoform"], name="source"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "MAPK14"],
            "site": ["Y182", "Y182"],
            "display_id": ["MAPK14;Y182;", "MAPK14;Y182;"],
            "site_sequence": ["SEQ_A", "SEQ_B"],
            "protein_id": ["P28482", "P28482"],
            "isoform_id": [pd.NA, "isoform-2"],
        },
        index=phospho.index.copy(),
    )
    site_keys = pd.Series(
        [
            encode_site_key(
                ProteinScopedPhosphositeKey(
                    organism="rat",
                    protein_namespace="protein_id",
                    protein_identifier="P28482",
                    residue="Y",
                    position=182,
                )
            ),
            encode_site_key(
                ProteinScopedPhosphositeKey(
                    organism="rat",
                    protein_namespace="protein_id",
                    protein_identifier="P28482",
                    residue="Y",
                    position=182,
                    isoform_id="isoform-2",
                )
            ),
        ],
        index=phospho.index.copy(),
        name="site_key",
    )
    display_ids = pd.Series(
        ["MAPK14;Y182;", "MAPK14;Y182;"],
        index=phospho.index.copy(),
        name="display_id",
    )

    with pytest.raises(PhosPyInputError, match="mixed isoform scope"):
        DuplicateSiteResolver().resolve(
            phospho=phospho,
            site_metadata=site_metadata,
            scientific_row_key=site_keys,
            constructed_display_id=display_ids,
            duplicate_site_policy=DATASET_SITE_MATRIX_DUPLICATE_POLICY_FIRST,
        )
