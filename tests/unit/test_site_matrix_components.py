from __future__ import annotations

import pandas as pd
import pandas.testing as pdt
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from phospy.api.configs import (
    DATASET_SITE_MATRIX_DUPLICATE_POLICY_AGGREGATE_MEAN,
    DATASET_SITE_MATRIX_DUPLICATE_POLICY_AGGREGATE_MEDIAN,
    DATASET_SITE_MATRIX_DUPLICATE_POLICY_ERROR,
    DATASET_SITE_MATRIX_DUPLICATE_POLICY_FIRST,
    DATASET_SITE_MATRIX_DUPLICATE_POLICY_MAX_MEAN_SIGNAL,
    DATASET_SITE_MATRIX_MISSING_DATA_POLICY_DROP_ANY_MISSING,
    DATASET_SITE_MATRIX_POLICY_BUILD_FROM_METADATA,
)
from phospy.errors.input import PhosPyInputError
from phospy.science.datasets.preprocessing.models import (
    DuplicateSiteResolutionResult,
    PreprocessingPlan,
    PreprocessingState,
)
from phospy.science.datasets.preprocessing.stages.site_matrix import SiteMatrixStage
from phospy.science.datasets.preprocessing.stages.site_matrix_components import (
    DuplicateSiteResolver,
    MetadataConflictDetector,
    MissingDataSiteFilter,
    SequenceSupportFilter,
    SiteMatrixAssembler,
    SiteMatrixProvenanceBuilder,
    SiteMatrixRowAuditBuilder,
)

_PROPERTY_SETTINGS = settings(
    max_examples=30,
    deadline=None,
    derandomize=True,
)


def _duplicate_policy_inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.Series]:
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
            "site_sequence": ["SEQ_A", "SEQ_B", "SEQ_C"],
            "protein_id": ["PROT_A", "PROT_A", "PROT_C"],
            "uid": ["A", "B", "C"],
        },
        index=phospho.index.copy(),
    )
    constructed_site_id = pd.Series(
        ["MAPK14;Y182;", "MAPK14;Y182;", "AKT1;T308;"],
        index=phospho.index.copy(),
        name="site_id",
    )
    return phospho, site_metadata, constructed_site_id


def _constructed_site_ids_strategy(
    *,
    min_size: int = 1,
    max_size: int = 8,
) -> st.SearchStrategy[list[str]]:
    return st.lists(
        st.integers(min_value=1, max_value=4),
        min_size=min_size,
        max_size=max_size,
    ).map(lambda ids: [f"G{site_id};S{site_id};" for site_id in ids])


def test_sequence_support_filter_preserves_supported_and_dropped_rows() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 2.0, 3.0, 4.0],
            "sample_b": [5.0, 6.0, 7.0, 8.0],
        },
        index=pd.Index(["row_a", "row_b", "row_c", "row_d"], name="source_row"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["A", "B", "C", "D"],
            "site": ["S1", "S2", "S3", "S4"],
            "site_sequence": ["SEQ_A", "   ", pd.NA, "SEQ_D"],
        },
        index=phospho.index.copy(),
    )
    constructed_site_id = pd.Series(
        ["A;S1;", "B;S2;", "C;S3;", "D;S4;"],
        index=phospho.index.copy(),
        name="site_id",
    )

    result = SequenceSupportFilter().filter(
        phospho=phospho,
        site_metadata=site_metadata,
        constructed_site_id=constructed_site_id,
    )

    assert result.phospho.index.tolist() == ["row_a", "row_d"]
    assert result.site_metadata.index.tolist() == ["row_a", "row_d"]
    assert result.constructed_site_id.tolist() == ["A;S1;", "D;S4;"]
    assert result.dropped_row_count == 2
    assert result.dropped_rows == (("row_b", "B;S2;"), ("row_c", "C;S3;"))


@pytest.mark.parametrize(
    (
        "policy",
        "minimum_observed_values",
        "expected_rows",
        "expected_required",
        "dropped",
    ),
    [
        (
            DATASET_SITE_MATRIX_MISSING_DATA_POLICY_DROP_ANY_MISSING,
            None,
            ["row_a"],
            2,
            (("row_b", "B;S2;", 1), ("row_c", "C;S3;", 0)),
        ),
        (
            "retain_missing",
            None,
            ["row_a", "row_b", "row_c"],
            0,
            (),
        ),
        (
            "require_min_observed_values",
            1,
            ["row_a", "row_b"],
            1,
            (("row_c", "C;S3;", 0),),
        ),
    ],
)
def test_missing_data_site_filter_preserves_policy_specific_row_selection(
    policy: str,
    minimum_observed_values: int | None,
    expected_rows: list[str],
    expected_required: int,
    dropped: tuple[tuple[str, str, int], ...],
) -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 2.0, float("nan")],
            "sample_b": [3.0, float("nan"), float("nan")],
        },
        index=pd.Index(["row_a", "row_b", "row_c"], name="source_row"),
    )
    constructed_site_id = pd.Series(
        ["A;S1;", "B;S2;", "C;S3;"],
        index=phospho.index.copy(),
        name="site_id",
    )

    result = MissingDataSiteFilter().filter(
        phospho=phospho,
        constructed_site_id=constructed_site_id,
        missing_data_policy=policy,
        minimum_observed_values=minimum_observed_values,
    )

    assert result.phospho.index.tolist() == expected_rows
    assert result.required_observed_count == expected_required
    assert result.dropped_rows == dropped
    assert result.dropped_row_count == len(dropped)


def test_site_matrix_assembler_preserves_index_order_and_dropped_row_ids() -> None:
    duplicate_site_result = DuplicateSiteResolutionResult(
        phospho=pd.DataFrame(
            {
                "sample_a": [8.0, 2.0],
                "sample_b": [8.5, 3.0],
            },
            index=pd.Index(["MAPK14;Y182;", "AKT1;T308;"], name="site_id"),
        ),
        site_metadata=pd.DataFrame(
            {
                "gene_symbol": ["MAPK14", "AKT1"],
                "site": ["Y182", "T308"],
                "site_sequence": ["SEQ_B", "SEQ_A"],
            },
            index=pd.Index(["MAPK14;Y182;", "AKT1;T308;"], name="site_id"),
        ),
        dropped_row_count=1,
        duplicate_site_resolution=pd.DataFrame(
            {
                "site_id": ["MAPK14;Y182;", "MAPK14;Y182;"],
                "source_row_id": ["row_a", "row_b"],
                "retained": [True, False],
            }
        ),
        metadata_conflicts=pd.DataFrame(),
    )

    result = SiteMatrixAssembler().assemble(
        duplicate_site_result=duplicate_site_result,
        output_index_name="input_row",
        dropped_missing_sequence_rows=(("row_missing", "AKT1;S9;"),),
        dropped_incomplete_rows=(("row_incomplete", "AKT1;T308;", 1),),
    )

    assert result.phospho.index.tolist() == ["AKT1;T308;", "MAPK14;Y182;"]
    assert result.phospho.index.name == "input_row"
    assert result.site_metadata.index.tolist() == ["AKT1;T308;", "MAPK14;Y182;"]
    assert result.site_metadata.index.name == "input_row"
    assert result.dropped_missing_sequence_row_ids == ("row_missing",)
    assert result.dropped_incomplete_row_ids == ("row_incomplete",)
    assert result.duplicate_dropped_row_ids == ("row_b",)
    assert result.dropped_row_ids == ("row_missing", "row_incomplete", "row_b")


def test_site_matrix_provenance_builder_preserves_fields_and_diagnostics() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [2.0, 8.0],
            "sample_b": [3.0, 8.5],
        },
        index=pd.Index(["AKT1;T308;", "MAPK14;Y182;"], name="input_row"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["AKT1", "MAPK14"],
            "site": ["T308", "Y182"],
            "site_sequence": ["SEQ_A", "SEQ_B"],
        },
        index=phospho.index.copy(),
    )
    duplicate_resolution = pd.DataFrame(
        {
            "site_id": ["MAPK14;Y182;"],
            "source_row_id": ["row_b"],
            "retained": [False],
            "source_rows": [("row_a", "row_b")],
            "dropped_reason": [pd.NA],
        }
    )

    result = SiteMatrixProvenanceBuilder().build(
        phospho=phospho,
        site_metadata=site_metadata,
        input_rows=4,
        dropped_missing_sequence=1,
        dropped_incomplete_values=1,
        missing_data_policy=DATASET_SITE_MATRIX_MISSING_DATA_POLICY_DROP_ANY_MISSING,
        required_observed_count=2,
        deduplicated_site_rows=1,
        duplicate_site_policy=DATASET_SITE_MATRIX_DUPLICATE_POLICY_FIRST,
        site_matrix_policy=DATASET_SITE_MATRIX_POLICY_BUILD_FROM_METADATA,
        dropped_missing_sequence_row_ids=("row_missing",),
        dropped_incomplete_row_ids=("row_incomplete",),
        dropped_row_ids=("row_missing", "row_incomplete", "row_b"),
        duplicate_site_resolution=duplicate_resolution,
        duplicate_aggregation_diagnostics={
            "aggregation_method": "first",
            "missing_value_policy": "not_applicable_row_selection",
            "duplicate_group_count": 1,
            "rows_collapsed_count": 1,
            "missing_cells_before_aggregation": 0,
            "missing_cells_after_aggregation": 0,
            "aggregation_reduced_missingness": False,
            "metadata_resolution_policy": "retain_earliest_input_row_per_site",
        },
    )

    assert result.row_drop_stats == {
        "input_rows": 4,
        "dropped_missing_sequence": 1,
        "dropped_incomplete_values": 1,
        "missing_data_policy": "drop_any_missing",
        "required_observed_count": 2,
        "deduplicated_site_rows": 1,
        "duplicate_site_policy": "first",
        "retained_rows": 2,
    }
    assert result.site_matrix_provenance == {
        "dropped_missing_sequence_row_ids": ("row_missing",),
        "dropped_incomplete_row_ids": ("row_incomplete",),
        "dropped_row_ids": ("row_missing", "row_incomplete", "row_b"),
        "duplicate_site_policy": "first",
        "missing_data_policy": "drop_any_missing",
        "required_observed_count": 2,
        "final_constructed_site_ids": ("AKT1;T308;", "MAPK14;Y182;"),
        "duplicate_aggregation": {
            "aggregation_method": "first",
            "missing_value_policy": "not_applicable_row_selection",
            "duplicate_group_count": 1,
            "rows_collapsed_count": 1,
            "missing_cells_before_aggregation": 0,
            "missing_cells_after_aggregation": 0,
            "aggregation_reduced_missingness": False,
            "metadata_resolution_policy": "retain_earliest_input_row_per_site",
        },
    }
    assert result.diagnostics["final_constructed_site_ids"] == [
        "AKT1;T308;",
        "MAPK14;Y182;",
    ]
    duplicate_decisions = result.diagnostics["duplicate_site_decisions"]
    assert duplicate_decisions[0]["source_rows"] == ["row_a", "row_b"]
    assert duplicate_decisions[0]["dropped_reason"] is None
    assert result.diagnostics["duplicate_aggregation"]["aggregation_method"] == "first"
    assert (
        result.diagnostics["duplicate_aggregation"]["missing_value_policy"]
        == "not_applicable_row_selection"
    )
    assert result.phospho.attrs["site_matrix_policy"] == "build_from_metadata"
    assert result.site_metadata.attrs["site_matrix_policy"] == "build_from_metadata"


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
                index=pd.Index(["MAPK14;Y182;", "AKT1;T308;"], name="site_id"),
            ),
            pd.DataFrame(
                {
                    "gene_symbol": ["MAPK14", "AKT1"],
                    "site": ["Y182", "T308"],
                    "site_sequence": ["SEQ_A", "SEQ_C"],
                    "protein_id": ["PROT_A", "PROT_C"],
                    "uid": ["A", "C"],
                },
                index=pd.Index(["MAPK14;Y182;", "AKT1;T308;"], name="site_id"),
            ),
        ),
        (
            DATASET_SITE_MATRIX_DUPLICATE_POLICY_MAX_MEAN_SIGNAL,
            pd.DataFrame(
                {
                    "sample_a": [5.0, 3.0],
                    "sample_b": [6.0, 4.0],
                },
                index=pd.Index(["AKT1;T308;", "MAPK14;Y182;"], name="site_id"),
            ),
            pd.DataFrame(
                {
                    "gene_symbol": ["AKT1", "MAPK14"],
                    "site": ["T308", "Y182"],
                    "site_sequence": ["SEQ_C", "SEQ_B"],
                    "protein_id": ["PROT_C", "PROT_A"],
                    "uid": ["C", "B"],
                },
                index=pd.Index(["AKT1;T308;", "MAPK14;Y182;"], name="site_id"),
            ),
        ),
        (
            DATASET_SITE_MATRIX_DUPLICATE_POLICY_AGGREGATE_MEAN,
            pd.DataFrame(
                {
                    "sample_a": [2.0, 5.0],
                    "sample_b": [3.0, 6.0],
                },
                index=pd.Index(["MAPK14;Y182;", "AKT1;T308;"], name="site_id"),
            ),
            pd.DataFrame(
                {
                    "gene_symbol": ["MAPK14", "AKT1"],
                    "site": ["Y182", "T308"],
                    "site_sequence": [pd.NA, "SEQ_C"],
                    "protein_id": ["PROT_A", "PROT_C"],
                    "uid": [pd.NA, "C"],
                },
                index=pd.Index(["MAPK14;Y182;", "AKT1;T308;"], name="site_id"),
            ),
        ),
        (
            DATASET_SITE_MATRIX_DUPLICATE_POLICY_AGGREGATE_MEDIAN,
            pd.DataFrame(
                {
                    "sample_a": [2.0, 5.0],
                    "sample_b": [3.0, 6.0],
                },
                index=pd.Index(["MAPK14;Y182;", "AKT1;T308;"], name="site_id"),
            ),
            pd.DataFrame(
                {
                    "gene_symbol": ["MAPK14", "AKT1"],
                    "site": ["Y182", "T308"],
                    "site_sequence": [pd.NA, "SEQ_C"],
                    "protein_id": ["PROT_A", "PROT_C"],
                    "uid": [pd.NA, "C"],
                },
                index=pd.Index(["MAPK14;Y182;", "AKT1;T308;"], name="site_id"),
            ),
        ),
    ],
)
def test_duplicate_site_resolver_preserves_outputs_for_supported_policies(
    policy: str,
    expected_phospho: pd.DataFrame,
    expected_site_metadata: pd.DataFrame,
) -> None:
    phospho, site_metadata, constructed_site_id = _duplicate_policy_inputs()
    result = DuplicateSiteResolver().resolve(
        phospho=phospho,
        site_metadata=site_metadata,
        constructed_site_id=constructed_site_id,
        duplicate_site_policy=policy,
    )

    pdt.assert_frame_equal(result.phospho, expected_phospho)
    pdt.assert_frame_equal(result.site_metadata, expected_site_metadata)


@given(
    constructed_site_ids=_constructed_site_ids_strategy(min_size=1, max_size=8),
    duplicate_site_policy=st.sampled_from(
        (
            DATASET_SITE_MATRIX_DUPLICATE_POLICY_FIRST,
            DATASET_SITE_MATRIX_DUPLICATE_POLICY_MAX_MEAN_SIGNAL,
            DATASET_SITE_MATRIX_DUPLICATE_POLICY_AGGREGATE_MEAN,
            DATASET_SITE_MATRIX_DUPLICATE_POLICY_AGGREGATE_MEDIAN,
        )
    ),
)
@_PROPERTY_SETTINGS
def test_duplicate_site_resolver_property_outputs_unique_site_ids(
    constructed_site_ids: list[str],
    duplicate_site_policy: str,
) -> None:
    site_protein_ids: dict[str, str] = {}
    for site_id in constructed_site_ids:
        if site_id in site_protein_ids:
            continue
        site_protein_ids[site_id] = (
            f"PROT_{site_id.split(';')[0]}_{site_id.split(';')[1]}"
        )

    row_ids = [f"row_{idx}" for idx in range(len(constructed_site_ids))]
    phospho = pd.DataFrame(
        {
            "sample_a": [float(idx + 1) for idx in range(len(constructed_site_ids))],
            "sample_b": [float(idx + 11) for idx in range(len(constructed_site_ids))],
        },
        index=pd.Index(row_ids, name="source_row"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": [site_id.split(";")[0] for site_id in constructed_site_ids],
            "site": [site_id.split(";")[1] for site_id in constructed_site_ids],
            "site_sequence": [
                ("A" * 15) + str(site).strip().upper()[0] + ("A" * 15)
                for site in [site_id.split(";")[1] for site_id in constructed_site_ids]
            ],
            "protein_id": [
                site_protein_ids[site_id] for site_id in constructed_site_ids
            ],
        },
        index=phospho.index.copy(),
    )
    constructed_site_id = pd.Series(
        constructed_site_ids,
        index=phospho.index.copy(),
        name="site_id",
    )

    result = DuplicateSiteResolver().resolve(
        phospho=phospho,
        site_metadata=site_metadata,
        constructed_site_id=constructed_site_id,
        duplicate_site_policy=duplicate_site_policy,
    )

    assert result.phospho.index.is_unique
    assert result.site_metadata.index.is_unique
    assert len(result.phospho.index) == len(set(constructed_site_ids))


def test_duplicate_site_resolver_reports_rows_and_metadata_conflicts() -> None:
    phospho, site_metadata, constructed_site_id = _duplicate_policy_inputs()
    result = DuplicateSiteResolver().resolve(
        phospho=phospho,
        site_metadata=site_metadata,
        constructed_site_id=constructed_site_id,
        duplicate_site_policy=DATASET_SITE_MATRIX_DUPLICATE_POLICY_FIRST,
    )

    resolution = result.duplicate_site_resolution.sort_values("source_row_id")
    assert resolution["source_row_id"].tolist() == ["row_a", "row_b"]
    assert resolution["retained"].tolist() == [True, False]
    assert "input order" in str(
        resolution.loc[resolution["source_row_id"] == "row_a", "retained_reason"].item()
    )
    assert "selected first" in str(
        resolution.loc[resolution["source_row_id"] == "row_b", "dropped_reason"].item()
    )
    assert bool(
        resolution.loc[
            resolution["source_row_id"] == "row_a", "metadata_conflict_detected"
        ].item()
    )

    conflicts = result.metadata_conflicts
    uid_conflicts = conflicts.loc[conflicts["field"] == "uid"]
    assert uid_conflicts.shape[0] == 1
    assert uid_conflicts.iloc[0]["site_id"] == "MAPK14;Y182;"
    assert uid_conflicts.iloc[0]["source_row_ids"] == ("row_a", "row_b")


def test_duplicate_site_resolver_aggregate_preserves_identical_metadata() -> None:
    phospho = pd.DataFrame(
        {"sample_a": [1.0, 3.0], "sample_b": [2.0, 4.0]},
        index=pd.Index(["row_a", "row_b"], name="source_row"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "MAPK14"],
            "site": ["Y182", "Y182"],
            "site_sequence": ["SEQ_A", "SEQ_A"],
            "protein_id": ["PROT_A", "PROT_A"],
            "uid": ["U1", "U1"],
        },
        index=phospho.index.copy(),
    )
    constructed_site_id = pd.Series(
        ["MAPK14;Y182;", "MAPK14;Y182;"],
        index=phospho.index.copy(),
        name="site_id",
    )

    result = DuplicateSiteResolver().resolve(
        phospho=phospho,
        site_metadata=site_metadata,
        constructed_site_id=constructed_site_id,
        duplicate_site_policy=DATASET_SITE_MATRIX_DUPLICATE_POLICY_AGGREGATE_MEAN,
    )

    assert result.site_metadata.loc["MAPK14;Y182;", "site_sequence"] == "SEQ_A"
    assert result.site_metadata.loc["MAPK14;Y182;", "protein_id"] == "PROT_A"
    assert result.site_metadata.loc["MAPK14;Y182;", "uid"] == "U1"
    assert result.metadata_conflicts.empty


def test_duplicate_site_resolver_error_policy_raises_for_duplicates() -> None:
    phospho, site_metadata, constructed_site_id = _duplicate_policy_inputs()
    with pytest.raises(
        PhosPyInputError,
        match="duplicate constructed site identifiers",
    ):
        DuplicateSiteResolver().resolve(
            phospho=phospho,
            site_metadata=site_metadata,
            constructed_site_id=constructed_site_id,
            duplicate_site_policy=DATASET_SITE_MATRIX_DUPLICATE_POLICY_ERROR,
        )


def test_duplicate_site_resolver_rejects_conflicting_protein_identity_collisions() -> (
    None
):
    phospho, site_metadata, constructed_site_id = _duplicate_policy_inputs()
    site_metadata.loc["row_b", "protein_accession"] = "P28482-2"
    site_metadata.loc["row_a", "protein_accession"] = "P28482-1"

    with pytest.raises(
        PhosPyInputError,
        match="conflicting scientific identities for duplicate display site IDs",
    ):
        DuplicateSiteResolver().resolve(
            phospho=phospho,
            site_metadata=site_metadata,
            constructed_site_id=constructed_site_id,
            duplicate_site_policy=DATASET_SITE_MATRIX_DUPLICATE_POLICY_FIRST,
        )


def test_metadata_conflict_detector_detects_distinct_values_only() -> None:
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "MAPK14", "MAPK14"],
            "site": ["Y182", "Y182", "Y182"],
            "site_sequence": ["SEQ_A", "SEQ_A", "SEQ_A"],
            "protein_id": ["P_A", "P_B", "P_B"],
        },
        index=pd.Index(["row_a", "row_b", "row_c"], name="source_row"),
    )
    constructed_site_id = pd.Series(
        ["MAPK14;Y182;", "MAPK14;Y182;", "MAPK14;Y182;"],
        index=site_metadata.index,
        name="site_id",
    )

    conflicts = MetadataConflictDetector().detect(
        site_metadata=site_metadata,
        constructed_site_id=constructed_site_id,
    )

    assert conflicts.shape[0] == 1
    assert conflicts.iloc[0]["site_id"] == "MAPK14;Y182;"
    assert conflicts.iloc[0]["field"] == "protein_id"
    assert conflicts.iloc[0]["values"] == ("P_A", "P_B")
    assert conflicts.iloc[0]["n_distinct_values"] == 2
    assert conflicts.iloc[0]["source_row_ids"] == ("row_a", "row_b", "row_c")


def test_site_matrix_row_audit_builder_tracks_dropped_and_duplicate_actions() -> None:
    phospho, site_metadata, constructed_site_id = _duplicate_policy_inputs()
    duplicate_result = DuplicateSiteResolver().resolve(
        phospho=phospho,
        site_metadata=site_metadata,
        constructed_site_id=constructed_site_id,
        duplicate_site_policy=DATASET_SITE_MATRIX_DUPLICATE_POLICY_FIRST,
    )
    records = SiteMatrixRowAuditBuilder().build(
        dropped_missing_sequence_rows=(("row_missing_sequence", "AKT1;S9;"),),
        dropped_incomplete_rows=(("row_incomplete", "AKT1;T308;", 1),),
        duplicate_site_resolution=duplicate_result.duplicate_site_resolution,
        site_matrix_policy=DATASET_SITE_MATRIX_POLICY_BUILD_FROM_METADATA,
        site_matrix_missing_data_policy=DATASET_SITE_MATRIX_MISSING_DATA_POLICY_DROP_ANY_MISSING,
        site_matrix_duplicate_site_policy=DATASET_SITE_MATRIX_DUPLICATE_POLICY_FIRST,
        required_observed_count=2,
    )

    by_source = {row.source_row_id: row for row in records}
    assert by_source["row_missing_sequence"].action == "dropped"
    assert (
        "site_sequence is missing or blank" in by_source["row_missing_sequence"].reason
    )
    assert by_source["row_incomplete"].action == "dropped"
    assert (
        by_source["row_incomplete"].parameter_snapshot["required_observed_count"] == 2
    )
    assert by_source["row_a"].action == "retained"
    assert by_source["row_b"].action == "collapsed"
    assert by_source["row_b"].retained_row_id == "row_a"
    assert by_source["row_b"].parameter_snapshot["metadata_conflict_detected"] is True


def test_site_matrix_stage_orchestration_remains_stable_for_representative_fixture() -> (
    None
):
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 3.0, 2.0, 4.0],
            "sample_b": [2.0, 4.0, 5.0, float("nan")],
        },
        index=pd.Index(["row_a", "row_b", "row_missing", "row_incomplete"], name="src"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "MAPK14", "AKT1", "GSK3B"],
            "site": ["Y182", "Y182", "T308", "S9"],
            "site_sequence": ["SEQ_A", "SEQ_B", "", "SEQ_D"],
            "protein_id": ["P_A", "P_A", "P_C", "P_D"],
        },
        index=phospho.index.copy(),
    )
    plan = PreprocessingPlan(
        stage_order=("site_matrix",),
        site_matrix_policy=DATASET_SITE_MATRIX_POLICY_BUILD_FROM_METADATA,
        site_matrix_duplicate_site_policy=DATASET_SITE_MATRIX_DUPLICATE_POLICY_FIRST,
        site_matrix_missing_data_policy=DATASET_SITE_MATRIX_MISSING_DATA_POLICY_DROP_ANY_MISSING,
    )
    state = PreprocessingState(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=None,
        total=None,
        plan=plan,
    )

    result = SiteMatrixStage().run(state)
    next_state = result.state

    assert next_state.phospho.index.tolist() == ["MAPK14;Y182;"]
    assert next_state.site_metadata.index.tolist() == ["MAPK14;Y182;"]
    assert next_state.site_metadata.loc["MAPK14;Y182;", "protein_id"] == "P_A"
    assert result.diagnostics["dropped_row_ids"] == (
        "row_missing",
        "row_incomplete",
        "row_b",
    )
    assert result.diagnostics["dropped_row_count"] == 3

    assert next_state.duplicate_site_resolution is not None
    assert next_state.metadata_conflicts is not None
    assert next_state.row_audit is not None

    duplicate_rows = next_state.duplicate_site_resolution.sort_values("source_row_id")
    assert duplicate_rows["source_row_id"].tolist() == ["row_a", "row_b"]
    assert duplicate_rows["retained"].tolist() == [True, False]
    assert duplicate_rows["resolution_policy"].tolist() == ["first", "first"]

    conflicts = next_state.metadata_conflicts
    sequence_conflicts = conflicts.loc[conflicts["field"] == "site_sequence"]
    assert sequence_conflicts.shape[0] == 1
    assert sequence_conflicts.iloc[0]["source_row_ids"] == ("row_a", "row_b")

    site_matrix_audit = next_state.row_audit.loc[
        next_state.row_audit["stage"] == "site_matrix"
    ]
    assert site_matrix_audit.shape[0] == 4
    assert set(site_matrix_audit["action"]) == {"dropped", "retained", "collapsed"}


def test_site_matrix_stage_records_explicit_duplicate_missing_value_policy() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, float("nan"), 4.0],
            "sample_b": [float("nan"), 3.0, 6.0],
        },
        index=pd.Index(["row_a", "row_b", "row_c"], name="src"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "MAPK14", "AKT1"],
            "site": ["Y182", "Y182", "T308"],
            "site_sequence": ["SEQ_A", "SEQ_A", "SEQ_C"],
            "protein_id": ["P_A", "P_A", "P_C"],
        },
        index=phospho.index.copy(),
    )
    plan = PreprocessingPlan(
        stage_order=("site_matrix",),
        site_matrix_policy=DATASET_SITE_MATRIX_POLICY_BUILD_FROM_METADATA,
        site_matrix_duplicate_site_policy=DATASET_SITE_MATRIX_DUPLICATE_POLICY_AGGREGATE_MEAN,
        site_matrix_missing_data_policy="retain_missing",
    )
    state = PreprocessingState(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=None,
        total=None,
        plan=plan,
    )

    result = SiteMatrixStage().run(state)
    diagnostics = result.diagnostics["diagnostics"]
    duplicate_aggregation = diagnostics["duplicate_aggregation"]
    assert duplicate_aggregation["aggregation_method"] == "aggregate_mean"
    assert duplicate_aggregation["missing_value_policy"] == "skip_missing_values"
    assert duplicate_aggregation["aggregation_reduced_missingness"] is True
    assert (
        result.state.phospho.attrs["site_matrix_provenance"]["duplicate_aggregation"][
            "missing_value_policy"
        ]
        == "skip_missing_values"
    )
