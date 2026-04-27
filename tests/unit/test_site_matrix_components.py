from __future__ import annotations

import pandas as pd
import pandas.testing as pdt
import pytest

from phospy.api.configs import (
    DATASET_SITE_MATRIX_DUPLICATE_POLICY_AGGREGATE_MEAN,
    DATASET_SITE_MATRIX_DUPLICATE_POLICY_AGGREGATE_MEDIAN,
    DATASET_SITE_MATRIX_DUPLICATE_POLICY_ERROR,
    DATASET_SITE_MATRIX_DUPLICATE_POLICY_FIRST,
    DATASET_SITE_MATRIX_DUPLICATE_POLICY_MAX_MEAN_SIGNAL,
    DATASET_SITE_MATRIX_MISSING_DATA_POLICY_DROP_ANY_MISSING,
    DATASET_SITE_MATRIX_POLICY_BUILD_FROM_METADATA,
)
from phospy.datasets.preprocessing.models import PreprocessingPlan, PreprocessingState
from phospy.datasets.preprocessing.stages.site_matrix import SiteMatrixStage
from phospy.datasets.preprocessing.stages.site_matrix_components import (
    DuplicateSiteResolver,
    MetadataConflictDetector,
    SiteMatrixRowAuditBuilder,
)
from phospy.errors.input import PhosPyInputError


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
            "protein_id": ["PROT_A", "PROT_B", "PROT_C"],
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
                    "protein_id": ["PROT_C", "PROT_B"],
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
                    "site_sequence": ["SEQ_A", "SEQ_C"],
                    "protein_id": ["PROT_A", "PROT_C"],
                    "uid": ["A", "C"],
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
                    "site_sequence": ["SEQ_A", "SEQ_C"],
                    "protein_id": ["PROT_A", "PROT_C"],
                    "uid": ["A", "C"],
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
    protein_conflicts = conflicts.loc[conflicts["field"] == "protein_id"]
    assert protein_conflicts.shape[0] == 1
    assert protein_conflicts.iloc[0]["site_id"] == "MAPK14;Y182;"
    assert protein_conflicts.iloc[0]["source_row_ids"] == ("row_a", "row_b")


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
            "protein_id": ["P_A", "P_B", "P_C", "P_D"],
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
    protein_conflicts = conflicts.loc[conflicts["field"] == "protein_id"]
    assert protein_conflicts.shape[0] == 1
    assert protein_conflicts.iloc[0]["source_row_ids"] == ("row_a", "row_b")

    site_matrix_audit = next_state.row_audit.loc[
        next_state.row_audit["stage"] == "site_matrix"
    ]
    assert site_matrix_audit.shape[0] == 4
    assert set(site_matrix_audit["action"]) == {"dropped", "retained", "collapsed"}
