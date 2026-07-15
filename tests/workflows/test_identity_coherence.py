from __future__ import annotations

import pandas as pd
import pytest

from phospy import DifferentialAnalysisWorkflow, KinaseWorkflow, SignalomeWorkflow
from phospy.errors import WorkflowValidationError
from phospy.science.signalomes.constants import (
    DISPLAY_ID_COLUMN,
    EXPANDED_SIGNALOME_ROW_KIND_COLUMN,
    EXPANDED_SIGNALOME_ROW_KIND_SITE,
    SITE_KEY_COLUMN,
)
from tests.support.site_keys import site_key_context_columns
from tests.support.unsafe_dataset_states import unsafe_set_dataset_site_metadata_columns
from tests.support.workflow_identity_coherence import (
    DUPLICATE_DISPLAY_ID,
    build_duplicate_display_differential_request,
    build_duplicate_display_kinase_dataset,
    build_duplicate_display_kinase_request,
    build_duplicate_display_signalome_request,
    drop_site_metadata_column,
)


def _assert_duplicate_display_rows(
    table: pd.DataFrame,
    *,
    expected_site_keys: list[str],
) -> None:
    assert table.index.astype(str).tolist() == expected_site_keys
    assert table.loc[:, SITE_KEY_COLUMN].astype(str).tolist() == expected_site_keys
    assert table.loc[:, DISPLAY_ID_COLUMN].astype(str).tolist() == [
        DUPLICATE_DISPLAY_ID,
        DUPLICATE_DISPLAY_ID,
    ]
    assert table.loc[:, SITE_KEY_COLUMN].is_unique
    assert int(table.loc[:, DISPLAY_ID_COLUMN].nunique()) == 1


def _assert_signalome_protein_grouping_metadata_message(message: str) -> None:
    assert "Missing signalome protein grouping metadata: protein_id" in message
    assert "dataset.site_metadata.protein_id" in message
    assert "not canonical protein identity" in message
    assert "protein_namespace" in message
    assert "protein_identifier" in message
    assert "does not infer protein_id from gene_symbol or display_id" in message
    assert "identity requirement failed" not in message


def test_differential_workflow_preserves_duplicate_display_ids_by_site_key() -> None:
    request = build_duplicate_display_differential_request()
    expected_site_keys = request.dataset.phospho.index.astype(str).tolist()
    expected_context = site_key_context_columns(expected_site_keys)

    result = DifferentialAnalysisWorkflow().run(request)
    table = result.table_for("B_vs_A")

    _assert_duplicate_display_rows(table, expected_site_keys=expected_site_keys)
    assert table.shape[0] == 2
    for column_name in ("organism", "protein_namespace", "protein_identifier"):
        assert column_name in table.columns
        assert (
            table.loc[:, column_name].astype(str).tolist()
            == (expected_context[column_name])
        )
    assert set(("logFC", "t", "P.Value", "adj.P.Val")).issubset(table.columns)
    observed_logfc = table.loc[expected_site_keys, "logFC"].astype(float).to_dict()
    assert observed_logfc == {
        expected_site_keys[0]: pytest.approx(1.0),
        expected_site_keys[1]: pytest.approx(0.05),
    }


def test_kinase_workflow_preserves_duplicate_display_ids_by_site_key() -> None:
    request = build_duplicate_display_kinase_request()
    expected_site_keys = request.dataset.phospho.index.astype(str).tolist()

    result = KinaseWorkflow().run(request)

    profile_scores = result.scoring_result.profile_scores
    fusion_scores = result.scoring_result.rank_weighted_fusion_scores
    prediction_matrix = result.prediction_result.pred_mat
    substrate_list = result.prediction_result.substrate_list

    assert profile_scores.index.astype(str).tolist() == expected_site_keys
    assert profile_scores.index.is_unique
    assert fusion_scores is not None
    assert fusion_scores.index.astype(str).tolist() == expected_site_keys
    assert fusion_scores.index.is_unique
    assert prediction_matrix.index.astype(str).tolist() == expected_site_keys
    assert prediction_matrix.index.is_unique
    assert prediction_matrix.shape[0] == 2
    assert substrate_list is not None
    assert substrate_list.loc[:, SITE_KEY_COLUMN].astype(str).tolist() == (
        expected_site_keys
    )
    assert substrate_list.loc[:, DISPLAY_ID_COLUMN].astype(str).tolist() == [
        DUPLICATE_DISPLAY_ID,
        DUPLICATE_DISPLAY_ID,
    ]
    assert substrate_list.loc[:, SITE_KEY_COLUMN].is_unique
    assert int(substrate_list.loc[:, DISPLAY_ID_COLUMN].nunique()) == 1


def test_signalome_workflow_preserves_duplicate_display_ids_by_site_key() -> None:
    request = build_duplicate_display_signalome_request()
    expected_site_keys = request.kinase_result.dataset.phospho.index.astype(
        str
    ).tolist()

    result = SignalomeWorkflow().run(request)
    module_assignments = result.module_assignments.table
    site_membership = result.site_membership
    expanded_signalome = result.expanded_signalome

    _assert_duplicate_display_rows(
        module_assignments,
        expected_site_keys=expected_site_keys,
    )
    assert module_assignments.shape[0] == 2
    assert site_membership is not None
    assert site_membership.loc[:, SITE_KEY_COLUMN].astype(str).tolist() == (
        expected_site_keys
    )
    assert site_membership.loc[:, DISPLAY_ID_COLUMN].astype(str).tolist() == [
        DUPLICATE_DISPLAY_ID,
        DUPLICATE_DISPLAY_ID,
    ]
    assert site_membership.loc[:, SITE_KEY_COLUMN].is_unique
    assert int(site_membership.loc[:, DISPLAY_ID_COLUMN].nunique()) == 1

    assert expanded_signalome is not None
    site_rows = expanded_signalome.loc[
        expanded_signalome.loc[:, EXPANDED_SIGNALOME_ROW_KIND_COLUMN].astype(str)
        == EXPANDED_SIGNALOME_ROW_KIND_SITE,
        :,
    ]
    assert site_rows.loc[:, SITE_KEY_COLUMN].nunique() == 2
    assert set(site_rows.loc[:, SITE_KEY_COLUMN].astype(str)) == set(expected_site_keys)
    assert set(site_rows.loc[:, DISPLAY_ID_COLUMN].astype(str)) == {
        DUPLICATE_DISPLAY_ID
    }


def test_signalome_workflow_does_not_repair_missing_site_key_identity() -> None:
    dataset = build_duplicate_display_kinase_dataset()
    request = build_duplicate_display_signalome_request(dataset=dataset)
    drop_site_metadata_column(dataset, SITE_KEY_COLUMN)

    with pytest.raises(
        WorkflowValidationError,
        match="missing required columns: site_key.*does not repair weak dataset identity",
    ):
        SignalomeWorkflow().run(request)


def test_signalome_workflow_rejects_missing_protein_id_grouping_metadata() -> None:
    dataset = build_duplicate_display_kinase_dataset()
    request = build_duplicate_display_signalome_request(dataset=dataset)
    drop_site_metadata_column(dataset, "protein_id")

    with pytest.raises(WorkflowValidationError) as exc_info:
        SignalomeWorkflow().run(request)

    message = str(exc_info.value)
    assert "is missing required columns: protein_id" in message
    _assert_signalome_protein_grouping_metadata_message(message)


def test_signalome_workflow_rejects_blank_protein_id_grouping_metadata() -> None:
    dataset = build_duplicate_display_kinase_dataset()
    request = build_duplicate_display_signalome_request(dataset=dataset)
    unsafe_set_dataset_site_metadata_columns(
        dataset,
        {"protein_id": ["", "   "]},
    )

    with pytest.raises(WorkflowValidationError) as exc_info:
        SignalomeWorkflow().run(request)

    message = str(exc_info.value)
    assert "to contain non-empty string values" in message
    _assert_signalome_protein_grouping_metadata_message(message)
