from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from phospy.errors import DatasetValidationError, PhosPyInputError
from phospy.errors.validation import WorkflowValidationError
from phospy.science.differential.models import (
    DifferentialAnalysisResult,
    EmpiricalBayesPriorDiagnostics,
)
from phospy.tables.datasets import SiteMetadataTable
from phospy.validation.workflows.identity import (
    DIFFERENTIAL_IDENTITY_CONTRACT,
    KINASE_IDENTITY_CONTRACT,
    enforce_workflow_site_identity_contract,
)
from tests.support.site_keys import protein_site_key, site_key_context_columns


def _site_index() -> pd.Index:
    return pd.Index(
        [
            protein_site_key(protein_identifier="P28482", site="Y182"),
            protein_site_key(protein_identifier="Q99999", site="Y182"),
        ],
        name="site_key",
    )


def _identity_metadata() -> pd.DataFrame:
    site_index = _site_index()
    return pd.DataFrame(
        {
            "site_key": site_index.astype(str).tolist(),
            "display_id": ["MAPK14;Y182;", "MAPK14;Y182;"],
            **site_key_context_columns(site_index),
            "gene_symbol": ["MAPK14", "MAPK14"],
            "site": ["Y182", "Y182"],
            "site_sequence": ["AAAAAAAYAAAAAAA", "AAAAAAAYAAAAAAA"],
            "protein_id": ["P28482", "Q99999"],
        },
        index=site_index.copy(),
    )


def _as_result_table(site_metadata: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "site_key": site_metadata.loc[:, "site_key"].astype(str).tolist(),
            "display_id": site_metadata.loc[:, "display_id"].astype(str).tolist(),
            "organism": site_metadata.loc[:, "organism"].astype(str).tolist(),
            "protein_namespace": site_metadata.loc[:, "protein_namespace"]
            .astype(str)
            .tolist(),
            "protein_identifier": site_metadata.loc[:, "protein_identifier"]
            .astype(str)
            .tolist(),
            "gene_symbol": site_metadata.loc[:, "gene_symbol"].astype(str).tolist(),
            "site": site_metadata.loc[:, "site"].astype(str).tolist(),
            "logFC": [1.0, -1.0],
            "t": [2.0, -2.0],
            "P.Value": [0.05, 0.10],
            "adj.P.Val": [0.10, 0.10],
        },
        index=site_metadata.index.copy(),
    )


def _prior_diagnostics(index: pd.Index) -> EmpiricalBayesPriorDiagnostics:
    return EmpiricalBayesPriorDiagnostics(
        method="standard",
        robust=False,
        trend=False,
        winsor_tail_p=(0.05, 0.1),
        base_prior_variance=1.0,
        base_prior_degrees_of_freedom=10.0,
        robust_outlier_count=0,
        robust_outlier_fraction=0.0,
        winsorized_low_count=0,
        winsorized_high_count=0,
        prior_variance=pd.Series(np.ones(index.size), index=index.copy()),
        prior_degrees_of_freedom=pd.Series(
            np.full(index.size, 10.0), index=index.copy()
        ),
    )


def _result_with_table(table: pd.DataFrame) -> DifferentialAnalysisResult:
    index = table.index.copy()
    return DifferentialAnalysisResult(
        residual_variance=pd.Series(np.ones(index.size), index=index.copy()),
        posterior_residual_variance=pd.Series(np.ones(index.size), index=index.copy()),
        prior_residual_variance=pd.Series(np.ones(index.size), index=index.copy()),
        prior_degrees_of_freedom_series_value=pd.Series(
            np.full(index.size, 10.0),
            index=index.copy(),
        ),
        prior_variance=1.0,
        prior_degrees_of_freedom=10.0,
        residual_degrees_of_freedom=4.0,
        empirical_bayes_method="standard",
        empirical_bayes_robust=False,
        empirical_bayes_trend=False,
        prior_diagnostics=_prior_diagnostics(index),
        mean_variance_trend_diagnostics=None,
        contrast_tables={"B_vs_A": table},
    )


def _validate_dataset_boundary(site_metadata: pd.DataFrame) -> None:
    SiteMetadataTable(frame=site_metadata, expected_index=site_metadata.index.copy())


def _validate_differential_workflow_identity(site_metadata: pd.DataFrame) -> None:
    enforce_workflow_site_identity_contract(
        site_metadata=site_metadata,
        expected_index=site_metadata.index.copy(),
        expected_index_field_name="test.dataset.phospho.index",
        field_name="test.dataset.site_metadata",
        contract=DIFFERENTIAL_IDENTITY_CONTRACT,
        error_type=WorkflowValidationError,
        allow_opaque_site_values=False,
    )


def _validate_kinase_workflow_identity(site_metadata: pd.DataFrame) -> None:
    enforce_workflow_site_identity_contract(
        site_metadata=site_metadata,
        expected_index=site_metadata.index.copy(),
        expected_index_field_name="test.dataset.phospho.index",
        field_name="test.dataset.site_metadata",
        contract=KINASE_IDENTITY_CONTRACT,
        error_type=WorkflowValidationError,
        allow_opaque_site_values=False,
    )


def _validate_result_boundary(site_metadata: pd.DataFrame) -> None:
    _result_with_table(_as_result_table(site_metadata))


def test_duplicate_display_id_with_distinct_site_keys_is_accepted_across_boundaries() -> (
    None
):
    site_metadata = _identity_metadata()

    _validate_dataset_boundary(site_metadata)
    _validate_differential_workflow_identity(site_metadata)
    _validate_result_boundary(site_metadata)


def test_duplicate_site_key_is_rejected_across_identity_boundaries() -> None:
    site_metadata = _identity_metadata()
    duplicate_index = pd.Index(
        [site_metadata.index[0], site_metadata.index[0]],
        name="site_key",
    )
    site_metadata = site_metadata.set_axis(duplicate_index, axis="index")
    site_metadata.loc[:, "site_key"] = duplicate_index.astype(str).tolist()

    with pytest.raises(DatasetValidationError, match="unique site_key"):
        _validate_dataset_boundary(site_metadata)
    with pytest.raises(WorkflowValidationError, match="unique site_key"):
        _validate_differential_workflow_identity(site_metadata)
    with pytest.raises(PhosPyInputError, match="must be unique"):
        _validate_result_boundary(site_metadata)


def test_site_key_column_index_mismatch_is_rejected_across_identity_boundaries() -> (
    None
):
    site_metadata = _identity_metadata()
    off_index_key = protein_site_key(protein_identifier="P31749", site="T308")
    site_metadata.loc[site_metadata.index[0], "site_key"] = off_index_key

    with pytest.raises(DatasetValidationError, match="site_key.*must exactly match"):
        _validate_dataset_boundary(site_metadata)
    with pytest.raises(WorkflowValidationError, match="site_key.*must exactly match"):
        _validate_differential_workflow_identity(site_metadata)
    with pytest.raises(PhosPyInputError, match="site_key.*must exactly match"):
        _validate_result_boundary(site_metadata)


def test_site_key_metadata_mismatch_is_rejected_across_identity_boundaries() -> None:
    site_metadata = _identity_metadata()
    site_metadata.loc[site_metadata.index[0], "protein_identifier"] = "P99999"

    with pytest.raises(DatasetValidationError, match="metadata-derived"):
        _validate_dataset_boundary(site_metadata)
    with pytest.raises(WorkflowValidationError, match="metadata-derived"):
        _validate_differential_workflow_identity(site_metadata)
    with pytest.raises(PhosPyInputError, match="protein_identifier is incoherent"):
        _validate_result_boundary(site_metadata)


def test_boundaries_require_different_identity_strictness_without_duplicate_logic() -> (
    None
):
    site_metadata = _identity_metadata().drop(columns=["site_sequence"])

    with pytest.raises(
        DatasetValidationError, match="missing required columns: site_sequence"
    ):
        _validate_dataset_boundary(site_metadata)
    _validate_differential_workflow_identity(site_metadata)
    _validate_result_boundary(site_metadata)
    with pytest.raises(
        WorkflowValidationError, match="requires centred sequence context"
    ):
        _validate_kinase_workflow_identity(site_metadata)
