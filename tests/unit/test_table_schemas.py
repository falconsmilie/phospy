from __future__ import annotations

import pandas as pd
import pytest

from phospy.errors.validation import (
    DatasetValidationError,
    PhosPyValidationError,
    ReferenceValidationError,
    WorkflowValidationError,
)
from phospy.references.models import Organism, ReferenceBundle
from phospy.tables.activity import ActivityMatrix, ActivityTargetTable
from phospy.tables.datasets import (
    PhosphoIntensityMatrix,
    SampleMetadataTable,
    SiteMetadataTable,
    TotalProteinMatrix,
)
from phospy.tables.kinase import KinasePredictionMatrix
from phospy.tables.references import KinaseSubstrateReference, SiteSequenceReference
from phospy.tables.signalome import SignalomeProteinSiteContext, SignalomeSiteContext


def _phospho_frame() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample_a": [1.0, 2.0],
            "sample_b": [1.5, 2.5],
        },
        index=pd.Index(["MAPK14;Y182;", "AKT1;T308;"], name="site_id"),
    )


def _site_metadata_frame(index: pd.Index) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "AKT1"],
            "site": ["Y182", "T308"],
            "site_sequence": ["A" * 31, "B" * 31],
            "protein_id": ["P28482", "P31749"],
        },
        index=index.copy(),
    )


def _site_membership_frame() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "site_id": ["MAPK14;Y182;", "AKT1;T308;"],
            "protein_id": ["P28482", "P31749"],
            "site_cluster": [1, pd.NA],
            "protein_module_id": [1, 0],
            "included_in_module_table": [True, False],
            "excluded_reason": ["", "dropped_all_missing_downstream_scores"],
            "gene_symbol": ["MAPK14", "AKT1"],
            "top_kinase": ["MAP2K6", ""],
            "top_kinase_score": [0.91, float("nan")],
            "top_kinase_weight": [0.83, float("nan")],
            "n_supported_kinases": [2, 0],
        }
    )


def _protein_site_context_frame() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "protein_id": ["P28482"],
            "n_sites": [2],
            "site_ids": ['["MAPK14;Y182;","MAPK14;T185;"]'],
            "site_clusters": ["[1,2]"],
            "n_distinct_site_clusters": [2],
            "protein_module_id": [1],
            "multi_site_protein": [True],
            "ambiguous_module_context": [True],
            "gene_symbol": ["MAPK14"],
            "top_kinases_by_site": ['{"MAPK14;Y182;":"MAP2K6","MAPK14;T185;":"MAPK1"}'],
            "module_ids_by_site": ['{"MAPK14;Y182;":1,"MAPK14;T185;":2}'],
        }
    )


def test_dataset_schema_valid_phospho_matrix_passes() -> None:
    wrapper = PhosphoIntensityMatrix(frame=_phospho_frame())
    assert wrapper.frame.shape == (2, 2)


def test_dataset_schema_non_numeric_phospho_column_fails() -> None:
    bad = _phospho_frame().astype(object)
    bad.loc[:, "sample_a"] = ["x", "y"]
    with pytest.raises(DatasetValidationError, match="numeric columns"):
        PhosphoIntensityMatrix(frame=bad)


def test_dataset_schema_missing_phospho_value_fails() -> None:
    bad = _phospho_frame().copy(deep=True)
    bad.loc["MAPK14;Y182;", "sample_a"] = float("nan")
    with pytest.raises(DatasetValidationError, match="must not contain missing values"):
        PhosphoIntensityMatrix(frame=bad)


def test_dataset_schema_duplicate_phospho_index_fails() -> None:
    bad = _phospho_frame().copy(deep=True)
    bad.index = pd.Index(["MAPK14;Y182;", "MAPK14;Y182;"])
    with pytest.raises(DatasetValidationError, match="index must be unique"):
        PhosphoIntensityMatrix(frame=bad)


def test_dataset_schema_duplicate_phospho_columns_fail() -> None:
    bad = _phospho_frame().copy(deep=True)
    bad.columns = pd.Index(["sample_a", "sample_a"])
    with pytest.raises(DatasetValidationError, match="columns must be unique"):
        PhosphoIntensityMatrix(frame=bad)


def test_dataset_schema_non_canonical_site_index_fails() -> None:
    bad = _phospho_frame().copy(deep=True)
    bad.index = pd.Index(["MAPK14;Y182; ", "AKT1;T308;"])
    with pytest.raises(DatasetValidationError, match="canonical site identifiers"):
        PhosphoIntensityMatrix(frame=bad)


def test_dataset_schema_site_metadata_missing_required_columns_fails() -> None:
    phospho = PhosphoIntensityMatrix(frame=_phospho_frame())
    bad = _site_metadata_frame(phospho.frame.index).drop(columns=["gene_symbol"])
    with pytest.raises(DatasetValidationError, match="missing required columns"):
        SiteMetadataTable(frame=bad, expected_index=phospho.frame.index)


def test_dataset_schema_site_metadata_index_mismatch_fails() -> None:
    phospho = PhosphoIntensityMatrix(frame=_phospho_frame())
    bad = _site_metadata_frame(pd.Index(["MAPK14;Y182;", "GSK3B;S9;"], name="site_id"))
    with pytest.raises(DatasetValidationError, match="must exactly match"):
        SiteMetadataTable(frame=bad, expected_index=phospho.frame.index)


def test_dataset_schema_site_metadata_identity_mismatch_fails() -> None:
    phospho = PhosphoIntensityMatrix(frame=_phospho_frame())
    bad = _site_metadata_frame(phospho.frame.index).copy(deep=True)
    bad.loc["MAPK14;Y182;", "gene_symbol"] = "MAPK1"
    with pytest.raises(DatasetValidationError, match="site-identity coherence failed"):
        SiteMetadataTable(frame=bad, expected_index=phospho.frame.index)


def test_dataset_schema_sample_metadata_index_mismatch_fails() -> None:
    phospho = PhosphoIntensityMatrix(frame=_phospho_frame())
    sample_metadata = pd.DataFrame(
        {"group": ["g1", "g2"]},
        index=pd.Index(["sample_a", "sample_x"], name="sample_id"),
    )
    with pytest.raises(DatasetValidationError, match="must exactly match"):
        SampleMetadataTable(frame=sample_metadata, expected_index=phospho.frame.columns)


def test_dataset_schema_total_matrix_column_mismatch_fails() -> None:
    phospho = PhosphoIntensityMatrix(frame=_phospho_frame())
    total = pd.DataFrame(
        {"sample_a": [1.0], "sample_x": [2.0]},
        index=pd.Index(["P28482"], name="protein_id"),
    )
    with pytest.raises(DatasetValidationError, match="must exactly match"):
        TotalProteinMatrix(frame=total, expected_sample_index=phospho.frame.columns)


def test_reference_schema_valid_kinase_substrate_reference_passes() -> None:
    frame = pd.DataFrame(
        {
            "kinase": ["MAP2K6", "AKT1"],
            "substrate_site": ["MAPK14;Y182;", "AKT1;T308;"],
        }
    )
    wrapper = KinaseSubstrateReference(frame=frame)
    assert wrapper.frame.shape == (2, 2)


def test_reference_schema_missing_required_column_fails() -> None:
    with pytest.raises(ReferenceValidationError, match="missing required columns"):
        KinaseSubstrateReference(
            frame=pd.DataFrame({"kinase": ["MAP2K6"]}),
        )


def test_reference_schema_duplicate_pairs_fail() -> None:
    frame = pd.DataFrame(
        {
            "kinase": ["MAP2K6", "MAP2K6"],
            "substrate_site": ["MAPK14;Y182;", "MAPK14;Y182;"],
        }
    )
    with pytest.raises(ReferenceValidationError, match="duplicate"):
        KinaseSubstrateReference(frame=frame)


def test_reference_schema_non_canonical_kinase_fails() -> None:
    frame = pd.DataFrame(
        {
            "kinase": [" MAP2K6 "],
            "substrate_site": ["MAPK14;Y182;"],
        }
    )
    with pytest.raises(ReferenceValidationError, match="canonical non-empty string"):
        KinaseSubstrateReference(frame=frame)


def test_reference_schema_non_canonical_substrate_site_fails() -> None:
    frame = pd.DataFrame(
        {
            "kinase": ["MAP2K6"],
            "substrate_site": ["MAPK14;Y182; "],
        }
    )
    with pytest.raises(ReferenceValidationError, match="canonical site identifiers"):
        KinaseSubstrateReference(frame=frame)


def test_reference_schema_missing_site_sequence_for_substrate_fails() -> None:
    with pytest.raises(ReferenceValidationError, match="missing sequence entries"):
        ReferenceBundle(
            organism=Organism.RAT,
            kinase_substrate_map=pd.DataFrame(
                {
                    "kinase": ["MAP2K6"],
                    "substrate_site": ["MAPK14;Y182;"],
                }
            ),
            site_sequences=pd.DataFrame(
                {"site_sequence": ["A" * 31]},
                index=pd.Index(["AKT1;T308;"], name="site_id"),
            ),
        )


def test_prediction_schema_valid_prediction_matrix_passes() -> None:
    pred_mat = pd.DataFrame(
        {"MAP2K6": [0.9, 0.8]},
        index=pd.Index(["MAPK14;Y182;", "AKT1;T308;"], name="site_id"),
    )
    wrapper = KinasePredictionMatrix(frame=pred_mat)
    assert wrapper.frame.shape == (2, 1)


def test_prediction_schema_duplicate_kinase_columns_fail() -> None:
    pred_mat = pd.DataFrame(
        [[0.9, 0.8]],
        index=pd.Index(["MAPK14;Y182;"], name="site_id"),
        columns=pd.Index(["MAP2K6", "MAP2K6"], name="kinase"),
    )
    with pytest.raises(PhosPyValidationError, match="columns must be unique"):
        KinasePredictionMatrix(frame=pred_mat)


def test_prediction_schema_non_canonical_site_index_fails() -> None:
    pred_mat = pd.DataFrame(
        {"MAP2K6": [0.9]},
        index=pd.Index(["MAPK14;Y182; "], name="site_id"),
    )
    with pytest.raises(PhosPyValidationError, match="canonical site identifiers"):
        KinasePredictionMatrix(frame=pred_mat)


def test_prediction_schema_out_of_range_score_fails() -> None:
    pred_mat = pd.DataFrame(
        {"MAP2K6": [1.2]},
        index=pd.Index(["MAPK14;Y182;"], name="site_id"),
    )
    with pytest.raises(PhosPyValidationError, match="between 0.0 and 1.0"):
        KinasePredictionMatrix(frame=pred_mat)


def test_activity_schema_non_numeric_matrix_fails() -> None:
    bad = pd.DataFrame(
        {"sample_a": ["x"]},
        index=pd.Index(["MAP2K6"], name="kinase"),
    )
    with pytest.raises(PhosPyValidationError, match="numeric columns"):
        ActivityMatrix(frame=bad)


def test_activity_schema_target_table_missing_required_columns_fails() -> None:
    with pytest.raises(PhosPyValidationError, match="missing required columns"):
        ActivityTargetTable(
            frame=pd.DataFrame({"site_id": ["MAPK14;Y182;"], "score": [0.8]}),
        )


def test_signalome_schema_valid_site_membership_table_passes() -> None:
    wrapper = SignalomeSiteContext(frame=_site_membership_frame())
    assert wrapper.frame.shape[0] == 2


def test_signalome_schema_missing_required_site_membership_column_fails() -> None:
    bad = _site_membership_frame().drop(columns=["top_kinase_weight"])
    with pytest.raises(WorkflowValidationError, match="missing required columns"):
        SignalomeSiteContext(frame=bad)


def test_signalome_schema_invalid_site_id_fails() -> None:
    bad = _site_membership_frame().copy(deep=True)
    bad.loc[0, "site_id"] = "MAPK14;Y182; "
    with pytest.raises(WorkflowValidationError, match="canonical site identifiers"):
        SignalomeSiteContext(frame=bad)


def test_signalome_schema_invalid_boolean_integer_numeric_columns_fail() -> None:
    bad = _site_membership_frame().copy(deep=True)
    bad = bad.astype({"included_in_module_table": object})
    bad.loc[0, "included_in_module_table"] = "yes"
    with pytest.raises(WorkflowValidationError, match="boolean"):
        SignalomeSiteContext(frame=bad)


def test_signalome_schema_valid_protein_site_context_passes() -> None:
    wrapper = SignalomeProteinSiteContext(frame=_protein_site_context_frame())
    assert wrapper.frame.shape[0] == 1


def test_signalome_schema_malformed_json_columns_fail() -> None:
    bad = _protein_site_context_frame().copy(deep=True)
    bad.loc[0, "site_ids"] = "[invalid"
    with pytest.raises(WorkflowValidationError, match="parseable JSON"):
        SignalomeProteinSiteContext(frame=bad)


def test_reference_site_sequence_wrapper_missing_required_column_fails() -> None:
    with pytest.raises(ReferenceValidationError, match="missing required columns"):
        SiteSequenceReference(
            frame=pd.DataFrame(
                {"sequence": ["A" * 31]},
                index=pd.Index(["MAPK14;Y182;"], name="site_id"),
            ),
        )
