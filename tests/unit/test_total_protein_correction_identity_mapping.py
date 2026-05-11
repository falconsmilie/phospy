from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from phospy.api.configs import (
    DATASET_TOTAL_PROTEIN_CORRECTION_IDENTITY_MATCHING_POLICY_GENE_SYMBOL_NORMALISED,
    DATASET_TOTAL_PROTEIN_CORRECTION_UNMATCHED_POLICY_ALLOW_UNCORRECTED,
    DatasetIntensityTransformConfig,
    DatasetPreprocessingConfig,
    DatasetTotalProteinCorrectionConfig,
    DatasetTotalProteinCorrectionIdentityConfig,
)
from phospy.datasets.builders.executor import DatasetBuildExecutor
from phospy.datasets.builders.interpreter import DatasetBuildRequestInterpreter
from phospy.datasets.builders.preprocessing import DatasetPreprocessor
from phospy.datasets.preprocessing.models import PreprocessingPlan
from phospy.errors.input import PhosPyInputError


def _phospho() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample_a": [10.0, 20.0, 30.0],
            "sample_b": [11.0, 21.0, 31.0],
        },
        index=pd.Index(
            ["MAPK14;Y182;", "MAPK14;T180;", "AKT1;T308;"],
            name="site_id",
        ),
    )


def _site_metadata() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "MAPK14", "AKT1"],
            "protein_accession": ["P53778", "P53778", "P31749"],
            "protein_group_id": ["PG_A", "PG_A", "PG_B"],
            "site": ["Y182", "T180", "T308"],
            "site_sequence": ["AAAAAYAAAAA", "AAAAATAAAAA", "AAAAATAAAAA"],
            "localisation_confidence": [0.95, 0.94, 0.93],
        },
        index=_phospho().index.copy(),
    )


def _total_with_index(index: list[str]) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample_a": [1.0 + i for i in range(len(index))],
            "sample_b": [2.0 + i for i in range(len(index))],
        },
        index=pd.Index(index, name="protein_id"),
    )


def _run_preprocessor(
    *,
    identity: DatasetTotalProteinCorrectionIdentityConfig,
    total: pd.DataFrame,
    site_metadata: pd.DataFrame | None = None,
) -> tuple[pd.DataFrame, dict[str, object]]:
    phospho = _phospho()
    metadata = _site_metadata() if site_metadata is None else site_metadata
    plan = PreprocessingPlan.from_config(
        DatasetPreprocessingConfig(
            intensity_transform=DatasetIntensityTransformConfig(policy="log2"),
            total_protein_correction=DatasetTotalProteinCorrectionConfig(
                policy="subtract_log_total",
                identity=identity,
            ),
        )
    )
    preprocessed = DatasetPreprocessor().run(
        phospho=phospho,
        site_metadata=metadata,
        sample_metadata=None,
        total=total,
        plan=plan,
    )
    trace = preprocessed.preprocessing_trace or ()
    diagnostics = {}
    for item in trace:
        if item.stage == "total_protein_correction":
            diagnostics = dict(item.diagnostics)
    return preprocessed.phospho, diagnostics


def test_total_correction_direct_accession_mapping_succeeds() -> None:
    total = _total_with_index(["P53778", "P31749"])
    corrected, diagnostics = _run_preprocessor(
        identity=DatasetTotalProteinCorrectionIdentityConfig(
            mode="direct",
            phosphosite_key="protein_accession",
            total_protein_key="__index__",
        ),
        total=total,
    )
    assert diagnostics["identity_mode"] == "direct"
    assert diagnostics["identity_matching_policy"] == "strict"
    assert diagnostics["corrected_row_count"] == 3
    assert diagnostics["quantitative_meaning"] == "phospho_total_log_ratio"
    assert diagnostics["total_rows_used_by_multiple_phosphosites"] == 1
    assert diagnostics["gene_symbol_matching_used"] is False
    assert corrected.shape == _phospho().shape


def test_total_correction_exact_match_values_unchanged() -> None:
    total = _total_with_index(["P53778", "P31749"])
    corrected, _ = _run_preprocessor(
        identity=DatasetTotalProteinCorrectionIdentityConfig(
            mode="direct",
            phosphosite_key="protein_accession",
            total_protein_key="__index__",
            matching_policy="strict",
        ),
        total=total,
    )
    phospho = _phospho()
    expected = phospho.copy(deep=True)
    expected.loc[:, :] = np.log2(expected.to_numpy(copy=False) + 1.0)
    expected.loc["MAPK14;Y182;", :] = (
        expected.loc["MAPK14;Y182;", :] - np.log2(total.loc["P53778", :] + 1.0)
    ).to_numpy(copy=False)
    expected.loc["MAPK14;T180;", :] = (
        expected.loc["MAPK14;T180;", :] - np.log2(total.loc["P53778", :] + 1.0)
    ).to_numpy(copy=False)
    expected.loc["AKT1;T308;", :] = (
        expected.loc["AKT1;T308;", :] - np.log2(total.loc["P31749", :] + 1.0)
    ).to_numpy(copy=False)
    pd.testing.assert_frame_equal(corrected, expected)


def test_total_correction_direct_protein_group_mapping_succeeds() -> None:
    total = _total_with_index(["PG_A", "PG_B"])
    _, diagnostics = _run_preprocessor(
        identity=DatasetTotalProteinCorrectionIdentityConfig(
            mode="direct",
            phosphosite_key="protein_group_id",
            total_protein_key="__index__",
        ),
        total=total,
    )
    assert diagnostics["corrected_row_count"] == 3
    assert diagnostics["identity_matching_policy"] == "strict"
    assert diagnostics["gene_symbol_matching_used"] is False


def test_total_correction_explicit_gene_symbol_mapping_succeeds() -> None:
    total = _total_with_index(["MAPK14", "AKT1"])
    _, diagnostics = _run_preprocessor(
        identity=DatasetTotalProteinCorrectionIdentityConfig(
            mode="direct",
            phosphosite_key="gene_symbol",
            total_protein_key="__index__",
            matching_policy=DATASET_TOTAL_PROTEIN_CORRECTION_IDENTITY_MATCHING_POLICY_GENE_SYMBOL_NORMALISED,
        ),
        total=total,
    )
    assert diagnostics["gene_symbol_matching_used"] is True
    assert diagnostics["gene_symbol_identity_warning"] is not None
    assert diagnostics["identity_matching_policy"] == "gene_symbol_normalised"


def test_total_correction_strict_gene_symbol_identity_match_succeeds() -> None:
    total = _total_with_index(["MAPK14", "AKT1"])
    _, diagnostics = _run_preprocessor(
        identity=DatasetTotalProteinCorrectionIdentityConfig(
            mode="direct",
            phosphosite_key="gene_symbol",
            total_protein_key="__index__",
            matching_policy="strict",
        ),
        total=total,
    )
    assert diagnostics["corrected_row_count"] == 3
    assert diagnostics["identity_matching_policy"] == "strict"
    assert diagnostics["gene_symbol_matching_used"] is False
    assert diagnostics["gene_symbol_identity_warning"] is None


def test_total_correction_strict_rejects_case_only_gene_symbol_matches() -> None:
    total = _total_with_index(["mapk14", "akt1"])
    with pytest.raises(
        PhosPyInputError,
        match="requires complete phosphosite-to-total mapping",
    ):
        _run_preprocessor(
            identity=DatasetTotalProteinCorrectionIdentityConfig(
                mode="direct",
                phosphosite_key="gene_symbol",
                total_protein_key="__index__",
                matching_policy="strict",
            ),
            total=total,
        )


def test_total_correction_lossy_gene_symbol_normalised_accepts_case_only_matches() -> (
    None
):
    total = _total_with_index(["mapk14", "akt1"])
    _, diagnostics = _run_preprocessor(
        identity=DatasetTotalProteinCorrectionIdentityConfig(
            mode="direct",
            phosphosite_key="gene_symbol",
            total_protein_key="__index__",
            matching_policy=DATASET_TOTAL_PROTEIN_CORRECTION_IDENTITY_MATCHING_POLICY_GENE_SYMBOL_NORMALISED,
        ),
        total=total,
    )
    assert diagnostics["corrected_row_count"] == 3
    assert diagnostics["identity_matching_policy"] == "gene_symbol_normalised"
    assert diagnostics["gene_symbol_matching_used"] is True


def test_total_correction_ambiguous_lossy_gene_symbol_matching_fails() -> None:
    total = _total_with_index(["MAPK14", "mapk14", "AKT1"])
    with pytest.raises(PhosPyInputError, match="duplicate keys"):
        _run_preprocessor(
            identity=DatasetTotalProteinCorrectionIdentityConfig(
                mode="direct",
                phosphosite_key="gene_symbol",
                total_protein_key="__index__",
                matching_policy=DATASET_TOTAL_PROTEIN_CORRECTION_IDENTITY_MATCHING_POLICY_GENE_SYMBOL_NORMALISED,
            ),
            total=total,
        )


def test_total_correction_does_not_fallback_to_gene_symbol_implicitly() -> None:
    total = _total_with_index(["MAPK14", "AKT1"])
    with pytest.raises(
        PhosPyInputError,
        match="requires complete phosphosite-to-total mapping",
    ):
        _run_preprocessor(
            identity=DatasetTotalProteinCorrectionIdentityConfig(
                mode="direct",
                phosphosite_key="protein_accession",
                total_protein_key="__index__",
            ),
            total=total,
        )


def test_total_correction_gene_symbol_normalised_requires_gene_symbol_keys() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="matching_policy='gene_symbol_normalised' requires at least one gene_symbol identity key",
    ):
        DatasetTotalProteinCorrectionIdentityConfig(
            mode="direct",
            phosphosite_key="protein_accession",
            total_protein_key="__index__",
            matching_policy=DATASET_TOTAL_PROTEIN_CORRECTION_IDENTITY_MATCHING_POLICY_GENE_SYMBOL_NORMALISED,
        )


def test_total_correction_mapping_table_mode_succeeds() -> None:
    total = _total_with_index(["TP_MAPK", "TP_AKT"])
    mapping_table = pd.DataFrame(
        {
            "site_ref": ["P53778", "P31749"],
            "total_ref": ["TP_MAPK", "TP_AKT"],
        }
    )
    _, diagnostics = _run_preprocessor(
        identity=DatasetTotalProteinCorrectionIdentityConfig(
            mode="mapping_table",
            phosphosite_key="protein_accession",
            total_protein_key="__index__",
            mapping_table=mapping_table,
            mapping_phosphosite_key="site_ref",
            mapping_total_protein_key="total_ref",
        ),
        total=total,
    )
    assert diagnostics["identity_mode"] == "mapping_table"
    assert diagnostics["identity_matching_policy"] == "strict"
    assert isinstance(diagnostics.get("mapping_table_fingerprint"), str)


def test_total_correction_mapping_table_unknown_phosphosite_fails() -> None:
    total = _total_with_index(["TP_MAPK", "TP_AKT"])
    mapping_table = pd.DataFrame(
        {
            "site_ref": ["UNKNOWN_SITE"],
            "total_ref": ["TP_MAPK"],
        }
    )
    with pytest.raises(PhosPyInputError, match="unknown phosphosite keys"):
        _run_preprocessor(
            identity=DatasetTotalProteinCorrectionIdentityConfig(
                mode="mapping_table",
                phosphosite_key="protein_accession",
                total_protein_key="__index__",
                mapping_table=mapping_table,
                mapping_phosphosite_key="site_ref",
                mapping_total_protein_key="total_ref",
            ),
            total=total,
        )


def test_total_correction_mapping_table_unknown_total_fails() -> None:
    total = _total_with_index(["TP_MAPK", "TP_AKT"])
    mapping_table = pd.DataFrame(
        {
            "site_ref": ["P53778"],
            "total_ref": ["UNKNOWN_TOTAL"],
        }
    )
    with pytest.raises(PhosPyInputError, match="unknown total-protein keys"):
        _run_preprocessor(
            identity=DatasetTotalProteinCorrectionIdentityConfig(
                mode="mapping_table",
                phosphosite_key="protein_accession",
                total_protein_key="__index__",
                mapping_table=mapping_table,
                mapping_phosphosite_key="site_ref",
                mapping_total_protein_key="total_ref",
            ),
            total=total,
        )


def test_total_correction_mapping_table_phosphosite_to_multiple_total_fails() -> None:
    total = _total_with_index(["TP_MAPK", "TP_AKT"])
    mapping_table = pd.DataFrame(
        {
            "site_ref": ["P53778", "P53778"],
            "total_ref": ["TP_MAPK", "TP_AKT"],
        }
    )
    with pytest.raises(
        PhosPyInputError,
        match="maps one phosphosite key to multiple total-protein keys",
    ):
        _run_preprocessor(
            identity=DatasetTotalProteinCorrectionIdentityConfig(
                mode="mapping_table",
                phosphosite_key="protein_accession",
                total_protein_key="__index__",
                mapping_table=mapping_table,
                mapping_phosphosite_key="site_ref",
                mapping_total_protein_key="total_ref",
            ),
            total=total,
        )


def test_total_correction_many_phosphosites_to_one_total_succeeds() -> None:
    total = _total_with_index(["TP_SHARED", "TP_AKT"])
    mapping_table = pd.DataFrame(
        {
            "site_ref": ["P53778", "P31749"],
            "total_ref": ["TP_SHARED", "TP_AKT"],
        }
    )
    _, diagnostics = _run_preprocessor(
        identity=DatasetTotalProteinCorrectionIdentityConfig(
            mode="mapping_table",
            phosphosite_key="protein_accession",
            total_protein_key="__index__",
            mapping_table=mapping_table,
            mapping_phosphosite_key="site_ref",
            mapping_total_protein_key="total_ref",
        ),
        total=total,
    )
    assert diagnostics["total_rows_used_by_multiple_phosphosites"] == 1


def test_total_correction_duplicate_total_identity_keys_fail() -> None:
    total = pd.DataFrame(
        {"sample_a": [1.0, 2.0], "sample_b": [2.0, 3.0]},
        index=pd.Index(["P53778", "P53778"], name="protein_id"),
    )
    with pytest.raises(PhosPyInputError, match="duplicate keys"):
        _run_preprocessor(
            identity=DatasetTotalProteinCorrectionIdentityConfig(
                mode="direct",
                phosphosite_key="protein_accession",
                total_protein_key="__index__",
            ),
            total=total,
        )


def test_total_correction_duplicate_mapping_rows_fail() -> None:
    total = _total_with_index(["TP_MAPK", "TP_AKT"])
    mapping_table = pd.DataFrame(
        {
            "site_ref": ["P53778", "P53778"],
            "total_ref": ["TP_MAPK", "TP_MAPK"],
        }
    )
    with pytest.raises(
        PhosPyInputError, match="duplicate phosphosite-to-total mapping rows"
    ):
        _run_preprocessor(
            identity=DatasetTotalProteinCorrectionIdentityConfig(
                mode="mapping_table",
                phosphosite_key="protein_accession",
                total_protein_key="__index__",
                mapping_table=mapping_table,
                mapping_phosphosite_key="site_ref",
                mapping_total_protein_key="total_ref",
            ),
            total=total,
        )


def test_total_correction_null_or_empty_identity_keys_fail() -> None:
    total = _total_with_index(["P53778", "P31749"])
    bad_metadata = _site_metadata()
    bad_metadata.loc["MAPK14;T180;", "protein_accession"] = ""
    with pytest.raises(PhosPyInputError, match="contains null/empty identifiers"):
        _run_preprocessor(
            identity=DatasetTotalProteinCorrectionIdentityConfig(
                mode="direct",
                phosphosite_key="protein_accession",
                total_protein_key="__index__",
            ),
            total=total,
            site_metadata=bad_metadata,
        )


def test_total_correction_unmatched_rows_can_be_retained_when_configured() -> None:
    total = _total_with_index(["P53778"])
    corrected, diagnostics = _run_preprocessor(
        identity=DatasetTotalProteinCorrectionIdentityConfig(
            mode="direct",
            phosphosite_key="protein_accession",
            total_protein_key="__index__",
            unmatched_policy=DATASET_TOTAL_PROTEIN_CORRECTION_UNMATCHED_POLICY_ALLOW_UNCORRECTED,
        ),
        total=total,
    )
    assert diagnostics["uncorrected_row_count"] == 1
    assert (
        diagnostics["quantitative_meaning"]
        == "mixed_phospho_total_log_ratio_and_phosphosite_log_abundance"
    )
    assert diagnostics["corrected_row_count"] == 2
    assert diagnostics["corrected_phosphosite_row_ids"] == [
        "MAPK14;T180;",
        "MAPK14;Y182;",
    ]
    assert diagnostics["corrected_phosphosite_to_total_protein_row_id"] == {
        "MAPK14;T180;": "P53778",
        "MAPK14;Y182;": "P53778",
    }
    assert diagnostics["unmatched_phosphosite_row_ids"] == ["AKT1;T308;"]
    assert diagnostics["uncorrected_phosphosite_row_reasons"] == {
        "AKT1;T308;": (
            "no_matching_total_protein_row_retained_by_"
            "unmatched_policy_allow_uncorrected"
        )
    }
    assert diagnostics["unused_total_protein_row_count"] == 0
    assert list(corrected.index) == list(_phospho().index)


def test_total_correction_reports_unused_total_rows() -> None:
    total = _total_with_index(["P53778", "P31749", "UNUSED_X"])
    _, diagnostics = _run_preprocessor(
        identity=DatasetTotalProteinCorrectionIdentityConfig(
            mode="direct",
            phosphosite_key="protein_accession",
            total_protein_key="__index__",
        ),
        total=total,
    )
    assert diagnostics["unused_total_protein_row_count"] == 1
    assert diagnostics["unused_total_protein_row_ids"] == ["UNUSED_X"]


def test_total_correction_mapping_table_fingerprint_changes_with_mapping_table() -> (
    None
):
    total = _total_with_index(["TP_MAPK", "TP_AKT"])
    mapping_a = pd.DataFrame({"site_ref": ["P53778"], "total_ref": ["TP_MAPK"]})
    mapping_b = pd.DataFrame({"site_ref": ["P53778"], "total_ref": ["TP_AKT"]})
    _, diagnostics_a = _run_preprocessor(
        identity=DatasetTotalProteinCorrectionIdentityConfig(
            mode="mapping_table",
            phosphosite_key="protein_accession",
            total_protein_key="__index__",
            mapping_table=mapping_a,
            mapping_phosphosite_key="site_ref",
            mapping_total_protein_key="total_ref",
            unmatched_policy=DATASET_TOTAL_PROTEIN_CORRECTION_UNMATCHED_POLICY_ALLOW_UNCORRECTED,
        ),
        total=total,
    )
    _, diagnostics_b = _run_preprocessor(
        identity=DatasetTotalProteinCorrectionIdentityConfig(
            mode="mapping_table",
            phosphosite_key="protein_accession",
            total_protein_key="__index__",
            mapping_table=mapping_b,
            mapping_phosphosite_key="site_ref",
            mapping_total_protein_key="total_ref",
            unmatched_policy=DATASET_TOTAL_PROTEIN_CORRECTION_UNMATCHED_POLICY_ALLOW_UNCORRECTED,
        ),
        total=total,
    )
    assert (
        diagnostics_a["mapping_table_fingerprint"]
        != diagnostics_b["mapping_table_fingerprint"]
    )


def test_total_correction_provenance_records_identity_policy() -> None:
    total = _total_with_index(["P53778", "P31749"])
    request_config = DatasetPreprocessingConfig(
        intensity_transform=DatasetIntensityTransformConfig(policy="log2"),
        total_protein_correction=DatasetTotalProteinCorrectionConfig(
            policy="subtract_log_total",
            identity=DatasetTotalProteinCorrectionIdentityConfig(
                mode="direct",
                phosphosite_key="protein_accession",
                total_protein_key="__index__",
            ),
        ),
    )
    from phospy.api.requests import DatasetBuildRequest

    interpreted = DatasetBuildRequestInterpreter().run(
        DatasetBuildRequest(
            phospho=_phospho(),
            site_metadata=_site_metadata(),
            total=total,
            preprocessing_config=request_config,
        )
    )
    built = DatasetBuildExecutor().run(interpreted)
    assert built.provenance is not None
    stage = next(
        item
        for item in built.provenance.preprocessing_stages
        if item.stage == "total_protein_correction"
    )
    assert stage.diagnostics is not None
    assert stage.diagnostics["identity_mode"] == "direct"
    assert stage.diagnostics["identity_matching_policy"] == "strict"
    assert stage.diagnostics["phosphosite_key"] == "protein_accession"
    assert stage.diagnostics["total_protein_key"] == "__index__"
    preprocessing_plan = built.provenance.workflow_parameters["preprocessing_plan"]
    assert isinstance(preprocessing_plan, dict)
    identity_payload = preprocessing_plan["total_protein_correction_identity_policy"]
    assert isinstance(identity_payload, dict)
    assert identity_payload["mode"] == "direct"
    assert identity_payload["matching_policy"] == "strict"
