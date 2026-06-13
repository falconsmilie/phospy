from __future__ import annotations

import pandas as pd
import pandas.testing as pdt

from phospy.api import DatasetProteinAwarePreparationConfig
from phospy.science.datasets.preprocessing.protein_aware_alignment import (
    PROTEIN_AWARE_REASON_AMBIGUOUS_PROTEIN_MAPPING,
    PROTEIN_AWARE_REASON_MISSING_TOTAL_PROTEIN_ROW,
    PROTEIN_AWARE_REASON_SAMPLE_MISMATCH,
    ProteinAwarePreparationEligibility,
)
from phospy.science.datasets.preprocessing.protein_aware_preparation import (
    ProteinAwarePreparationStage,
)
from phospy.science.datasets.preprocessing.protein_mapping import (
    ProteinMappingConfig,
)
from phospy.science.transformations.models import (
    IntensityScaleState,
    MatrixIntensityScaleState,
)


def _state() -> IntensityScaleState:
    return IntensityScaleState(
        phospho=MatrixIntensityScaleState.log2(established_by="test.phospho"),
        total=MatrixIntensityScaleState.log2(established_by="test.total"),
    )


def _config(
    protein_mapping_policy: str = "require_unambiguous",
) -> DatasetProteinAwarePreparationConfig:
    return DatasetProteinAwarePreparationConfig(
        policy="prepare_model_inputs",
        protein_mapping_policy=protein_mapping_policy,
    )


def _mapping_config() -> ProteinMappingConfig:
    return ProteinMappingConfig(protein_identifier_columns=("protein_id",))


def test_protein_aware_preparation_exact_one_site_to_one_protein() -> None:
    phospho = pd.DataFrame(
        {"sample_a": [1.0], "sample_b": [2.0]},
        index=pd.Index(["MAPK14;Y182;"], name="site_id"),
    )
    original = phospho.copy(deep=True)
    site_metadata = pd.DataFrame(
        {"protein_id": ["P53778"], "gene_symbol": ["MAPK14"], "site": ["Y182"]},
        index=phospho.index.copy(),
    )
    total = pd.DataFrame(
        {"sample_a": [10.0], "sample_b": [11.0]},
        index=pd.Index(["P53778"], name="protein_id"),
    )

    result = ProteinAwarePreparationStage().run(
        phospho=phospho,
        site_metadata=site_metadata,
        total=total,
        transformation_state=_state(),
        config=_config(),
        mapping_config=_mapping_config(),
    )

    assert result is not None
    assert result.report.eligible_site_keys == ("MAPK14;Y182;",)
    assert result.report.fallback_site_keys == ()
    assert result.report.excluded_site_keys == ()
    assert result.matched_pairs.to_dict(orient="records") == [
        {
            "site_key": "MAPK14;Y182;",
            "protein_identifier": "P53778",
            "total_protein_row_key": "P53778",
        }
    ]
    pdt.assert_frame_equal(result.protein_covariate_matrix, total)
    pdt.assert_frame_equal(phospho, original)


def test_protein_aware_preparation_missing_total_row_fallback() -> None:
    phospho = pd.DataFrame(
        {"sample_a": [1.0, 2.0], "sample_b": [3.0, 4.0]},
        index=pd.Index(["MAPK14;Y182;", "AKT1;T308;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "protein_id": ["P53778", "P31749"],
            "gene_symbol": ["MAPK14", "AKT1"],
            "site": ["Y182", "T308"],
        },
        index=phospho.index.copy(),
    )
    total = pd.DataFrame(
        {"sample_a": [10.0], "sample_b": [11.0]},
        index=pd.Index(["P53778"], name="protein_id"),
    )

    result = ProteinAwarePreparationStage().run(
        phospho=phospho,
        site_metadata=site_metadata,
        total=total,
        transformation_state=_state(),
        config=_config("allow_missing_with_report"),
        mapping_config=_mapping_config(),
    )

    assert result is not None
    assert result.report.eligible_site_keys == ("MAPK14;Y182;",)
    assert result.report.fallback_site_keys == ("AKT1;T308;",)
    assert result.report.excluded_site_keys == ()
    assert result.protein_covariate_matrix.index.tolist() == ["P53778"]
    missing = result.missing_protein_abundance_diagnostics
    assert missing.to_dict(orient="records") == [
        {
            "site_key": "AKT1;T308;",
            "protein_identifier": "P31749",
            "mapping_status": "missing_total_protein_row",
            "reason": PROTEIN_AWARE_REASON_MISSING_TOTAL_PROTEIN_ROW,
        }
    ]


def test_protein_aware_preparation_ambiguous_mapping_is_excluded() -> None:
    phospho = pd.DataFrame(
        {"sample_a": [1.0], "sample_b": [2.0]},
        index=pd.Index(["MAPK14;Y182;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "protein_id": ["P53778", "Q16539"],
            "gene_symbol": ["MAPK14", "MAPK14"],
            "site": ["Y182", "Y182"],
        },
        index=pd.Index(["MAPK14;Y182;", "MAPK14;Y182;"], name="site_id"),
    )
    total = pd.DataFrame(
        {"sample_a": [10.0, 12.0], "sample_b": [11.0, 13.0]},
        index=pd.Index(["P53778", "Q16539"], name="protein_id"),
    )

    result = ProteinAwarePreparationStage().run(
        phospho=phospho,
        site_metadata=site_metadata,
        total=total,
        transformation_state=_state(),
        config=_config("allow_missing_with_report"),
        mapping_config=_mapping_config(),
    )

    assert result is not None
    assert result.report.eligible_site_keys == ()
    assert result.report.fallback_site_keys == ()
    assert result.report.excluded_site_keys == ("MAPK14;Y182;",)
    assert result.matched_pairs.empty
    assert result.protein_covariate_matrix.empty
    ambiguous = result.ambiguous_mapping_diagnostics
    assert ambiguous.loc[0, "candidate_protein_identifiers"] == (
        "P53778",
        "Q16539",
    )
    assert ambiguous.loc[0, "reason"] == (
        PROTEIN_AWARE_REASON_AMBIGUOUS_PROTEIN_MAPPING
    )


def test_protein_aware_preparation_sample_mismatch_is_diagnostic_exclusion() -> None:
    phospho = pd.DataFrame(
        {"sample_a": [1.0], "sample_b": [2.0]},
        index=pd.Index(["MAPK14;Y182;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {"protein_id": ["P53778"], "gene_symbol": ["MAPK14"], "site": ["Y182"]},
        index=phospho.index.copy(),
    )
    total = pd.DataFrame(
        {"sample_a": [10.0]},
        index=pd.Index(["P53778"], name="protein_id"),
    )

    result = ProteinAwarePreparationStage().run(
        phospho=phospho,
        site_metadata=site_metadata,
        total=total,
        transformation_state=_state(),
        config=_config(),
        mapping_config=_mapping_config(),
    )

    assert result is not None
    assert result.report.eligible_site_keys == ()
    assert result.report.excluded_site_keys == ("MAPK14;Y182;",)
    assert result.sample_alignment_diagnostics.missing_total_protein_samples == (
        "sample_b",
    )
    table = result.site_eligibility_table
    assert table.loc[0, "eligibility"] == (
        ProteinAwarePreparationEligibility.EXCLUDED_FROM_PREPARATION.value
    )
    assert PROTEIN_AWARE_REASON_SAMPLE_MISMATCH in table.loc[0, "reasons"]
    assert result.protein_covariate_matrix.columns.tolist() == [
        "sample_a",
        "sample_b",
    ]
    assert result.protein_covariate_matrix.empty
