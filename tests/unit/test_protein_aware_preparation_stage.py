from __future__ import annotations

import pandas as pd
import pandas.testing as pdt

from phospy.advanced import DatasetProteinAwarePreparationConfig
from phospy.preparation_records.protein_aware import (
    ProteinAwarePreparationEligibility,
)
from phospy.preprocessing.protein_aware.alignment import (
    PROTEIN_AWARE_REASON_AMBIGUOUS_PROTEIN_MAPPING,
    PROTEIN_AWARE_REASON_INCOMPATIBLE_TRANSFORMATION_STATE,
    PROTEIN_AWARE_REASON_MISSING_TOTAL_PROTEIN_ROW,
    PROTEIN_AWARE_REASON_SAMPLE_MISMATCH,
)
from phospy.preprocessing.protein_aware.mapping import (
    ProteinMappingConfig,
)
from phospy.preprocessing.results.protein_aware import (
    ProteinAwarePreparationStageOutput,
)
from phospy.preprocessing.stages.protein_aware import (
    ProteinAwarePreparationStage,
)
from phospy.quantitative.models import (
    IntensityScaleState,
    MatrixIntensityScaleState,
)


def _state() -> IntensityScaleState:
    return IntensityScaleState(
        phospho=MatrixIntensityScaleState.log2(established_by="test.phospho"),
        total=MatrixIntensityScaleState.log2(established_by="test.total"),
    )


def _incompatible_state() -> IntensityScaleState:
    return IntensityScaleState(
        phospho=MatrixIntensityScaleState.log2(established_by="test.phospho"),
        total=MatrixIntensityScaleState.linear(established_by="test.total"),
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


def _site_keys(
    output: ProteinAwarePreparationStageOutput,
    eligibility: ProteinAwarePreparationEligibility,
) -> tuple[str, ...]:
    return tuple(
        row.site_key
        for row in output.site_eligibility
        if row.eligibility is eligibility
    )


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
    assert _site_keys(
        result,
        ProteinAwarePreparationEligibility.ELIGIBLE_FOR_PROTEIN_AWARE_PREPARATION,
    ) == ("MAPK14;Y182;",)
    assert (
        _site_keys(result, ProteinAwarePreparationEligibility.FALLBACK_TO_PHOSPHO_ONLY)
        == ()
    )
    assert (
        _site_keys(result, ProteinAwarePreparationEligibility.EXCLUDED_FROM_PREPARATION)
        == ()
    )
    assert result.policy_parameters["preparation_mode"] == (
        "aligned_model_input_preparation_only"
    )
    assert result.policy_parameters["modifies_phospho_matrix"] is False
    assert result.policy_parameters["performs_total_protein_subtraction"] is (False)
    assert result.policy_parameters["performs_normalisation"] is False
    assert result.policy_parameters["performs_differential_modelling"] is False
    assert result.policy_parameters["claims_msstatsptm_equivalence"] is False
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
    assert _site_keys(
        result,
        ProteinAwarePreparationEligibility.ELIGIBLE_FOR_PROTEIN_AWARE_PREPARATION,
    ) == ("MAPK14;Y182;",)
    assert _site_keys(
        result, ProteinAwarePreparationEligibility.FALLBACK_TO_PHOSPHO_ONLY
    ) == ("AKT1;T308;",)
    assert (
        _site_keys(result, ProteinAwarePreparationEligibility.EXCLUDED_FROM_PREPARATION)
        == ()
    )
    assert result.protein_covariate_matrix.index.tolist() == ["P53778"]
    missing = result.mapping_diagnostics.missing_protein_abundance
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
    assert (
        _site_keys(
            result,
            ProteinAwarePreparationEligibility.ELIGIBLE_FOR_PROTEIN_AWARE_PREPARATION,
        )
        == ()
    )
    assert (
        _site_keys(result, ProteinAwarePreparationEligibility.FALLBACK_TO_PHOSPHO_ONLY)
        == ()
    )
    assert _site_keys(
        result, ProteinAwarePreparationEligibility.EXCLUDED_FROM_PREPARATION
    ) == ("MAPK14;Y182;",)
    assert result.matched_pairs.empty
    assert result.protein_covariate_matrix.empty
    ambiguous = result.mapping_diagnostics.ambiguous_mapping
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
    assert (
        _site_keys(
            result,
            ProteinAwarePreparationEligibility.ELIGIBLE_FOR_PROTEIN_AWARE_PREPARATION,
        )
        == ()
    )
    assert _site_keys(
        result, ProteinAwarePreparationEligibility.EXCLUDED_FROM_PREPARATION
    ) == ("MAPK14;Y182;",)
    assert result.sample_alignment.missing_total_protein_samples == ("sample_b",)
    row = result.site_eligibility[0]
    assert (
        row.eligibility is ProteinAwarePreparationEligibility.EXCLUDED_FROM_PREPARATION
    )
    assert PROTEIN_AWARE_REASON_SAMPLE_MISMATCH in row.reasons
    assert result.protein_covariate_matrix.columns.tolist() == [
        "sample_a",
        "sample_b",
    ]
    assert result.protein_covariate_matrix.empty


def test_protein_aware_preparation_transformation_mismatch_is_diagnostic_exclusion() -> (
    None
):
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
        transformation_state=_incompatible_state(),
        config=_config(),
        mapping_config=_mapping_config(),
    )

    assert result is not None
    assert (
        _site_keys(
            result,
            ProteinAwarePreparationEligibility.ELIGIBLE_FOR_PROTEIN_AWARE_PREPARATION,
        )
        == ()
    )
    assert _site_keys(
        result, ProteinAwarePreparationEligibility.EXCLUDED_FROM_PREPARATION
    ) == ("MAPK14;Y182;",)
    assert result.transformation_state is not None
    assert result.transformation_state.compatible is False
    row = result.site_eligibility[0]
    assert PROTEIN_AWARE_REASON_INCOMPATIBLE_TRANSFORMATION_STATE in row.reasons
    assert result.protein_covariate_matrix.empty
    pdt.assert_frame_equal(phospho, original)
