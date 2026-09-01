from __future__ import annotations

import json
from dataclasses import replace
from typing import cast

import numpy as np
import pandas as pd
import pytest

from phospy.advanced import (
    DIFFERENTIAL_PROTEIN_AWARE_METHOD_PROTEIN_COVARIATE_ADJUSTED_MODERATED_LINEAR_MODEL_V1,
    PAIRED_DESIGN_POLICY_DUPLICATE_CORRELATION,
    DifferentialAnalysisConfig,
)
from phospy.api import (
    Contrast,
    DifferentialAnalysisRequest,
    DifferentialAnalysisWorkflow,
    ExperimentalDesign,
    Organism,
    SampleDesignRecord,
)
from phospy.errors.input import PhosPyInputError
from phospy.errors.workflows import WorkflowBoundaryError
from phospy.science.differential.models import (
    ContrastMatrix,
    DesignMatrix,
    DifferentialAnalysisResult,
    EmpiricalBayesPriorDiagnostics,
    ProteinAwareDifferentialDiagnostics,
)
from phospy.science.differential.models.protein_aware import (
    ProteinAwareDifferentialComputationRequest,
    ProteinAwareDifferentialComputationResult,
)
from phospy.science.differential.models.tables import (
    DIFFERENTIAL_RESULT_REASON_PROTEIN_AUGMENTED_DESIGN_RANK_DEFICIENT,
    DIFFERENTIAL_RESULT_REASON_PROTEIN_CONTRAST_NON_ESTIMABLE,
    DIFFERENTIAL_RESULT_REASON_PROTEIN_COVARIATE_ZERO_VARIANCE,
    DIFFERENTIAL_RESULT_REASON_PROTEIN_PREPARATION_EXCLUDED,
    DIFFERENTIAL_RESULT_REASON_PROTEIN_PREPARATION_FALLBACK,
    DIFFERENTIAL_RESULT_STATUS_COLUMN,
    DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN,
    DIFFERENTIAL_RESULT_STATUS_TESTED,
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_AUGMENTED_DESIGN_INVALID,
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_CONTRAST_NON_ESTIMABLE,
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_COVARIATE_INVALID,
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_PREPARATION_INELIGIBLE,
)
from phospy.workflows.differential.interpreter import DifferentialAnalysisInterpreter
from phospy.workflows.differential.models import (
    DifferentialFeatureEligibilityInputs,
    InterpretedDifferentialAnalysisRequest,
    ProteinAwareDifferentialResolvedInputs,
)
from phospy.workflows.differential.result_assembly import DifferentialResultAssembler
from phospy.workflows.differential.validator import DifferentialAnalysisValidator
from tests.support.analysis_ready_dataset_factories import (
    trusted_analysis_ready_dataset_from_tables,
)
from tests.support.intensity_scale_states import (
    supported_log2_intensity_scale_state,
    supported_log2_processing_state,
)
from tests.support.site_keys import protein_site_key_index, site_key_context_columns

_SAMPLE_IDS = ("A_1", "A_2", "A_3", "B_1", "B_2", "B_3")
_PROTEINS = ("MAPK14", "GSK3B", "AKT1", "PRKACA", "RPS6", "EIF4EBP1")
_SITES = ("Y182", "S9", "T308", "S339", "S235", "T37")
_TOTAL_ROW_KEYS = (
    "protein:MAPK14",
    None,
    None,
    "protein:AKT1",
    "protein:PRKACA",
    "protein:RPS6",
)


def test_protein_aware_assembly_expands_to_full_index_with_withholding_statuses() -> (
    None
):
    result = _protein_aware_result()
    table = result.table_for("B_vs_A")

    assert table.index.tolist() == _site_index().tolist()
    assert table["site_key"].tolist() == _site_index().tolist()
    assert table[DIFFERENTIAL_RESULT_STATUS_COLUMN].tolist() == [
        DIFFERENTIAL_RESULT_STATUS_TESTED,
        DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_PREPARATION_INELIGIBLE,
        DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_PREPARATION_INELIGIBLE,
        DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_COVARIATE_INVALID,
        DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_AUGMENTED_DESIGN_INVALID,
        DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_CONTRAST_NON_ESTIMABLE,
    ]
    assert table[DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN].tolist() == [
        "",
        DIFFERENTIAL_RESULT_REASON_PROTEIN_PREPARATION_FALLBACK,
        DIFFERENTIAL_RESULT_REASON_PROTEIN_PREPARATION_EXCLUDED,
        DIFFERENTIAL_RESULT_REASON_PROTEIN_COVARIATE_ZERO_VARIANCE,
        DIFFERENTIAL_RESULT_REASON_PROTEIN_AUGMENTED_DESIGN_RANK_DEFICIENT,
        DIFFERENTIAL_RESULT_REASON_PROTEIN_CONTRAST_NON_ESTIMABLE,
    ]
    assert np.isfinite(
        table.loc[table.index[0], ["logFC", "t", "P.Value", "adj.P.Val"]]
    ).all()
    assert table.iloc[1:][["logFC", "t", "P.Value", "adj.P.Val"]].isna().all().all()


def test_protein_aware_diagnostics_are_defensive_and_site_scoped() -> None:
    result = _protein_aware_result()
    diagnostics = result.protein_aware_diagnostics
    assert diagnostics is not None

    assert diagnostics.method_id == (
        DIFFERENTIAL_PROTEIN_AWARE_METHOD_PROTEIN_COVARIATE_ADJUSTED_MODERATED_LINEAR_MODEL_V1
    )
    assert diagnostics.claim_status == "experimental"
    assert diagnostics.preparation_policy == "prepare_model_inputs"
    assert diagnostics.protein_mapping_policy == "require_unambiguous"
    assert diagnostics.protein_covariate_centered is True
    assert diagnostics.protein_covariate_standardized is False
    assert diagnostics.total_site_count == 6
    assert diagnostics.ordinary_eligible_site_count == 6
    assert diagnostics.protein_preparation_eligible_site_count == 4
    assert diagnostics.tested_site_count == 1
    assert diagnostics.withheld_site_count == 5
    assert diagnostics.execution_sample_order == _SAMPLE_IDS
    assert diagnostics.base_design_rank == 2
    assert diagnostics.expected_augmented_rank == 3
    assert diagnostics.common_augmented_rank == 3
    assert diagnostics.common_augmented_residual_degrees_of_freedom == pytest.approx(
        3.0
    )
    assert diagnostics.distinct_matched_protein_row_count == 4
    assert diagnostics.fitted_protein_row_count == 1
    assert diagnostics.median_augmented_condition_number == pytest.approx(11.0)
    assert result.diagnostics.singular_values == ()
    assert result.diagnostics.condition_number == pytest.approx(11.0)
    assert result.diagnostics.model_type == "protein_covariate_adjusted_moderated_ols"

    per_site = diagnostics.per_site_diagnostics_dataframe()
    assert per_site.index.tolist() == _site_index().tolist()
    assert per_site.loc[
        per_site.index[0], "protein_covariate_coefficient"
    ] == pytest.approx(0.25)
    assert per_site.iloc[1:]["protein_covariate_coefficient"].isna().all()
    assert (
        per_site.loc[
            per_site.index[-1],
            "protein_contrast_estimability_status",
        ]
        == "non_estimable"
    )

    per_site.iloc[0, per_site.columns.get_loc("protein_covariate_coefficient")] = -99.0
    exported_again = diagnostics.per_site_diagnostics_dataframe()
    assert exported_again.loc[
        exported_again.index[0],
        "protein_covariate_coefficient",
    ] == pytest.approx(0.25)

    via_result = result.protein_aware_site_diagnostics_dataframe()
    assert via_result is not None
    via_result.iloc[
        0, via_result.columns.get_loc("protein_covariate_coefficient")
    ] = -99.0
    exported_third = result.protein_aware_site_diagnostics_dataframe()
    assert exported_third is not None
    assert exported_third.loc[
        exported_third.index[0],
        "protein_covariate_coefficient",
    ] == pytest.approx(0.25)


def test_protein_aware_scientific_equality_uses_diagnostics_content() -> None:
    left = _protein_aware_result()
    same = _protein_aware_result()
    different = _protein_aware_result(protein_coefficient=0.40)

    assert left.scientifically_equals(same)
    assert not left.scientifically_equals(different)


def test_protein_aware_payload_is_json_compatible_and_deterministic() -> None:
    result = _protein_aware_result()

    first_payload = result.to_payload()
    second_payload = result.to_payload()

    assert first_payload == second_payload
    assert "protein_aware_diagnostics" in first_payload
    json.dumps(first_payload, sort_keys=True)
    serialized_keys = json.dumps(first_payload, sort_keys=True)
    assert "protein_covariate_p_value" not in serialized_keys
    assert "protein_covariate_q_value" not in serialized_keys


def test_ordinary_payload_structure_has_no_protein_aware_diagnostics() -> None:
    result = DifferentialAnalysisWorkflow().run(_request())
    payload = result.to_payload()

    assert result.protein_aware_diagnostics is None
    assert result.protein_aware_site_diagnostics_dataframe() is None
    assert set(payload) == {
        "caveats",
        "diagnostics",
        "workflow_provenance",
        "policy_provenance",
        "empirical_bayes",
        "contrast_tables",
        "feature_eligibility",
    }
    assert "protein_aware_diagnostics" not in payload


def test_duplicate_correlation_result_assembly_has_no_protein_aware_diagnostics() -> (
    None
):
    result = DifferentialAnalysisWorkflow().run(
        _request(
            config=DifferentialAnalysisConfig(
                paired_design_policy=PAIRED_DESIGN_POLICY_DUPLICATE_CORRELATION
            ),
            include_blocks=True,
        )
    )

    assert result.protein_aware_diagnostics is None
    assert result.protein_aware_site_diagnostics_dataframe() is None
    assert "protein_aware_diagnostics" not in result.to_payload()
    assert result.diagnostics.model_type == "moderated_gls_duplicate_correlation"


def test_protein_aware_result_rejects_misaligned_diagnostics_index() -> None:
    result = _protein_aware_result()
    diagnostics = _required_protein_aware_diagnostics(result)
    per_site = diagnostics.per_site_diagnostics_dataframe()
    original_first_site = str(per_site.index[0])
    mismatched_first_site = f"{original_first_site}|different"
    per_site = per_site.rename(index={original_first_site: mismatched_first_site})
    per_site.loc[mismatched_first_site, "site_key"] = mismatched_first_site
    misaligned_diagnostics = _protein_aware_diagnostics_copy(
        diagnostics,
        per_site_diagnostics=per_site,
    )

    with pytest.raises(PhosPyInputError, match="per-site diagnostics index"):
        _result_with_protein_aware_diagnostics(
            result,
            protein_aware_diagnostics=misaligned_diagnostics,
        )


@pytest.mark.parametrize(
    ("field_name", "value", "message"),
    (
        ("claim_status", "validated", "claim_status"),
        ("method_id", "different_method", "method_id"),
        ("protein_covariate_centered", False, "protein_covariate_centered"),
        ("protein_covariate_standardized", True, "protein_covariate_standardized"),
        (
            "protein_covariate_centering_policy",
            "z_scored",
            "protein_covariate_centering_policy",
        ),
        ("protein_covariate_imputation_policy", "mean", "imputation_policy"),
        ("fallback_policy", "fallback_to_ordinary", "fallback_policy"),
    ),
)
def test_protein_aware_diagnostics_rejects_adr_inconsistent_contract_values(
    field_name: str,
    value: object,
    message: str,
) -> None:
    diagnostics = _required_protein_aware_diagnostics(_protein_aware_result())

    with pytest.raises(PhosPyInputError, match=message):
        _protein_aware_diagnostics_copy_with_override(
            diagnostics,
            field_name=field_name,
            value=value,
        )


@pytest.mark.parametrize(
    ("mismatch", "seam"),
    (
        ("sample_order", "request_sample_order"),
        ("base_design", "request_base_design"),
        ("base_contrasts", "request_base_contrasts"),
    ),
)
def test_protein_aware_assembly_rejects_resolved_inputs_from_different_request(
    mismatch: str,
    seam: str,
) -> None:
    interpreted = _interpreted_request()
    resolved_inputs = _resolved_inputs(interpreted, protein_coefficient=0.25)
    mismatched_resolved_inputs = _resolved_inputs_with_request_mismatch(
        resolved_inputs,
        mismatch=mismatch,
    )
    computation_result = _computation_result(
        interpreted,
        mismatched_resolved_inputs,
        protein_coefficient=0.25,
    )

    with pytest.raises(WorkflowBoundaryError) as exc_info:
        DifferentialResultAssembler().run_protein_aware(
            request=interpreted,
            resolved_inputs=mismatched_resolved_inputs,
            computation_result=computation_result,
            workflow_provenance={"stage": "pa006-test"},
        )

    assert exc_info.value.seam is not None
    assert exc_info.value.seam.endswith(seam)


def test_protein_aware_assembly_carries_computation_time_failure_diagnostics() -> None:
    interpreted = _interpreted_request()
    resolved_inputs = _resolved_inputs(
        interpreted,
        protein_coefficient=0.25,
        tested_positions=(0, 5),
    )
    failed_site_id = resolved_inputs.tested_site_ids[1]
    computation_result = _computation_result(
        interpreted,
        resolved_inputs,
        protein_coefficient=0.25,
        failure_site_id=failed_site_id,
        failure_total_protein_row_key="protein:RPS6",
    )

    result = DifferentialResultAssembler().run_protein_aware(
        request=interpreted,
        resolved_inputs=resolved_inputs,
        computation_result=computation_result,
        workflow_provenance={"stage": "pa006-test"},
    )

    table = result.table_for("B_vs_A")
    assert table.loc[failed_site_id, DIFFERENTIAL_RESULT_STATUS_COLUMN] == (
        DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_AUGMENTED_DESIGN_INVALID
    )
    assert table.loc[failed_site_id, DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN] == (
        DIFFERENTIAL_RESULT_REASON_PROTEIN_AUGMENTED_DESIGN_RANK_DEFICIENT
    )
    assert (
        table.loc[failed_site_id, ["logFC", "t", "P.Value", "adj.P.Val"]].isna().all()
    )

    per_site = result.protein_aware_site_diagnostics_dataframe()
    assert per_site is not None
    assert per_site.loc[
        failed_site_id, "protein_augmented_design_rank"
    ] == pytest.approx(2.0)
    assert per_site.loc[
        failed_site_id,
        "protein_augmented_design_condition_number",
    ] == pytest.approx(1.0e11)
    assert pd.isna(per_site.loc[failed_site_id, "protein_covariate_coefficient"])
    assert (
        per_site.loc[failed_site_id, "protein_aware_failure_message"]
        == "protein augmented design is rank deficient"
    )


def _protein_aware_result(
    *,
    protein_coefficient: float = 0.25,
) -> DifferentialAnalysisResult:
    interpreted = _interpreted_request()
    resolved_inputs = _resolved_inputs(
        interpreted,
        protein_coefficient=protein_coefficient,
    )
    computation_result = _computation_result(
        interpreted,
        resolved_inputs,
        protein_coefficient=protein_coefficient,
    )
    return DifferentialResultAssembler().run_protein_aware(
        request=interpreted,
        resolved_inputs=resolved_inputs,
        computation_result=computation_result,
        workflow_provenance={"stage": "pa006-test"},
    )


def _resolved_inputs(
    interpreted: InterpretedDifferentialAnalysisRequest,
    *,
    protein_coefficient: float,
    tested_positions: tuple[int, ...] = (0,),
) -> ProteinAwareDifferentialResolvedInputs:
    site_ids = tuple(str(value) for value in interpreted.result_identity_metadata.index)
    tested_site_ids = tuple(site_ids[position] for position in tested_positions)
    base_design = cast(DesignMatrix, interpreted.computation_request.design)
    base_contrasts = cast(ContrastMatrix, interpreted.computation_request.contrasts)
    matched_pair_rows = [
        {
            "site_key": site_ids[position],
            "protein_identifier": str(_PROTEINS[position]),
            "total_protein_row_key": str(_TOTAL_ROW_KEYS[position]),
        }
        for position in tested_positions
    ]
    matched_pairs = pd.DataFrame(matched_pair_rows)
    covariate_values_by_row_key = {
        "protein:MAPK14": [10.0, 10.2, 10.4, 11.2, 11.4, 11.6],
        "protein:RPS6": [12.0, 12.3, 12.6, 13.2, 13.5, 13.8],
    }
    tested_total_row_keys = tuple(
        dict.fromkeys(str(row["total_protein_row_key"]) for row in matched_pair_rows)
    )
    covariates = pd.DataFrame(
        [covariate_values_by_row_key[row_key] for row_key in tested_total_row_keys],
        index=pd.Index(tested_total_row_keys, name="total_protein_row_key"),
        columns=pd.Index(_SAMPLE_IDS, name="sample_id"),
    )
    computation_request = ProteinAwareDifferentialComputationRequest(
        phosphosite_matrix=interpreted.computation_request.matrix.loc[
            list(tested_site_ids),
            list(_SAMPLE_IDS),
        ],
        base_design=base_design,
        base_contrasts=base_contrasts,
        sample_order=_SAMPLE_IDS,
        matched_pairs=matched_pairs,
        resolved_protein_covariates=covariates,
        empirical_bayes=interpreted.computation_request.empirical_bayes,
        multiple_testing_method=interpreted.computation_request.multiple_testing_method,
        method_id=(
            DIFFERENTIAL_PROTEIN_AWARE_METHOD_PROTEIN_COVARIATE_ADJUSTED_MODERATED_LINEAR_MODEL_V1
        ),
    )
    feature_metadata = _protein_aware_feature_metadata(
        site_ids,
        protein_coefficient=protein_coefficient,
    )
    feature_metadata.loc[
        list(tested_site_ids),
        DIFFERENTIAL_RESULT_STATUS_COLUMN,
    ] = DIFFERENTIAL_RESULT_STATUS_TESTED
    feature_metadata.loc[
        list(tested_site_ids),
        DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN,
    ] = ""
    feature_metadata.loc[list(tested_site_ids), "protein_aware_tested"] = True
    feature_eligibility = DifferentialFeatureEligibilityInputs(
        feature_metadata=feature_metadata,
        result_status=feature_metadata[DIFFERENTIAL_RESULT_STATUS_COLUMN],
        testable_feature_ids=tested_site_ids,
        attach_to_result_tables=True,
    )
    return ProteinAwareDifferentialResolvedInputs(
        computation_request=computation_request,
        feature_eligibility_inputs=feature_eligibility,
        matched_pairs=matched_pairs,
        candidate_matched_pairs=pd.DataFrame(
            [
                {
                    "site_key": site_ids[0],
                    "protein_identifier": "MAPK14",
                    "total_protein_row_key": "protein:MAPK14",
                },
                {
                    "site_key": site_ids[3],
                    "protein_identifier": "AKT1",
                    "total_protein_row_key": "protein:AKT1",
                },
                {
                    "site_key": site_ids[4],
                    "protein_identifier": "PRKACA",
                    "total_protein_row_key": "protein:PRKACA",
                },
                {
                    "site_key": site_ids[5],
                    "protein_identifier": "RPS6",
                    "total_protein_row_key": "protein:RPS6",
                },
            ]
        ),
        resolved_protein_covariates=covariates,
        site_eligibility_metadata=feature_metadata,
        full_site_ids=site_ids,
        tested_site_ids=tested_site_ids,
        sample_order=_SAMPLE_IDS,
        base_design=base_design,
        base_contrasts=base_contrasts,
        method_id=(
            DIFFERENTIAL_PROTEIN_AWARE_METHOD_PROTEIN_COVARIATE_ADJUSTED_MODERATED_LINEAR_MODEL_V1
        ),
        preparation_policy="prepare_model_inputs",
        protein_mapping_policy="require_unambiguous",
        eligibility_counts=(
            ("total_site_count", 6),
            ("ordinary_testable_site_count", 6),
            ("protein_preparation_candidate_site_count", 4),
            ("protein_aware_tested_site_count", len(tested_site_ids)),
        ),
        status_counts=_status_counts(feature_metadata),
        reason_counts=_reason_counts(feature_metadata),
    )


def _computation_result(
    interpreted: InterpretedDifferentialAnalysisRequest,
    resolved_inputs: ProteinAwareDifferentialResolvedInputs,
    *,
    protein_coefficient: float,
    failure_site_id: str | None = None,
    failure_total_protein_row_key: str | None = None,
) -> ProteinAwareDifferentialComputationResult:
    failure_site_ids = () if failure_site_id is None else (failure_site_id,)
    tested_index = pd.Index(
        tuple(
            site_id
            for site_id in resolved_inputs.tested_site_ids
            if site_id not in failure_site_ids
        ),
        name="site_key",
    )
    prior_diagnostics = EmpiricalBayesPriorDiagnostics(
        method="standard",
        robust=False,
        trend=False,
        winsor_tail_p=(0.05, 0.1),
        base_prior_variance=0.5,
        base_prior_degrees_of_freedom=4.0,
        robust_outlier_count=0,
        robust_outlier_fraction=0.0,
        winsorized_low_count=0,
        winsorized_high_count=0,
        prior_variance=pd.Series([0.5], index=tested_index, name="prior_variance"),
        prior_degrees_of_freedom=pd.Series(
            [4.0],
            index=tested_index,
            name="prior_degrees_of_freedom",
        ),
    )
    return ProteinAwareDifferentialComputationResult(
        residual_variance=pd.Series([0.42], index=tested_index, name="s2"),
        posterior_residual_variance=pd.Series(
            [0.45],
            index=tested_index,
            name="s2_post",
        ),
        prior_residual_variance=pd.Series(
            [0.5],
            index=tested_index,
            name="s2_prior",
        ),
        prior_degrees_of_freedom_series_value=pd.Series(
            [4.0],
            index=tested_index,
            name="df_prior",
        ),
        prior_variance=0.5,
        prior_degrees_of_freedom=4.0,
        residual_degrees_of_freedom=3.0,
        empirical_bayes_method="standard",
        empirical_bayes_robust=False,
        empirical_bayes_trend=False,
        prior_diagnostics=prior_diagnostics,
        mean_variance_trend_diagnostics=None,
        contrast_tables={
            "B_vs_A": pd.DataFrame(
                {
                    "logFC": [1.25],
                    "t": [2.5],
                    "P.Value": [0.04],
                    "adj.P.Val": [0.04],
                },
                index=tested_index,
            )
        },
        protein_coefficient=pd.Series(
            [protein_coefficient],
            index=tested_index,
            name="protein_covariate_coefficient",
        ),
        site_diagnostics=pd.DataFrame(
            {
                "site_key": tested_index.tolist(),
                "protein_identifier": ["MAPK14"],
                "total_protein_row_key": ["protein:MAPK14"],
                "protein_coefficient": [protein_coefficient],
                "mean_intensity": [12.0],
            },
            index=tested_index,
        ),
        augmented_design_diagnostics=pd.DataFrame(
            {
                "total_protein_row_key": ["protein:MAPK14"],
                "sample_count": [6],
                "coefficient_count": [3],
                "rank": [3],
                "residual_degrees_of_freedom": [3.0],
                "condition_number": [11.0],
                "max_condition_number": [1.0e10],
                "decomposition_method": [
                    interpreted.design_decomposition.decomposition_method
                ],
                "solver": [interpreted.design_decomposition.solver],
                "column_scale_method": [
                    interpreted.design_decomposition.column_scale_method
                ],
                "rank_tolerance_policy": [
                    interpreted.design_decomposition.rank_tolerance_policy
                ],
                "rank_tolerance": [interpreted.design_decomposition.rank_tolerance],
            },
            index=pd.Index(["protein:MAPK14"], name="total_protein_row_key"),
        ),
        site_failure_diagnostics=_site_failure_diagnostics(
            failure_site_id=failure_site_id,
            total_protein_row_key=failure_total_protein_row_key,
        ),
        augmented_design_failure_diagnostics=_augmented_design_failure_diagnostics(
            total_protein_row_key=failure_total_protein_row_key,
        ),
        tested_site_ids=tested_index.tolist(),
        method_id=(
            DIFFERENTIAL_PROTEIN_AWARE_METHOD_PROTEIN_COVARIATE_ADJUSTED_MODERATED_LINEAR_MODEL_V1
        ),
    )


def _site_failure_diagnostics(
    *,
    failure_site_id: str | None,
    total_protein_row_key: str | None,
) -> pd.DataFrame | None:
    if failure_site_id is None or total_protein_row_key is None:
        return None
    return pd.DataFrame(
        {
            "site_key": [failure_site_id],
            "total_protein_row_key": [total_protein_row_key],
            DIFFERENTIAL_RESULT_STATUS_COLUMN: [
                DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_AUGMENTED_DESIGN_INVALID
            ],
            DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN: [
                DIFFERENTIAL_RESULT_REASON_PROTEIN_AUGMENTED_DESIGN_RANK_DEFICIENT
            ],
            "failure_message": ["protein augmented design is rank deficient"],
        },
        index=pd.Index([failure_site_id], name="site_key"),
    )


def _augmented_design_failure_diagnostics(
    *,
    total_protein_row_key: str | None,
) -> pd.DataFrame | None:
    if total_protein_row_key is None:
        return None
    return pd.DataFrame(
        {
            "total_protein_row_key": [total_protein_row_key],
            DIFFERENTIAL_RESULT_STATUS_COLUMN: [
                DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_AUGMENTED_DESIGN_INVALID
            ],
            DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN: [
                DIFFERENTIAL_RESULT_REASON_PROTEIN_AUGMENTED_DESIGN_RANK_DEFICIENT
            ],
            "sample_count": [6],
            "coefficient_count": [3],
            "rank": [2],
            "residual_degrees_of_freedom": [4.0],
            "condition_number": [1.0e11],
            "max_condition_number": [1.0e10],
            "protein_raw_mean": [12.9],
            "protein_raw_standard_deviation": [0.72],
            "protein_centered_variance": [0.52],
            "failure_message": ["protein augmented design is rank deficient"],
        },
        index=pd.Index([total_protein_row_key], name="total_protein_row_key"),
    )


def _protein_aware_feature_metadata(
    site_ids: tuple[str, ...],
    *,
    protein_coefficient: float,
) -> pd.DataFrame:
    statuses = [
        DIFFERENTIAL_RESULT_STATUS_TESTED,
        DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_PREPARATION_INELIGIBLE,
        DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_PREPARATION_INELIGIBLE,
        DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_COVARIATE_INVALID,
        DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_AUGMENTED_DESIGN_INVALID,
        DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_CONTRAST_NON_ESTIMABLE,
    ]
    reasons = [
        "",
        DIFFERENTIAL_RESULT_REASON_PROTEIN_PREPARATION_FALLBACK,
        DIFFERENTIAL_RESULT_REASON_PROTEIN_PREPARATION_EXCLUDED,
        DIFFERENTIAL_RESULT_REASON_PROTEIN_COVARIATE_ZERO_VARIANCE,
        DIFFERENTIAL_RESULT_REASON_PROTEIN_AUGMENTED_DESIGN_RANK_DEFICIENT,
        DIFFERENTIAL_RESULT_REASON_PROTEIN_CONTRAST_NON_ESTIMABLE,
    ]
    centered_variance = [0.53, np.nan, np.nan, 0.0, 0.70, 0.65]
    return pd.DataFrame(
        {
            "site_key": list(site_ids),
            "protein_aware_method_id": [
                DIFFERENTIAL_PROTEIN_AWARE_METHOD_PROTEIN_COVARIATE_ADJUSTED_MODERATED_LINEAR_MODEL_V1
            ]
            * len(site_ids),
            "protein_aware_centering_policy": ["mean_centered_no_standardization"]
            * len(site_ids),
            "protein_aware_candidate": [True, False, False, True, True, True],
            "protein_aware_tested": [True, False, False, False, False, False],
            "protein_aware_preparation_eligibility": [
                "eligible_for_protein_aware_preparation",
                "fallback_to_phospho_only",
                "excluded_from_protein_aware_preparation",
                "eligible_for_protein_aware_preparation",
                "eligible_for_protein_aware_preparation",
                "eligible_for_protein_aware_preparation",
            ],
            "protein_aware_preparation_reasons": [
                ("matched_protein_available",),
                ("missing_total_protein_row",),
                ("ambiguous_protein_mapping",),
                ("matched_protein_available",),
                ("matched_protein_available",),
                ("matched_protein_available",),
            ],
            "protein_aware_mapping_status": [
                "matched",
                "missing_total_protein_row",
                "ambiguous",
                "matched",
                "matched",
                "matched",
            ],
            "protein_identifier": list(_PROTEINS),
            "total_protein_row_key": list(_TOTAL_ROW_KEYS),
            "protein_covariate_raw_mean": [10.8, np.nan, np.nan, 10.0, 11.0, 12.0],
            "protein_covariate_raw_standard_deviation": [
                0.72,
                np.nan,
                np.nan,
                0.0,
                0.84,
                0.81,
            ],
            "protein_covariate_centered_variance": centered_variance,
            "protein_augmented_design_rank": [3, np.nan, np.nan, 3, 2, 3],
            "protein_augmented_design_residual_degrees_of_freedom": [
                3.0,
                np.nan,
                np.nan,
                3.0,
                4.0,
                3.0,
            ],
            "protein_augmented_design_condition_number": [
                11.0,
                np.nan,
                np.nan,
                8.0,
                1.0e11,
                9.5,
            ],
            "protein_augmented_design_max_condition_number": [1.0e10] * len(site_ids),
            "protein_aware_failure_message": [
                "",
                "protein-aware preparation reported fallback",
                "protein-aware preparation excluded this phosphosite",
                "protein covariate has zero centered variance",
                "protein augmented design is rank deficient",
                "protein-aware contrast is not estimable",
            ],
            DIFFERENTIAL_RESULT_STATUS_COLUMN: statuses,
            DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN: reasons,
            "protein_covariate_coefficient": [
                protein_coefficient,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
            ],
        },
        index=pd.Index(site_ids, name="site_key"),
    )


def _required_protein_aware_diagnostics(
    result: DifferentialAnalysisResult,
) -> ProteinAwareDifferentialDiagnostics:
    diagnostics = result.protein_aware_diagnostics
    assert diagnostics is not None
    return diagnostics


def _protein_aware_diagnostics_copy(
    diagnostics: ProteinAwareDifferentialDiagnostics,
    *,
    method_id: str | None = None,
    claim_status: str | None = None,
    protein_covariate_centered: bool | None = None,
    protein_covariate_standardized: bool | None = None,
    protein_covariate_centering_policy: str | None = None,
    protein_covariate_imputation_policy: str | None = None,
    fallback_policy: str | None = None,
    per_site_diagnostics: pd.DataFrame | None = None,
) -> ProteinAwareDifferentialDiagnostics:
    return ProteinAwareDifferentialDiagnostics(
        method_id=diagnostics.method_id if method_id is None else method_id,
        claim_status=diagnostics.claim_status if claim_status is None else claim_status,
        preparation_policy=diagnostics.preparation_policy,
        protein_mapping_policy=diagnostics.protein_mapping_policy,
        protein_covariate_centered=(
            diagnostics.protein_covariate_centered
            if protein_covariate_centered is None
            else protein_covariate_centered
        ),
        protein_covariate_standardized=(
            diagnostics.protein_covariate_standardized
            if protein_covariate_standardized is None
            else protein_covariate_standardized
        ),
        protein_covariate_centering_policy=(
            diagnostics.protein_covariate_centering_policy
            if protein_covariate_centering_policy is None
            else protein_covariate_centering_policy
        ),
        protein_covariate_imputation_policy=(
            diagnostics.protein_covariate_imputation_policy
            if protein_covariate_imputation_policy is None
            else protein_covariate_imputation_policy
        ),
        fallback_policy=diagnostics.fallback_policy
        if fallback_policy is None
        else fallback_policy,
        total_site_count=diagnostics.total_site_count,
        ordinary_eligible_site_count=diagnostics.ordinary_eligible_site_count,
        protein_preparation_eligible_site_count=(
            diagnostics.protein_preparation_eligible_site_count
        ),
        tested_site_count=diagnostics.tested_site_count,
        withheld_site_count=diagnostics.withheld_site_count,
        status_counts=diagnostics.status_counts,
        reason_counts=diagnostics.reason_counts,
        execution_sample_order=diagnostics.execution_sample_order,
        base_design_rank=diagnostics.base_design_rank,
        base_residual_degrees_of_freedom=(diagnostics.base_residual_degrees_of_freedom),
        expected_augmented_rank=diagnostics.expected_augmented_rank,
        common_augmented_rank=diagnostics.common_augmented_rank,
        common_augmented_residual_degrees_of_freedom=(
            diagnostics.common_augmented_residual_degrees_of_freedom
        ),
        distinct_matched_protein_row_count=(
            diagnostics.distinct_matched_protein_row_count
        ),
        fitted_protein_row_count=diagnostics.fitted_protein_row_count,
        condition_number_summary_scope=diagnostics.condition_number_summary_scope,
        condition_number_summary_statistic=(
            diagnostics.condition_number_summary_statistic
        ),
        min_augmented_condition_number=diagnostics.min_augmented_condition_number,
        median_augmented_condition_number=(
            diagnostics.median_augmented_condition_number
        ),
        max_augmented_condition_number=diagnostics.max_augmented_condition_number,
        per_site_diagnostics=(
            diagnostics.per_site_diagnostics_dataframe()
            if per_site_diagnostics is None
            else per_site_diagnostics
        ),
    )


def _protein_aware_diagnostics_copy_with_override(
    diagnostics: ProteinAwareDifferentialDiagnostics,
    *,
    field_name: str,
    value: object,
) -> ProteinAwareDifferentialDiagnostics:
    if field_name == "claim_status":
        return _protein_aware_diagnostics_copy(
            diagnostics,
            claim_status=cast(str, value),
        )
    if field_name == "method_id":
        return _protein_aware_diagnostics_copy(
            diagnostics,
            method_id=cast(str, value),
        )
    if field_name == "protein_covariate_centered":
        return _protein_aware_diagnostics_copy(
            diagnostics,
            protein_covariate_centered=cast(bool, value),
        )
    if field_name == "protein_covariate_standardized":
        return _protein_aware_diagnostics_copy(
            diagnostics,
            protein_covariate_standardized=cast(bool, value),
        )
    if field_name == "protein_covariate_centering_policy":
        return _protein_aware_diagnostics_copy(
            diagnostics,
            protein_covariate_centering_policy=cast(str, value),
        )
    if field_name == "protein_covariate_imputation_policy":
        return _protein_aware_diagnostics_copy(
            diagnostics,
            protein_covariate_imputation_policy=cast(str, value),
        )
    if field_name == "fallback_policy":
        return _protein_aware_diagnostics_copy(
            diagnostics,
            fallback_policy=cast(str, value),
        )
    raise AssertionError(f"unsupported diagnostics override: {field_name}")


def _result_with_protein_aware_diagnostics(
    result: DifferentialAnalysisResult,
    *,
    protein_aware_diagnostics: ProteinAwareDifferentialDiagnostics,
) -> DifferentialAnalysisResult:
    return DifferentialAnalysisResult.from_trusted_owned(
        residual_variance=result.residual_variance_series(),
        posterior_residual_variance=result.posterior_residual_variance_series(),
        prior_residual_variance=result.prior_residual_variance_series(),
        prior_degrees_of_freedom_series_value=result.prior_degrees_of_freedom_series(),
        prior_variance=result.prior_variance,
        prior_degrees_of_freedom=result.prior_degrees_of_freedom,
        residual_degrees_of_freedom=result.residual_degrees_of_freedom,
        empirical_bayes_method=result.empirical_bayes_method,
        empirical_bayes_robust=result.empirical_bayes_robust,
        empirical_bayes_trend=result.empirical_bayes_trend,
        prior_diagnostics=result.prior_diagnostics,
        mean_variance_trend_diagnostics=result.mean_variance_trend_diagnostics,
        contrast_tables=result.contrast_tables,
        diagnostics=result.diagnostics,
        policy_provenance=result.policy_provenance,
        workflow_provenance=result.workflow_provenance,
        caveats=result.caveats,
        input_dataset_preprocessing_report=result.input_dataset_preprocessing_report,
        feature_eligibility=result.feature_eligibility,
        protein_aware_diagnostics=protein_aware_diagnostics,
    )


def _resolved_inputs_with_request_mismatch(
    resolved_inputs: ProteinAwareDifferentialResolvedInputs,
    *,
    mismatch: str,
) -> ProteinAwareDifferentialResolvedInputs:
    sample_order = resolved_inputs.sample_order
    base_design = resolved_inputs.base_design
    base_contrasts = resolved_inputs.base_contrasts
    phosphosite_matrix = resolved_inputs.computation_request.phosphosite_matrix
    resolved_protein_covariates = resolved_inputs.resolved_protein_covariates
    if mismatch == "sample_order":
        sample_order = tuple(reversed(sample_order))
        base_design = DesignMatrix(base_design.frame.loc[list(sample_order)])
        phosphosite_matrix = phosphosite_matrix.loc[:, list(sample_order)]
        resolved_protein_covariates = resolved_protein_covariates.loc[
            :,
            list(sample_order),
        ]
    elif mismatch == "base_design":
        frame = base_design.to_dataframe()
        frame.iloc[0, 0] = float(frame.iloc[0, 0]) + 0.125
        base_design = DesignMatrix(frame)
    elif mismatch == "base_contrasts":
        frame = base_contrasts.to_dataframe()
        frame.iloc[0, 0] = float(frame.iloc[0, 0]) + 0.125
        base_contrasts = ContrastMatrix(frame)
    else:
        raise AssertionError(f"unsupported mismatch case: {mismatch}")
    computation_request = ProteinAwareDifferentialComputationRequest(
        phosphosite_matrix=phosphosite_matrix,
        base_design=base_design,
        base_contrasts=base_contrasts,
        sample_order=sample_order,
        matched_pairs=resolved_inputs.matched_pairs,
        resolved_protein_covariates=resolved_protein_covariates,
        empirical_bayes=resolved_inputs.computation_request.empirical_bayes,
        multiple_testing_method=(
            resolved_inputs.computation_request.multiple_testing_method
        ),
        method_id=resolved_inputs.method_id,
    )
    return replace(
        resolved_inputs,
        computation_request=computation_request,
        resolved_protein_covariates=resolved_protein_covariates,
        sample_order=sample_order,
        base_design=base_design,
        base_contrasts=base_contrasts,
    )


def _status_counts(feature_metadata: pd.DataFrame) -> tuple[tuple[str, int], ...]:
    values = [
        str(value) for value in feature_metadata[DIFFERENTIAL_RESULT_STATUS_COLUMN]
    ]
    return tuple((value, values.count(value)) for value in dict.fromkeys(values))


def _reason_counts(feature_metadata: pd.DataFrame) -> tuple[tuple[str, int], ...]:
    values = [
        str(value).strip()
        for value in feature_metadata[DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN]
        if str(value).strip()
    ]
    return tuple((value, values.count(value)) for value in dict.fromkeys(values))


def _interpreted_request(
    *,
    config: DifferentialAnalysisConfig | None = None,
    include_blocks: bool = False,
) -> InterpretedDifferentialAnalysisRequest:
    return DifferentialAnalysisInterpreter().run(
        DifferentialAnalysisValidator().run(
            _request(config=config, include_blocks=include_blocks)
        )
    )


def _request(
    *,
    config: DifferentialAnalysisConfig | None = None,
    include_blocks: bool = False,
) -> DifferentialAnalysisRequest:
    return DifferentialAnalysisRequest(
        dataset=_dataset(),
        design=ExperimentalDesign(
            samples=_sample_design(include_blocks=include_blocks)
        ),
        contrasts=(
            Contrast(
                name="B_vs_A",
                numerator_condition="B",
                denominator_condition="A",
            ),
        ),
        config=DifferentialAnalysisConfig() if config is None else config,
    )


def _sample_design(*, include_blocks: bool) -> tuple[SampleDesignRecord, ...]:
    block_ids = ("donor_1", "donor_2", "donor_3", "donor_1", "donor_2", "donor_3")
    conditions = ("A", "A", "A", "B", "B", "B")
    records: list[SampleDesignRecord] = []
    for sample_id, condition, block_id in zip(
        _SAMPLE_IDS,
        conditions,
        block_ids,
        strict=True,
    ):
        records.append(
            SampleDesignRecord(
                sample_id=sample_id,
                condition=condition,
                biological_replicate_id=block_id,
                block_id=block_id if include_blocks else None,
            )
        )
    return tuple(records)


def _dataset():
    site_index = _site_index()
    phospho = pd.DataFrame(
        [
            [10.0, 10.4, 10.8, 11.2, 11.4, 11.8],
            [9.5, 9.7, 9.9, 10.1, 10.4, 10.5],
            [8.0, 8.2, 8.5, 9.0, 9.2, 9.5],
            [12.0, 12.1, 12.3, 11.9, 12.2, 12.4],
            [7.0, 7.2, 7.4, 8.1, 8.3, 8.6],
            [13.0, 12.8, 13.2, 13.5, 13.7, 13.9],
        ],
        index=site_index,
        columns=pd.Index(_SAMPLE_IDS, name="sample_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "site_key": site_index.tolist(),
            "display_id": [
                f"{protein};{site};"
                for protein, site in zip(_PROTEINS, _SITES, strict=True)
            ],
            **site_key_context_columns(site_index),
            "gene_symbol": list(_PROTEINS),
            "site": list(_SITES),
            "site_sequence": [("A" * 15) + site[0] + ("A" * 15) for site in _SITES],
            "protein_id": list(_PROTEINS),
        },
        index=site_index.copy(),
    )
    return trusted_analysis_ready_dataset_from_tables(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        intensity_scale_state=supported_log2_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_log2_processing_state(has_total_matrix=False),
    )


def _site_index() -> pd.Index:
    return protein_site_key_index(
        protein_identifiers=list(_PROTEINS),
        sites=list(_SITES),
    )
