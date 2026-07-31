from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy import stats

import phospy.science.differential as differential_public
import phospy.science.differential.aggregation as aggregation_public
from phospy import AnalysisReadyDatasetBuilder, DifferentialAnalysisWorkflow
from phospy.api import (
    Contrast,
    DatasetBuildRequest,
    DatasetIntensityTransformConfig,
    DatasetPreprocessingConfig,
    DifferentialAnalysisRequest,
    ExperimentalDesign,
    Organism,
    SampleDesignRecord,
)
from phospy.api.requests import (
    DATASET_MULTI_SITE_POLICY_KEEP_JOINT,
    DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE,
)
from phospy.errors import PhosPyInputError
from phospy.science.differential.aggregation import (
    PEPTIDE_DIFFERENTIAL_CONSISTENCY_POLICY,
    PEPTIDE_DIFFERENTIAL_MAPPING_WEIGHT_POLICY_REJECT,
    PEPTIDE_DIFFERENTIAL_STATISTIC_DISTRIBUTION_MODERATED_T,
    PEPTIDE_TO_SITE_AGGREGATION_LEVEL_SAMPLE_INTENSITY,
    PEPTIDE_TO_SITE_AGGREGATION_LEVEL_SINGLE_ESTIMATE,
    PEPTIDE_TO_SITE_AGGREGATION_SUPPORT_STATUS,
    PEPTIDE_TO_SITE_DEPENDENCE_POLICY_INDEPENDENT_SOURCES,
    PEPTIDE_TO_SITE_DEPENDENCE_POLICY_SAME_EXPERIMENT_CORRELATED,
    PEPTIDE_TO_SITE_FIXED_EFFECT_MIN_ASYMPTOTIC_MODERATED_DF,
    PEPTIDE_TO_SITE_MAPPING_POLICY_EXPLICIT_SITE_ID,
    PEPTIDE_TO_SITE_MAPPING_POLICY_KEEP_JOINT,
    PEPTIDE_TO_SITE_UNCERTAINTY_METHOD_FIXED_EFFECT_INVERSE_VARIANCE,
    PEPTIDE_TO_SITE_UNCERTAINTY_METHOD_SINGLE_ESTIMATE_PASSTHROUGH,
    PEPTIDE_TO_SITE_UNCERTAINTY_METHOD_STOUFFER_SIGNED_P,
    PeptideDifferentialEstimateTable,
    PeptideToSiteAggregationConfig,
    PeptideToSiteAggregator,
    signed_z_from_t_statistic,
)
from phospy.science.differential.aggregation.experimental import (
    EXPERIMENTAL_INTERNAL_API,
    EXPERIMENTAL_INTERNAL_REASON,
)
from phospy.science.statistics.multiple_testing import adjust_p_values

ROOT = Path(__file__).resolve().parents[2]


def _estimate_frame(
    rows: tuple[dict[str, object], ...] | None = None,
) -> pd.DataFrame:
    if rows is None:
        rows = (
            _estimate_row(
                peptide_id="pep_1",
                site_id="MAPK1;S10;",
                effect=1.25,
                standard_error=0.5,
                statistic=2.5,
                moderated_degrees_of_freedom=4.0,
                source_experiment_id="run_1",
            ),
        )
    return pd.DataFrame(rows)


def _estimate_row(
    *,
    peptide_id: str,
    site_id: str,
    effect: float,
    standard_error: float,
    statistic: float,
    moderated_degrees_of_freedom: float,
    source_experiment_id: str,
    p_value: float | None = None,
    residual_degrees_of_freedom: float = 3.0,
    contrast_id: str = "B_vs_A",
    contrast_orientation: str = "B_minus_A",
    effect_scale: str = "log2_fold_change",
    effect_unit: str = "log2_ratio",
    model_estimator_id: str = "limma_moderated_ols",
    statistic_distribution: str = PEPTIDE_DIFFERENTIAL_STATISTIC_DISTRIBUTION_MODERATED_T,
    uncertainty_method_version: str = "limma_ebayes_moderated_t_v1",
    dependence_policy: str = PEPTIDE_TO_SITE_DEPENDENCE_POLICY_INDEPENDENT_SOURCES,
    mapping_policy: str = PEPTIDE_TO_SITE_MAPPING_POLICY_EXPLICIT_SITE_ID,
    mapping_uncertainty: bool = False,
) -> dict[str, object]:
    return {
        "site_id": site_id,
        "peptide_id": peptide_id,
        "contrast_id": contrast_id,
        "contrast_orientation": contrast_orientation,
        "effect_scale": effect_scale,
        "effect_unit": effect_unit,
        "model_estimator_id": model_estimator_id,
        "statistic_distribution": statistic_distribution,
        "uncertainty_method_version": uncertainty_method_version,
        "effect": effect,
        "standard_error": standard_error,
        "statistic": statistic,
        "p_value": (
            float(p_value)
            if p_value is not None
            else float(
                2.0
                * stats.t.sf(
                    abs(statistic),
                    df=moderated_degrees_of_freedom,
                )
            )
        ),
        "residual_degrees_of_freedom": residual_degrees_of_freedom,
        "moderated_degrees_of_freedom": moderated_degrees_of_freedom,
        "source_experiment_id": source_experiment_id,
        "dependence_policy": dependence_policy,
        "peptide_to_site_mapping_policy": mapping_policy,
        "mapping_uncertainty": mapping_uncertainty,
    }


def _aggregate(
    frame: pd.DataFrame,
    *,
    config: PeptideToSiteAggregationConfig | None = None,
):
    return PeptideToSiteAggregator().run(
        PeptideDifferentialEstimateTable(frame),
        config=config,
        contrast_name="B_vs_A",
    )


def test_peptide_to_site_typed_route_is_supported_public_api() -> None:
    expected = {
        "PeptideDifferentialEstimateTable",
        "PeptideToSiteAggregationConfig",
        "PeptideToSiteAggregationResult",
        "PeptideToSiteAggregator",
        "PEPTIDE_DIFFERENTIAL_CONSISTENCY_POLICY",
        "PEPTIDE_DIFFERENTIAL_MAPPING_WEIGHT_POLICY_REJECT",
        "PEPTIDE_DIFFERENTIAL_STATISTIC_DISTRIBUTION_MODERATED_T",
        "PEPTIDE_TO_SITE_MAPPING_POLICY_EXCLUDE_FROM_STATISTICAL_MODEL",
        "SUPPORTED_PEPTIDE_TO_SITE_MAPPING_POLICIES",
        "PEPTIDE_TO_SITE_UNCERTAINTY_METHOD_STOUFFER_SIGNED_P",
        "signed_z_from_two_sided_p_value",
    }

    assert expected <= set(aggregation_public.__all__)
    assert expected <= set(differential_public.__all__)
    assert PEPTIDE_TO_SITE_AGGREGATION_SUPPORT_STATUS == (
        "supported_typed_estimate_combination_v2"
    )
    assert EXPERIMENTAL_INTERNAL_API is False
    assert "replaced by a supported typed estimate-combination contract" in (
        EXPERIMENTAL_INTERNAL_REASON
    )


def test_typed_estimate_model_rejects_unknown_mapping_policy() -> None:
    frame = _estimate_frame()
    frame.loc[:, "peptide_to_site_mapping_policy"] = "implicit_best_guess"

    with pytest.raises(PhosPyInputError, match="peptide_to_site_mapping_policy"):
        PeptideDifferentialEstimateTable(frame)


def test_typed_estimate_model_rejects_effect_statistic_sign_mismatch() -> None:
    frame = _estimate_frame(
        (
            _estimate_row(
                peptide_id="pep_bad",
                site_id="MAPK1;S10;",
                effect=2.0,
                standard_error=0.5,
                statistic=-4.0,
                moderated_degrees_of_freedom=4.0,
                source_experiment_id="run_1",
            ),
        )
    )

    with pytest.raises(PhosPyInputError, match="effect/statistic signs"):
        PeptideDifferentialEstimateTable(frame)


def test_typed_estimate_model_rejects_effect_se_statistic_mismatch() -> None:
    frame = _estimate_frame(
        (
            _estimate_row(
                peptide_id="pep_bad",
                site_id="MAPK1;S10;",
                effect=2.0,
                standard_error=0.5,
                statistic=3.0,
                moderated_degrees_of_freedom=4.0,
                source_experiment_id="run_1",
            ),
        )
    )

    with pytest.raises(PhosPyInputError, match="effect / standard_error"):
        PeptideDifferentialEstimateTable(frame)


def test_typed_estimate_model_rejects_p_value_statistic_df_mismatch() -> None:
    frame = _estimate_frame(
        (
            _estimate_row(
                peptide_id="pep_bad",
                site_id="MAPK1;S10;",
                effect=2.0,
                standard_error=0.5,
                statistic=4.0,
                p_value=0.9,
                moderated_degrees_of_freedom=4.0,
                source_experiment_id="run_1",
            ),
        )
    )

    with pytest.raises(PhosPyInputError, match="two-sided moderated_t"):
        PeptideDifferentialEstimateTable(frame)


@pytest.mark.parametrize(
    ("field_name", "replacement", "match"),
    (
        ("contrast_id", "A_vs_B", "contrast_id"),
        ("contrast_orientation", "A_minus_B", "contrast_orientation"),
        ("effect_scale", "natural_log_fold_change", "effect_scale"),
        ("effect_unit", "natural_log_ratio", "effect_unit"),
        ("model_estimator_id", "wald_glm", "model_estimator_id"),
        ("statistic_distribution", "wald_z", "statistic_distribution"),
        ("uncertainty_method_version", "other_moderation_v2", "uncertainty"),
    ),
)
def test_typed_estimate_model_rejects_mixed_estimate_identity(
    field_name: str,
    replacement: str,
    match: str,
) -> None:
    first = _estimate_row(
        peptide_id="pep_1",
        site_id="MAPK1;S10;",
        effect=1.0,
        standard_error=0.5,
        statistic=2.0,
        moderated_degrees_of_freedom=4.0,
        source_experiment_id="run_1",
    )
    second = _estimate_row(
        peptide_id="pep_2",
        site_id="MAPK1;S10;",
        effect=1.0,
        standard_error=0.5,
        statistic=2.0,
        moderated_degrees_of_freedom=4.0,
        source_experiment_id="run_2",
    )
    second[field_name] = replacement
    frame = _estimate_frame((first, second))

    with pytest.raises(PhosPyInputError, match=match):
        PeptideDifferentialEstimateTable(frame)


def test_mapping_weight_is_rejected_in_posthoc_estimate_lane() -> None:
    frame = _estimate_frame()
    frame.loc[:, "mapping_weight"] = 0.5

    with pytest.raises(PhosPyInputError, match="mapping_weight"):
        PeptideDifferentialEstimateTable(frame)


def test_aggregation_contrast_name_must_match_input_contrast_identity() -> None:
    estimates = PeptideDifferentialEstimateTable(_estimate_frame())

    with pytest.raises(PhosPyInputError, match="contrast_name"):
        PeptideToSiteAggregator().run(
            estimates,
            config=PeptideToSiteAggregationConfig(),
            contrast_name="A_vs_B",
        )


def test_aggregation_default_contrast_name_uses_input_contrast_identity() -> None:
    estimates = PeptideDifferentialEstimateTable(_estimate_frame())

    result = PeptideToSiteAggregator().run(
        estimates,
        config=PeptideToSiteAggregationConfig(),
    )

    assert result.contrast_name == "B_vs_A"


@pytest.mark.parametrize("degrees_of_freedom", (2.0, 3.0, 4.0, 10.0))
def test_finite_df_t_to_z_conversion_uses_signed_p_value(
    degrees_of_freedom: float,
) -> None:
    t_statistic = 2.0
    expected_p = float(2.0 * stats.t.sf(abs(t_statistic), df=degrees_of_freedom))
    expected_z = float(stats.norm.isf(expected_p / 2.0))

    observed = signed_z_from_t_statistic(t_statistic, degrees_of_freedom)

    assert observed == pytest.approx(expected_z, rel=1e-12, abs=0.0)
    assert observed != pytest.approx(t_statistic, rel=0.0, abs=1e-6)


def test_large_df_t_to_z_conversion_approaches_t_limit() -> None:
    observed = signed_z_from_t_statistic(2.0, 1_000_000.0)

    assert observed == pytest.approx(2.0, rel=1e-5, abs=1e-5)


def test_signed_t_to_z_conversion_preserves_direction() -> None:
    positive = signed_z_from_t_statistic(2.0, 4.0)
    negative = signed_z_from_t_statistic(-2.0, 4.0)

    assert positive > 0
    assert negative == pytest.approx(-positive, rel=1e-12, abs=0.0)


def test_single_estimate_reproduces_original_effect_statistic_and_p_value() -> None:
    input_frame = _estimate_frame()
    result = _aggregate(input_frame)
    table = result.to_dataframe()
    row = table.loc["MAPK1;S10;"]
    input_row = input_frame.iloc[0]
    normal_p = float(2.0 * stats.norm.sf(abs(float(input_row["statistic"]))))

    assert float(row["logFC"]) == float(input_row["effect"])
    assert float(row["standard_error"]) == float(input_row["standard_error"])
    assert float(row["uncertainty_statistic"]) == float(input_row["statistic"])
    assert float(row["P.Value"]) == float(input_row["p_value"])
    assert float(row["P.Value"]) != pytest.approx(normal_p)
    assert float(row["adj.P.Val"]) == float(input_row["p_value"])
    assert row["contrast_id"] == "B_vs_A"
    assert row["contrast_orientation"] == "B_minus_A"
    assert row["effect_scale"] == "log2_fold_change"
    assert row["effect_unit"] == "log2_ratio"
    assert row["model_estimator_id"] == "limma_moderated_ols"
    assert row["input_statistic_distribution"] == "moderated_t"
    assert row["input_uncertainty_method_version"] == "limma_ebayes_moderated_t_v1"
    assert row["aggregation_level"] == PEPTIDE_TO_SITE_AGGREGATION_LEVEL_SINGLE_ESTIMATE
    assert (
        row["uncertainty_method"]
        == PEPTIDE_TO_SITE_UNCERTAINTY_METHOD_SINGLE_ESTIMATE_PASSTHROUGH
    )
    assert row["statistic_distribution"] == "moderated_t"


def test_zero_effect_and_zero_statistic_pass_through_explicitly() -> None:
    frame = _estimate_frame(
        (
            _estimate_row(
                peptide_id="pep_zero",
                site_id="MAPK1;S10;",
                effect=0.0,
                standard_error=0.5,
                statistic=0.0,
                moderated_degrees_of_freedom=3.0,
                source_experiment_id="run_1",
            ),
        )
    )
    result = _aggregate(frame)
    row = result.to_dataframe().loc["MAPK1;S10;"]

    assert float(row["logFC"]) == 0.0
    assert float(row["uncertainty_statistic"]) == 0.0
    assert float(row["P.Value"]) == 1.0


def test_zero_effect_and_zero_statistic_require_unit_p_value() -> None:
    frame = _estimate_frame(
        (
            _estimate_row(
                peptide_id="pep_zero_bad",
                site_id="MAPK1;S10;",
                effect=0.0,
                standard_error=0.5,
                statistic=0.0,
                p_value=0.9,
                moderated_degrees_of_freedom=3.0,
                source_experiment_id="run_1",
            ),
        )
    )

    with pytest.raises(PhosPyInputError, match="zero-effect moderated_t"):
        PeptideDifferentialEstimateTable(frame)


def test_zero_effect_and_zero_statistic_combine_to_unit_p_value() -> None:
    frame = _estimate_frame(
        (
            _estimate_row(
                peptide_id="pep_zero_1",
                site_id="MAPK1;S10;",
                effect=0.0,
                standard_error=0.5,
                statistic=0.0,
                moderated_degrees_of_freedom=3.0,
                source_experiment_id="run_1",
            ),
            _estimate_row(
                peptide_id="pep_zero_2",
                site_id="MAPK1;S10;",
                effect=0.0,
                standard_error=0.7,
                statistic=0.0,
                moderated_degrees_of_freedom=10.0,
                source_experiment_id="run_2",
            ),
        )
    )
    result = _aggregate(frame)
    row = result.to_dataframe().loc["MAPK1;S10;"]

    assert float(row["uncertainty_statistic"]) == 0.0
    assert float(row["P.Value"]) == 1.0
    assert row["statistic_distribution"] == "standard_normal_z"


def test_conflicting_independent_effects_reduce_site_level_significance() -> None:
    frame = _estimate_frame(
        (
            _estimate_row(
                peptide_id="pep_up",
                site_id="MAPK1;S10;",
                effect=1.0,
                standard_error=0.25,
                statistic=4.0,
                moderated_degrees_of_freedom=10.0,
                source_experiment_id="run_1",
            ),
            _estimate_row(
                peptide_id="pep_down",
                site_id="MAPK1;S10;",
                effect=-1.0,
                standard_error=0.25,
                statistic=-4.0,
                moderated_degrees_of_freedom=10.0,
                source_experiment_id="run_2",
            ),
        )
    )

    result = _aggregate(frame)
    row = result.to_dataframe().loc["MAPK1;S10;"]

    assert abs(float(row["logFC"])) < 1e-12
    assert abs(float(row["uncertainty_statistic"])) < 1e-12
    assert float(row["P.Value"]) == pytest.approx(1.0)


def test_same_experiment_peptide_estimates_are_rejected() -> None:
    frame = _estimate_frame(
        (
            _estimate_row(
                peptide_id="pep_1",
                site_id="MAPK1;S10;",
                effect=1.0,
                standard_error=0.5,
                statistic=2.0,
                moderated_degrees_of_freedom=4.0,
                source_experiment_id="run_1",
            ),
            _estimate_row(
                peptide_id="pep_2",
                site_id="MAPK1;S10;",
                effect=1.1,
                standard_error=0.5,
                statistic=2.2,
                moderated_degrees_of_freedom=4.0,
                source_experiment_id="run_1",
            ),
        )
    )

    with pytest.raises(PhosPyInputError, match="same-experiment peptide estimates"):
        _aggregate(frame)


def test_same_experiment_rejection_precedes_minimum_evidence_attrition() -> None:
    frame = _estimate_frame(
        (
            _estimate_row(
                peptide_id="pep_1",
                site_id="MAPK1;S10;",
                effect=1.0,
                standard_error=0.5,
                statistic=2.0,
                moderated_degrees_of_freedom=4.0,
                source_experiment_id="run_1",
            ),
            _estimate_row(
                peptide_id="pep_2",
                site_id="MAPK1;S10;",
                effect=1.1,
                standard_error=0.5,
                statistic=2.2,
                moderated_degrees_of_freedom=4.0,
                source_experiment_id="run_1",
            ),
        )
    )

    with pytest.raises(PhosPyInputError, match="same-experiment peptide estimates"):
        _aggregate(
            frame, config=PeptideToSiteAggregationConfig(min_estimates_per_site=3)
        )


def test_declared_correlated_same_sample_estimates_are_rejected() -> None:
    frame = _estimate_frame(
        (
            _estimate_row(
                peptide_id="pep_1",
                site_id="MAPK1;S10;",
                effect=1.0,
                standard_error=0.5,
                statistic=2.0,
                moderated_degrees_of_freedom=4.0,
                source_experiment_id="run_1",
                dependence_policy=(
                    PEPTIDE_TO_SITE_DEPENDENCE_POLICY_SAME_EXPERIMENT_CORRELATED
                ),
            ),
            _estimate_row(
                peptide_id="pep_2",
                site_id="MAPK1;S10;",
                effect=1.1,
                standard_error=0.5,
                statistic=2.2,
                moderated_degrees_of_freedom=4.0,
                source_experiment_id="run_2",
                dependence_policy=(
                    PEPTIDE_TO_SITE_DEPENDENCE_POLICY_SAME_EXPERIMENT_CORRELATED
                ),
            ),
        )
    )

    with pytest.raises(PhosPyInputError, match="correlated peptide estimates"):
        _aggregate(frame)


def test_independent_estimates_can_be_combined_with_signed_p_stouffer() -> None:
    frame = _estimate_frame(
        (
            _estimate_row(
                peptide_id="pep_run_1",
                site_id="MAPK1;S10;",
                effect=1.0,
                standard_error=0.5,
                statistic=2.0,
                moderated_degrees_of_freedom=4.0,
                source_experiment_id="run_1",
            ),
            _estimate_row(
                peptide_id="pep_run_2",
                site_id="MAPK1;S10;",
                effect=1.0,
                standard_error=0.5,
                statistic=2.0,
                moderated_degrees_of_freedom=4.0,
                source_experiment_id="run_2",
            ),
        )
    )
    input_p = float(2.0 * stats.t.sf(abs(2.0), df=4.0))
    expected_single_z = float(stats.norm.isf(input_p / 2.0))
    expected_z = np.sqrt(2.0) * expected_single_z

    result = _aggregate(frame)
    row = result.to_dataframe().loc["MAPK1;S10;"]

    assert int(row["n_peptides_used"]) == 2
    assert float(row["uncertainty_statistic"]) == pytest.approx(expected_z)
    assert float(row["P.Value"]) == pytest.approx(2.0 * stats.norm.sf(abs(expected_z)))
    assert (
        row["uncertainty_method"]
        == PEPTIDE_TO_SITE_UNCERTAINTY_METHOD_STOUFFER_SIGNED_P
    )
    assert "run_1|run_2" == row["source_experiment_ids"]


def test_fixed_effect_combination_consumes_typed_standard_error_for_large_df() -> None:
    large_df = PEPTIDE_TO_SITE_FIXED_EFFECT_MIN_ASYMPTOTIC_MODERATED_DF * 10.0
    frame = _estimate_frame(
        (
            _estimate_row(
                peptide_id="pep_1",
                site_id="MAPK1;S10;",
                effect=1.0,
                standard_error=0.5,
                statistic=2.0,
                residual_degrees_of_freedom=large_df,
                moderated_degrees_of_freedom=large_df,
                source_experiment_id="run_1",
            ),
            _estimate_row(
                peptide_id="pep_2",
                site_id="MAPK1;S10;",
                effect=3.0,
                standard_error=1.0,
                statistic=3.0,
                residual_degrees_of_freedom=large_df,
                moderated_degrees_of_freedom=large_df,
                source_experiment_id="run_2",
            ),
        )
    )
    result = _aggregate(
        frame,
        config=PeptideToSiteAggregationConfig(
            uncertainty_method=(
                PEPTIDE_TO_SITE_UNCERTAINTY_METHOD_FIXED_EFFECT_INVERSE_VARIANCE
            )
        ),
    )
    row = result.to_dataframe().loc["MAPK1;S10;"]
    expected_effect = (4.0 * 1.0 + 1.0 * 3.0) / 5.0
    expected_standard_error = np.sqrt(1.0 / 5.0)
    expected_z = expected_effect / expected_standard_error

    assert float(row["logFC"]) == pytest.approx(expected_effect)
    assert float(row["standard_error"]) == pytest.approx(expected_standard_error)
    assert float(row["uncertainty_statistic"]) == pytest.approx(expected_z)
    assert float(row["P.Value"]) == pytest.approx(2.0 * stats.norm.sf(abs(expected_z)))


def test_fixed_effect_combination_rejects_small_finite_df_inputs() -> None:
    frame = _estimate_frame(
        (
            _estimate_row(
                peptide_id="pep_1",
                site_id="MAPK1;S10;",
                effect=1.0,
                standard_error=0.5,
                statistic=2.0,
                moderated_degrees_of_freedom=4.0,
                source_experiment_id="run_1",
            ),
            _estimate_row(
                peptide_id="pep_2",
                site_id="MAPK1;S10;",
                effect=1.0,
                standard_error=0.5,
                statistic=2.0,
                moderated_degrees_of_freedom=4.0,
                source_experiment_id="run_2",
            ),
        )
    )

    with pytest.raises(PhosPyInputError, match="asymptotic-normal input envelope"):
        _aggregate(
            frame,
            config=PeptideToSiteAggregationConfig(
                uncertainty_method=(
                    PEPTIDE_TO_SITE_UNCERTAINTY_METHOD_FIXED_EFFECT_INVERSE_VARIANCE
                )
            ),
        )


@pytest.mark.parametrize("method", ("none", "bonferroni", "holm"))
def test_multiple_testing_correction_is_configurable(method: str) -> None:
    frame = _estimate_frame(
        (
            _estimate_row(
                peptide_id="pep_site_1",
                site_id="MAPK1;S10;",
                effect=1.0,
                standard_error=0.5,
                statistic=2.0,
                moderated_degrees_of_freedom=4.0,
                source_experiment_id="run_1",
            ),
            _estimate_row(
                peptide_id="pep_site_2",
                site_id="MAPK1;T12;",
                effect=0.5,
                standard_error=0.5,
                statistic=1.0,
                moderated_degrees_of_freedom=4.0,
                source_experiment_id="run_1",
            ),
        )
    )
    result = _aggregate(
        frame,
        config=PeptideToSiteAggregationConfig(multiple_testing_method=method),
    )
    table = result.to_dataframe()

    np.testing.assert_allclose(
        table.loc[:, "adj.P.Val"].to_numpy(dtype=float),
        adjust_p_values(table.loc[:, "P.Value"].to_numpy(dtype=float), method=method),
    )
    assert set(table.loc[:, "correction_method"].astype(str)) == {method}


def test_multi_site_mapping_policy_is_visible_in_output_and_provenance() -> None:
    frame = _estimate_frame(
        (
            _estimate_row(
                peptide_id="pep_joint",
                site_id="MAPK1;S10,T12;",
                effect=1.0,
                standard_error=0.5,
                statistic=2.0,
                moderated_degrees_of_freedom=4.0,
                source_experiment_id="run_1",
                mapping_policy=PEPTIDE_TO_SITE_MAPPING_POLICY_KEEP_JOINT,
                mapping_uncertainty=True,
            ),
        )
    )
    result = _aggregate(frame)
    table = result.to_dataframe()
    row = table.loc["MAPK1;S10,T12;"]

    assert row["peptide_to_site_mapping_policy"] == (
        PEPTIDE_TO_SITE_MAPPING_POLICY_KEEP_JOINT
    )
    assert int(row["multi_site_estimate_count"]) == 1
    assert result.provenance["multi_site_estimate_count"] == 1
    assert result.provenance["peptide_to_site_mapping_policies"] == (
        PEPTIDE_TO_SITE_MAPPING_POLICY_KEEP_JOINT,
    )


def test_minimum_evidence_attrition_is_explicit_and_not_corrected_as_finite() -> None:
    result = _aggregate(
        _estimate_frame(),
        config=PeptideToSiteAggregationConfig(min_estimates_per_site=2),
    )
    row = result.to_dataframe().loc["MAPK1;S10;"]

    assert np.isnan(float(row["P.Value"]))
    assert np.isnan(float(row["adj.P.Val"]))
    assert int(row["n_peptides_used"]) == 1
    assert result.warnings
    assert result.provenance["withheld_below_minimum_site_count"] == 1


def test_result_records_aggregation_dependence_uncertainty_and_correction_methods() -> (
    None
):
    result = _aggregate(_estimate_frame())
    row = result.to_dataframe().iloc[0]

    assert row["aggregation_level"] == PEPTIDE_TO_SITE_AGGREGATION_LEVEL_SINGLE_ESTIMATE
    assert "single_estimate" in str(row["dependence_assumption"])
    assert (
        row["uncertainty_method"]
        == PEPTIDE_TO_SITE_UNCERTAINTY_METHOD_SINGLE_ESTIMATE_PASSTHROUGH
    )
    assert row["correction_method"] == "benjamini_hochberg"
    assert (
        result.provenance["preferred_phospy_lane"]
        == "resolve peptide evidence at sample-intensity level before differential model fitting"
    )
    assert result.provenance["input_contrast_id"] == "B_vs_A"
    assert result.provenance["input_contrast_orientation"] == "B_minus_A"
    assert result.provenance["input_effect_scale"] == "log2_fold_change"
    assert result.provenance["input_effect_unit"] == "log2_ratio"
    assert result.provenance["input_model_estimator_id"] == "limma_moderated_ols"
    assert result.provenance["input_statistic_distribution"] == "moderated_t"
    assert (
        result.provenance["input_uncertainty_method_version"]
        == "limma_ebayes_moderated_t_v1"
    )
    assert result.provenance["consistency_policy"] == (
        PEPTIDE_DIFFERENTIAL_CONSISTENCY_POLICY
    )
    assert result.provenance["mapping_weight_policy"] == (
        PEPTIDE_DIFFERENTIAL_MAPPING_WEIGHT_POLICY_REJECT
    )
    assert result.scientific_policies[0].parameters["multiple_testing_method"] == (
        "benjamini_hochberg"
    )


def test_provenance_serializes_new_estimate_identity_fields() -> None:
    result = _aggregate(_estimate_frame())
    policy_payload = result.scientific_policies[0].to_payload()
    policy_parameters = policy_payload["parameters"]

    assert isinstance(policy_parameters, dict)
    assert policy_parameters["support_status"] == (
        "supported_typed_estimate_combination_v2"
    )
    assert policy_parameters["input_contrast_id"] == "B_vs_A"
    assert policy_parameters["input_contrast_orientation"] == "B_minus_A"
    assert policy_parameters["input_effect_scale"] == "log2_fold_change"
    assert policy_parameters["input_effect_unit"] == "log2_ratio"
    assert policy_parameters["input_model_estimator_id"] == "limma_moderated_ols"
    assert policy_parameters["input_statistic_distribution"] == "moderated_t"
    assert (
        policy_parameters["input_uncertainty_method_version"]
        == "limma_ebayes_moderated_t_v1"
    )
    assert policy_parameters["consistency_policy"] == (
        PEPTIDE_DIFFERENTIAL_CONSISTENCY_POLICY
    )
    assert policy_parameters["mapping_weight_policy"] == (
        PEPTIDE_DIFFERENTIAL_MAPPING_WEIGHT_POLICY_REJECT
    )


def test_provenance_is_deterministic_for_input_order_changes() -> None:
    rows = (
        _estimate_row(
            peptide_id="pep_2",
            site_id="MAPK1;T12;",
            effect=0.5,
            standard_error=0.5,
            statistic=1.0,
            moderated_degrees_of_freedom=4.0,
            source_experiment_id="run_2",
        ),
        _estimate_row(
            peptide_id="pep_1",
            site_id="MAPK1;S10;",
            effect=1.0,
            standard_error=0.5,
            statistic=2.0,
            moderated_degrees_of_freedom=4.0,
            source_experiment_id="run_1",
        ),
    )
    forward = _aggregate(_estimate_frame(rows)).provenance
    reverse = _aggregate(_estimate_frame(tuple(reversed(rows)))).provenance

    assert dict(forward) == dict(reverse)


def test_old_raw_table_evidence_route_is_rejected() -> None:
    with pytest.raises(PhosPyInputError, match="raw peptide_differential_table"):
        PeptideToSiteAggregator().run_table(
            peptide_differential_table=pd.DataFrame(
                {"logFC": [1.0], "t": [2.0], "P.Value": [0.1]}
            ),
            evidence=object(),
            config=PeptideToSiteAggregationConfig(),
        )


def test_uncertainty_handling_stays_out_of_validators_interpreters_and_assemblers() -> (
    None
):
    forbidden_files = (
        ROOT / "src" / "phospy" / "validation" / "workflows" / "differential.py",
        ROOT / "src" / "phospy" / "workflows" / "differential" / "interpreter.py",
        ROOT / "src" / "phospy" / "workflows" / "differential" / "result_assembly.py",
    )
    forbidden_fragments = (
        "PeptideDifferentialEstimateTable",
        "PeptideToSiteAggregator",
        "signed_z_from_t_statistic",
        "signed_z_from_two_sided_p_value",
        "stouffer",
    )

    for path in forbidden_files:
        source = path.read_text(encoding="utf-8")
        for fragment in forbidden_fragments:
            assert fragment not in source

    executor_source = (
        ROOT
        / "src"
        / "phospy"
        / "science"
        / "differential"
        / "aggregation"
        / "executor.py"
    ).read_text(encoding="utf-8")
    assert "signed_z_from_t_statistic" in executor_source
    assert "adjust_p_values" in executor_source


def test_preferred_lane_resolves_peptide_intensities_before_core_differential_model() -> (
    None
):
    peptide_evidence = pd.DataFrame(
        {
            "peptide_row_id": ["pep_1", "pep_2"],
            "site_id": ["MAPK1;S10;", "GSK3B;S9;"],
            "unique_feature_id": ["feat_1", "feat_2"],
            "gene_symbol": ["MAPK1", "GSK3B"],
            "protein_accession": ["P28482", "P49841"],
            "site_string": ["S10", "S9"],
            "A_1": [100.0, 80.0],
            "A_2": [105.0, 84.0],
            "B_1": [200.0, 82.0],
            "B_2": [210.0, 83.0],
            "peptide_sequence": ["AAASAAA", "AAASAAA"],
            "modified_peptide_sequence": ["AAA[pS]AAA", "AAA[pS]AAA"],
            "multi_site": [False, False],
            "provenance_source": ["fixture", "fixture"],
            "localisation_confidence": [0.99, 0.98],
            "site_sequence": [
                "AAAAAAAAAAAAAAASAAAAAAAAAAAAAAA",
                "AAAAAAAAAAAAAAASAAAAAAAAAAAAAAA",
            ],
        }
    )
    dataset = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            site_resolution_mode=DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE,
            peptide_evidence=peptide_evidence,
            peptide_evidence_sample_intensity_columns=("A_1", "A_2", "B_1", "B_2"),
            multi_site_policy=DATASET_MULTI_SITE_POLICY_KEEP_JOINT,
            organism=Organism.HUMAN,
            input_intensity_scale="linear",
            preprocessing_config=DatasetPreprocessingConfig(
                intensity_transform=DatasetIntensityTransformConfig(
                    policy="log2",
                    pseudocount=1.0,
                )
            ),
        )
    )
    assert dataset.provenance is not None
    peptide_resolution = dataset.provenance.workflow_parameters[
        "peptide_evidence_resolution"
    ]
    assert isinstance(peptide_resolution, Mapping)
    assert peptide_resolution["aggregation_policy"] is not None
    assert (
        PEPTIDE_TO_SITE_AGGREGATION_LEVEL_SAMPLE_INTENSITY
        == "sample_intensity_resolution_before_differential_model"
    )

    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(
                sample_id="A_1",
                condition="A",
                biological_replicate_id="A_r1",
            ),
            SampleDesignRecord(
                sample_id="A_2",
                condition="A",
                biological_replicate_id="A_r2",
            ),
            SampleDesignRecord(
                sample_id="B_1",
                condition="B",
                biological_replicate_id="B_r1",
            ),
            SampleDesignRecord(
                sample_id="B_2",
                condition="B",
                biological_replicate_id="B_r2",
            ),
        )
    )
    result = DifferentialAnalysisWorkflow().run(
        DifferentialAnalysisRequest(
            dataset=dataset,
            design=design,
            contrasts=(
                Contrast(
                    name="B_vs_A",
                    numerator_condition="B",
                    denominator_condition="A",
                ),
            ),
        )
    )
    table = result.table_for("B_vs_A")

    assert "B_vs_A" in result.contrast_tables
    assert float(table.loc[table["display_id"] == "MAPK1;S10;", "logFC"].iloc[0]) > 0
