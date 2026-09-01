from __future__ import annotations

from collections.abc import Mapping
from typing import Any, cast

import numpy as np
import numpy.typing as npt
import pandas as pd
import pytest

import phospy.science.differential.protein_covariate_adjusted as kernel_module
from phospy import AnalysisReadyDatasetBuilder, DifferentialAnalysisWorkflow
from phospy.advanced import (
    PAIRED_DESIGN_POLICY_DUPLICATE_CORRELATION,
    PAIRED_DESIGN_POLICY_FIXED_BLOCK,
    DatasetProteinAwarePreparationConfig,
    DifferentialAnalysisConfig,
    DifferentialProteinAwareModelConfig,
    EmpiricalBayesConfig,
    TechnicalReplicatePolicy,
)
from phospy.api import (
    Contrast,
    DatasetBuildRequest,
    DatasetPreprocessingConfig,
    DifferentialAnalysisRequest,
    ExperimentalDesign,
    FixedEffectCovariate,
    Organism,
    SampleDesignRecord,
)
from phospy.errors import WorkflowBoundaryError
from phospy.provenance.hashing import fingerprint_table_strict
from phospy.science.differential.linear_model import (
    DIFFERENTIAL_LINEAR_MODEL_MAX_CONDITION_NUMBER,
)
from phospy.science.differential.models.protein_aware import (
    ProteinAwareDifferentialComputationRequest,
)
from phospy.science.differential.models.provenance import (
    DIFFERENTIAL_PROTEIN_AWARE_CLAIM_STATUS_EXPERIMENTAL,
    DIFFERENTIAL_PROTEIN_AWARE_CONTRAST_MATRIX_FINGERPRINT_NAME,
    DIFFERENTIAL_PROTEIN_AWARE_COVARIATE_MATRIX_FINGERPRINT_NAME,
    DIFFERENTIAL_PROTEIN_AWARE_DESIGN_MATRIX_FINGERPRINT_NAME,
    DIFFERENTIAL_PROTEIN_AWARE_MATCHED_PAIRS_FINGERPRINT_NAME,
    DIFFERENTIAL_PROTEIN_AWARE_PHOSPHO_MATRIX_FINGERPRINT_NAME,
    DIFFERENTIAL_PROTEIN_AWARE_SITE_ELIGIBILITY_FINGERPRINT_NAME,
    DifferentialProteinAwareInputFingerprints,
)
from phospy.science.differential.models.tables import (
    DIFFERENTIAL_RESULT_REASON_PROTEIN_AUGMENTED_DESIGN_ILL_CONDITIONED,
    DIFFERENTIAL_RESULT_REASON_PROTEIN_AUGMENTED_DESIGN_RANK_DEFICIENT,
    DIFFERENTIAL_RESULT_REASON_PROTEIN_COVARIATE_NON_FINITE,
    DIFFERENTIAL_RESULT_REASON_PROTEIN_COVARIATE_ZERO_VARIANCE,
    DIFFERENTIAL_RESULT_REASON_PROTEIN_PREPARATION_EXCLUDED,
    DIFFERENTIAL_RESULT_REASON_PROTEIN_PREPARATION_FALLBACK,
    DIFFERENTIAL_RESULT_STATUS_COLUMN,
    DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN,
    DIFFERENTIAL_RESULT_STATUS_TESTED,
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_AUGMENTED_DESIGN_INVALID,
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_COVARIATE_INVALID,
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_PREPARATION_INELIGIBLE,
)
from phospy.science.differential.protein_covariate_adjusted import (
    PROTEIN_AWARE_CENTERING_POLICY,
    PROTEIN_AWARE_COVARIATE_COEFFICIENT_NAME,
    ProteinCovariateAdjustedDifferentialKernel,
)
from phospy.science.statistics.multiple_testing import (
    MULTIPLE_TESTING_CORRECTION_BONFERRONI,
    adjust_p_values,
)
from tests.support.unsafe_dataset_states import (
    unsafe_mark_dataset_total_protein_correction_applied,
)

_FloatArray = npt.NDArray[np.float64]

SAMPLES = ("A_1", "A_2", "A_3", "A_4", "B_1", "B_2", "B_3", "B_4")
CONDITION_A = np.array([1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0])
CONDITION_B = np.array([0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0])
BASE_DESIGN = pd.DataFrame(
    {"A": CONDITION_A, "B": CONDITION_B},
    index=pd.Index(SAMPLES, name="sample_id"),
)
BASE_CONTRASTS = pd.DataFrame(
    {"B_vs_A": [-1.0, 1.0]},
    index=pd.Index(["A", "B"], name="coefficient"),
)
NOISE_SEEDS = (
    np.array([0.8, -1.3, 0.4, -0.7, 1.1, -0.2, 0.6, -0.9]),
    np.array([-0.6, 0.2, 1.4, -1.1, 0.5, -0.8, 1.2, -0.4]),
    np.array([1.5, -0.3, -1.0, 0.7, -1.2, 0.9, 0.1, -0.7]),
)


def test_protein_differential_known_effects_match_independent_oracle() -> None:
    shifted_protein = np.array([0.0, 1.0, 2.0, 3.0, 3.0, 4.0, 5.0, 6.0])
    independent_protein = np.array([3.0, 1.0, 4.0, 2.0, 5.0, 2.0, 6.0, 3.0])
    exact_protein = np.array([2.0, 3.0, 5.0, 7.0, 1.0, 4.0, 6.0, 8.0])
    stabilizer_protein = np.array([1.2, 0.4, 2.3, 1.7, 2.6, 1.1, 3.0, 2.2])
    site_specs = {
        "protein_only_apparent": (10.0, 10.0, 1.5, shifted_protein, 0.02),
        "phospho_and_protein": (7.0, 8.25, -0.7, independent_protein, 0.03),
        "exact_coefficient_identity": (3.0, 4.2, 0.6, exact_protein, 0.0),
        "moderation_stabilizer": (5.0, 5.35, 0.25, stabilizer_protein, 0.04),
    }
    matrix = _matrix(
        {
            site: _site_values(
                beta={"A": beta_a, "B": beta_b},
                gamma=gamma,
                protein=protein,
                noise_scale=noise_scale,
                noise_seed=NOISE_SEEDS[position % len(NOISE_SEEDS)],
            )
            for position, (
                site,
                (beta_a, beta_b, gamma, protein, noise_scale),
            ) in enumerate(site_specs.items())
        }
    )
    request = _kernel_request(
        matrix=matrix,
        proteins=_proteins(
            {f"protein_{site}": spec[3] for site, spec in site_specs.items()}
        ),
        pairs=_pairs(
            {
                site: (f"ID_{position}", f"protein_{site}")
                for position, site in enumerate(site_specs, start=1)
            }
        ),
    )

    result = ProteinCovariateAdjustedDifferentialKernel().run(request)

    assert result.tested_site_ids == tuple(site_specs)
    ordinary_apparent = _ordinary_contrast_effect(
        matrix.loc["protein_only_apparent"].to_numpy(dtype=float),
    )
    assert ordinary_apparent == pytest.approx(4.5, abs=0.05)
    assert result.table_for("B_vs_A").loc[
        "protein_only_apparent", "logFC"
    ] == pytest.approx(0.0, abs=1.0e-11)
    assert result.protein_coefficient_series().loc[
        "protein_only_apparent"
    ] == pytest.approx(1.5, abs=1.0e-11)
    assert result.table_for("B_vs_A").loc[
        "phospho_and_protein", "logFC"
    ] == pytest.approx(1.25, abs=1.0e-11)
    assert result.protein_coefficient_series().loc[
        "phospho_and_protein"
    ] == pytest.approx(-0.7, abs=1.0e-11)
    assert result.table_for("B_vs_A").loc[
        "exact_coefficient_identity", "logFC"
    ] == pytest.approx(1.2, abs=1.0e-12)

    for site, (_, _, _, protein, _) in site_specs.items():
        oracle = _raw_oracle(
            values=matrix.loc[site].to_numpy(dtype=float),
            protein=protein,
        )
        np.testing.assert_allclose(
            result.coefficient_dataframe().loc[site].to_numpy(dtype=float),
            oracle["coefficients"],
            rtol=1.0e-11,
            atol=1.0e-11,
        )
        np.testing.assert_allclose(
            result.residuals_dataframe().loc[site].to_numpy(dtype=float),
            oracle["residuals"],
            rtol=1.0e-11,
            atol=1.0e-11,
        )
        np.testing.assert_allclose(
            result.contrast_standard_error_scale_dataframe()
            .loc[site]
            .to_numpy(dtype=float),
            oracle["contrast_scales"],
            rtol=1.0e-11,
            atol=1.0e-11,
        )


def test_protein_differential_unrelated_covariate_preserves_condition_effect() -> None:
    protein = np.array([0.0, 1.0, 3.0, 4.0, 0.0, 1.0, 3.0, 4.0])
    values = _site_values(
        beta={"A": 6.0, "B": 7.4},
        gamma=0.0,
        protein=protein,
        noise_scale=0.025,
    )
    request = _kernel_request(
        matrix=_matrix({"unrelated_site": values}),
        proteins=_proteins({"unrelated_protein": protein}),
        pairs=_pairs({"unrelated_site": ("ID_1", "unrelated_protein")}),
    )

    result = ProteinCovariateAdjustedDifferentialKernel().run(request)

    adjusted = float(result.table_for("B_vs_A").loc["unrelated_site", "logFC"])
    ordinary = _ordinary_contrast_effect(values)
    assert adjusted == pytest.approx(ordinary, rel=1.0e-12, abs=1.0e-12)
    assert adjusted == pytest.approx(1.4, rel=1.0e-12, abs=1.0e-12)
    assert result.protein_coefficient_series().loc["unrelated_site"] == pytest.approx(
        0.0,
        abs=1.0e-12,
    )


def test_protein_differential_shared_protein_row_keeps_independent_site_coefficients() -> (
    None
):
    shared_protein = np.array([0.0, 1.0, 2.0, 4.0, 1.0, 3.0, 5.0, 6.0])
    matrix = _matrix(
        {
            "shared_site_a": _site_values(
                beta={"A": 4.0, "B": 5.0},
                gamma=0.4,
                protein=shared_protein,
                noise_scale=0.03,
            ),
            "shared_site_b": _site_values(
                beta={"A": 8.0, "B": 7.5},
                gamma=-0.9,
                protein=shared_protein,
                noise_scale=0.04,
                noise_seed=NOISE_SEEDS[1],
            ),
        }
    )
    request = _kernel_request(
        matrix=matrix,
        proteins=_proteins({"shared_protein": shared_protein}),
        pairs=_pairs(
            {
                "shared_site_a": ("ID_1", "shared_protein"),
                "shared_site_b": ("ID_2", "shared_protein"),
            }
        ),
    )

    result = ProteinCovariateAdjustedDifferentialKernel().run(request)

    assert result.augmented_design_diagnostics_dataframe().index.tolist() == [
        "shared_protein"
    ]
    assert result.coefficient_dataframe().columns.tolist() == [
        "A",
        "B",
        PROTEIN_AWARE_COVARIATE_COEFFICIENT_NAME,
    ]
    assert result.augmented_design_diagnostics_dataframe().loc[
        "shared_protein", "protein_covariate_contrast_weights"
    ] == (0.0,)
    assert result.protein_coefficient_series().loc["shared_site_a"] == pytest.approx(
        0.4,
        abs=1.0e-11,
    )
    assert result.protein_coefficient_series().loc["shared_site_b"] == pytest.approx(
        -0.9,
        abs=1.0e-11,
    )
    scales = result.contrast_standard_error_scale_dataframe()["B_vs_A"]
    assert scales.loc["shared_site_a"] == pytest.approx(
        scales.loc["shared_site_b"],
        rel=0.0,
        abs=0.0,
    )


def test_protein_differential_invalid_covariate_and_confounding_statuses_are_typed() -> (
    None
):
    valid = np.array([0.0, 1.0, 2.0, 3.0, 0.5, 1.5, 2.5, 3.5])
    constant = np.full(len(SAMPLES), 4.2)
    non_finite = np.array([0.2, 0.4, np.nan, 0.6, 0.8, 1.0, 1.2, 1.4])
    collinear = CONDITION_B.copy()
    near_collinear = CONDITION_B + np.array(
        [-1.0e-12, 0.0, 1.0e-12, 2.0e-12, -1.0e-12, 0.0, 1.0e-12, 2.0e-12]
    )
    request = _kernel_request(
        matrix=_matrix(
            {
                "valid_site": _site_values(
                    beta={"A": 5.0, "B": 6.0},
                    gamma=0.3,
                    protein=valid,
                ),
                "constant_site": _site_values(
                    beta={"A": 5.0, "B": 6.0},
                    gamma=0.3,
                    protein=valid,
                ),
                "nonfinite_site": _site_values(
                    beta={"A": 5.0, "B": 6.0},
                    gamma=0.3,
                    protein=valid,
                ),
                "collinear_site": _site_values(
                    beta={"A": 5.0, "B": 6.0},
                    gamma=0.3,
                    protein=valid,
                ),
                "near_collinear_site": _site_values(
                    beta={"A": 5.0, "B": 6.0},
                    gamma=0.3,
                    protein=valid,
                ),
            }
        ),
        proteins=_proteins(
            {
                "valid_protein": valid,
                "constant_protein": constant,
                "nonfinite_protein": non_finite,
                "collinear_protein": collinear,
                "near_collinear_protein": near_collinear,
            }
        ),
        pairs=_pairs(
            {
                "valid_site": ("ID_valid", "valid_protein"),
                "constant_site": ("ID_constant", "constant_protein"),
                "nonfinite_site": ("ID_nonfinite", "nonfinite_protein"),
                "collinear_site": ("ID_collinear", "collinear_protein"),
                "near_collinear_site": ("ID_near_collinear", "near_collinear_protein"),
            }
        ),
    )

    result = ProteinCovariateAdjustedDifferentialKernel().run(request)

    assert result.tested_site_ids == ("valid_site",)
    assert result.residuals_dataframe().columns.tolist() == list(SAMPLES)
    failures = result.site_failure_diagnostics_dataframe()
    assert failures.loc["constant_site", DIFFERENTIAL_RESULT_STATUS_COLUMN] == (
        DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_COVARIATE_INVALID
    )
    assert failures.loc["constant_site", DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN] == (
        DIFFERENTIAL_RESULT_REASON_PROTEIN_COVARIATE_ZERO_VARIANCE
    )
    assert failures.loc["nonfinite_site", DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN] == (
        DIFFERENTIAL_RESULT_REASON_PROTEIN_COVARIATE_NON_FINITE
    )
    assert failures.loc["collinear_site", DIFFERENTIAL_RESULT_STATUS_COLUMN] == (
        DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_AUGMENTED_DESIGN_INVALID
    )
    assert failures.loc["near_collinear_site", DIFFERENTIAL_RESULT_STATUS_COLUMN] == (
        DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_AUGMENTED_DESIGN_INVALID
    )
    group_failures = result.augmented_design_failure_diagnostics_dataframe()
    assert (
        group_failures.loc[
            "collinear_protein", DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN
        ]
        == DIFFERENTIAL_RESULT_REASON_PROTEIN_AUGMENTED_DESIGN_RANK_DEFICIENT
    )
    assert group_failures.loc[
        "near_collinear_protein", DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN
    ] in {
        DIFFERENTIAL_RESULT_REASON_PROTEIN_AUGMENTED_DESIGN_RANK_DEFICIENT,
        DIFFERENTIAL_RESULT_REASON_PROTEIN_AUGMENTED_DESIGN_ILL_CONDITIONED,
    }
    assert "nonfinite_site" not in result.table_for("B_vs_A").index


def test_protein_differential_fixed_covariate_and_block_designs_recover_effects() -> (
    None
):
    dose_design = BASE_DESIGN.assign(
        dose=np.array([-1.5, -0.5, 0.5, 1.5, -1.2, -0.2, 0.8, 1.8])
    )
    block_design = BASE_DESIGN.assign(
        block_2=np.array([0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0]),
        block_3=np.array([0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0]),
        block_4=np.array([0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]),
    )
    protein = np.array([0.0, 2.0, 1.0, 3.0, 1.5, 4.0, 2.5, 5.0])

    for name, design, beta in (
        (
            "dose_site",
            dose_design,
            {"A": 6.0, "B": 7.1, "dose": -0.35},
        ),
        (
            "block_site",
            block_design,
            {
                "A": 4.0,
                "B": 5.4,
                "block_2": 0.2,
                "block_3": -0.15,
                "block_4": 0.4,
            },
        ),
    ):
        contrasts = pd.DataFrame(
            {"B_vs_A": [-1.0, 1.0, *([0.0] * (design.shape[1] - 2))]},
            index=pd.Index(design.columns, name="coefficient"),
        )
        values = _site_values(
            beta=beta,
            gamma=0.55,
            protein=protein,
            base_design=design,
            noise_scale=0.025,
        )
        result = ProteinCovariateAdjustedDifferentialKernel().run(
            _kernel_request(
                matrix=_matrix({name: values}),
                proteins=_proteins({f"protein_{name}": protein}),
                pairs=_pairs({name: ("ID_1", f"protein_{name}")}),
                design=design,
                contrasts=contrasts,
            )
        )
        oracle = _raw_oracle(
            values=values,
            protein=protein,
            base_design=design,
            base_contrasts=contrasts,
        )

        assert result.coefficient_dataframe().columns.tolist() == [
            *list(design.columns),
            PROTEIN_AWARE_COVARIATE_COEFFICIENT_NAME,
        ]
        assert result.augmented_design_diagnostics_dataframe().iloc[0][
            "protein_covariate_contrast_weights"
        ] == (0.0,)
        np.testing.assert_allclose(
            result.coefficient_dataframe().loc[name].to_numpy(dtype=float),
            oracle["coefficients"],
            rtol=1.0e-11,
            atol=1.0e-11,
        )
        assert result.table_for("B_vs_A").loc[name, "logFC"] == pytest.approx(
            float(oracle["contrast_effects"][0]),
            rel=1.0e-11,
            abs=1.0e-11,
        )


def test_protein_differential_public_workflow_fixed_covariate_and_fixed_block_recover_known_effects() -> (
    None
):
    protein = np.array([0.0, 2.0, 1.0, 3.0, 1.5, 4.0, 2.5, 5.0])
    dose = np.array([-1.5, -0.5, 0.5, 1.5, -1.2, -0.2, 0.8, 1.8])
    dose_design = _expected_public_design_matrix(
        _workflow_design_with_continuous_covariate("dose", dose),
        paired_design_policy="reject",
    )
    block_design = _expected_public_design_matrix(
        _workflow_design(paired_design_policy=PAIRED_DESIGN_POLICY_FIXED_BLOCK),
        paired_design_policy=PAIRED_DESIGN_POLICY_FIXED_BLOCK,
    )

    for design, paired_design_policy, base_design, beta, gamma, noise_seed in (
        (
            _workflow_design_with_continuous_covariate("dose", dose),
            "reject",
            dose_design,
            {"A": 6.0, "B": 7.1, "dose": -0.35},
            0.55,
            NOISE_SEEDS[0],
        ),
        (
            _workflow_design(paired_design_policy=PAIRED_DESIGN_POLICY_FIXED_BLOCK),
            PAIRED_DESIGN_POLICY_FIXED_BLOCK,
            block_design,
            {
                "A": 4.0,
                "B": 5.4,
                "block[pair_2]": 0.2,
                "block[pair_3]": -0.15,
                "block[pair_4]": 0.4,
            },
            -0.3,
            NOISE_SEEDS[1],
        ),
    ):
        values = _site_values(
            beta=beta,
            gamma=gamma,
            protein=protein,
            base_design=base_design,
            noise_seed=noise_seed,
            noise_scale=0.02,
        )
        dataset = _single_site_workflow_dataset(
            protein_id="P1",
            total_protein_id="P1",
            phosphosite_values=values,
            total_protein_values=protein,
        )
        request = DifferentialAnalysisRequest(
            dataset=dataset,
            design=design,
            contrasts=(
                Contrast(
                    name="B_vs_A",
                    numerator_condition="B",
                    denominator_condition="A",
                ),
            ),
            config=DifferentialAnalysisConfig(
                paired_design_policy=paired_design_policy,  # type: ignore[arg-type]
                protein_aware_model=DifferentialProteinAwareModelConfig(),
            ),
        )

        result = DifferentialAnalysisWorkflow().run(request)

        site_key = _site_key_for_metadata_value(dataset, "protein_id", "P1")
        contrast_matrix = _expected_public_contrast_matrix(
            coefficient_labels=tuple(str(column) for column in base_design.columns),
            contrasts=request.contrasts,
        )
        oracle = _raw_oracle(
            values=values,
            protein=protein,
            base_design=base_design,
            base_contrasts=contrast_matrix,
        )
        diagnostics = result.protein_aware_diagnostics
        assert diagnostics is not None
        per_site = diagnostics.per_site_diagnostics_dataframe()

        assert result.table_for("B_vs_A").loc[site_key, "logFC"] == pytest.approx(
            float(oracle["contrast_effects"][0]),
            rel=1.0e-11,
            abs=1.0e-11,
        )
        assert per_site.loc[site_key, "protein_covariate_coefficient"] == pytest.approx(
            gamma, rel=1.0e-11, abs=1.0e-11
        )
        assert diagnostics.base_design_rank == int(
            np.linalg.matrix_rank(base_design.to_numpy(dtype=float))
        )
        assert diagnostics.expected_augmented_rank == int(
            np.linalg.matrix_rank(_augmented_design(base_design, protein))
        )
        assert per_site.loc[
            site_key, "protein_augmented_design_condition_number"
        ] == pytest.approx(
            _scaled_condition_number(_augmented_design(base_design, protein)),
            rel=1.0e-12,
            abs=1.0e-12,
        )
        policy = _protein_aware_policy(result)
        assert policy.input_fingerprints.design_matrix == fingerprint_table_strict(
            base_design,
            name=DIFFERENTIAL_PROTEIN_AWARE_DESIGN_MATRIX_FINGERPRINT_NAME,
        )
        assert policy.input_fingerprints.contrast_matrix == fingerprint_table_strict(
            contrast_matrix,
            name=DIFFERENTIAL_PROTEIN_AWARE_CONTRAST_MATRIX_FINGERPRINT_NAME,
        )


def test_protein_differential_row_order_changes_presentation_not_site_values() -> None:
    protein_a = np.array([0.0, 1.0, 2.0, 3.0, 0.5, 1.5, 2.5, 3.5])
    protein_b = np.array([2.0, 0.5, 3.0, 1.0, 2.8, 1.2, 3.6, 1.9])
    matrix = _matrix(
        {
            "site_a": _site_values(
                beta={"A": 4.0, "B": 5.0},
                gamma=0.2,
                protein=protein_a,
            ),
            "site_b": _site_values(
                beta={"A": 3.0, "B": 2.4},
                gamma=-0.6,
                protein=protein_b,
                noise_seed=NOISE_SEEDS[1],
            ),
            "site_c": _site_values(
                beta={"A": 6.0, "B": 6.7},
                gamma=0.9,
                protein=protein_a,
                noise_seed=NOISE_SEEDS[2],
            ),
        }
    )
    pairs = _pairs(
        {
            "site_a": ("ID_A", "protein_a"),
            "site_b": ("ID_B", "protein_b"),
            "site_c": ("ID_C", "protein_a"),
        }
    )
    proteins = _proteins({"protein_a": protein_a, "protein_b": protein_b})
    baseline = ProteinCovariateAdjustedDifferentialKernel().run(
        _kernel_request(matrix=matrix, proteins=proteins, pairs=pairs)
    )
    reordered_sites = ["site_c", "site_a", "site_b"]
    reordered = ProteinCovariateAdjustedDifferentialKernel().run(
        _kernel_request(
            matrix=matrix.loc[reordered_sites],
            proteins=proteins.loc[["protein_b", "protein_a"]],
            pairs=pairs.set_index("site_key").loc[reordered_sites].reset_index(),
        )
    )

    assert reordered.tested_site_ids == tuple(reordered_sites)
    pd.testing.assert_frame_equal(
        baseline.table_for("B_vs_A").sort_index(),
        reordered.table_for("B_vs_A").sort_index(),
        check_exact=False,
        rtol=1.0e-12,
        atol=1.0e-12,
    )
    pd.testing.assert_frame_equal(
        baseline.coefficient_dataframe().sort_index(),
        reordered.coefficient_dataframe().sort_index(),
        check_exact=False,
        rtol=1.0e-12,
        atol=1.0e-12,
    )


@pytest.mark.parametrize(
    "empirical_bayes",
    [
        EmpiricalBayesConfig(method="standard", trend=False),
        EmpiricalBayesConfig(method="robust", trend=False),
        EmpiricalBayesConfig(method="standard", trend=True),
        EmpiricalBayesConfig(method="robust", trend=True),
    ],
)
def test_protein_differential_empirical_bayes_uses_only_tested_augmented_rows(
    monkeypatch: pytest.MonkeyPatch,
    empirical_bayes: EmpiricalBayesConfig,
) -> None:
    protein_a = np.array([0.0, 1.0, 2.0, 3.0, 0.5, 1.5, 2.5, 3.5])
    protein_b = np.array([2.0, 0.5, 3.0, 1.0, 2.8, 1.2, 3.6, 1.9])
    constant = np.full(len(SAMPLES), 4.0)
    matrix = _matrix(
        {
            "tested_a": _site_values(
                beta={"A": 4.0, "B": 5.0},
                gamma=0.2,
                protein=protein_a,
            ),
            "tested_b": _site_values(
                beta={"A": 3.0, "B": 2.4},
                gamma=-0.6,
                protein=protein_b,
                noise_seed=NOISE_SEEDS[1],
            ),
            "tested_c": _site_values(
                beta={"A": 6.0, "B": 6.7},
                gamma=0.9,
                protein=protein_a,
                noise_seed=NOISE_SEEDS[2],
            ),
            "withheld_constant": _site_values(
                beta={"A": 5.0, "B": 5.2},
                gamma=0.4,
                protein=protein_a,
            ),
        }
    )
    request = _kernel_request(
        matrix=matrix,
        proteins=_proteins(
            {"protein_a": protein_a, "protein_b": protein_b, "constant": constant}
        ),
        pairs=_pairs(
            {
                "tested_a": ("ID_A", "protein_a"),
                "tested_b": ("ID_B", "protein_b"),
                "tested_c": ("ID_C", "protein_a"),
                "withheld_constant": ("ID_D", "constant"),
            }
        ),
        contrasts=_multiple_contrasts(),
        empirical_bayes=empirical_bayes,
        multiple_testing_method=MULTIPLE_TESTING_CORRECTION_BONFERRONI,
    )
    captured: dict[str, Any] = {}
    real_fit = kernel_module.fit_empirical_bayes

    def _spy_fit_empirical_bayes(**kwargs: Any) -> Any:
        captured.update(kwargs)
        return real_fit(**kwargs)

    monkeypatch.setattr(kernel_module, "fit_empirical_bayes", _spy_fit_empirical_bayes)

    result = ProteinCovariateAdjustedDifferentialKernel().run(request)

    expected_sites = ("tested_a", "tested_b", "tested_c")
    assert result.tested_site_ids == expected_sites
    expected_variances = np.array(
        [
            _raw_oracle(
                values=matrix.loc[site].to_numpy(dtype=float),
                protein=protein_a if site != "tested_b" else protein_b,
            )["residual_variance"]
            for site in expected_sites
        ],
        dtype=float,
    )
    np.testing.assert_allclose(
        cast(_FloatArray, captured["variances"]),
        expected_variances,
        rtol=1.0e-11,
        atol=1.0e-11,
    )
    assert captured["residual_dof"] == pytest.approx(5.0, abs=1.0e-12)
    assert captured["method"] == empirical_bayes.method
    assert captured["trend"] is empirical_bayes.trend
    assert captured["winsor_tail_p"] == empirical_bayes.winsor_tail_p
    np.testing.assert_allclose(
        cast(_FloatArray, captured["mean_intensity"]),
        matrix.loc[list(expected_sites)].mean(axis=1).to_numpy(dtype=float),
        rtol=1.0e-12,
        atol=1.0e-12,
    )
    assert "withheld_constant" not in result.prior_residual_variance_series().index
    for contrast_name in ("B_vs_A", "A_level"):
        table = result.table_for(contrast_name)
        expected_adjusted = adjust_p_values(
            table["P.Value"].to_numpy(dtype=float),
            method=MULTIPLE_TESTING_CORRECTION_BONFERRONI,
        )
        np.testing.assert_allclose(table["adj.P.Val"], expected_adjusted)
        np.testing.assert_allclose(
            table["adj.P.Val"],
            np.clip(table["P.Value"].to_numpy(dtype=float) * len(expected_sites), 0, 1),
        )


def test_protein_differential_public_workflow_reports_no_fallback_mapping_attrition() -> (
    None
):
    dataset = _workflow_dataset(
        protein_ids=("P1", "P_missing", "P3;P4"),
        total_protein_ids=("P1", "P3", "P4"),
        protein_mapping_policy="allow_missing_with_report",
    )
    request = _workflow_request(dataset)

    adjusted = DifferentialAnalysisWorkflow().run(request)
    ordinary = DifferentialAnalysisWorkflow().run(_ordinary_workflow_request(dataset))

    adjusted_table = adjusted.table_for("B_vs_A")
    ordinary_table = ordinary.table_for("B_vs_A")
    missing_site = _site_key_for_metadata_value(dataset, "protein_id", "P_missing")
    ambiguous_site = _site_key_for_metadata_value(dataset, "protein_id", "P3;P4")
    tested_site = _site_key_for_metadata_value(dataset, "protein_id", "P1")

    assert adjusted_table.loc[tested_site, DIFFERENTIAL_RESULT_STATUS_COLUMN] == (
        DIFFERENTIAL_RESULT_STATUS_TESTED
    )
    assert adjusted_table.loc[missing_site, DIFFERENTIAL_RESULT_STATUS_COLUMN] == (
        DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_PREPARATION_INELIGIBLE
    )
    assert adjusted_table.loc[
        missing_site, DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN
    ] == (DIFFERENTIAL_RESULT_REASON_PROTEIN_PREPARATION_FALLBACK)
    assert adjusted_table.loc[ambiguous_site, DIFFERENTIAL_RESULT_STATUS_COLUMN] == (
        DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_PREPARATION_INELIGIBLE
    )
    assert (
        adjusted_table.loc[ambiguous_site, DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN]
        == DIFFERENTIAL_RESULT_REASON_PROTEIN_PREPARATION_EXCLUDED
    )
    assert (
        adjusted_table.loc[
            [missing_site, ambiguous_site], ["logFC", "t", "P.Value", "adj.P.Val"]
        ]
        .isna()
        .all()
        .all()
    )
    assert np.isfinite(
        ordinary_table.loc[
            [missing_site, ambiguous_site], ["logFC", "t", "P.Value", "adj.P.Val"]
        ].to_numpy(dtype=float)
    ).all()


def test_protein_differential_sample_reordering_is_equivalent_and_fingerprinted() -> (
    None
):
    dataset = _workflow_dataset(
        protein_ids=("P1", "P2", "P3"),
        total_protein_ids=("P1", "P2", "P3"),
    )
    first_order = ("B_4", "A_2", "B_2", "A_4", "B_1", "A_1")
    second_order = ("A_1", "B_1", "A_2", "B_2", "A_4", "B_4")
    first_request = _workflow_request(
        dataset,
        sample_ids=first_order,
        allow_design_subset=True,
    )
    repeated_first = _workflow_request(
        dataset,
        sample_ids=first_order,
        allow_design_subset=True,
    )
    second_request = _workflow_request(
        dataset,
        sample_ids=second_order,
        allow_design_subset=True,
    )

    first = DifferentialAnalysisWorkflow().run(first_request)
    repeated = DifferentialAnalysisWorkflow().run(repeated_first)
    second = DifferentialAnalysisWorkflow().run(second_request)

    pd.testing.assert_frame_equal(
        first.table_for("B_vs_A"), repeated.table_for("B_vs_A")
    )
    pd.testing.assert_frame_equal(
        first.table_for("B_vs_A"),
        second.table_for("B_vs_A"),
        check_exact=False,
        rtol=1.0e-12,
        atol=1.0e-12,
    )
    first_policy = _protein_aware_policy(first)
    second_policy = _protein_aware_policy(second)
    assert (
        first_policy.claim_status
        == DIFFERENTIAL_PROTEIN_AWARE_CLAIM_STATUS_EXPERIMENTAL
    )
    assert first_policy.execution_sample_order == first_order
    assert second_policy.execution_sample_order == second_order
    assert first_policy.input_fingerprints == _expected_input_fingerprints(
        first_request
    )
    assert _protein_aware_policy(repeated).input_fingerprints == (
        first_policy.input_fingerprints
    )
    assert (
        first_policy.input_fingerprints.phospho_matrix.exact_hash_value
        != second_policy.input_fingerprints.phospho_matrix.exact_hash_value
    )
    assert (
        first_policy.input_fingerprints.protein_covariate_matrix.exact_hash_value
        != second_policy.input_fingerprints.protein_covariate_matrix.exact_hash_value
    )
    assert (
        first_policy.input_fingerprints.design_matrix.exact_hash_value
        != second_policy.input_fingerprints.design_matrix.exact_hash_value
    )


def test_protein_differential_public_diagnostics_match_model_facts() -> None:
    protein = np.array([0.0, 1.0, 2.0, 3.0, 0.5, 1.5, 2.5, 3.5])
    dataset = _workflow_dataset(
        protein_ids=("P1", "P2", "P3"),
        total_protein_ids=("P1", "P2", "P3"),
        protein_values={
            "P1": protein,
            "P2": protein + 0.7,
            "P3": np.array([2.0, 0.5, 3.0, 1.0, 2.8, 1.2, 3.6, 1.9]),
        },
    )
    result = DifferentialAnalysisWorkflow().run(_workflow_request(dataset))

    diagnostics = result.protein_aware_diagnostics
    assert diagnostics is not None
    assert (
        diagnostics.claim_status == DIFFERENTIAL_PROTEIN_AWARE_CLAIM_STATUS_EXPERIMENTAL
    )
    assert diagnostics.protein_covariate_centered is True
    assert diagnostics.protein_covariate_standardized is False
    assert diagnostics.protein_covariate_imputation_policy == "none"
    assert diagnostics.fallback_policy == "no_fallback_to_ordinary_differential_lane"
    assert diagnostics.execution_sample_order == SAMPLES
    assert diagnostics.base_design_rank == 2
    assert diagnostics.base_residual_degrees_of_freedom == pytest.approx(6.0)
    assert diagnostics.expected_augmented_rank == 3
    assert diagnostics.common_augmented_rank == 3
    assert diagnostics.common_augmented_residual_degrees_of_freedom == pytest.approx(
        5.0
    )
    assert result.diagnostics.singular_values == ()
    assert any(
        "median condition number across independently fitted total-protein-row"
        in warning
        for warning in result.diagnostics.warnings
    )

    site_key = _site_key_for_metadata_value(dataset, "protein_id", "P1")
    per_site = diagnostics.per_site_diagnostics_dataframe()
    expected_centered_std = float(np.std(protein - float(np.mean(protein)), ddof=1))
    assert per_site.loc[site_key, "protein_covariate_raw_mean"] == pytest.approx(
        float(np.mean(protein)),
        abs=1.0e-12,
    )
    assert per_site.loc[
        site_key, "protein_covariate_raw_standard_deviation"
    ] == pytest.approx(float(np.std(protein, ddof=1)), abs=1.0e-12)
    assert per_site.loc[
        site_key, "protein_covariate_centered_standard_deviation"
    ] == pytest.approx(expected_centered_std, abs=1.0e-12)
    assert (
        per_site.loc[site_key, "protein_augmented_design_rank"]
        == diagnostics.common_augmented_rank
    )
    assert per_site.loc[
        site_key, "protein_augmented_design_residual_degrees_of_freedom"
    ] == pytest.approx(5.0)
    assert per_site.loc[site_key, "protein_covariate_coefficient"] == pytest.approx(
        0.45,
        abs=1.0e-10,
    )
    expected_condition_number = _scaled_condition_number(
        _augmented_design(BASE_DESIGN, protein)
    )
    assert per_site.loc[
        site_key, "protein_augmented_design_condition_number"
    ] == pytest.approx(expected_condition_number, rel=1.0e-12, abs=1.0e-12)


@pytest.mark.parametrize(
    ("case", "expected_seam"),
    [
        ("all_withheld", "differential.protein_aware_inputs.all_sites_withheld"),
        (
            "prior_subtraction",
            "differential.protein_aware_inputs.prior_total_protein_subtraction",
        ),
        (
            "technical_aggregation",
            "differential.protein_aware_inputs.technical_replicate_aggregation",
        ),
        (
            "duplicate_correlation",
            "differential.protein_aware_inputs.duplicate_correlation",
        ),
    ],
)
def test_protein_differential_global_failure_boundaries_are_fail_closed(
    case: str,
    expected_seam: str,
) -> None:
    constant_proteins = (
        frozenset({"P1", "P2", "P3"}) if case == "all_withheld" else frozenset()
    )
    dataset = _workflow_dataset(
        protein_ids=("P1", "P2", "P3"),
        total_protein_ids=("P1", "P2", "P3"),
        constant_total_proteins=constant_proteins,
    )
    if case == "prior_subtraction":
        unsafe_mark_dataset_total_protein_correction_applied(dataset)
    request = _workflow_request(
        dataset,
        paired_design_policy=(
            PAIRED_DESIGN_POLICY_DUPLICATE_CORRELATION
            if case == "duplicate_correlation"
            else "reject"
        ),
        technical_replicate_policy=(
            TechnicalReplicatePolicy.MEAN
            if case == "technical_aggregation"
            else TechnicalReplicatePolicy.REJECT
        ),
        technical_replicates=case == "technical_aggregation",
    )

    with pytest.raises(WorkflowBoundaryError, match=expected_seam) as exc_info:
        DifferentialAnalysisWorkflow().run(request)

    if case == "all_withheld":
        assert exc_info.value.details["eligibility_counts"] == {
            "total_site_count": 3,
            "ordinary_testable_site_count": 3,
            "protein_preparation_candidate_site_count": 3,
            "protein_aware_tested_site_count": 0,
        }
        assert exc_info.value.details["status_counts"] == {
            DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_COVARIATE_INVALID: 3
        }
        assert exc_info.value.details["reason_counts"] == {
            DIFFERENTIAL_RESULT_REASON_PROTEIN_COVARIATE_ZERO_VARIANCE: 3
        }


def _kernel_request(
    *,
    matrix: pd.DataFrame,
    proteins: pd.DataFrame,
    pairs: pd.DataFrame,
    design: pd.DataFrame = BASE_DESIGN,
    contrasts: pd.DataFrame = BASE_CONTRASTS,
    empirical_bayes: EmpiricalBayesConfig | None = None,
    multiple_testing_method: str = "benjamini_hochberg",
) -> ProteinAwareDifferentialComputationRequest:
    return ProteinAwareDifferentialComputationRequest(
        phosphosite_matrix=matrix,
        base_design=design,
        base_contrasts=contrasts,
        sample_order=tuple(str(value) for value in design.index.tolist()),
        matched_pairs=pairs,
        resolved_protein_covariates=proteins.loc[:, list(design.index)],
        empirical_bayes=empirical_bayes if empirical_bayes else EmpiricalBayesConfig(),
        multiple_testing_method=multiple_testing_method,
    )


def _site_values(
    *,
    beta: Mapping[str, float],
    gamma: float,
    protein: _FloatArray,
    base_design: pd.DataFrame = BASE_DESIGN,
    noise_seed: _FloatArray = NOISE_SEEDS[0],
    noise_scale: float = 0.02,
) -> _FloatArray:
    centered = protein - float(np.mean(protein))
    augmented = _augmented_design(base_design, protein)
    beta_vector = np.array([float(beta[column]) for column in base_design.columns])
    signal = base_design.to_numpy(dtype=float) @ beta_vector + gamma * centered
    if noise_scale == 0.0:
        return np.asarray(signal, dtype=np.float64)
    noise = _orthogonal_noise(augmented, noise_seed[: len(centered)])
    return np.asarray(signal + noise_scale * noise, dtype=np.float64)


def _orthogonal_noise(design: _FloatArray, seed: _FloatArray) -> _FloatArray:
    projected = design @ np.linalg.lstsq(design, seed, rcond=None)[0]
    residual = seed - projected
    norm = float(np.linalg.norm(residual))
    if norm <= 0.0:
        raise AssertionError("test fixture noise seed is not independent of design")
    return np.asarray(residual / norm, dtype=np.float64)


def _matrix(values_by_site: Mapping[str, _FloatArray]) -> pd.DataFrame:
    return pd.DataFrame(
        {
            sample: [float(values[position]) for values in values_by_site.values()]
            for position, sample in enumerate(SAMPLES)
        },
        index=pd.Index(tuple(values_by_site), name="site_key"),
    )


def _proteins(values_by_protein: Mapping[str, _FloatArray]) -> pd.DataFrame:
    return pd.DataFrame(
        {
            sample: [float(values[position]) for values in values_by_protein.values()]
            for position, sample in enumerate(SAMPLES)
        },
        index=pd.Index(tuple(values_by_protein), name="total_protein_row_key"),
    )


def _pairs(mapping: Mapping[str, tuple[str, str]]) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "site_key": list(mapping),
            "protein_identifier": [value[0] for value in mapping.values()],
            "total_protein_row_key": [value[1] for value in mapping.values()],
        }
    )


def _multiple_contrasts() -> pd.DataFrame:
    return pd.DataFrame(
        {"B_vs_A": [-1.0, 1.0], "A_level": [1.0, 0.0]},
        index=pd.Index(["A", "B"], name="coefficient"),
    )


def _augmented_design(
    base_design: pd.DataFrame,
    protein: _FloatArray,
) -> _FloatArray:
    centered = protein - float(np.mean(protein))
    return np.column_stack([base_design.to_numpy(dtype=float), centered])


def _augmented_contrasts(base_contrasts: pd.DataFrame) -> _FloatArray:
    return np.vstack(
        [
            base_contrasts.to_numpy(dtype=float),
            np.zeros((1, int(base_contrasts.shape[1])), dtype=float),
        ]
    )


def _raw_oracle(
    *,
    values: _FloatArray,
    protein: _FloatArray,
    base_design: pd.DataFrame = BASE_DESIGN,
    base_contrasts: pd.DataFrame = BASE_CONTRASTS,
) -> dict[str, _FloatArray | float]:
    design = _augmented_design(base_design, protein)
    contrasts = _augmented_contrasts(base_contrasts)
    coefficients = np.linalg.lstsq(design, values, rcond=None)[0]
    fitted = design @ coefficients
    residuals = np.asarray(values - fitted, dtype=np.float64)
    rank = int(np.linalg.matrix_rank(design))
    residual_dof = float(design.shape[0] - rank)
    covariance = np.linalg.inv(design.T @ design)
    contrast_covariance = contrasts.T @ covariance @ contrasts
    return {
        "coefficients": np.asarray(coefficients, dtype=np.float64),
        "residuals": residuals,
        "residual_variance": float(residuals @ residuals / residual_dof),
        "contrast_effects": np.asarray(coefficients @ contrasts, dtype=np.float64),
        "contrast_scales": np.asarray(np.sqrt(np.diag(contrast_covariance))),
    }


def _ordinary_contrast_effect(values: _FloatArray) -> float:
    coefficients = np.linalg.lstsq(
        BASE_DESIGN.to_numpy(dtype=float),
        values,
        rcond=None,
    )[0]
    return float(coefficients @ BASE_CONTRASTS.to_numpy(dtype=float).ravel())


def _scaled_condition_number(design: _FloatArray) -> float:
    column_scales = np.linalg.norm(design, axis=0)
    scaled = design[:, column_scales > 0.0] / column_scales[column_scales > 0.0]
    _, singular_values, _ = np.linalg.svd(scaled, full_matrices=False)
    return float(singular_values[0] / singular_values[-1])


def _workflow_dataset(
    *,
    protein_ids: tuple[str, ...],
    total_protein_ids: tuple[str, ...],
    protein_values: Mapping[str, _FloatArray] | None = None,
    constant_total_proteins: frozenset[str] = frozenset(),
    protein_mapping_policy: str = "require_unambiguous",
):
    base_values = {
        "P1": np.array([0.0, 1.0, 2.0, 3.0, 0.5, 1.5, 2.5, 3.5]),
        "P2": np.array([2.0, 0.5, 3.0, 1.0, 2.8, 1.2, 3.6, 1.9]),
        "P3": np.array([1.2, 0.4, 2.3, 1.7, 2.6, 1.1, 3.0, 2.2]),
        "P4": np.array([3.0, 2.5, 4.0, 3.2, 4.5, 3.7, 5.0, 4.2]),
    }
    selected_values = {
        **base_values,
        **(dict(protein_values) if protein_values else {}),
    }
    total = pd.DataFrame(
        {
            sample: [
                (
                    7.5
                    if protein_id in constant_total_proteins
                    else float(selected_values[protein_id][position])
                )
                for protein_id in total_protein_ids
            ]
            for position, sample in enumerate(SAMPLES)
        },
        index=pd.Index(total_protein_ids, name="protein_id"),
    )
    phospho_rows = []
    for position, protein_id in enumerate(protein_ids):
        source_id = protein_id.split(";", maxsplit=1)[0]
        protein = selected_values.get(source_id, selected_values["P1"])
        phospho_rows.append(
            _site_values(
                beta={"A": 5.0 + position, "B": 5.8 + position},
                gamma=0.45,
                protein=protein,
                noise_seed=NOISE_SEEDS[position % len(NOISE_SEEDS)],
                noise_scale=0.025,
            )
        )
    phospho = pd.DataFrame(
        {
            sample: [float(values[position]) for values in phospho_rows]
            for position, sample in enumerate(SAMPLES)
        },
        index=pd.Index(
            [f"GENE{idx};S{idx};" for idx, _ in enumerate(protein_ids, start=1)],
            name="site_id",
        ),
    )
    site_metadata = pd.DataFrame(
        {
            "protein_id": list(protein_ids),
            "gene_symbol": [f"GENE{idx}" for idx in range(1, len(protein_ids) + 1)],
            "site": [f"S{idx}" for idx in range(1, len(protein_ids) + 1)],
            "site_sequence": ["A" * 15 + "S" + "A" * 15 for _ in protein_ids],
            "localisation_confidence": [0.99 for _ in protein_ids],
        },
        index=phospho.index.copy(),
    )
    return AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            total=total,
            organism=Organism.RAT,
            input_intensity_scale="log2",
            preprocessing_config=DatasetPreprocessingConfig(
                protein_aware_preparation=DatasetProteinAwarePreparationConfig(
                    policy="prepare_model_inputs",
                    protein_mapping_policy=protein_mapping_policy,
                )
            ),
        )
    )


def _single_site_workflow_dataset(
    *,
    protein_id: str,
    total_protein_id: str,
    phosphosite_values: _FloatArray,
    total_protein_values: _FloatArray,
):
    phospho = pd.DataFrame(
        {
            sample: [float(phosphosite_values[position])]
            for position, sample in enumerate(SAMPLES)
        },
        index=pd.Index(["GENE1;S1;"], name="site_id"),
    )
    total = pd.DataFrame(
        {
            sample: [float(total_protein_values[position])]
            for position, sample in enumerate(SAMPLES)
        },
        index=pd.Index([total_protein_id], name="protein_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "protein_id": [protein_id],
            "gene_symbol": ["GENE1"],
            "site": ["S1"],
            "site_sequence": ["A" * 15 + "S" + "A" * 15],
            "localisation_confidence": [0.99],
        },
        index=phospho.index.copy(),
    )
    return AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            total=total,
            organism=Organism.RAT,
            input_intensity_scale="log2",
            preprocessing_config=DatasetPreprocessingConfig(
                protein_aware_preparation=DatasetProteinAwarePreparationConfig(
                    policy="prepare_model_inputs",
                    protein_mapping_policy="require_unambiguous",
                )
            ),
        )
    )


def _workflow_request(
    dataset: object,
    *,
    sample_ids: tuple[str, ...] = SAMPLES,
    allow_design_subset: bool = False,
    paired_design_policy: str = "reject",
    technical_replicate_policy: TechnicalReplicatePolicy = TechnicalReplicatePolicy.REJECT,
    technical_replicates: bool = False,
) -> DifferentialAnalysisRequest:
    return DifferentialAnalysisRequest(
        dataset=dataset,
        design=_workflow_design(
            sample_ids=sample_ids,
            paired_design_policy=paired_design_policy,
            technical_replicates=technical_replicates,
        ),
        contrasts=(
            Contrast(
                name="B_vs_A",
                numerator_condition="B",
                denominator_condition="A",
            ),
        ),
        config=DifferentialAnalysisConfig(
            allow_design_subset=allow_design_subset,
            paired_design_policy=paired_design_policy,  # type: ignore[arg-type]
            technical_replicate_policy=technical_replicate_policy,
            protein_aware_model=DifferentialProteinAwareModelConfig(),
        ),
    )


def _ordinary_workflow_request(dataset: object) -> DifferentialAnalysisRequest:
    return DifferentialAnalysisRequest(
        dataset=dataset,
        design=_workflow_design(),
        contrasts=(
            Contrast(
                name="B_vs_A",
                numerator_condition="B",
                denominator_condition="A",
            ),
        ),
    )


def _workflow_design(
    *,
    sample_ids: tuple[str, ...] = SAMPLES,
    paired_design_policy: str = "reject",
    technical_replicates: bool = False,
) -> ExperimentalDesign:
    return ExperimentalDesign(
        samples=tuple(
            SampleDesignRecord(
                sample_id=sample_id,
                condition=sample_id.split("_", maxsplit=1)[0],
                biological_replicate_id=_biological_replicate_id(
                    sample_id,
                    technical_replicates=technical_replicates,
                ),
                technical_replicate_id=(
                    _technical_replicate_id(sample_id) if technical_replicates else None
                ),
                block_id=(
                    f"pair_{sample_id.split('_', maxsplit=1)[1]}"
                    if paired_design_policy
                    in {
                        PAIRED_DESIGN_POLICY_DUPLICATE_CORRELATION,
                        PAIRED_DESIGN_POLICY_FIXED_BLOCK,
                    }
                    else None
                ),
            )
            for sample_id in sample_ids
        )
    )


def _workflow_design_with_continuous_covariate(
    covariate_name: str,
    values: _FloatArray,
) -> ExperimentalDesign:
    return ExperimentalDesign(
        samples=tuple(
            SampleDesignRecord(
                sample_id=sample_id,
                condition=sample_id.split("_", maxsplit=1)[0],
                biological_replicate_id=f"{sample_id}_bio",
                covariates={covariate_name: float(values[position])},
            )
            for position, sample_id in enumerate(SAMPLES)
        ),
        fixed_effects=(FixedEffectCovariate.continuous(covariate_name),),
    )


def _biological_replicate_id(
    sample_id: str,
    *,
    technical_replicates: bool,
) -> str:
    if not technical_replicates:
        return f"{sample_id}_bio"
    condition, replicate = sample_id.split("_", maxsplit=1)
    if replicate in {"1", "2"}:
        return f"{condition}_bio_technical_pair"
    return f"{sample_id}_bio"


def _technical_replicate_id(sample_id: str) -> str:
    return f"tech_{sample_id.split('_', maxsplit=1)[1]}"


def _site_key_for_metadata_value(
    dataset: object,
    column_name: str,
    value: str,
) -> str:
    metadata = cast(Any, dataset).site_metadata
    matches = metadata.index[metadata.loc[:, column_name].astype(str) == value]
    if len(matches) != 1:
        raise AssertionError(f"expected one site with {column_name}={value!r}")
    return str(matches[0])


def _protein_aware_policy(result: object) -> Any:
    policy = cast(Any, result).policy_provenance
    if policy is None or policy.protein_aware is None:
        raise AssertionError("protein-aware policy provenance is required")
    return policy.protein_aware


def _expected_input_fingerprints(
    request: DifferentialAnalysisRequest,
) -> DifferentialProteinAwareInputFingerprints:
    dataset = cast(Any, request.dataset)
    sidecar = dataset.protein_aware_preparation
    if sidecar is None:
        raise AssertionError("expected dataset-owned protein-aware preparation")
    sample_order = request.design.sample_ids()
    phospho_matrix = dataset.phospho.loc[:, list(sample_order)]
    matched_pairs = _expected_matched_pairs(
        sidecar.matched_pairs_dataframe(),
        site_order=tuple(str(value) for value in phospho_matrix.index.tolist()),
    )
    protein_covariates = sidecar.protein_covariate_matrix_dataframe().loc[
        _distinct_values(matched_pairs["total_protein_row_key"]),
        list(sample_order),
    ]
    design_matrix = _expected_public_design_matrix(
        request.design,
        paired_design_policy=request.config.paired_design_policy,
    )
    contrast_matrix = _expected_public_contrast_matrix(
        coefficient_labels=tuple(str(value) for value in design_matrix.columns),
        contrasts=request.contrasts,
    )
    site_eligibility_metadata = _expected_all_tested_site_eligibility_metadata(
        phospho_matrix=phospho_matrix,
        matched_pairs=matched_pairs,
        protein_covariates=protein_covariates,
        sidecar_site_eligibility=sidecar.site_eligibility_table,
        design_matrix=design_matrix,
        method_id=str(request.config.protein_aware_model.method),
    )
    return DifferentialProteinAwareInputFingerprints(
        phospho_matrix=fingerprint_table_strict(
            phospho_matrix,
            name=DIFFERENTIAL_PROTEIN_AWARE_PHOSPHO_MATRIX_FINGERPRINT_NAME,
        ),
        protein_matched_pairs=fingerprint_table_strict(
            matched_pairs,
            name=DIFFERENTIAL_PROTEIN_AWARE_MATCHED_PAIRS_FINGERPRINT_NAME,
        ),
        protein_covariate_matrix=fingerprint_table_strict(
            protein_covariates,
            name=DIFFERENTIAL_PROTEIN_AWARE_COVARIATE_MATRIX_FINGERPRINT_NAME,
        ),
        protein_site_eligibility=fingerprint_table_strict(
            site_eligibility_metadata,
            name=DIFFERENTIAL_PROTEIN_AWARE_SITE_ELIGIBILITY_FINGERPRINT_NAME,
        ),
        design_matrix=fingerprint_table_strict(
            design_matrix,
            name=DIFFERENTIAL_PROTEIN_AWARE_DESIGN_MATRIX_FINGERPRINT_NAME,
        ),
        contrast_matrix=fingerprint_table_strict(
            contrast_matrix,
            name=DIFFERENTIAL_PROTEIN_AWARE_CONTRAST_MATRIX_FINGERPRINT_NAME,
        ),
    )


def _expected_matched_pairs(
    matched_pairs: pd.DataFrame,
    *,
    site_order: tuple[str, ...],
) -> pd.DataFrame:
    by_site = matched_pairs.set_index("site_key", drop=False)
    expected = by_site.loc[
        list(site_order),
        ["site_key", "protein_identifier", "total_protein_row_key"],
    ].reset_index(drop=True)
    return pd.DataFrame(
        expected.to_numpy(dtype=object),
        columns=["site_key", "protein_identifier", "total_protein_row_key"],
    )


def _distinct_values(values: pd.Series) -> list[str]:
    distinct: list[str] = []
    seen: set[str] = set()
    for value in values.astype(str).tolist():
        if value in seen:
            continue
        seen.add(value)
        distinct.append(value)
    return distinct


def _expected_public_design_matrix(
    design: ExperimentalDesign,
    *,
    paired_design_policy: str,
) -> pd.DataFrame:
    records = design.samples
    condition_labels = _condition_labels(records)
    data: dict[str, list[float]] = {
        condition: [1.0 if record.condition == condition else 0.0 for record in records]
        for condition in condition_labels
    }
    for covariate in design.fixed_effects:
        if not covariate.include_in_model:
            continue
        if str(covariate.kind) != "continuous":
            raise AssertionError(
                "test helper only constructs expected continuous covariates"
            )
        data[covariate.name] = [
            float(record.covariates[covariate.name]) for record in records
        ]
    if paired_design_policy == PAIRED_DESIGN_POLICY_FIXED_BLOCK:
        block_values = [
            str(record.block_id) for record in records if record.block_id is not None
        ]
        if len(block_values) != len(records):
            raise AssertionError("fixed-block test design must include every block_id")
        for level in sorted(set(block_values))[1:]:
            data[f"block[{level}]"] = [
                1.0 if record.block_id == level else 0.0 for record in records
            ]
    frame = pd.DataFrame(
        data,
        index=pd.Index([record.sample_id for record in records], name="sample"),
        dtype=float,
    )
    frame.columns = pd.Index(tuple(data), name="coefficient")
    return frame


def _condition_labels(records: tuple[SampleDesignRecord, ...]) -> tuple[str, ...]:
    labels: list[str] = []
    seen: set[str] = set()
    for record in records:
        if record.condition in seen:
            continue
        seen.add(record.condition)
        labels.append(record.condition)
    return tuple(labels)


def _expected_public_contrast_matrix(
    *,
    coefficient_labels: tuple[str, ...],
    contrasts: tuple[Contrast, ...],
) -> pd.DataFrame:
    data: dict[str, list[float]] = {}
    for contrast in contrasts:
        weights = []
        for coefficient in coefficient_labels:
            if coefficient == contrast.numerator_condition:
                weights.append(1.0)
            elif coefficient == contrast.denominator_condition:
                weights.append(-1.0)
            else:
                weights.append(0.0)
        data[contrast.name] = weights
    frame = pd.DataFrame(
        data,
        index=pd.Index(coefficient_labels, name="coefficient"),
        dtype=float,
    )
    frame.columns = pd.Index(tuple(data), name="contrast")
    return frame


def _expected_all_tested_site_eligibility_metadata(
    *,
    phospho_matrix: pd.DataFrame,
    matched_pairs: pd.DataFrame,
    protein_covariates: pd.DataFrame,
    sidecar_site_eligibility: pd.DataFrame,
    design_matrix: pd.DataFrame,
    method_id: str,
) -> pd.DataFrame:
    sidecar_by_site = sidecar_site_eligibility.set_index("site_key", drop=False)
    pairs_by_site = matched_pairs.set_index("site_key", drop=False)
    columns: dict[str, list[Any]] = {
        "site_key": [],
        "analysed_value_count": [],
        "observed_value_count": [],
        "invalid_numeric_value_count": [],
        "unique_observed_value_count": [],
        "protein_aware_method_id": [],
        "protein_aware_centering_policy": [],
        "protein_aware_candidate": [],
        "protein_aware_tested": [],
        "protein_aware_preparation_eligibility": [],
        "protein_aware_preparation_reasons": [],
        "protein_aware_mapping_status": [],
        "protein_identifier": [],
        "total_protein_row_key": [],
        "protein_covariate_raw_mean": [],
        "protein_covariate_raw_standard_deviation": [],
        "protein_covariate_centered_variance": [],
        "protein_augmented_design_rank": [],
        "protein_augmented_design_residual_degrees_of_freedom": [],
        "protein_augmented_design_condition_number": [],
        "protein_augmented_design_max_condition_number": [],
        "protein_aware_failure_message": [],
        DIFFERENTIAL_RESULT_STATUS_COLUMN: [],
        DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN: [],
    }
    for site_key in phospho_matrix.index.astype(str).tolist():
        values = phospho_matrix.loc[site_key].to_numpy(dtype=float)
        finite = np.isfinite(values)
        if not finite.all() or int(np.unique(values[finite]).size) <= 1:
            raise AssertionError("fingerprint helper expects all sites to be tested")
        pair = pairs_by_site.loc[site_key]
        total_row_key = str(pair["total_protein_row_key"])
        protein = protein_covariates.loc[total_row_key].to_numpy(dtype=float)
        centered = protein - float(np.mean(protein))
        augmented_frame = design_matrix.copy(deep=True)
        augmented_frame.loc[:, PROTEIN_AWARE_COVARIATE_COEFFICIENT_NAME] = centered
        augmented = augmented_frame.to_numpy(dtype=float)
        rank = int(np.linalg.matrix_rank(augmented))
        sidecar_row = sidecar_by_site.loc[site_key]

        columns["site_key"].append(site_key)
        columns["analysed_value_count"].append(int(values.size))
        columns["observed_value_count"].append(int(finite.sum()))
        columns["invalid_numeric_value_count"].append(int((~finite).sum()))
        columns["unique_observed_value_count"].append(int(np.unique(values).size))
        columns["protein_aware_method_id"].append(method_id)
        columns["protein_aware_centering_policy"].append(PROTEIN_AWARE_CENTERING_POLICY)
        columns["protein_aware_candidate"].append(True)
        columns["protein_aware_tested"].append(True)
        columns["protein_aware_preparation_eligibility"].append(
            str(sidecar_row["eligibility"])
        )
        columns["protein_aware_preparation_reasons"].append(
            tuple(sidecar_row["reasons"])
        )
        columns["protein_aware_mapping_status"].append(
            str(sidecar_row["mapping_status"])
        )
        columns["protein_identifier"].append(str(pair["protein_identifier"]))
        columns["total_protein_row_key"].append(total_row_key)
        columns["protein_covariate_raw_mean"].append(float(np.mean(protein)))
        columns["protein_covariate_raw_standard_deviation"].append(
            float(np.std(protein, ddof=1))
        )
        columns["protein_covariate_centered_variance"].append(
            float(np.var(centered, ddof=1))
        )
        columns["protein_augmented_design_rank"].append(rank)
        columns["protein_augmented_design_residual_degrees_of_freedom"].append(
            float(len(protein) - rank)
        )
        columns["protein_augmented_design_condition_number"].append(
            _scaled_condition_number(augmented)
        )
        columns["protein_augmented_design_max_condition_number"].append(
            DIFFERENTIAL_LINEAR_MODEL_MAX_CONDITION_NUMBER
        )
        columns["protein_aware_failure_message"].append("")
        columns[DIFFERENTIAL_RESULT_STATUS_COLUMN].append(
            DIFFERENTIAL_RESULT_STATUS_TESTED
        )
        columns[DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN].append(
            "Feature has finite, non-constant values for the differential design "
            "samples."
        )
    return pd.DataFrame(columns, index=phospho_matrix.index.copy())
