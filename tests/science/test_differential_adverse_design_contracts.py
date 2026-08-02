from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from phospy.api import (
    BatchCovariate,
    ContinuousCovariate,
    Contrast,
    DifferentialAnalysisConfig,
    DifferentialAnalysisRequest,
    DifferentialAnalysisWorkflow,
    EmpiricalBayesConfig,
    ExperimentalDesign,
    Organism,
    SampleDesignRecord,
)
from phospy.api.configs import (
    DIFFERENTIAL_RELIABILITY_PROFILE_EXPLORATORY_SINGLE_REPLICATE,
)
from phospy.errors import WorkflowValidationError
from phospy.science.differential.empirical_bayes import fit_empirical_bayes
from tests.support.analysis_ready_dataset_factories import (
    trusted_analysis_ready_dataset_from_tables,
)
from tests.support.intensity_scale_states import (
    supported_log2_intensity_scale_state,
    supported_log2_processing_state,
)
from tests.support.site_keys import protein_site_key_index, site_key_context_columns

pytestmark = pytest.mark.release_gate


class _ExecutorSpy:
    def __init__(self) -> None:
        self.calls = 0

    def run(self, request: object) -> object:
        del request
        self.calls += 1
        raise AssertionError("executor should not be called")


def _dataset(
    matrix: pd.DataFrame,
    *,
    genes: tuple[str, ...] = ("MAPK14", "GSK3B", "AKT1", "RPS6KB1"),
    sites: tuple[str, ...] = ("Y182", "S9", "T308", "T389"),
):
    selected_genes = genes[: matrix.shape[0]]
    selected_sites = sites[: matrix.shape[0]]
    site_index = protein_site_key_index(
        protein_identifiers=selected_genes,
        sites=selected_sites,
    )
    phospho = matrix.copy(deep=True)
    phospho.index = site_index
    site_metadata = pd.DataFrame(
        {
            "site_key": site_index.tolist(),
            "display_id": [
                f"{gene};{site};"
                for gene, site in zip(selected_genes, selected_sites, strict=True)
            ],
            **site_key_context_columns(site_index),
            "gene_symbol": list(selected_genes),
            "site": list(selected_sites),
            "site_sequence": [
                ("A" * 15) + site[0].upper() + ("A" * 15) for site in selected_sites
            ],
            "protein_id": list(selected_genes),
            "localisation_confidence": [0.95] * int(matrix.shape[0]),
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


def _contrast() -> tuple[Contrast, ...]:
    return (
        Contrast(
            name="B_vs_A",
            numerator_condition="B",
            denominator_condition="A",
        ),
    )


def _run(
    matrix: pd.DataFrame,
    design: ExperimentalDesign,
    *,
    minimum_condition_replicates: int = 2,
    reliability_profile: str = "production",
    empirical_bayes: EmpiricalBayesConfig | None = None,
):
    return DifferentialAnalysisWorkflow().run(
        DifferentialAnalysisRequest(
            dataset=_dataset(matrix),
            design=design,
            contrasts=_contrast(),
            config=DifferentialAnalysisConfig(
                reliability_profile=reliability_profile,  # type: ignore[arg-type]
                minimum_condition_replicates=minimum_condition_replicates,
                empirical_bayes=(
                    EmpiricalBayesConfig()
                    if empirical_bayes is None
                    else empirical_bayes
                ),
            ),
        )
    )


def test_unbalanced_group_design_runs_with_finite_statistics() -> None:
    samples = [f"A_{idx}" for idx in range(1, 4)] + [f"B_{idx}" for idx in range(1, 6)]
    matrix = pd.DataFrame(
        {
            sample: [
                1.0 + 0.05 * position + (0.55 if sample.startswith("B") else 0.0),
                2.0 - 0.03 * position + (0.10 if sample.startswith("B") else 0.0),
                0.5 + 0.02 * position - (0.30 if sample.startswith("B") else 0.0),
                1.2 + 0.04 * position + (0.20 if sample.startswith("B") else 0.0),
            ]
            for position, sample in enumerate(samples)
        },
        index=pd.Index(["site1", "site2", "site3", "site4"], name="site_id"),
    )
    design = ExperimentalDesign(
        samples=tuple(
            SampleDesignRecord(
                sample_id=sample,
                condition=sample.split("_", maxsplit=1)[0],
                biological_replicate_id=sample,
            )
            for sample in samples
        )
    )

    result = _run(matrix, design)
    table = result.table_for("B_vs_A")

    assert result.residual_degrees_of_freedom == pytest.approx(6.0)
    assert result.policy_provenance is not None
    assert result.policy_provenance.replicates.condition_replicate_counts == (
        ("A", 3),
        ("B", 5),
    )
    assert np.isfinite(table.loc[:, ["logFC", "t", "P.Value", "adj.P.Val"]]).all().all()


def test_continuous_covariate_adjustment_is_executable_and_recorded() -> None:
    samples = ("A_1", "A_2", "A_3", "B_1", "B_2", "B_3")
    dose = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0], dtype=float)
    condition_b = np.array([0.0, 0.0, 0.0, 1.0, 1.0, 1.0], dtype=float)
    residual = np.array([0.02, -0.04, 0.03, -0.01, 0.04, -0.02], dtype=float)
    matrix = pd.DataFrame(
        {
            sample: [
                10.0 + 2.0 * dose[position] + residual[position],
                5.0
                + condition_b[position]
                + 0.5 * dose[position]
                + residual[::-1][position],
                1.0 + 0.2 * position + 0.01 * (-1) ** position,
                3.0
                + 0.3 * condition_b[position]
                - 0.2 * dose[position]
                + residual[position] * 0.5,
            ]
            for position, sample in enumerate(samples)
        },
        index=pd.Index(["site1", "site2", "site3", "site4"], name="site_id"),
    )
    unadjusted_design = ExperimentalDesign(
        samples=tuple(
            SampleDesignRecord(
                sample_id=sample,
                condition=sample.split("_", maxsplit=1)[0],
                biological_replicate_id=sample,
            )
            for sample in samples
        )
    )
    adjusted_design = ExperimentalDesign(
        samples=tuple(
            SampleDesignRecord(
                sample_id=sample,
                condition=sample.split("_", maxsplit=1)[0],
                biological_replicate_id=sample,
                covariates={"dose": float(dose[position])},
            )
            for position, sample in enumerate(samples)
        ),
        fixed_effects=(ContinuousCovariate("dose"),),
    )

    unadjusted = _run(matrix, unadjusted_design)
    adjusted = _run(matrix, adjusted_design)

    assert adjusted.policy_provenance is not None
    assert adjusted.policy_provenance.design.covariates[0].name == "dose"
    assert adjusted.policy_provenance.design.coefficient_labels == ("A", "B", "dose")
    assert abs(unadjusted.table_for("B_vs_A").iloc[0]["logFC"]) > 5.0
    assert abs(adjusted.table_for("B_vs_A").iloc[0]["logFC"]) < 0.1


def test_near_rank_deficient_estimable_contrast_runs_without_normal_equations() -> None:
    samples = ("A_1", "A_2", "A_3", "B_1", "B_2", "B_3")
    condition_b = np.array([0.0, 0.0, 0.0, 1.0, 1.0, 1.0], dtype=float)
    almost_b = condition_b + np.array(
        [-1.0e-6, 0.0, 1.0e-6, -1.0e-6, 0.0, 1.0e-6],
        dtype=float,
    )
    residual = np.array([0.02, -0.01, 0.03, -0.02, 0.01, -0.03], dtype=float)
    matrix = pd.DataFrame(
        {
            sample: [
                2.0
                + 0.4 * condition_b[position]
                + 0.7 * almost_b[position]
                + residual[position],
                -1.0
                + 0.2 * condition_b[position]
                - 0.3 * almost_b[position]
                + residual[::-1][position],
                0.5
                + 0.1 * position
                + 0.2 * almost_b[position]
                + residual[position] * 0.5,
                1.5
                - 0.1 * position
                + 0.1 * condition_b[position]
                + residual[::-1][position] * 0.25,
            ]
            for position, sample in enumerate(samples)
        },
        index=pd.Index(["site1", "site2", "site3", "site4"], name="site_id"),
    )
    design = ExperimentalDesign(
        samples=tuple(
            SampleDesignRecord(
                sample_id=sample,
                condition=sample.split("_", maxsplit=1)[0],
                biological_replicate_id=sample,
                covariates={"almost_B": float(almost_b[position])},
            )
            for position, sample in enumerate(samples)
        ),
        fixed_effects=(ContinuousCovariate("almost_B"),),
    )

    result = _run(matrix, design)
    table = result.table_for("B_vs_A")

    assert result.policy_provenance is not None
    assert result.policy_provenance.design.rank == 3
    assert result.policy_provenance.design.condition_number > 1.0e5
    assert np.isfinite(table.loc[:, ["logFC", "t", "P.Value", "adj.P.Val"]]).all().all()


def test_non_estimable_contrast_design_fails_before_numerical_fitting() -> None:
    matrix = pd.DataFrame(
        {
            "A_1": [1.0, 2.0, 3.0, 4.0],
            "A_2": [1.1, 2.1, 3.1, 4.1],
            "B_1": [1.7, 2.2, 3.3, 4.4],
            "B_2": [1.8, 2.3, 3.4, 4.5],
        },
        index=pd.Index(["site1", "site2", "site3", "site4"], name="site_id"),
    )
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(sample_id="A_1", condition="A", batch="batch_1"),
            SampleDesignRecord(sample_id="A_2", condition="A", batch="batch_1"),
            SampleDesignRecord(sample_id="B_1", condition="B", batch="batch_2"),
            SampleDesignRecord(sample_id="B_2", condition="B", batch="batch_2"),
        ),
        fixed_effects=(
            # Batch is perfectly confounded with condition, so B - A is not
            # estimable in this explicit design.
            BatchCovariate(),
        ),
    )
    executor = _ExecutorSpy()

    with pytest.raises(WorkflowValidationError, match="rank deficient.*confounded"):
        DifferentialAnalysisWorkflow._with_components(executor=executor).run(
            DifferentialAnalysisRequest(
                dataset=_dataset(matrix),
                design=design,
                contrasts=_contrast(),
            )
        )
    assert executor.calls == 0


def test_very_small_positive_residual_degrees_of_freedom_runs_with_guarded_output() -> (
    None
):
    matrix = pd.DataFrame(
        {
            "A_1": [1.00, 2.00, 3.00, 4.00],
            "A_2": [1.10, 2.05, 3.10, 4.10],
            "B_1": [1.65, 2.30, 3.20, 4.45],
        },
        index=pd.Index(["site1", "site2", "site3", "site4"], name="site_id"),
    )
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(sample_id="A_1", condition="A"),
            SampleDesignRecord(sample_id="A_2", condition="A"),
            SampleDesignRecord(sample_id="B_1", condition="B"),
        )
    )

    result = _run(
        matrix,
        design,
        minimum_condition_replicates=1,
        reliability_profile=(
            DIFFERENTIAL_RELIABILITY_PROFILE_EXPLORATORY_SINGLE_REPLICATE
        ),
        empirical_bayes=EmpiricalBayesConfig(method="standard", trend=False),
    )

    assert result.residual_degrees_of_freedom == pytest.approx(1.0)
    assert (
        np.isfinite(
            result.table_for("B_vs_A").loc[:, ["logFC", "t", "P.Value", "adj.P.Val"]]
        )
        .all()
        .all()
    )


def test_empirical_bayes_trend_above_1024_features_uses_large_branch_contract() -> None:
    n_features = 1500
    mean_intensity = np.linspace(7.5, 13.5, n_features, dtype=float)
    log_variances = (
        -1.4
        + 0.18 * np.sin(mean_intensity * 1.4)
        - 0.04 * (mean_intensity - float(np.mean(mean_intensity)))
    )

    result = fit_empirical_bayes(
        variances=np.exp(log_variances),
        residual_dof=8.0,
        method="standard",
        trend=True,
        winsor_tail_p=(0.05, 0.10),
        mean_intensity=mean_intensity,
    )

    assert result.prior_variance.shape == (n_features,)
    assert result.fitted_log_prior_variance is not None
    assert result.fitted_log_prior_variance.shape == (n_features,)
    assert np.isfinite(result.prior_variance).all()
    assert np.ptp(result.fitted_log_prior_variance) > 0.1
