from __future__ import annotations

import json

import pandas as pd

from phospy import AnalysisReadyPhosphoDataset, DifferentialAnalysisWorkflow
from phospy.api import (
    BatchCovariate,
    Contrast,
    DifferentialAnalysisRequest,
    ExperimentalDesign,
    Organism,
    SampleDesignRecord,
)
from phospy.science.datasets.models import DatasetPreprocessingReport
from phospy.science.datasets.preprocessing.batch_correction import (
    BatchCorrectionDiagnostics,
    BatchCorrectionPolicy,
    BatchCorrectionReport,
)
from tests.support.intensity_scale_states import (
    supported_log2_intensity_scale_state,
    supported_log2_processing_state,
)
from tests.support.site_keys import protein_site_key_index, site_key_context_columns


def _dataset(*, with_batch_correction: bool = False) -> AnalysisReadyPhosphoDataset:
    genes = ["MAPK14", "GSK3B", "AKT1"]
    sites = ["Y182", "S9", "T308"]
    site_index = protein_site_key_index(protein_identifiers=genes, sites=sites)
    phospho = pd.DataFrame(
        {
            "A_1": [1.0, 2.0, 0.9],
            "A_2": [1.2, 2.1, 1.1],
            "B_1": [2.0, 1.8, 0.8],
            "B_2": [2.2, 2.0, 1.0],
        },
        index=site_index.copy(),
    )
    site_metadata = pd.DataFrame(
        {
            "site_key": site_index.tolist(),
            "display_id": ["MAPK14;Y182;", "GSK3B;S9;", "AKT1;T308;"],
            **site_key_context_columns(site_index),
            "gene_symbol": genes,
            "site": sites,
            "site_sequence": [
                ("A" * 15) + str(site).strip().upper()[0] + ("A" * 15) for site in sites
            ],
            "protein_id": genes,
        },
        index=site_index.copy(),
    )
    preprocessing_report = None
    if with_batch_correction:
        preprocessing_report = DatasetPreprocessingReport.from_rows(
            batch_correction=BatchCorrectionReport(
                status="applied",
                policy=BatchCorrectionPolicy(
                    method="linear_residualize_batch",
                    batch_column="batch",
                    condition_column="condition",
                ),
                diagnostics=BatchCorrectionDiagnostics(
                    number_of_batches=2,
                    batch_levels=("run_1", "run_2"),
                    condition_levels=("A", "B"),
                    confounding_check_status="passed",
                    matrix_shape_before=phospho.shape,
                    matrix_shape_after=phospho.shape,
                ),
            )
        )
    return AnalysisReadyPhosphoDataset(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        intensity_scale_state=supported_log2_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_log2_processing_state(has_total_matrix=False),
        preprocessing_report=preprocessing_report,
    )


def _design(*, with_batch_covariate: bool = False) -> ExperimentalDesign:
    fixed_effects = (BatchCovariate(),) if with_batch_covariate else ()
    batches = (
        ("run_1", "run_2", "run_1", "run_2") if with_batch_covariate else (None,) * 4
    )
    return ExperimentalDesign(
        samples=(
            SampleDesignRecord(
                sample_id="A_1",
                condition="A",
                biological_replicate_id="A_r1",
                batch=batches[0],
            ),
            SampleDesignRecord(
                sample_id="A_2",
                condition="A",
                biological_replicate_id="A_r2",
                batch=batches[1],
            ),
            SampleDesignRecord(
                sample_id="B_1",
                condition="B",
                biological_replicate_id="B_r1",
                batch=batches[2],
            ),
            SampleDesignRecord(
                sample_id="B_2",
                condition="B",
                biological_replicate_id="B_r2",
                batch=batches[3],
            ),
        ),
        fixed_effects=fixed_effects,
    )


def _contrast() -> tuple[Contrast, ...]:
    return (
        Contrast(
            name="B_vs_A",
            numerator_condition="B",
            denominator_condition="A",
        ),
    )


def test_differential_result_diagnostics_are_present_and_serialized() -> None:
    result = DifferentialAnalysisWorkflow().run(
        DifferentialAnalysisRequest(
            dataset=_dataset(),
            design=_design(),
            contrasts=_contrast(),
        )
    )

    diagnostics = result.diagnostics
    assert diagnostics.model_type == "moderated_ols_fixed_effect"
    assert diagnostics.design_columns == ("A", "B")
    assert [definition.name for definition in diagnostics.contrast_definitions] == [
        "B_vs_A"
    ]
    assert diagnostics.rank == 2
    assert diagnostics.n_samples == 4
    assert diagnostics.n_sites == 3
    assert diagnostics.residual_degrees_of_freedom == 2.0
    assert diagnostics.variance_method == "ordinary_least_squares_residual_variance"
    assert diagnostics.moderation_method == "empirical_bayes_standard"
    assert diagnostics.multiple_testing_method == "benjamini_hochberg"
    assert diagnostics.imputation_policy == "reject"
    assert (
        diagnostics.missing_value_policy
        == "reject_missing_values_before_differential_execution"
    )
    assert diagnostics.intensity_scale == "log2"
    assert diagnostics.normalisation_state == "none"
    assert "mixed-effects differential modelling" in " ".join(
        diagnostics.unsupported_assumptions
    )

    payload = result.to_payload()
    assert payload["diagnostics"] == diagnostics.to_payload()
    json.dumps(payload, allow_nan=False)


def test_differential_result_diagnostics_warn_on_batch_assumption_boundaries() -> None:
    result = DifferentialAnalysisWorkflow().run(
        DifferentialAnalysisRequest(
            dataset=_dataset(with_batch_correction=True),
            design=_design(with_batch_covariate=True),
            contrasts=_contrast(),
        )
    )

    diagnostics = result.diagnostics
    assert diagnostics.batch_or_covariate_terms == ("batch[run_2]",)
    warning_text = " ".join(diagnostics.warnings).lower()
    unsupported_text = " ".join(diagnostics.unsupported_assumptions).lower()
    assert "fixed-effect covariates" in warning_text
    assert "upstream batch correction" in warning_text
    assert "not full batch correction" in unsupported_text
    assert "not revalidated or rerun" in unsupported_text
