from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from phospy import AnalysisReadyPhosphoDataset, DifferentialAnalysisWorkflow
from phospy.api import (
    Contrast,
    DifferentialAnalysisRequest,
    ExperimentalDesign,
    Organism,
    SampleDesignRecord,
)
from tests.support.intensity_scale_states import (
    supported_log2_intensity_scale_state,
    supported_log2_processing_state,
)
from tests.support.site_keys import (
    site_key_context_columns,
    site_key_index_from_display_ids,
)

ROOT = Path(__file__).resolve().parents[2]
PARITY_FIXTURE_ROOT = ROOT / "tests" / "fixtures" / "rewrite_parity"

pytestmark = [pytest.mark.golden, pytest.mark.reproducibility, pytest.mark.release_gate]


def _require_common_provenance_fields(path: Path) -> None:
    text = path.read_text(encoding="utf-8")
    assert "Generated with R version" in text
    assert "limma version:" in text
    assert "Design:" in text
    assert "Contrasts" in text


def _require_expected_schema(path: Path) -> None:
    frame = pd.read_csv(path)
    assert list(frame.columns) == ["site_id", "logFC", "t", "P.Value", "adj.P.Val"]


def test_differential_limma_fixture_provenance_is_source_labelled() -> None:
    base_fixture = PARITY_FIXTURE_ROOT / "differential_r_reference" / "PROVENANCE.md"
    envelope_fixture = (
        PARITY_FIXTURE_ROOT / "differential_limma_envelope" / "PROVENANCE.md"
    )

    assert base_fixture.is_file()
    assert envelope_fixture.is_file()
    _require_common_provenance_fields(base_fixture)
    _require_common_provenance_fields(envelope_fixture)


def test_differential_limma_expected_tables_keep_stable_schema() -> None:
    _require_expected_schema(
        PARITY_FIXTURE_ROOT / "differential_r_reference" / "limma_B_vs_A.csv"
    )
    _require_expected_schema(
        PARITY_FIXTURE_ROOT / "differential_r_reference" / "limma_C_vs_A.csv"
    )
    _require_expected_schema(
        PARITY_FIXTURE_ROOT / "differential_limma_envelope" / "limma_B_vs_A.csv"
    )
    _require_expected_schema(
        PARITY_FIXTURE_ROOT / "differential_limma_envelope" / "limma_A_vs_B.csv"
    )


def _dataset() -> AnalysisReadyPhosphoDataset:
    display_ids = ["MAPK14;Y182;", "GSK3B;S9;"]
    site_index = site_key_index_from_display_ids(
        display_ids,
        protein_namespace="gene_symbol",
    )
    phospho = pd.DataFrame(
        {
            "A_1": [1.0, 2.0],
            "A_2": [1.1, 2.1],
            "B_1": [2.0, 1.9],
            "B_2": [2.2, 2.2],
        },
        index=site_index,
    )
    site_metadata = pd.DataFrame(
        {
            "site_key": site_index.tolist(),
            "display_id": display_ids,
            **site_key_context_columns(site_index),
            "gene_symbol": ["MAPK14", "GSK3B"],
            "site": ["Y182", "S9"],
            "site_sequence": [
                ("A" * 15) + "Y" + ("A" * 15),
                ("A" * 15) + "S" + ("A" * 15),
            ],
            "protein_id": ["MAPK14", "GSK3B"],
        },
        index=site_index.copy(),
    )
    return AnalysisReadyPhosphoDataset(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        intensity_scale_state=supported_log2_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_log2_processing_state(has_total_matrix=False),
    )


def _request() -> DifferentialAnalysisRequest:
    return DifferentialAnalysisRequest(
        dataset=_dataset(),
        design=ExperimentalDesign(
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
        ),
        contrasts=(
            Contrast(
                name="B_vs_A",
                numerator_condition="B",
                denominator_condition="A",
            ),
        ),
    )


def test_differential_policy_provenance_snapshot_is_stable() -> None:
    result = DifferentialAnalysisWorkflow().run(_request())
    assert result.policy_provenance is not None
    policy = result.policy_provenance
    snapshot = {
        "design": {
            "formula": policy.design.formula,
            "description": policy.design.description,
            "sample_labels": list(policy.design.sample_labels),
            "coefficient_labels": list(policy.design.coefficient_labels),
            "condition_columns": list(policy.design.condition_columns),
            "covariates": [
                {
                    "name": covariate.name,
                    "kind": covariate.kind,
                    "columns": list(covariate.columns),
                    "levels": list(covariate.levels),
                    "reference_level": covariate.reference_level,
                    "unused_levels": list(covariate.unused_levels),
                }
                for covariate in policy.design.covariates
            ],
            "sample_count": policy.design.sample_count,
            "coefficient_count": policy.design.coefficient_count,
            "rank": policy.design.rank,
            "residual_degrees_of_freedom": policy.design.residual_degrees_of_freedom,
            "paired_design_policy": policy.design.paired_design_policy,
            "block_id_field_name": policy.design.block_id_field_name,
            "block_count": policy.design.block_count,
            "block_levels": list(policy.design.block_levels),
            "block_levels_included": list(policy.design.block_levels_included),
            "block_reference_level": policy.design.block_reference_level,
            "block_columns": list(policy.design.block_columns),
            "block_column_names": list(policy.design.block_column_names),
            "condition_coverage_rule": policy.design.condition_coverage_rule,
            "limitations": list(policy.design.limitations),
            "rank_validation_status": policy.design.rank_validation_status,
            "estimability_validation_status": (
                policy.design.estimability_validation_status
            ),
        },
        "contrasts": [
            {
                "name": contrast.name,
                "numerator_condition": contrast.numerator_condition,
                "denominator_condition": contrast.denominator_condition,
                "coefficients": list(contrast.coefficients),
                "description": contrast.description,
            }
            for contrast in policy.contrasts
        ],
        "replicates": {
            "minimum_condition_replicates": (
                policy.replicates.minimum_condition_replicates
            ),
            "technical_replicate_policy": policy.replicates.technical_replicate_policy,
            "condition_replicate_counts": list(
                policy.replicates.condition_replicate_counts
            ),
            "technical_replicate_groups": [
                {
                    "condition": group.condition,
                    "biological_replicate_id": group.biological_replicate_id,
                    "output_sample_id": group.output_sample_id,
                    "input_sample_ids": list(group.input_sample_ids),
                    "technical_replicate_ids": list(group.technical_replicate_ids),
                    "n_technical_replicates": group.n_technical_replicates,
                }
                for group in policy.replicates.technical_replicate_groups
            ],
        },
        "empirical_bayes": {
            "method": policy.empirical_bayes.method,
            "robust": policy.empirical_bayes.robust,
            "trend": policy.empirical_bayes.trend,
            "winsor_tail_p": list(policy.empirical_bayes.winsor_tail_p),
        },
        "statistical_testing": {
            "test_statistic": policy.statistical_testing.test_statistic,
            "p_value_method": policy.statistical_testing.p_value_method,
            "adjusted_p_value_method": (
                policy.statistical_testing.adjusted_p_value_method
            ),
            "input_intensity_scale": (policy.statistical_testing.input_intensity_scale),
            "logfc_interpretation": (policy.statistical_testing.logfc_interpretation),
            "allow_suspicious_declared_input_scale": (
                policy.statistical_testing.allow_suspicious_declared_input_scale
            ),
        },
        "missing_values": {
            "policy": policy.missing_values.policy,
            "stage": policy.missing_values.stage,
        },
        "unsupported_design": {
            "intentionally_rejected_features": list(
                policy.unsupported_design.intentionally_rejected_features
            ),
            "enforcement_stage": policy.unsupported_design.enforcement_stage,
            "policy": policy.unsupported_design.policy,
        },
    }
    assert snapshot == {
        "design": {
            "formula": "~0 + condition",
            "description": "condition-only fixed-effect design",
            "sample_labels": ["A_1", "A_2", "B_1", "B_2"],
            "coefficient_labels": ["A", "B"],
            "condition_columns": ["A", "B"],
            "covariates": [],
            "sample_count": 4,
            "coefficient_count": 2,
            "rank": 2,
            "residual_degrees_of_freedom": 2.0,
            "paired_design_policy": "reject",
            "block_id_field_name": "block_id",
            "block_count": 0,
            "block_levels": [],
            "block_levels_included": [],
            "block_reference_level": None,
            "block_columns": [],
            "block_column_names": [],
            "condition_coverage_rule": (
                "block terms are not constructed under "
                "paired_design_policy='reject'; explicit block_id values are "
                "rejected before design-matrix construction"
            ),
            "limitations": [
                "paired_design_policy='reject' does not construct fixed-block terms",
                (
                    "explicit block_id metadata is rejected unless "
                    "paired_design_policy='fixed_block'"
                ),
                (
                    "unpaired condition and covariate workflows do not fit "
                    "duplicateCorrelation, mixed-effects, or random subject-effect "
                    "models"
                ),
            ],
            "rank_validation_status": "validated_full_rank",
            "estimability_validation_status": "validated_estimable",
        },
        "contrasts": [
            {
                "name": "B_vs_A",
                "numerator_condition": "B",
                "denominator_condition": "A",
                "coefficients": [("A", -1.0), ("B", 1.0)],
                "description": (
                    "condition contrast B - A; non-condition coefficients fixed at 0"
                ),
            }
        ],
        "replicates": {
            "minimum_condition_replicates": 2,
            "technical_replicate_policy": "reject",
            "condition_replicate_counts": [("A", 2), ("B", 2)],
            "technical_replicate_groups": [],
        },
        "empirical_bayes": {
            "method": "standard",
            "robust": False,
            "trend": False,
            "winsor_tail_p": [0.05, 0.1],
        },
        "statistical_testing": {
            "test_statistic": "moderated_t",
            "p_value_method": "two_sided_t_distribution_survival_function",
            "adjusted_p_value_method": "benjamini_hochberg",
            "input_intensity_scale": "log2",
            "logfc_interpretation": (
                "fitted condition contrast on the established log2 phosphosite "
                "intensity scale"
            ),
            "allow_suspicious_declared_input_scale": False,
        },
        "missing_values": {
            "policy": "reject_missing_values_before_differential_execution",
            "stage": "analysis_ready_dataset_boundary",
        },
        "unsupported_design": {
            "intentionally_rejected_features": [
                (
                    "correlated repeated-measure differential modelling beyond "
                    "explicit fixed blocks"
                ),
                "duplicateCorrelation-style correlated-replicate modelling",
                "mixed-effects differential modelling",
                "random subject-effect differential modelling",
            ],
            "enforcement_stage": (
                "validation.workflows.differential.ExperimentalDesignContractValidator"
            ),
            "policy": "reject_unsupported_design_features_before_execution",
        },
    }
