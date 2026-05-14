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
    phospho = pd.DataFrame(
        {
            "A_1": [1.0, 2.0],
            "A_2": [1.1, 2.1],
            "B_1": [2.0, 1.9],
            "B_2": [2.2, 2.2],
        },
        index=pd.Index(["MAPK14;Y182;", "GSK3B;S9;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "GSK3B"],
            "site": ["Y182", "S9"],
            "site_sequence": [
                ("A" * 15) + "Y" + ("A" * 15),
                ("A" * 15) + "S" + ("A" * 15),
            ],
            "protein_id": ["MAPK14", "GSK3B"],
        },
        index=phospho.index.copy(),
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
            "sample_labels": list(policy.design.sample_labels),
            "coefficient_labels": list(policy.design.coefficient_labels),
            "sample_count": policy.design.sample_count,
            "coefficient_count": policy.design.coefficient_count,
            "rank": policy.design.rank,
            "residual_degrees_of_freedom": policy.design.residual_degrees_of_freedom,
        },
        "contrasts": [
            {
                "name": contrast.name,
                "numerator_condition": contrast.numerator_condition,
                "denominator_condition": contrast.denominator_condition,
                "coefficients": list(contrast.coefficients),
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
        },
    }
    assert snapshot == {
        "design": {
            "formula": "~0 + condition",
            "sample_labels": ["A_1", "A_2", "B_1", "B_2"],
            "coefficient_labels": ["A", "B"],
            "sample_count": 4,
            "coefficient_count": 2,
            "rank": 2,
            "residual_degrees_of_freedom": 2.0,
        },
        "contrasts": [
            {
                "name": "B_vs_A",
                "numerator_condition": "B",
                "denominator_condition": "A",
                "coefficients": [("A", -1.0), ("B", 1.0)],
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
        },
        "missing_values": {
            "policy": "reject_missing_values_before_differential_execution",
            "stage": "analysis_ready_dataset_boundary",
        },
        "unsupported_design": {
            "intentionally_rejected_features": [
                "batch-aware differential modelling",
                "blocking/paired/repeated-measure differential modelling",
            ],
            "enforcement_stage": (
                "validation.workflows.differential.ExperimentalDesignContractValidator"
            ),
        },
    }
