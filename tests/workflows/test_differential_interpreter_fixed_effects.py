from __future__ import annotations

import numpy as np
import pandas as pd
import pandas.testing as pdt

from phospy import AnalysisReadyPhosphoDataset
from phospy.advanced import DifferentialAnalysisConfig
from phospy.api import (
    CategoricalCovariate,
    ContinuousCovariate,
    Contrast,
    DifferentialAnalysisRequest,
    DifferentialAnalysisWorkflow,
    ExperimentalDesign,
    Organism,
    SampleDesignRecord,
)
from phospy.contracts.configs import (
    PAIRED_DESIGN_POLICY_FIXED_BLOCK,
    PAIRED_DESIGN_POLICY_REJECT,
)
from phospy.workflows.differential.models import (
    InterpretedDifferentialAnalysisRequest,
)
from tests.support.analysis_ready_dataset_factories import (
    trusted_analysis_ready_dataset_from_tables,
)
from tests.support.intensity_scale_states import (
    supported_log2_intensity_scale_state,
    supported_log2_processing_state,
)
from tests.support.site_keys import protein_site_key_index, site_key_context_columns


def _dataset_from_phospho(
    *,
    phospho: pd.DataFrame,
    genes: tuple[str, ...],
    sites: tuple[str, ...],
) -> AnalysisReadyPhosphoDataset:
    site_metadata = pd.DataFrame(
        {
            "site_key": phospho.index.astype(str).tolist(),
            "display_id": [
                f"{gene};{site};" for gene, site in zip(genes, sites, strict=True)
            ],
            **site_key_context_columns(phospho.index),
            "gene_symbol": genes,
            "site": sites,
            "site_sequence": [
                ("A" * 15) + str(site).strip().upper()[0] + ("A" * 15) for site in sites
            ],
            "protein_id": genes,
        },
        index=phospho.index.copy(),
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


def _dataset() -> AnalysisReadyPhosphoDataset:
    genes = ("MAPK14", "GSK3B")
    sites = ("Y182", "S9")
    site_index = protein_site_key_index(protein_identifiers=genes, sites=sites)
    phospho = pd.DataFrame(
        {
            "A_1": [1.0, 2.0],
            "A_2": [1.1, 2.1],
            "B_1": [2.0, 1.9],
            "B_2": [2.2, 2.2],
        },
        index=site_index,
    )
    return _dataset_from_phospho(phospho=phospho, genes=genes, sites=sites)


def _continuous_adjustment_dataset() -> AnalysisReadyPhosphoDataset:
    genes = ["MAPK14", "GSK3B", "AKT1"]
    sites = ["Y182", "S9", "T308"]
    site_index = protein_site_key_index(protein_identifiers=genes, sites=sites)
    samples = ("A_1", "A_2", "A_3", "B_1", "B_2", "B_3")
    dose = (0.0, 1.0, 2.0, 3.0, 4.0, 5.0)
    noise = (0.10, -0.12, 0.06, -0.08, 0.09, -0.04)
    values = {
        "site1": [10.0 + 2.0 * dose[idx] + noise[idx] for idx in range(6)],
        "site2": [
            5.0
            + (1.0 if sample.startswith("B") else 0.0)
            + 0.5 * dose[idx]
            + tuple(reversed(noise))[idx]
            for idx, sample in enumerate(samples)
        ],
        "site3": [1.0 + 0.2 * float(idx) + noise[idx] * 0.3 for idx in range(6)],
    }
    phospho = pd.DataFrame(
        {
            sample: [feature_values[idx] for feature_values in values.values()]
            for idx, sample in enumerate(samples)
        },
        index=site_index,
    )
    site_metadata = pd.DataFrame(
        {
            "site_key": site_index.tolist(),
            "display_id": [
                "MAPK14;Y182;",
                "GSK3B;S9;",
                "AKT1;T308;",
            ],
            **site_key_context_columns(site_index),
            "gene_symbol": genes,
            "site": sites,
            "site_sequence": [
                ("A" * 15) + str(site).strip().upper()[0] + ("A" * 15) for site in sites
            ],
            "protein_id": genes,
        },
        index=phospho.index.copy(),
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


def _block_effect_dataset() -> AnalysisReadyPhosphoDataset:
    genes = ("MAPK14", "GSK3B", "AKT1")
    sites = ("Y182", "S9", "T308")
    site_index = protein_site_key_index(protein_identifiers=genes, sites=sites)
    phospho = pd.DataFrame(
        {
            "A_pair_1_rep_1": [10.0, 5.0, 1.0],
            "A_pair_1_rep_2": [10.4, 5.2, 1.2],
            "B_pair_1": [12.2, 5.1, 1.1],
            "A_pair_2": [0.2, -4.8, 0.8],
            "B_pair_2_rep_1": [2.0, -4.9, 1.0],
            "B_pair_2_rep_2": [2.4, -5.1, 0.9],
        },
        index=site_index,
    )
    return _dataset_from_phospho(phospho=phospho, genes=genes, sites=sites)


def _block_effect_design(*, include_blocks: bool) -> ExperimentalDesign:
    records = (
        ("A_pair_1_rep_1", "A", "pair_1"),
        ("A_pair_1_rep_2", "A", "pair_1"),
        ("B_pair_1", "B", "pair_1"),
        ("A_pair_2", "A", "pair_2"),
        ("B_pair_2_rep_1", "B", "pair_2"),
        ("B_pair_2_rep_2", "B", "pair_2"),
    )
    return ExperimentalDesign(
        samples=tuple(
            SampleDesignRecord(
                sample_id=sample_id,
                condition=condition,
                biological_replicate_id=f"{condition}_{idx}",
                block_id=block_id if include_blocks else None,
            )
            for idx, (sample_id, condition, block_id) in enumerate(records, start=1)
        )
    )


def _b_vs_a_contrast() -> tuple[Contrast, ...]:
    return (
        Contrast(
            name="B_vs_A",
            numerator_condition="B",
            denominator_condition="A",
        ),
    )


def test_differential_interpreter_passes_fixed_effect_inputs_to_executor() -> None:
    sentinel = object()
    captured: dict[str, InterpretedDifferentialAnalysisRequest] = {}

    class _ExecutorSpy:
        def run(self, request: InterpretedDifferentialAnalysisRequest) -> object:
            captured["request"] = request
            return sentinel

    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(
                sample_id="A_1",
                condition="A",
                biological_replicate_id="A_r1",
                covariates={"sex": "F"},
            ),
            SampleDesignRecord(
                sample_id="A_2",
                condition="A",
                biological_replicate_id="A_r2",
                covariates={"sex": "M"},
            ),
            SampleDesignRecord(
                sample_id="B_1",
                condition="B",
                biological_replicate_id="B_r1",
                covariates={"sex": "F"},
            ),
            SampleDesignRecord(
                sample_id="B_2",
                condition="B",
                biological_replicate_id="B_r2",
                covariates={"sex": "M"},
            ),
        ),
        fixed_effects=(CategoricalCovariate("sex"),),
    )
    result = DifferentialAnalysisWorkflow._with_components(
        executor=_ExecutorSpy(),  # type: ignore[arg-type]
    ).run(
        DifferentialAnalysisRequest(
            dataset=_dataset(),
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

    assert result is sentinel
    interpreted = captured["request"]
    execution_design = interpreted.execution_design
    assert execution_design is not None
    assert execution_design.formula == "~0 + condition + sex"
    assert execution_design.covariate_columns[0].columns == ("sex[M]",)
    assert interpreted.computation_request.design.to_dataframe().columns.tolist() == [
        "A",
        "B",
        "sex[M]",
    ]


def test_differential_block_fixed_block_inputs_are_passed_to_executor() -> None:
    sentinel = object()
    captured: dict[str, InterpretedDifferentialAnalysisRequest] = {}

    class _ExecutorSpy:
        def run(self, request: InterpretedDifferentialAnalysisRequest) -> object:
            captured["request"] = request
            return sentinel

    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(sample_id="A_1", condition="A", block_id="pair_1"),
            SampleDesignRecord(sample_id="A_2", condition="A", block_id="pair_2"),
            SampleDesignRecord(sample_id="B_1", condition="B", block_id="pair_1"),
            SampleDesignRecord(sample_id="B_2", condition="B", block_id="pair_2"),
        )
    )

    result = DifferentialAnalysisWorkflow._with_components(
        executor=_ExecutorSpy(),  # type: ignore[arg-type]
    ).run(
        DifferentialAnalysisRequest(
            dataset=_dataset(),
            design=design,
            contrasts=(
                Contrast(
                    name="B_vs_A",
                    numerator_condition="B",
                    denominator_condition="A",
                ),
            ),
            config=DifferentialAnalysisConfig(
                paired_design_policy=PAIRED_DESIGN_POLICY_FIXED_BLOCK,
            ),
        )
    )

    assert result is sentinel
    interpreted = captured["request"]
    execution_design = interpreted.execution_design
    assert execution_design is not None
    assert execution_design.formula == "~0 + condition + block"
    assert execution_design.description == "fixed-effect design: ~0 + condition + block"
    assert execution_design.paired_design_policy == PAIRED_DESIGN_POLICY_FIXED_BLOCK
    assert execution_design.block_column_metadata is not None
    assert execution_design.block_column_metadata.levels == ("pair_1", "pair_2")
    assert execution_design.block_column_metadata.reference_level == "pair_1"
    assert execution_design.block_column_metadata.columns == (
        ("pair_2", "block[pair_2]"),
    )
    assert interpreted.computation_request.design.to_dataframe().columns.tolist() == [
        "A",
        "B",
        "block[pair_2]",
    ]


def test_differential_block_valid_paired_two_condition_design_executes() -> None:
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(sample_id="A_1", condition="A", block_id="pair_1"),
            SampleDesignRecord(sample_id="A_2", condition="A", block_id="pair_2"),
            SampleDesignRecord(sample_id="B_1", condition="B", block_id="pair_1"),
            SampleDesignRecord(sample_id="B_2", condition="B", block_id="pair_2"),
        )
    )

    result = DifferentialAnalysisWorkflow().run(
        DifferentialAnalysisRequest(
            dataset=_dataset(),
            design=design,
            contrasts=_b_vs_a_contrast(),
            config=DifferentialAnalysisConfig(
                paired_design_policy=PAIRED_DESIGN_POLICY_FIXED_BLOCK,
            ),
        )
    )

    table = result.table_for("B_vs_A")
    assert np.isfinite(table.loc[:, "logFC"]).all()
    assert np.isfinite(table.loc[:, "t"]).all()
    assert result.policy_provenance is not None
    assert result.policy_provenance.design.formula == "~0 + condition + block"
    assert result.policy_provenance.design.paired_design_policy == (
        PAIRED_DESIGN_POLICY_FIXED_BLOCK
    )
    assert result.policy_provenance.design.block_levels == ("pair_1", "pair_2")
    assert result.policy_provenance.design.block_reference_level == "pair_1"
    assert result.policy_provenance.design.block_columns == (
        ("pair_2", "block[pair_2]"),
    )


def test_differential_fixed_block_provenance_records_block_count_and_columns() -> None:
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(sample_id="A_1", condition="A", block_id="pair_1"),
            SampleDesignRecord(sample_id="A_2", condition="A", block_id="pair_2"),
            SampleDesignRecord(sample_id="B_1", condition="B", block_id="pair_1"),
            SampleDesignRecord(sample_id="B_2", condition="B", block_id="pair_2"),
        )
    )

    result = DifferentialAnalysisWorkflow().run(
        DifferentialAnalysisRequest(
            dataset=_dataset(),
            design=design,
            contrasts=_b_vs_a_contrast(),
            config=DifferentialAnalysisConfig(
                paired_design_policy=PAIRED_DESIGN_POLICY_FIXED_BLOCK,
            ),
        )
    )

    assert result.policy_provenance is not None
    design_provenance = result.policy_provenance.design
    assert design_provenance.block_id_field_name == "block_id"
    assert design_provenance.block_count == 2
    assert design_provenance.block_levels_included == ("pair_1", "pair_2")
    assert design_provenance.block_column_names == ("block[pair_2]",)
    assert "incomplete or partially covered blocks are rejected" in (
        design_provenance.condition_coverage_rule
    )
    assert "limma duplicateCorrelation" in design_provenance.limitations[1]


def test_differential_block_adjusted_contrast_differs_when_block_effect_exists() -> (
    None
):
    dataset = _block_effect_dataset()
    unblocked = DifferentialAnalysisWorkflow().run(
        DifferentialAnalysisRequest(
            dataset=dataset,
            design=_block_effect_design(include_blocks=False),
            contrasts=_b_vs_a_contrast(),
        )
    )
    adjusted = DifferentialAnalysisWorkflow().run(
        DifferentialAnalysisRequest(
            dataset=dataset,
            design=_block_effect_design(include_blocks=True),
            contrasts=_b_vs_a_contrast(),
            config=DifferentialAnalysisConfig(
                paired_design_policy=PAIRED_DESIGN_POLICY_FIXED_BLOCK,
            ),
        )
    )

    site_key = str(dataset.phospho.index[0])
    unblocked_log_fc = float(unblocked.table_for("B_vs_A").at[site_key, "logFC"])
    adjusted_log_fc = float(adjusted.table_for("B_vs_A").at[site_key, "logFC"])

    assert unblocked_log_fc < 0.0
    assert adjusted_log_fc > 1.5
    assert adjusted_log_fc != unblocked_log_fc


def test_differential_block_empirical_bayes_still_runs_with_fixed_block_design() -> (
    None
):
    result = DifferentialAnalysisWorkflow().run(
        DifferentialAnalysisRequest(
            dataset=_block_effect_dataset(),
            design=_block_effect_design(include_blocks=True),
            contrasts=_b_vs_a_contrast(),
            config=DifferentialAnalysisConfig(
                paired_design_policy=PAIRED_DESIGN_POLICY_FIXED_BLOCK,
            ),
        )
    )

    assert result.empirical_bayes_method == "standard"
    assert result.prior_diagnostics.method == "standard"
    assert np.isfinite(result.posterior_residual_variance_series()).all()
    assert result.residual_degrees_of_freedom == 3.0


def test_differential_block_reject_policy_condition_only_output_is_unchanged() -> None:
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(sample_id="A_1", condition="A"),
            SampleDesignRecord(sample_id="A_2", condition="A"),
            SampleDesignRecord(sample_id="B_1", condition="B"),
            SampleDesignRecord(sample_id="B_2", condition="B"),
        )
    )
    default_result = DifferentialAnalysisWorkflow().run(
        DifferentialAnalysisRequest(
            dataset=_dataset(),
            design=design,
            contrasts=_b_vs_a_contrast(),
        )
    )
    explicit_reject_result = DifferentialAnalysisWorkflow().run(
        DifferentialAnalysisRequest(
            dataset=_dataset(),
            design=design,
            contrasts=_b_vs_a_contrast(),
            config=DifferentialAnalysisConfig(
                paired_design_policy=PAIRED_DESIGN_POLICY_REJECT,
            ),
        )
    )

    pdt.assert_frame_equal(
        default_result.table_for("B_vs_A"),
        explicit_reject_result.table_for("B_vs_A"),
    )


def test_workflow_executes_continuous_covariate_adjusted_contrast() -> None:
    dataset = _continuous_adjustment_dataset()
    samples = ("A_1", "A_2", "A_3", "B_1", "B_2", "B_3")
    conditions = ("A", "A", "A", "B", "B", "B")
    doses = (0.0, 1.0, 2.0, 3.0, 4.0, 5.0)
    unadjusted_design = ExperimentalDesign(
        samples=tuple(
            SampleDesignRecord(
                sample_id=sample,
                condition=condition,
                biological_replicate_id=f"{condition}_{idx}",
            )
            for idx, (sample, condition) in enumerate(
                zip(samples, conditions, strict=True),
                start=1,
            )
        )
    )
    adjusted_design = ExperimentalDesign(
        samples=tuple(
            SampleDesignRecord(
                sample_id=sample,
                condition=condition,
                biological_replicate_id=f"{condition}_{idx}",
                covariates={"dose": dose},
            )
            for idx, (sample, condition, dose) in enumerate(
                zip(samples, conditions, doses, strict=True),
                start=1,
            )
        ),
        fixed_effects=(ContinuousCovariate("dose"),),
    )
    contrasts = (
        Contrast(
            name="B_vs_A",
            numerator_condition="B",
            denominator_condition="A",
        ),
    )

    unadjusted = DifferentialAnalysisWorkflow().run(
        DifferentialAnalysisRequest(
            dataset=dataset,
            design=unadjusted_design,
            contrasts=contrasts,
        )
    )
    adjusted = DifferentialAnalysisWorkflow().run(
        DifferentialAnalysisRequest(
            dataset=dataset,
            design=adjusted_design,
            contrasts=contrasts,
        )
    )

    site_key = str(dataset.phospho.index[0])
    assert abs(float(unadjusted.table_for("B_vs_A").at[site_key, "logFC"])) > 5.0
    assert abs(float(adjusted.table_for("B_vs_A").at[site_key, "logFC"])) < 0.1
    assert adjusted.policy_provenance is not None
    assert adjusted.policy_provenance.design.formula == "~0 + condition + dose"
    assert adjusted.policy_provenance.design.coefficient_labels == ("A", "B", "dose")
