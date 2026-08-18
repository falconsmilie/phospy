from __future__ import annotations

import pandas as pd

from phospy import AnalysisReadyPhosphoDataset
from phospy.advanced import DifferentialAnalysisConfig
from phospy.api import (
    BatchCovariate,
    CategoricalCovariate,
    ContinuousCovariate,
    Contrast,
    DifferentialAnalysisRequest,
    DifferentialAnalysisWorkflow,
    ExperimentalDesign,
    Organism,
    SampleDesignRecord,
)
from phospy.contracts.configs import PAIRED_DESIGN_POLICY_FIXED_BLOCK
from phospy.contracts.results import DifferentialFixedEffectCovariateProvenance
from phospy.science.differential.models import DifferentialPolicyProvenance
from tests.support.analysis_ready_dataset_factories import (
    trusted_analysis_ready_dataset_from_tables,
)
from tests.support.intensity_scale_states import (
    supported_log2_intensity_scale_state,
    supported_log2_processing_state,
)
from tests.support.site_keys import protein_site_key_index, site_key_context_columns


def _dataset() -> AnalysisReadyPhosphoDataset:
    genes = ["MAPK14", "GSK3B", "AKT1"]
    sites = ["Y182", "S9", "T308"]
    site_index = protein_site_key_index(protein_identifiers=genes, sites=sites)
    phospho = pd.DataFrame(
        {
            "A_1": [1.0, 2.0, 1.0],
            "A_2": [1.1, 2.1, 1.1],
            "B_1": [2.1, 2.0, 1.0],
            "B_2": [2.0, 2.2, 0.9],
        },
        index=site_index,
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
    return trusted_analysis_ready_dataset_from_tables(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        intensity_scale_state=supported_log2_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_log2_processing_state(has_total_matrix=False),
    )


def _policy_for(
    design: ExperimentalDesign,
    *,
    config: DifferentialAnalysisConfig | None = None,
) -> DifferentialPolicyProvenance:
    result = DifferentialAnalysisWorkflow().run(
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
            config=DifferentialAnalysisConfig() if config is None else config,
        )
    )
    assert result.policy_provenance is not None
    return result.policy_provenance


def _condition_only_design() -> ExperimentalDesign:
    return ExperimentalDesign(
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


def test_condition_only_design_provenance_records_design_and_validation_status() -> (
    None
):
    policy = _policy_for(_condition_only_design())

    assert policy.design.formula == "~0 + condition"
    assert policy.design.description == "condition-only fixed-effect design"
    assert policy.design.condition_columns == ("A", "B")
    assert policy.design.covariates == ()
    assert policy.design.paired_design_policy == "reject"
    assert policy.design.block_id_field_name == "block_id"
    assert policy.design.block_count == 0
    assert policy.design.block_levels == ()
    assert policy.design.block_levels_included == ()
    assert policy.design.block_reference_level is None
    assert policy.design.block_columns == ()
    assert policy.design.block_column_names == ()
    assert policy.design.condition_coverage_rule == (
        "block terms are not constructed under paired_design_policy='reject'; "
        "explicit block_id values are rejected before design-matrix construction"
    )
    assert policy.design.limitations == (
        "paired_design_policy='reject' does not construct fixed-block terms",
        (
            "explicit block_id metadata is rejected unless "
            "paired_design_policy='fixed_block' or "
            "paired_design_policy='duplicate_correlation'"
        ),
        (
            "unpaired condition and covariate workflows do not fit "
            "duplicate_correlation, mixed-effects, or random subject-effect models"
        ),
    )
    assert policy.design.rank_validation_status == "validated_full_rank"
    assert policy.design.estimability_validation_status == "validated_estimable"
    assert policy.contrasts[0].coefficients == (("A", -1.0), ("B", 1.0))
    assert policy.contrasts[0].description == (
        "condition contrast B - A; non-condition coefficients fixed at 0"
    )
    assert policy.unsupported_design.policy == (
        "reject_unsupported_design_features_before_execution"
    )
    assert (
        "correlated repeated-measure differential modelling beyond explicit "
        "fixed_block and duplicate_correlation policies"
    ) in policy.unsupported_design.intentionally_rejected_features
    assert "duplicateCorrelation-style correlated-replicate modelling" not in (
        policy.unsupported_design.intentionally_rejected_features
    )
    assert "mixed-effects differential modelling" in (
        policy.unsupported_design.intentionally_rejected_features
    )


def test_categorical_covariate_provenance_records_columns_and_kind() -> None:
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

    policy = _policy_for(design)
    covariate = policy.design.covariates[0]

    assert isinstance(covariate, DifferentialFixedEffectCovariateProvenance)
    assert policy.design.formula == "~0 + condition + sex"
    assert policy.design.condition_columns == ("A", "B")
    assert covariate.name == "sex"
    assert covariate.kind == "categorical"
    assert covariate.columns == ("sex[M]",)
    assert covariate.levels == ("F", "M")
    assert covariate.reference_level == "F"
    assert policy.contrasts[0].coefficients == (
        ("A", -1.0),
        ("B", 1.0),
        ("sex[M]", 0.0),
    )


def test_continuous_covariate_provenance_records_columns_and_kind() -> None:
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(
                sample_id="A_1",
                condition="A",
                biological_replicate_id="A_r1",
                covariates={"dose": 0.0},
            ),
            SampleDesignRecord(
                sample_id="A_2",
                condition="A",
                biological_replicate_id="A_r2",
                covariates={"dose": 1.0},
            ),
            SampleDesignRecord(
                sample_id="B_1",
                condition="B",
                biological_replicate_id="B_r1",
                covariates={"dose": 0.0},
            ),
            SampleDesignRecord(
                sample_id="B_2",
                condition="B",
                biological_replicate_id="B_r2",
                covariates={"dose": 1.0},
            ),
        ),
        fixed_effects=(ContinuousCovariate("dose"),),
    )

    policy = _policy_for(design)
    covariate = policy.design.covariates[0]

    assert policy.design.formula == "~0 + condition + dose"
    assert covariate.name == "dose"
    assert covariate.kind == "continuous"
    assert covariate.columns == ("dose",)
    assert covariate.levels == ()
    assert covariate.reference_level is None
    assert policy.contrasts[0].coefficients == (
        ("A", -1.0),
        ("B", 1.0),
        ("dose", 0.0),
    )


def test_batch_fixed_effect_provenance_records_batch_as_model_covariate() -> None:
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(
                sample_id="A_1",
                condition="A",
                biological_replicate_id="A_r1",
                batch="batch_1",
            ),
            SampleDesignRecord(
                sample_id="A_2",
                condition="A",
                biological_replicate_id="A_r2",
                batch="batch_2",
            ),
            SampleDesignRecord(
                sample_id="B_1",
                condition="B",
                biological_replicate_id="B_r1",
                batch="batch_1",
            ),
            SampleDesignRecord(
                sample_id="B_2",
                condition="B",
                biological_replicate_id="B_r2",
                batch="batch_2",
            ),
        ),
        fixed_effects=(BatchCovariate(),),
    )

    policy = _policy_for(design)
    covariate = policy.design.covariates[0]

    assert policy.design.formula == "~0 + condition + batch"
    assert policy.design.description == "fixed-effect design: ~0 + condition + batch"
    assert covariate.name == "batch"
    assert covariate.kind == "batch"
    assert covariate.columns == ("batch[batch_2]",)
    assert covariate.levels == ("batch_1", "batch_2")
    assert covariate.reference_level == "batch_1"


def test_fixed_block_provenance_records_block_columns_and_policy() -> None:
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(
                sample_id="A_1",
                condition="A",
                biological_replicate_id="A_r1",
                block_id="pair_1",
            ),
            SampleDesignRecord(
                sample_id="A_2",
                condition="A",
                biological_replicate_id="A_r2",
                block_id="pair_2",
            ),
            SampleDesignRecord(
                sample_id="B_1",
                condition="B",
                biological_replicate_id="B_r1",
                block_id="pair_1",
            ),
            SampleDesignRecord(
                sample_id="B_2",
                condition="B",
                biological_replicate_id="B_r2",
                block_id="pair_2",
            ),
        )
    )

    policy = _policy_for(
        design,
        config=DifferentialAnalysisConfig(
            paired_design_policy=PAIRED_DESIGN_POLICY_FIXED_BLOCK,
        ),
    )

    assert policy.design.formula == "~0 + condition + block"
    assert policy.design.description == "fixed-effect design: ~0 + condition + block"
    assert policy.design.paired_design_policy == PAIRED_DESIGN_POLICY_FIXED_BLOCK
    assert policy.design.block_id_field_name == "block_id"
    assert policy.design.block_count == 2
    assert policy.design.block_levels == ("pair_1", "pair_2")
    assert policy.design.block_levels_included == ("pair_1", "pair_2")
    assert policy.design.block_reference_level == "pair_1"
    assert policy.design.block_columns == (("pair_2", "block[pair_2]"),)
    assert policy.design.block_column_names == ("block[pair_2]",)
    assert policy.design.condition_coverage_rule == (
        "for every requested condition contrast, every block must contain both "
        "numerator and denominator conditions; incomplete or partially covered "
        "blocks are rejected before execution"
    )
    assert policy.design.limitations == (
        "fixed_block adds block_id levels as ordinary fixed-effect design columns",
        (
            "fixed_block does not estimate within-block correlation and is not "
            "limma duplicateCorrelation"
        ),
        "fixed_block does not fit mixed-effects or random subject-effect models",
        (
            "incomplete blocks are rejected before execution; samples are not "
            "silently dropped"
        ),
    )
    assert policy.contrasts[0].coefficients == (
        ("A", -1.0),
        ("B", 1.0),
        ("block[pair_2]", 0.0),
    )
