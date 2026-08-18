from __future__ import annotations

import pandas as pd
import pytest


def test_duplicate_correlation_scientific_contracts_are_not_public_exports() -> None:
    import phospy
    import phospy.advanced as advanced_api
    import phospy.api as public_api
    import phospy.api.configs as api_configs
    import phospy.contracts.configs as contract_configs
    import phospy.contracts.results as contract_results
    import phospy.science.differential.models as differential_models

    duplicate_correlation_symbols = {
        "DUPLICATE_CORRELATION_METHOD_REML_FISHER_TRIMMED_MEAN",
        "DUPLICATE_CORRELATION_TRIM_FRACTION",
        "DuplicateCorrelationBlockStructureSummary",
        "DuplicateCorrelationConsensusResult",
        "DuplicateCorrelationConsensusSummary",
        "DuplicateCorrelationFailureReason",
        "DuplicateCorrelationFeatureEstimate",
        "DuplicateCorrelationFeatureStatus",
        "DuplicateCorrelationWorkflowProvenance",
    }
    public_modules = (
        phospy,
        public_api,
        advanced_api,
        api_configs,
        contract_configs,
        contract_results,
        differential_models,
    )

    assert contract_configs.PAIRED_DESIGN_POLICY_DUPLICATE_CORRELATION == (
        "duplicate_correlation"
    )
    assert "duplicate_correlation" in contract_configs.SUPPORTED_PAIRED_DESIGN_POLICIES
    for module in public_modules:
        exported_names = set(getattr(module, "__all__", ()))
        assert duplicate_correlation_symbols.isdisjoint(exported_names)
        assert all(
            not hasattr(module, symbol_name)
            for symbol_name in duplicate_correlation_symbols
        )


def test_duplicate_correlation_execution_design_cannot_mix_fixed_block_columns() -> (
    None
):
    from phospy.contracts.configs import PAIRED_DESIGN_POLICY_DUPLICATE_CORRELATION
    from phospy.errors import WorkflowBoundaryError
    from phospy.science.differential.linear_model import decompose_differential_design
    from phospy.science.differential.models import ContrastMatrix, DesignMatrix
    from phospy.workflows.differential.models import (
        DifferentialBlockColumnMetadata,
        DifferentialConditionContrastVector,
        DifferentialExecutionDesignInputs,
    )

    sample_order = ("A_1", "B_1", "A_2", "B_2")
    design_frame = pd.DataFrame(
        {
            "A": [1.0, 0.0, 1.0, 0.0],
            "B": [0.0, 1.0, 0.0, 1.0],
        },
        index=pd.Index(sample_order, name="sample_id"),
    )
    contrast_frame = pd.DataFrame(
        {"B_vs_A": [-1.0, 1.0]},
        index=pd.Index(("A", "B"), name="coefficient"),
    )
    decomposition = decompose_differential_design(design_frame.to_numpy(dtype=float))

    valid_duplicate_design = DifferentialExecutionDesignInputs(
        design_matrix=DesignMatrix(design_frame),
        contrast_matrix=ContrastMatrix(contrast_frame),
        condition_contrast_vectors=(
            DifferentialConditionContrastVector(
                name="B_vs_A",
                numerator_condition="B",
                denominator_condition="A",
                coefficients=(("A", -1.0), ("B", 1.0)),
            ),
        ),
        covariate_columns=(),
        formula="~0 + condition",
        description="condition-only fixed-effect design",
        sample_order=sample_order,
        paired_design_policy=PAIRED_DESIGN_POLICY_DUPLICATE_CORRELATION,
        block_column_metadata=None,
        block_ids=("pair_1", "pair_1", "pair_2", "pair_2"),
        condition_labels=("A", "B"),
        coefficient_labels=("A", "B"),
        design_decomposition=decomposition,
    )

    assert valid_duplicate_design.block_column_metadata is None
    assert valid_duplicate_design.block_ids == (
        "pair_1",
        "pair_1",
        "pair_2",
        "pair_2",
    )

    with pytest.raises(
        WorkflowBoundaryError,
        match="do not combine fixed block coefficients",
    ):
        DifferentialExecutionDesignInputs(
            design_matrix=DesignMatrix(design_frame),
            contrast_matrix=ContrastMatrix(contrast_frame),
            condition_contrast_vectors=(
                DifferentialConditionContrastVector(
                    name="B_vs_A",
                    numerator_condition="B",
                    denominator_condition="A",
                    coefficients=(("A", -1.0), ("B", 1.0)),
                ),
            ),
            covariate_columns=(),
            formula="~0 + condition + block",
            description="invalid duplicate-correlation design",
            sample_order=sample_order,
            paired_design_policy=PAIRED_DESIGN_POLICY_DUPLICATE_CORRELATION,
            block_column_metadata=DifferentialBlockColumnMetadata(
                levels=("pair_1", "pair_2"),
                reference_level="pair_1",
                columns=(("pair_2", "block[pair_2]"),),
            ),
            block_ids=("pair_1", "pair_1", "pair_2", "pair_2"),
            condition_labels=("A", "B"),
            coefficient_labels=("A", "B"),
            design_decomposition=decomposition,
        )
