from __future__ import annotations

import pandas as pd
import pytest

from phospy import AnalysisReadyPhosphoDataset
from phospy.api import (
    PAIRED_DESIGN_POLICY_FIXED_BLOCK,
    BatchCovariate,
    CategoricalCovariate,
    ContinuousCovariate,
    Contrast,
    ExperimentalDesign,
    FixedEffectCovariate,
    Organism,
    SampleDesignRecord,
)
from phospy.errors import WorkflowValidationError
from phospy.validation.workflows.differential import (
    ExperimentalDesignContractValidator,
)
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
)
from tests.support.site_keys import protein_site_key_index, site_key_context_columns


def _dataset(
    *,
    sample_order: tuple[str, ...] = ("A_1", "A_2", "B_1", "B_2"),
) -> AnalysisReadyPhosphoDataset:
    gene_symbols = ["MAPK14", "GSK3B"]
    sites = ["Y182", "S9"]
    display_ids = ["MAPK14;Y182;", "GSK3B;S9;"]
    site_index = protein_site_key_index(
        protein_identifiers=gene_symbols,
        sites=sites,
    )
    phospho = pd.DataFrame(
        {
            "A_1": [1.0, 2.0],
            "A_2": [1.1, 2.2],
            "B_1": [2.0, 1.9],
            "B_2": [2.2, 2.1],
        },
        index=site_index,
    ).loc[:, list(sample_order)]
    site_metadata = pd.DataFrame(
        {
            "site_key": site_index.astype(str).tolist(),
            "display_id": display_ids,
            **site_key_context_columns(site_index),
            "gene_symbol": gene_symbols,
            "site": sites,
            "site_sequence": [
                ("A" * 15) + str(site).strip().upper()[0] + ("A" * 15) for site in sites
            ],
            "protein_id": gene_symbols,
        },
        index=phospho.index.copy(),
    )
    return AnalysisReadyPhosphoDataset(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
    )


def _design() -> ExperimentalDesign:
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


def _contrasts() -> tuple[Contrast, ...]:
    return (
        Contrast(
            name="B_vs_A",
            numerator_condition="B",
            denominator_condition="A",
        ),
    )


def _paired_block_design() -> ExperimentalDesign:
    return ExperimentalDesign(
        samples=(
            SampleDesignRecord(sample_id="A_1", condition="A", block_id="pair_1"),
            SampleDesignRecord(sample_id="A_2", condition="A", block_id="pair_2"),
            SampleDesignRecord(sample_id="B_1", condition="B", block_id="pair_1"),
            SampleDesignRecord(sample_id="B_2", condition="B", block_id="pair_2"),
        )
    )


def test_valid_simple_two_condition_design() -> None:
    validated = ExperimentalDesignContractValidator().run(
        dataset=_dataset(),
        design=_design(),
        contrasts=_contrasts(),
        allow_design_subset=False,
        minimum_condition_replicates=2,
    )
    assert validated.analysis_sample_ids == ("A_1", "A_2", "B_1", "B_2")
    assert validated.design_frame.columns.tolist() == ["A", "B"]
    assert validated.contrast_frame.columns.tolist() == ["B_vs_A"]


def test_missing_sample_in_design_fails() -> None:
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(sample_id="A_1", condition="A"),
            SampleDesignRecord(sample_id="A_2", condition="A"),
            SampleDesignRecord(sample_id="B_1", condition="B"),
        )
    )
    with pytest.raises(WorkflowValidationError, match="missing required dataset"):
        ExperimentalDesignContractValidator().run(
            dataset=_dataset(),
            design=design,
            contrasts=_contrasts(),
            allow_design_subset=False,
            minimum_condition_replicates=1,
        )


def test_extra_sample_in_design_fails() -> None:
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(sample_id="A_1", condition="A"),
            SampleDesignRecord(sample_id="A_2", condition="A"),
            SampleDesignRecord(sample_id="B_1", condition="B"),
            SampleDesignRecord(sample_id="B_2", condition="B"),
            SampleDesignRecord(sample_id="C_1", condition="C"),
        )
    )
    with pytest.raises(WorkflowValidationError, match="not present in dataset"):
        ExperimentalDesignContractValidator().run(
            dataset=_dataset(),
            design=design,
            contrasts=_contrasts(),
            allow_design_subset=False,
            minimum_condition_replicates=1,
        )


def test_duplicate_sample_ids_rejected() -> None:
    with pytest.raises(WorkflowValidationError, match="duplicate sample IDs"):
        ExperimentalDesign(
            samples=(
                SampleDesignRecord(sample_id="A_1", condition="A"),
                SampleDesignRecord(sample_id="A_1", condition="A"),
                SampleDesignRecord(sample_id="B_1", condition="B"),
                SampleDesignRecord(sample_id="B_2", condition="B"),
            )
        )


def test_empty_condition_labels_rejected() -> None:
    with pytest.raises(WorkflowValidationError, match="condition"):
        SampleDesignRecord(sample_id="A_1", condition="")


def test_contrast_unknown_condition_fails() -> None:
    with pytest.raises(WorkflowValidationError, match="unknown denominator condition"):
        ExperimentalDesignContractValidator().run(
            dataset=_dataset(),
            design=_design(),
            contrasts=(
                Contrast(
                    name="B_vs_A",
                    numerator_condition="B",
                    denominator_condition="A_missing",
                ),
            ),
            allow_design_subset=False,
            minimum_condition_replicates=2,
        )


def test_insufficient_replicates_fails() -> None:
    one_rep_per_condition = ExperimentalDesign(
        samples=(
            SampleDesignRecord(sample_id="A_1", condition="A"),
            SampleDesignRecord(sample_id="B_1", condition="B"),
        )
    )
    with pytest.raises(WorkflowValidationError, match="insufficient replicate counts"):
        ExperimentalDesignContractValidator().run(
            dataset=_dataset(sample_order=("A_1", "B_1")),
            design=one_rep_per_condition,
            contrasts=_contrasts(),
            allow_design_subset=False,
            minimum_condition_replicates=2,
        )


def test_optional_batch_field_validation_fails_when_partially_defined() -> None:
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(sample_id="A_1", condition="A", batch="batch_1"),
            SampleDesignRecord(sample_id="A_2", condition="A"),
            SampleDesignRecord(sample_id="B_1", condition="B", batch="batch_2"),
            SampleDesignRecord(sample_id="B_2", condition="B", batch="batch_2"),
        )
    )
    with pytest.raises(WorkflowValidationError, match="optional field 'batch'"):
        ExperimentalDesignContractValidator().run(
            dataset=_dataset(),
            design=design,
            contrasts=_contrasts(),
            allow_design_subset=False,
            minimum_condition_replicates=1,
        )


def test_differential_block_fixed_block_policy_rejects_partial_block_identifiers() -> (
    None
):
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(sample_id="A_1", condition="A", block_id="pair_1"),
            SampleDesignRecord(sample_id="A_2", condition="A"),
            SampleDesignRecord(sample_id="B_1", condition="B", block_id="pair_1"),
            SampleDesignRecord(sample_id="B_2", condition="B", block_id="pair_2"),
        )
    )
    with pytest.raises(
        WorkflowValidationError,
        match="fixed_block.*requires block_id.*A_2",
    ):
        ExperimentalDesignContractValidator().run(
            dataset=_dataset(),
            design=design,
            contrasts=_contrasts(),
            allow_design_subset=False,
            minimum_condition_replicates=1,
            paired_design_policy=PAIRED_DESIGN_POLICY_FIXED_BLOCK,
        )


def test_differential_block_metadata_rejected_under_default_policy() -> None:
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(sample_id="A_1", condition="A", block_id="pair_1"),
            SampleDesignRecord(sample_id="A_2", condition="A", block_id="pair_2"),
            SampleDesignRecord(sample_id="B_1", condition="B", block_id="pair_1"),
            SampleDesignRecord(sample_id="B_2", condition="B", block_id="pair_2"),
        )
    )

    with pytest.raises(
        WorkflowValidationError,
        match="paired_design_policy='fixed_block'",
    ):
        ExperimentalDesignContractValidator().run(
            dataset=_dataset(),
            design=design,
            contrasts=_contrasts(),
            allow_design_subset=False,
            minimum_condition_replicates=1,
        )


def test_differential_block_fixed_block_policy_requires_block_identifiers() -> None:
    with pytest.raises(
        WorkflowValidationError,
        match="fixed_block.*requires block_id.*A_1",
    ):
        ExperimentalDesignContractValidator().run(
            dataset=_dataset(),
            design=_design(),
            contrasts=_contrasts(),
            allow_design_subset=False,
            minimum_condition_replicates=1,
            paired_design_policy=PAIRED_DESIGN_POLICY_FIXED_BLOCK,
        )


def test_differential_block_fixed_block_valid_paired_two_condition_design() -> None:
    validated = ExperimentalDesignContractValidator().run(
        dataset=_dataset(),
        design=_paired_block_design(),
        contrasts=_contrasts(),
        allow_design_subset=False,
        minimum_condition_replicates=2,
        paired_design_policy=PAIRED_DESIGN_POLICY_FIXED_BLOCK,
    )

    assert validated.analysis_sample_ids == ("A_1", "A_2", "B_1", "B_2")
    assert validated.design_frame.columns.tolist() == ["A", "B", "block[pair_2]"]
    assert validated.contrast_frame.index.tolist() == ["A", "B", "block[pair_2]"]
    assert validated.contrast_frame.loc[:, "B_vs_A"].tolist() == [-1.0, 1.0, 0.0]
    assert validated.design_build_result is not None
    assert validated.design_build_result.block_levels == ("pair_1", "pair_2")


def test_differential_block_fixed_block_policy_accepts_metadata_for_later_validation() -> (
    None
):
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(sample_id="A_1", condition="A", block_id="pair_1"),
            SampleDesignRecord(sample_id="A_2", condition="A", block_id="pair_2"),
            SampleDesignRecord(sample_id="B_1", condition="B", block_id="pair_1"),
            SampleDesignRecord(sample_id="B_2", condition="B", block_id="pair_2"),
        )
    )

    with pytest.raises(
        WorkflowValidationError,
        match="unknown denominator condition",
    ):
        ExperimentalDesignContractValidator().run(
            dataset=_dataset(),
            design=design,
            contrasts=(
                Contrast(
                    name="B_vs_missing",
                    numerator_condition="B",
                    denominator_condition="missing",
                ),
            ),
            allow_design_subset=False,
            minimum_condition_replicates=1,
            paired_design_policy=PAIRED_DESIGN_POLICY_FIXED_BLOCK,
        )


def test_differential_block_fixed_block_rejects_incomplete_pair() -> None:
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(sample_id="A_1", condition="A", block_id="pair_1"),
            SampleDesignRecord(sample_id="A_2", condition="A", block_id="pair_2"),
            SampleDesignRecord(sample_id="B_1", condition="B", block_id="pair_1"),
            SampleDesignRecord(sample_id="B_2", condition="B", block_id="pair_3"),
        )
    )

    with pytest.raises(
        WorkflowValidationError,
        match="fixed_block.*at least 2 samples.*incomplete blocks.*pair_2",
    ):
        ExperimentalDesignContractValidator().run(
            dataset=_dataset(),
            design=design,
            contrasts=_contrasts(),
            allow_design_subset=False,
            minimum_condition_replicates=1,
            paired_design_policy=PAIRED_DESIGN_POLICY_FIXED_BLOCK,
        )


def test_differential_block_fixed_block_rejects_block_with_one_condition_only() -> None:
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(sample_id="A_1", condition="A", block_id="pair_1"),
            SampleDesignRecord(sample_id="A_2", condition="A", block_id="pair_2"),
            SampleDesignRecord(sample_id="B_1", condition="B", block_id="pair_1"),
            SampleDesignRecord(sample_id="B_2", condition="A", block_id="pair_2"),
        )
    )

    with pytest.raises(
        WorkflowValidationError,
        match="invalid block condition coverage.*pair_2 missing numerator condition 'B'",
    ):
        ExperimentalDesignContractValidator().run(
            dataset=_dataset(),
            design=design,
            contrasts=_contrasts(),
            allow_design_subset=False,
            minimum_condition_replicates=1,
            paired_design_policy=PAIRED_DESIGN_POLICY_FIXED_BLOCK,
        )


def test_differential_block_fixed_block_rejects_condition_confounded_with_block() -> (
    None
):
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(sample_id="A_1", condition="A", block_id="pair_1"),
            SampleDesignRecord(sample_id="A_2", condition="A", block_id="pair_1"),
            SampleDesignRecord(sample_id="B_1", condition="B", block_id="pair_2"),
            SampleDesignRecord(sample_id="B_2", condition="B", block_id="pair_2"),
        )
    )

    with pytest.raises(
        WorkflowValidationError,
        match="condition perfectly confounded with block.*A, B",
    ):
        ExperimentalDesignContractValidator().run(
            dataset=_dataset(),
            design=design,
            contrasts=_contrasts(),
            allow_design_subset=False,
            minimum_condition_replicates=1,
            paired_design_policy=PAIRED_DESIGN_POLICY_FIXED_BLOCK,
        )


def test_differential_block_fixed_block_rejects_rank_deficient_design() -> None:
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(
                sample_id="A_1",
                condition="A",
                block_id="pair_1",
                covariates={"sex": "F"},
            ),
            SampleDesignRecord(
                sample_id="A_2",
                condition="A",
                block_id="pair_2",
                covariates={"sex": "F"},
            ),
            SampleDesignRecord(
                sample_id="B_1",
                condition="B",
                block_id="pair_1",
                covariates={"sex": "M"},
            ),
            SampleDesignRecord(
                sample_id="B_2",
                condition="B",
                block_id="pair_2",
                covariates={"sex": "M"},
            ),
        ),
        fixed_effects=(CategoricalCovariate("sex"),),
    )

    with pytest.raises(
        WorkflowValidationError,
        match="rank deficient.*confounded",
    ):
        ExperimentalDesignContractValidator().run(
            dataset=_dataset(),
            design=design,
            contrasts=_contrasts(),
            allow_design_subset=False,
            minimum_condition_replicates=1,
            paired_design_policy=PAIRED_DESIGN_POLICY_FIXED_BLOCK,
        )


def test_differential_block_fixed_block_rejects_non_estimable_contrast() -> None:
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(sample_id="A_1", condition="A", block_id="pair_1"),
            SampleDesignRecord(sample_id="A_2", condition="C", block_id="pair_1"),
            SampleDesignRecord(sample_id="B_1", condition="B", block_id="pair_2"),
            SampleDesignRecord(sample_id="B_2", condition="C", block_id="pair_2"),
        )
    )

    with pytest.raises(
        WorkflowValidationError,
        match="contrast 'B_vs_A' is non-estimable.*no block contains both",
    ):
        ExperimentalDesignContractValidator().run(
            dataset=_dataset(),
            design=design,
            contrasts=_contrasts(),
            allow_design_subset=False,
            minimum_condition_replicates=1,
            paired_design_policy=PAIRED_DESIGN_POLICY_FIXED_BLOCK,
        )


def test_differential_block_fixed_block_rejects_silent_sample_drop() -> None:
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(sample_id="A_1", condition="A", block_id="pair_1"),
            SampleDesignRecord(sample_id="B_1", condition="B", block_id="pair_1"),
        )
    )

    with pytest.raises(
        WorkflowValidationError,
        match="no samples are silently dropped.*A_2, B_2",
    ):
        ExperimentalDesignContractValidator().run(
            dataset=_dataset(),
            design=design,
            contrasts=_contrasts(),
            allow_design_subset=True,
            minimum_condition_replicates=1,
            paired_design_policy=PAIRED_DESIGN_POLICY_FIXED_BLOCK,
        )


def test_contract_rejects_undeclared_covariate_values() -> None:
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(
                sample_id="A_1",
                condition="A",
                covariates={"dose": 1.0},
            ),
            SampleDesignRecord(sample_id="A_2", condition="A"),
            SampleDesignRecord(sample_id="B_1", condition="B"),
            SampleDesignRecord(sample_id="B_2", condition="B"),
        )
    )
    with pytest.raises(WorkflowValidationError, match="undeclared=dose"):
        ExperimentalDesignContractValidator().run(
            dataset=_dataset(),
            design=design,
            contrasts=_contrasts(),
            allow_design_subset=False,
            minimum_condition_replicates=1,
        )


def test_categorical_covariate_declaration_records_contract() -> None:
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(
                sample_id="A_1",
                condition="A",
                covariates={"sex": "F"},
            ),
            SampleDesignRecord(
                sample_id="A_2",
                condition="A",
                covariates={"sex": "M"},
            ),
            SampleDesignRecord(
                sample_id="B_1",
                condition="B",
                covariates={"sex": "F"},
            ),
            SampleDesignRecord(
                sample_id="B_2",
                condition="B",
                covariates={"sex": "M"},
            ),
        ),
        fixed_effects=(
            CategoricalCovariate(
                name="sex",
                required=True,
                include_in_model=False,
            ),
        ),
    )

    validated = ExperimentalDesignContractValidator().run(
        dataset=_dataset(),
        design=design,
        contrasts=_contrasts(),
        allow_design_subset=False,
        minimum_condition_replicates=1,
    )

    assert validated.design is design
    assert design.fixed_effects[0].name == "sex"
    assert design.fixed_effects[0].kind == "categorical"
    assert design.fixed_effects[0].required is True
    assert design.fixed_effects[0].include_in_model is False
    assert validated.design_frame.columns.tolist() == ["A", "B"]


def test_continuous_covariate_declaration_records_contract() -> None:
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(
                sample_id="A_1",
                condition="A",
                covariates={"dose": 0.0},
            ),
            SampleDesignRecord(
                sample_id="A_2",
                condition="A",
                covariates={"dose": 1.0},
            ),
            SampleDesignRecord(
                sample_id="B_1",
                condition="B",
                covariates={"dose": 2.0},
            ),
            SampleDesignRecord(
                sample_id="B_2",
                condition="B",
                covariates={"dose": 3.0},
            ),
        ),
        fixed_effects=(
            ContinuousCovariate(
                name="dose",
                required=True,
                include_in_model=False,
            ),
        ),
    )

    validated = ExperimentalDesignContractValidator().run(
        dataset=_dataset(),
        design=design,
        contrasts=_contrasts(),
        allow_design_subset=False,
        minimum_condition_replicates=1,
    )

    assert validated.design is design
    assert design.fixed_effects[0].name == "dose"
    assert design.fixed_effects[0].kind == "continuous"
    assert design.fixed_effects[0].required is True
    assert design.fixed_effects[0].include_in_model is False
    assert validated.design_frame.columns.tolist() == ["A", "B"]


def test_differential_validation_accepts_condition_and_categorical_covariate() -> None:
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(
                sample_id="A_1",
                condition="A",
                covariates={"sex": "F"},
            ),
            SampleDesignRecord(
                sample_id="A_2",
                condition="A",
                covariates={"sex": "M"},
            ),
            SampleDesignRecord(
                sample_id="B_1",
                condition="B",
                covariates={"sex": "F"},
            ),
            SampleDesignRecord(
                sample_id="B_2",
                condition="B",
                covariates={"sex": "M"},
            ),
        ),
        fixed_effects=(CategoricalCovariate("sex"),),
    )

    validated = ExperimentalDesignContractValidator().run(
        dataset=_dataset(),
        design=design,
        contrasts=_contrasts(),
        allow_design_subset=False,
        minimum_condition_replicates=1,
    )

    assert validated.design_frame.columns.tolist() == ["A", "B", "sex[M]"]
    assert validated.contrast_frame.index.tolist() == ["A", "B", "sex[M]"]
    assert validated.contrast_frame.loc[:, "B_vs_A"].tolist() == [-1.0, 1.0, 0.0]


def test_differential_validation_accepts_condition_and_continuous_covariate() -> None:
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(
                sample_id="A_1",
                condition="A",
                covariates={"dose": 0.0},
            ),
            SampleDesignRecord(
                sample_id="A_2",
                condition="A",
                covariates={"dose": 1.0},
            ),
            SampleDesignRecord(
                sample_id="B_1",
                condition="B",
                covariates={"dose": 0.0},
            ),
            SampleDesignRecord(
                sample_id="B_2",
                condition="B",
                covariates={"dose": 1.0},
            ),
        ),
        fixed_effects=(ContinuousCovariate("dose"),),
    )

    validated = ExperimentalDesignContractValidator().run(
        dataset=_dataset(),
        design=design,
        contrasts=_contrasts(),
        allow_design_subset=False,
        minimum_condition_replicates=1,
    )

    assert validated.design_frame.columns.tolist() == ["A", "B", "dose"]
    assert validated.contrast_frame.index.tolist() == ["A", "B", "dose"]
    assert validated.contrast_frame.loc[:, "B_vs_A"].tolist() == [-1.0, 1.0, 0.0]


def test_batch_covariate_declaration_records_contract() -> None:
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(sample_id="A_1", condition="A", batch="batch_1"),
            SampleDesignRecord(sample_id="A_2", condition="A", batch="batch_1"),
            SampleDesignRecord(sample_id="B_1", condition="B", batch="batch_2"),
            SampleDesignRecord(sample_id="B_2", condition="B", batch="batch_2"),
        ),
        fixed_effects=(BatchCovariate(include_in_model=False),),
    )
    validated = ExperimentalDesignContractValidator().run(
        dataset=_dataset(),
        design=design,
        contrasts=_contrasts(),
        allow_design_subset=False,
        minimum_condition_replicates=1,
    )
    assert validated.design is design
    assert design.fixed_effects[0].name == "batch"
    assert design.fixed_effects[0].kind == "batch"
    assert design.fixed_effects[0].required is True
    assert design.fixed_effects[0].include_in_model is False
    assert validated.design_frame.columns.tolist() == ["A", "B"]


def test_duplicate_covariate_names_rejected() -> None:
    with pytest.raises(WorkflowValidationError, match="duplicate covariate names"):
        ExperimentalDesign(
            samples=_design().samples,
            fixed_effects=(
                CategoricalCovariate("sex"),
                ContinuousCovariate("sex"),
            ),
        )


def test_unsupported_covariate_kind_rejected() -> None:
    with pytest.raises(WorkflowValidationError, match="unsupported covariate kind"):
        FixedEffectCovariate(
            name="age",
            kind="ordinal",  # type: ignore[arg-type]
        )


def test_differential_validation_rejects_missing_modelled_covariate() -> None:
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(
                sample_id="A_1",
                condition="A",
                covariates={"sex": "F"},
            ),
            SampleDesignRecord(sample_id="A_2", condition="A"),
            SampleDesignRecord(
                sample_id="B_1",
                condition="B",
                covariates={"sex": "F"},
            ),
            SampleDesignRecord(
                sample_id="B_2",
                condition="B",
                covariates={"sex": "M"},
            ),
        ),
        fixed_effects=(CategoricalCovariate("sex"),),
    )
    with pytest.raises(
        WorkflowValidationError,
        match="covariate 'sex' is required for modelling but missing.*A_2",
    ):
        ExperimentalDesignContractValidator().run(
            dataset=_dataset(),
            design=design,
            contrasts=_contrasts(),
            allow_design_subset=False,
            minimum_condition_replicates=1,
        )


def test_differential_validation_rejects_non_numeric_continuous_covariate() -> None:
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(
                sample_id="A_1",
                condition="A",
                covariates={"dose": 0.0},
            ),
            SampleDesignRecord(
                sample_id="A_2",
                condition="A",
                covariates={"dose": "2.5"},
            ),
            SampleDesignRecord(
                sample_id="B_1",
                condition="B",
                covariates={"dose": 0.0},
            ),
            SampleDesignRecord(
                sample_id="B_2",
                condition="B",
                covariates={"dose": 1.0},
            ),
        ),
        fixed_effects=(ContinuousCovariate("dose"),),
    )
    with pytest.raises(
        WorkflowValidationError,
        match="continuous covariate 'dose' must be numeric.*A_2",
    ):
        ExperimentalDesignContractValidator().run(
            dataset=_dataset(),
            design=design,
            contrasts=_contrasts(),
            allow_design_subset=False,
            minimum_condition_replicates=1,
        )


def test_differential_validation_rejects_single_level_categorical_covariate() -> None:
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(
                sample_id="A_1",
                condition="A",
                covariates={"sex": "F"},
            ),
            SampleDesignRecord(
                sample_id="A_2",
                condition="A",
                covariates={"sex": "F"},
            ),
            SampleDesignRecord(
                sample_id="B_1",
                condition="B",
                covariates={"sex": "F"},
            ),
            SampleDesignRecord(
                sample_id="B_2",
                condition="B",
                covariates={"sex": "F"},
            ),
        ),
        fixed_effects=(CategoricalCovariate("sex"),),
    )
    with pytest.raises(
        WorkflowValidationError,
        match="categorical fixed-effect covariate 'sex' must have at least two",
    ):
        ExperimentalDesignContractValidator().run(
            dataset=_dataset(),
            design=design,
            contrasts=_contrasts(),
            allow_design_subset=False,
            minimum_condition_replicates=1,
        )


def test_differential_validation_rejects_rank_deficient_fixed_effect_design() -> None:
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(
                sample_id="A_1",
                condition="A",
                covariates={"sex": "F"},
            ),
            SampleDesignRecord(
                sample_id="A_2",
                condition="A",
                covariates={"sex": "F"},
            ),
            SampleDesignRecord(
                sample_id="B_1",
                condition="B",
                covariates={"sex": "M"},
            ),
            SampleDesignRecord(
                sample_id="B_2",
                condition="B",
                covariates={"sex": "M"},
            ),
        ),
        fixed_effects=(CategoricalCovariate("sex"),),
    )
    with pytest.raises(
        WorkflowValidationError,
        match="rank deficient.*confounded",
    ):
        ExperimentalDesignContractValidator().run(
            dataset=_dataset(),
            design=design,
            contrasts=_contrasts(),
            allow_design_subset=False,
            minimum_condition_replicates=1,
        )


def test_differential_validation_rejects_confounded_condition_batch_design() -> None:
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(sample_id="A_1", condition="A", batch="batch_1"),
            SampleDesignRecord(sample_id="A_2", condition="A", batch="batch_1"),
            SampleDesignRecord(sample_id="B_1", condition="B", batch="batch_2"),
            SampleDesignRecord(sample_id="B_2", condition="B", batch="batch_2"),
        ),
        fixed_effects=(BatchCovariate(),),
    )
    with pytest.raises(
        WorkflowValidationError,
        match="rank deficient.*confounded",
    ):
        ExperimentalDesignContractValidator().run(
            dataset=_dataset(),
            design=design,
            contrasts=_contrasts(),
            allow_design_subset=False,
            minimum_condition_replicates=1,
        )


def test_differential_validation_rejects_invalid_contrast_vector_alignment() -> None:
    class _MissingConditionCoefficientBuilder:
        def run(
            self,
            *,
            design: ExperimentalDesign,
            condition_labels: tuple[str, ...],
        ) -> object:
            frame = pd.DataFrame(
                {"A": [1.0, 1.0, 0.0, 0.0]},
                index=pd.Index(design.sample_ids(), name="sample"),
            )
            frame.columns = pd.Index(["A"], name="coefficient")
            return type("_BuildResult", (), {"frame": frame})()

    with pytest.raises(
        WorkflowValidationError,
        match="contrast vectors are invalid.*missing condition coefficients: B",
    ):
        ExperimentalDesignContractValidator(
            design_matrix_builder=(  # type: ignore[arg-type]
                _MissingConditionCoefficientBuilder()
            ),
        ).run(
            dataset=_dataset(),
            design=_design(),
            contrasts=_contrasts(),
            allow_design_subset=False,
            minimum_condition_replicates=1,
        )


def test_deterministic_ordering_of_sample_ids() -> None:
    design = _design()
    validated = ExperimentalDesignContractValidator().run(
        dataset=_dataset(sample_order=("B_2", "A_1", "B_1", "A_2")),
        design=design,
        contrasts=_contrasts(),
        allow_design_subset=False,
        minimum_condition_replicates=2,
    )
    assert design.sample_ids() == ("A_1", "A_2", "B_1", "B_2")
    assert validated.analysis_sample_ids == ("A_1", "A_2", "B_1", "B_2")
