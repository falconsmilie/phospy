from __future__ import annotations

import pandas as pd
import pytest

from phospy import AnalysisReadyPhosphoDataset
from phospy.api import (
    Contrast,
    ExperimentalDesign,
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


def test_optional_block_field_validation_fails_when_partially_defined() -> None:
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(sample_id="A_1", condition="A", block="pair_1"),
            SampleDesignRecord(sample_id="A_2", condition="A"),
            SampleDesignRecord(sample_id="B_1", condition="B", block="pair_1"),
            SampleDesignRecord(sample_id="B_2", condition="B", block="pair_2"),
        )
    )
    with pytest.raises(WorkflowValidationError, match="optional field 'block'"):
        ExperimentalDesignContractValidator().run(
            dataset=_dataset(),
            design=design,
            contrasts=_contrasts(),
            allow_design_subset=False,
            minimum_condition_replicates=1,
        )


def test_supported_contract_reports_unsupported_batch_modelling() -> None:
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(sample_id="A_1", condition="A", batch="batch_1"),
            SampleDesignRecord(sample_id="A_2", condition="A", batch="batch_1"),
            SampleDesignRecord(sample_id="B_1", condition="B", batch="batch_2"),
            SampleDesignRecord(sample_id="B_2", condition="B", batch="batch_2"),
        )
    )
    with pytest.raises(WorkflowValidationError, match="unsupported design features"):
        ExperimentalDesignContractValidator().run(
            dataset=_dataset(),
            design=design,
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
