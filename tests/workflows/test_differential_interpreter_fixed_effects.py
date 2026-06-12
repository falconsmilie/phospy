from __future__ import annotations

import pandas as pd

from phospy import AnalysisReadyPhosphoDataset
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
from phospy.workflows.differential.models import (
    InterpretedDifferentialAnalysisRequest,
)
from tests.support.intensity_scale_states import (
    supported_log2_intensity_scale_state,
    supported_log2_processing_state,
)
from tests.support.site_keys import protein_site_key_index, site_key_context_columns


def _dataset() -> AnalysisReadyPhosphoDataset:
    genes = ["MAPK14", "GSK3B"]
    sites = ["Y182", "S9"]
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
    site_metadata = pd.DataFrame(
        {
            "site_key": site_index.tolist(),
            "display_id": ["MAPK14;Y182;", "GSK3B;S9;"],
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
    return AnalysisReadyPhosphoDataset(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        intensity_scale_state=supported_log2_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_log2_processing_state(has_total_matrix=False),
    )


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
    return AnalysisReadyPhosphoDataset(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        intensity_scale_state=supported_log2_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_log2_processing_state(has_total_matrix=False),
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
    result = DifferentialAnalysisWorkflow(
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
