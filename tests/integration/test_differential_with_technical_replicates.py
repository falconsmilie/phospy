from __future__ import annotations

import pandas as pd
import pytest

from phospy import AnalysisReadyPhosphoDataset, DifferentialAnalysisWorkflow
from phospy.advanced import (
    DifferentialAnalysisConfig,
    TechnicalReplicatePolicy,
)
from phospy.api import (
    Contrast,
    DifferentialAnalysisRequest,
    ExperimentalDesign,
    Organism,
    SampleDesignRecord,
)
from phospy.errors import WorkflowValidationError
from tests.support.analysis_ready_dataset_factories import (
    trusted_analysis_ready_dataset_from_tables,
)
from tests.support.intensity_scale_states import (
    supported_log2_intensity_scale_state,
    supported_log2_processing_state,
)
from tests.support.site_keys import (
    site_key_context_columns,
    site_key_index_from_display_ids,
)


def _dataset() -> AnalysisReadyPhosphoDataset:
    display_ids = ["MAPK14;Y182;", "GSK3B;S9;", "AKT1;T308;"]
    site_index = site_key_index_from_display_ids(
        display_ids,
        protein_namespace="gene_symbol",
    )
    phospho = pd.DataFrame(
        {
            "A1_T1": [1.0, 2.0, 1.0],
            "A1_T2": [1.2, 2.2, 1.1],
            "A2_T1": [1.1, 2.1, 1.0],
            "A2_T2": [1.3, 2.0, 1.2],
            "B1_T1": [2.0, 2.3, 0.9],
            "B1_T2": [2.1, 2.4, 0.8],
            "B2_T1": [2.2, 2.5, 0.7],
            "B2_T2": [2.3, 2.6, 0.6],
        },
        index=site_index,
    )
    site_metadata = pd.DataFrame(
        {
            "site_key": site_index.tolist(),
            "display_id": display_ids,
            **site_key_context_columns(site_index),
            "gene_symbol": ["MAPK14", "GSK3B", "AKT1"],
            "site": ["Y182", "S9", "T308"],
            "site_sequence": [
                ("A" * 15) + str(site).strip().upper()[0] + ("A" * 15)
                for site in ["Y182", "S9", "T308"]
            ],
            "protein_id": ["MAPK14", "GSK3B", "AKT1"],
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


def _design() -> ExperimentalDesign:
    return ExperimentalDesign(
        samples=(
            SampleDesignRecord(
                sample_id="A1_T1",
                condition="A",
                biological_replicate_id="A1",
                technical_replicate_id="T1",
            ),
            SampleDesignRecord(
                sample_id="A1_T2",
                condition="A",
                biological_replicate_id="A1",
                technical_replicate_id="T2",
            ),
            SampleDesignRecord(
                sample_id="A2_T1",
                condition="A",
                biological_replicate_id="A2",
                technical_replicate_id="T1",
            ),
            SampleDesignRecord(
                sample_id="A2_T2",
                condition="A",
                biological_replicate_id="A2",
                technical_replicate_id="T2",
            ),
            SampleDesignRecord(
                sample_id="B1_T1",
                condition="B",
                biological_replicate_id="B1",
                technical_replicate_id="T1",
            ),
            SampleDesignRecord(
                sample_id="B1_T2",
                condition="B",
                biological_replicate_id="B1",
                technical_replicate_id="T2",
            ),
            SampleDesignRecord(
                sample_id="B2_T1",
                condition="B",
                biological_replicate_id="B2",
                technical_replicate_id="T1",
            ),
            SampleDesignRecord(
                sample_id="B2_T2",
                condition="B",
                biological_replicate_id="B2",
                technical_replicate_id="T2",
            ),
        )
    )


def _contrast() -> tuple[Contrast, ...]:
    return (
        Contrast(
            name="B_vs_A",
            numerator_condition="B",
            denominator_condition="A",
        ),
    )


def test_differential_workflow_rejects_repeated_biological_replicates_by_default() -> (
    None
):
    with pytest.raises(
        WorkflowValidationError,
        match="Technical replicates require explicit aggregation",
    ):
        DifferentialAnalysisWorkflow().run(
            DifferentialAnalysisRequest(
                dataset=_dataset(),
                design=_design(),
                contrasts=_contrast(),
            )
        )


def test_differential_workflow_runs_after_mean_technical_replicate_aggregation() -> (
    None
):
    result = DifferentialAnalysisWorkflow().run(
        DifferentialAnalysisRequest(
            dataset=_dataset(),
            design=_design(),
            contrasts=_contrast(),
            config=DifferentialAnalysisConfig(
                technical_replicate_policy=TechnicalReplicatePolicy.MEAN
            ),
        )
    )

    table = result.table_for("B_vs_A")
    assert {"logFC", "t", "P.Value", "adj.P.Val"} <= set(table.columns)
    assert table.index.tolist() == _dataset().phospho.index.tolist()
    assert result.workflow_provenance is not None
    assert result.workflow_provenance["technical_replicate_policy"] == "mean"
    assert result.workflow_provenance["aggregation_policy"] == "mean"
    assert result.workflow_provenance["aggregation_method"] == "mean"
    assert result.workflow_provenance["matrices_aggregated"] == {
        "phospho": True,
        "total_protein": False,
    }
    assert result.workflow_provenance["both_phospho_and_total_aggregated"] is False
    groups = result.workflow_provenance["groups"]
    assert isinstance(groups, list)
    assert len(groups) == 4
    assert result.policy_provenance is not None
    assert result.policy_provenance.replicates.technical_replicate_policy == "mean"
    assert len(result.policy_provenance.replicates.technical_replicate_groups) == 4


def test_differential_policy_provenance_is_deterministic_for_identical_inputs() -> None:
    request = DifferentialAnalysisRequest(
        dataset=_dataset(),
        design=_design(),
        contrasts=_contrast(),
        config=DifferentialAnalysisConfig(
            technical_replicate_policy=TechnicalReplicatePolicy.MEAN
        ),
    )
    first = DifferentialAnalysisWorkflow().run(request)
    second = DifferentialAnalysisWorkflow().run(request)
    assert first.policy_provenance == second.policy_provenance
