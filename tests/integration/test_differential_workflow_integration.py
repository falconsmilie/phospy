from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from phospy import (
    AnalysisReadyDatasetBuilder,
    AnalysisReadyPhosphoDataset,
    DifferentialAnalysisWorkflow,
)
from phospy.advanced import DatasetIntensityTransformConfig
from phospy.api import (
    Contrast,
    DatasetBuildRequest,
    DatasetPreprocessingConfig,
    DifferentialAnalysisRequest,
    ExperimentalDesign,
    Organism,
    SampleDesignRecord,
)
from phospy.errors import WorkflowValidationError
from phospy.science.differential.models import (
    DIFFERENTIAL_RESULT_STATUS_COLUMN,
    DIFFERENTIAL_RESULT_STATUS_TESTED,
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_ALL_CONSTANT,
)
from phospy.science.sites.site_keys import (
    build_protein_scoped_site_key,
    encode_site_key,
)
from phospy.science.statistics.multiple_testing import adjust_p_values
from phospy.science.transformations.models import (
    IntensityScaleState,
    MatrixIntensityScaleState,
)
from tests.support.analysis_ready_dataset_factories import (
    trusted_analysis_ready_dataset_from_tables,
)
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
    supported_log2_intensity_scale_state,
    supported_log2_processing_state,
)
from tests.support.site_keys import site_key_context_columns
from tests.support.unsafe_dataset_states import (
    unsafe_replace_dataset_intensity_scale_state,
)


def _dataset() -> AnalysisReadyPhosphoDataset:
    phospho = pd.DataFrame(
        {
            "A_1": [1.0, 2.0, 1.0, 1.2],
            "A_2": [1.1, 2.1, 1.1, 1.1],
            "B_1": [2.1, 2.0, 1.0, 1.3],
            "B_2": [2.0, 2.2, 0.9, 1.2],
        },
        index=pd.Index(
            ["MAPK14;Y182;", "GSK3B;S9;", "AKT1;T308;", "RPS6KB1;T389;"],
            name="display_id",
        ),
    )
    display_ids = phospho.index.astype(str).tolist()
    parsed = [site_id.split(";") for site_id in display_ids]
    protein_ids = [parts[0] for parts in parsed]
    site_keys = [
        _site_key_for_display_id(display_id, protein_id=protein_id)
        for display_id, protein_id in zip(display_ids, protein_ids, strict=True)
    ]
    phospho.index = pd.Index(site_keys, name="site_key")
    site_metadata = pd.DataFrame(
        {
            "site_key": site_keys,
            "display_id": display_ids,
            **site_key_context_columns(site_keys),
            "gene_symbol": ["MAPK14", "GSK3B", "AKT1", "RPS6KB1"],
            "site": ["Y182", "S9", "T308", "T389"],
            "site_sequence": [
                ("A" * 15) + str(site).strip().upper()[0] + ("A" * 15)
                for site in ["Y182", "S9", "T308", "T389"]
            ],
            "localisation_confidence": [0.95] * phospho.shape[0],
            "protein_id": ["MAPK14", "GSK3B", "AKT1", "RPS6KB1"],
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


def _site_key_for_display_id(
    display_id: str,
    *,
    protein_id: str,
) -> str:
    parts = [token.strip() for token in display_id.split(";") if token.strip()]
    site = parts[1]
    key = build_protein_scoped_site_key(
        organism="rat",
        protein_namespace="protein_id",
        protein_identifier=protein_id,
        residue=site.upper()[0],
        position=int(site[1:]),
        field_name="tests.integration.test_differential_workflow_integration.site_key",
        error_type=ValueError,
    )
    return encode_site_key(key)


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


def test_differential_workflow_runs_on_analysis_ready_dataset() -> None:
    result = DifferentialAnalysisWorkflow().run(
        DifferentialAnalysisRequest(
            dataset=_dataset(),
            design=_design(),
            contrasts=_contrasts(),
        )
    )

    table = result.table_for("B_vs_A")
    assert list(table.columns) == [
        "site_key",
        "display_id",
        "organism",
        "protein_namespace",
        "protein_identifier",
        "gene_symbol",
        "site",
        "protein_id",
        "logFC",
        "t",
        "P.Value",
        "adj.P.Val",
    ]
    assert table.index.tolist() == _dataset().phospho.index.tolist()
    assert (table.loc[:, "adj.P.Val"] >= 0.0).all()
    assert (table.loc[:, "adj.P.Val"] <= 1.0).all()


def test_differential_workflow_excludes_withheld_rows_from_multiple_testing() -> None:
    base_dataset = _dataset()
    matrix = base_dataset.phospho
    matrix.iloc[0, :] = 5.0
    dataset = trusted_analysis_ready_dataset_from_tables(
        phospho=matrix,
        site_metadata=base_dataset.site_metadata,
        organism=base_dataset.organism,
        intensity_scale_state=base_dataset.intensity_scale_state,
        processing_state=base_dataset.processing_state,
    )

    result = DifferentialAnalysisWorkflow().run(
        DifferentialAnalysisRequest(
            dataset=dataset,
            design=_design(),
            contrasts=_contrasts(),
        )
    )
    table = result.table_for("B_vs_A")
    tested_rows = (
        table[DIFFERENTIAL_RESULT_STATUS_COLUMN] == DIFFERENTIAL_RESULT_STATUS_TESTED
    )
    withheld_rows = (
        table[DIFFERENTIAL_RESULT_STATUS_COLUMN]
        == DIFFERENTIAL_RESULT_STATUS_WITHHELD_ALL_CONSTANT
    )

    assert int(tested_rows.sum()) == 3
    assert int(withheld_rows.sum()) == 1
    assert (
        table.loc[withheld_rows, ["logFC", "t", "P.Value", "adj.P.Val"]]
        .isna()
        .all()
        .all()
    )
    expected_adjusted = adjust_p_values(
        table.loc[tested_rows, "P.Value"].to_numpy(dtype=float),
        method="benjamini_hochberg",
    )
    np.testing.assert_allclose(
        table.loc[tested_rows, "adj.P.Val"].to_numpy(dtype=float),
        expected_adjusted,
        rtol=1e-12,
        atol=1e-12,
    )


def test_differential_workflow_invalid_contrast_fails_before_execution() -> None:
    with pytest.raises(WorkflowValidationError):
        DifferentialAnalysisWorkflow().run(
            DifferentialAnalysisRequest(
                dataset=_dataset(),
                design=_design(),
                contrasts=(
                    Contrast(
                        name="B_vs_A",
                        numerator_condition="B",
                        denominator_condition="A_bad",
                    ),
                ),
            )
        )


def test_differential_workflow_misaligned_design_fails_before_execution() -> None:
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(sample_id="X1", condition="A"),
            SampleDesignRecord(sample_id="X2", condition="A"),
            SampleDesignRecord(sample_id="X3", condition="B"),
            SampleDesignRecord(sample_id="X4", condition="B"),
        )
    )
    with pytest.raises(WorkflowValidationError):
        DifferentialAnalysisWorkflow().run(
            DifferentialAnalysisRequest(
                dataset=_dataset(),
                design=design,
                contrasts=_contrasts(),
            )
        )


def test_differential_workflow_runs_on_builder_log2_dataset() -> None:
    phospho = pd.DataFrame(
        {
            "A_1": [10.0, 8.0, 6.0, 4.0],
            "A_2": [11.0, 8.5, 5.5, 4.5],
            "B_1": [15.0, 9.0, 6.5, 4.8],
            "B_2": [16.0, 9.2, 6.6, 5.0],
        },
        index=pd.Index(
            ["MAPK14;Y182;", "GSK3B;S9;", "AKT1;T308;", "RPS6KB1;T389;"],
            name="site_id",
        ),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "GSK3B", "AKT1", "RPS6KB1"],
            "site": ["Y182", "S9", "T308", "T389"],
            "site_sequence": [
                ("A" * 15) + str(site).strip().upper()[0] + ("A" * 15)
                for site in ["Y182", "S9", "T308", "T389"]
            ],
            "localisation_confidence": [0.95] * phospho.shape[0],
            "protein_id": ["MAPK14", "GSK3B", "AKT1", "RPS6KB1"],
        },
        index=phospho.index.copy(),
    )
    dataset = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
            preprocessing_config=DatasetPreprocessingConfig(
                intensity_transform=DatasetIntensityTransformConfig(
                    policy="log2",
                    pseudocount=1.0,
                )
            ),
        )
    )

    result = DifferentialAnalysisWorkflow().run(
        DifferentialAnalysisRequest(
            dataset=dataset,
            design=_design(),
            contrasts=_contrasts(),
        )
    )
    assert "B_vs_A" in result.contrast_tables


def test_documented_two_vs_two_differential_example_contract() -> None:
    phospho = pd.DataFrame(
        {
            "control_rep1": [8200.0, 9100.0],
            "control_rep2": [8000.0, 9000.0],
            "treatment_rep1": [16200.0, 9150.0],
            "treatment_rep2": [15800.0, 9050.0],
        },
        index=pd.Index(
            ["MAPK14;Y182;", "GSK3B;S9;"],
            name="site_id",
        ),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "GSK3B"],
            "site": ["Y182", "S9"],
            "site_sequence": [
                "MPRKSLVGTPYWMNQYAVNQKQTLRDLKQEN",
                "ATMSGRPRTTSFAESSKPVQQPSAFGQAAAL",
            ],
            "protein_id": ["MAPK14", "GSK3B"],
            "localisation_confidence": [0.95, 0.95],
        },
        index=phospho.index.copy(),
    )
    dataset = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
            preprocessing_config=DatasetPreprocessingConfig(
                intensity_transform=DatasetIntensityTransformConfig(
                    policy="log2",
                    pseudocount=1.0,
                )
            ),
        )
    )

    assert dataset.intensity_scale_state.kind.value == "log2"
    assert dataset.intensity_scale_state.is_established

    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(
                sample_id="control_rep1",
                condition="control",
                biological_replicate_id="control_r1",
            ),
            SampleDesignRecord(
                sample_id="control_rep2",
                condition="control",
                biological_replicate_id="control_r2",
            ),
            SampleDesignRecord(
                sample_id="treatment_rep1",
                condition="treatment",
                biological_replicate_id="treatment_r1",
            ),
            SampleDesignRecord(
                sample_id="treatment_rep2",
                condition="treatment",
                biological_replicate_id="treatment_r2",
            ),
        )
    )
    result = DifferentialAnalysisWorkflow().run(
        DifferentialAnalysisRequest(
            dataset=dataset,
            design=design,
            contrasts=(
                Contrast(
                    name="treatment_vs_control",
                    numerator_condition="treatment",
                    denominator_condition="control",
                ),
            ),
        )
    )

    table = result.table_for("treatment_vs_control")
    mapk14_site_key = _site_key_for_display_id("MAPK14;Y182;", protein_id="MAPK14")
    gsk3b_site_key = _site_key_for_display_id("GSK3B;S9;", protein_id="GSK3B")
    assert float(table.loc[mapk14_site_key, "logFC"]) > 0.7
    assert abs(float(table.loc[gsk3b_site_key, "logFC"])) < 0.1


def test_differential_workflow_rejects_linear_scale_before_execution() -> None:
    valid_dataset = _dataset()
    linear_dataset = trusted_analysis_ready_dataset_from_tables(
        phospho=valid_dataset.phospho,
        site_metadata=valid_dataset.site_metadata,
        organism=valid_dataset.organism,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
    )

    with pytest.raises(
        WorkflowValidationError,
        match="requires established log2-scale phospho intensities",
    ):
        DifferentialAnalysisWorkflow().run(
            DifferentialAnalysisRequest(
                dataset=linear_dataset,
                design=_design(),
                contrasts=_contrasts(),
            )
        )


def test_differential_workflow_rejects_unestablished_log2_scale_before_execution() -> (
    None
):
    dataset = _dataset()
    unsafe_replace_dataset_intensity_scale_state(
        dataset,
        IntensityScaleState(
            phospho=MatrixIntensityScaleState.log2(established_by="test.declaration"),
            total=None,
        ),
    )

    with pytest.raises(
        WorkflowValidationError,
        match="requires established log2-scale phospho intensities",
    ):
        DifferentialAnalysisWorkflow().run(
            DifferentialAnalysisRequest(
                dataset=dataset,
                design=_design(),
                contrasts=_contrasts(),
            )
        )


def test_differential_workflow_rejects_unknown_scale_before_execution() -> None:
    dataset = _dataset()
    unknown_state = IntensityScaleState.raw(has_total_matrix=False)
    unsafe_replace_dataset_intensity_scale_state(dataset, unknown_state)

    with pytest.raises(
        WorkflowValidationError,
        match="requires established log2-scale phospho intensities",
    ):
        DifferentialAnalysisWorkflow().run(
            DifferentialAnalysisRequest(
                dataset=dataset,
                design=_design(),
                contrasts=_contrasts(),
            )
        )
