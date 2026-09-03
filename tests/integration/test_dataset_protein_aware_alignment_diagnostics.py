from __future__ import annotations

import pandas as pd
import pytest

from phospy.advanced import (
    DatasetProteinAwarePreparationConfig,
    DifferentialAnalysisConfig,
    DifferentialProteinAwareModelConfig,
)
from phospy.api import (
    Contrast,
    DatasetBuildRequest,
    DatasetPreprocessingConfig,
    DifferentialAnalysisRequest,
    ExperimentalDesign,
    Organism,
    SampleDesignRecord,
)
from phospy.science.datasets.builders.executor import DatasetBuildExecutor
from phospy.science.datasets.builders.interpreter import DatasetBuildRequestInterpreter
from phospy.science.datasets.preprocessing.protein_aware_alignment import (
    ProteinAwareAlignmentConfig,
    ProteinAwareAlignmentEligibilityResolver,
)
from phospy.science.datasets.preprocessing.protein_mapping import (
    ProteinMappingConfig,
    ProteinMappingResolver,
)
from phospy.workflows.differential.protein_aware_inputs import (
    ProteinAwareDifferentialInputResolver,
)
from phospy.workflows.differential.validator import DifferentialAnalysisValidator

pytestmark = pytest.mark.integration


def test_dataset_builder_outputs_support_protein_aware_alignment_diagnostics() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 2.0],
            "sample_b": [3.0, 4.0],
        },
        index=pd.Index(["MAPK14;Y182;", "AKT1;T308;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "AKT1"],
            "site": ["Y182", "T308"],
            "protein_id": ["P53778", "P31749"],
            "site_sequence": ["AAAAAYAAAAA", "AAAAATAAAAA"],
            "localisation_confidence": [0.95, 0.96],
        },
        index=phospho.index.copy(),
    )
    total = pd.DataFrame(
        {
            "sample_a": [5.0, 6.0],
            "sample_b": [7.0, 8.0],
        },
        index=pd.Index(["P53778", "P31749"], name="protein_id"),
    )

    built = DatasetBuildExecutor().run(
        DatasetBuildRequestInterpreter().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=site_metadata,
                total=total,
                organism=Organism.RAT,
                input_intensity_scale="log2",
            )
        )
    )

    assert built.total is not None
    mapping_result = ProteinMappingResolver().run(
        site_metadata=built.site_metadata,
        phospho_matrix_index=built.phospho.index,
        total_protein_matrix_index=built.total.index,
        config=ProteinMappingConfig(protein_identifier_columns=("protein_id",)),
    )
    diagnostics = ProteinAwareAlignmentEligibilityResolver().run(
        phospho=built.phospho,
        total=built.total,
        mapping_result=mapping_result,
        intensity_scale_state=built.intensity_scale_state,
        config=ProteinAwareAlignmentConfig(
            protein_mapping_policy="allow_missing_with_report"
        ),
    )

    assert diagnostics.sample_alignment.sample_order_compatible is True
    assert diagnostics.transformation_state.compatible is True
    assert set(diagnostics.eligible_for_protein_aware_preparation) == set(
        built.phospho.index.astype(str).tolist()
    )
    assert diagnostics.fallback_to_phospho_only == ()
    assert diagnostics.excluded_from_preparation == ()
    assert built.processing_state.total_protein_correction.policy == "none"
    assert built.provenance is not None
    assert "protein_aware_preparation" not in built.provenance.workflow_parameters
    assert {stage.stage for stage in built.provenance.preprocessing_stages}.isdisjoint(
        {"protein_aware_preparation"}
    )


def test_dataset_builder_integrates_protein_aware_preparation_report() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 2.0],
            "sample_b": [3.0, 4.0],
        },
        index=pd.Index(["MAPK14;Y182;", "AKT1;T308;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "AKT1"],
            "site": ["Y182", "T308"],
            "protein_id": ["P53778", "P31749"],
            "site_sequence": ["AAAAAYAAAAA", "AAAAATAAAAA"],
            "localisation_confidence": [0.95, 0.96],
        },
        index=phospho.index.copy(),
    )
    total = pd.DataFrame(
        {
            "sample_a": [5.0, 6.0],
            "sample_b": [7.0, 8.0],
        },
        index=pd.Index(["P53778", "P31749"], name="protein_id"),
    )

    built = DatasetBuildExecutor().run(
        DatasetBuildRequestInterpreter().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=site_metadata,
                total=total,
                organism=Organism.RAT,
                input_intensity_scale="log2",
                preprocessing_config=DatasetPreprocessingConfig(
                    protein_aware_preparation=(
                        DatasetProteinAwarePreparationConfig(
                            policy="prepare_model_inputs"
                        )
                    )
                ),
            )
        )
    )

    assert built.protein_aware_preparation is not None
    preparation = built.protein_aware_preparation
    expected_site_keys = tuple(built.phospho.index.astype(str).tolist())
    assert preparation.report.eligible_site_keys == expected_site_keys
    assert preparation.report.fallback_site_keys == ()
    assert preparation.report.excluded_site_keys == ()
    assert preparation.protein_covariate_matrix.to_dict(orient="index") == {
        "P53778": {"sample_a": 5.0, "sample_b": 7.0},
        "P31749": {"sample_a": 6.0, "sample_b": 8.0},
    }
    assert built.preprocessing_report is not None
    assert built.preprocessing_report.protein_aware_preparation is preparation.report
    operation_stages = (
        built.preprocessing_report.operations.loc[:, "stage"].astype(str).tolist()
    )
    assert "protein_aware_preparation" in operation_stages
    operation = built.preprocessing_report.operations.loc[
        built.preprocessing_report.operations.loc[:, "stage"]
        == "protein_aware_preparation"
    ].iloc[0]
    assert operation["operation"] == "prepare_model_inputs"
    assert operation["parameters"]["eligible_site_count"] == 2
    assert operation["parameters"]["performs_model_adjustment"] is False
    assert operation["parameters"]["performs_differential_modelling"] is False
    assert operation["parameters"]["claims_msstatsptm_equivalence"] is False
    assert built.provenance is not None
    preparation_provenance = built.provenance.workflow_parameters[
        "protein_aware_preparation"
    ]
    assert preparation_provenance["status"] == "prepared"
    assert preparation_provenance["preparation_policy"] == "prepare_model_inputs"
    assert preparation_provenance["protein_mapping_policy"] == "require_unambiguous"
    assert preparation_provenance["eligible_site_count"] == 2
    assert preparation_provenance["performs_total_protein_subtraction"] is False
    assert preparation_provenance["performs_normalisation"] is False
    assert preparation_provenance["performs_differential_modelling"] is False
    assert preparation_provenance["claims_msstatsptm_equivalence"] is False


def test_dataset_builder_protein_aware_sidecar_feeds_private_differential_inputs() -> (
    None
):
    phospho = pd.DataFrame(
        {
            "A_1": [1.0, 2.0],
            "A_2": [1.1, 2.1],
            "B_1": [2.1, 2.0],
            "B_2": [2.0, 2.2],
        },
        index=pd.Index(["MAPK14;Y182;", "AKT1;T308;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "AKT1"],
            "site": ["Y182", "T308"],
            "protein_id": ["P53778", "P31749"],
            "site_sequence": ["AAAAAYAAAAA", "AAAAATAAAAA"],
            "localisation_confidence": [0.95, 0.96],
        },
        index=phospho.index.copy(),
    )
    total = pd.DataFrame(
        {
            "A_1": [10.0, 20.0],
            "A_2": [11.0, 21.0],
            "B_1": [13.0, 24.0],
            "B_2": [12.0, 23.0],
        },
        index=pd.Index(["P53778", "P31749"], name="protein_id"),
    )
    built = DatasetBuildExecutor().run(
        DatasetBuildRequestInterpreter().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=site_metadata,
                total=total,
                organism=Organism.RAT,
                input_intensity_scale="log2",
                preprocessing_config=DatasetPreprocessingConfig(
                    protein_aware_preparation=DatasetProteinAwarePreparationConfig(
                        policy="prepare_model_inputs"
                    )
                ),
            )
        )
    )
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(sample_id="A_1", condition="A"),
            SampleDesignRecord(sample_id="A_2", condition="A"),
            SampleDesignRecord(sample_id="B_1", condition="B"),
            SampleDesignRecord(sample_id="B_2", condition="B"),
        )
    )

    resolved = ProteinAwareDifferentialInputResolver().run(
        DifferentialAnalysisValidator().run(
            DifferentialAnalysisRequest(
                dataset=built,
                design=design,
                contrasts=(
                    Contrast(
                        name="B_vs_A",
                        numerator_condition="B",
                        denominator_condition="A",
                    ),
                ),
                config=DifferentialAnalysisConfig(
                    protein_aware_model=DifferentialProteinAwareModelConfig()
                ),
            )
        ),
        minimum_condition_replicates=2,
    )

    assert resolved.tested_site_ids == tuple(built.phospho.index.astype(str).tolist())
    assert tuple(resolved.resolved_protein_covariates.columns.astype(str)) == (
        "A_1",
        "A_2",
        "B_1",
        "B_2",
    )
    assert tuple(
        resolved.computation_request.matched_pairs.loc[:, "total_protein_row_key"]
        .astype(str)
        .tolist()
    ) == ("P53778", "P31749")
