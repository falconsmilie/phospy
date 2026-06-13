from __future__ import annotations

import pandas as pd
import pytest

from phospy.api import DatasetBuildRequest, Organism
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
    assert {stage.stage for stage in built.provenance.preprocessing_stages}.isdisjoint(
        {"protein_aware_preparation"}
    )
