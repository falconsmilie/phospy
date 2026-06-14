from __future__ import annotations

from dataclasses import replace

import pandas as pd
import pytest

from phospy import AnalysisReadyPhosphoDataset
from phospy.api import Organism
from phospy.errors import WorkflowValidationError
from phospy.science.sites.site_keys import (
    build_protein_scoped_site_key,
    encode_site_key,
)
from phospy.science.transformations.models import (
    IntensityScaleState,
    MatrixIntensityScaleState,
    QuantitativeMeaning,
)
from phospy.validation.workflows.differential import (
    DifferentialDatasetEligibilityValidator,
)
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
    supported_log2_intensity_scale_state,
    supported_log2_processing_state,
)
from tests.support.site_keys import site_key_context_columns


def _frames() -> tuple[pd.DataFrame, pd.DataFrame]:
    site_keys = [
        encode_site_key(
            build_protein_scoped_site_key(
                organism="rat",
                protein_namespace="protein_id",
                protein_identifier="MAPK14",
                residue="Y",
                position=182,
                field_name="test.site_key.mapk14",
                error_type=ValueError,
            )
        ),
        encode_site_key(
            build_protein_scoped_site_key(
                organism="rat",
                protein_namespace="protein_id",
                protein_identifier="AKT1",
                residue="T",
                position=308,
                field_name="test.site_key.akt1",
                error_type=ValueError,
            )
        ),
    ]
    display_ids = ["MAPK14;Y182;", "AKT1;T308;"]
    phospho = pd.DataFrame(
        {
            "A_1": [1.0, 2.0],
            "A_2": [1.1, 2.1],
            "B_1": [2.1, 2.0],
            "B_2": [2.0, 2.2],
        },
        index=pd.Index(site_keys, name="site_key"),
    )
    site_metadata = pd.DataFrame(
        {
            "site_key": site_keys,
            "display_id": display_ids,
            **site_key_context_columns(site_keys),
            "gene_symbol": ["MAPK14", "AKT1"],
            "site": ["Y182", "T308"],
            "site_sequence": [
                ("A" * 15) + str(site).strip().upper()[0] + ("A" * 15)
                for site in ["Y182", "T308"]
            ],
            "protein_id": ["MAPK14", "AKT1"],
        },
        index=phospho.index.copy(),
    )
    return phospho, site_metadata


def _dataset_with_log2_scale() -> AnalysisReadyPhosphoDataset:
    phospho, site_metadata = _frames()
    return AnalysisReadyPhosphoDataset(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        intensity_scale_state=supported_log2_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_log2_processing_state(has_total_matrix=False),
    )


def _dataset_with_linear_scale() -> AnalysisReadyPhosphoDataset:
    phospho, site_metadata = _frames()
    return AnalysisReadyPhosphoDataset(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
    )


def _dataset_with_imputed_missing_data_state() -> AnalysisReadyPhosphoDataset:
    dataset = _dataset_with_log2_scale()
    processing_state = dataset.processing_state
    return AnalysisReadyPhosphoDataset(
        phospho=dataset.phospho,
        site_metadata=dataset.site_metadata,
        organism=dataset.organism,
        intensity_scale_state=dataset.intensity_scale_state,
        processing_state=replace(
            processing_state,
            missing_data=replace(processing_state.missing_data, imputed=True),
        ),
    )


def test_differential_accepts_non_imputed_dataset() -> None:
    DifferentialDatasetEligibilityValidator().run(dataset=_dataset_with_log2_scale())


def test_differential_rejects_imputed_dataset_by_default() -> None:
    with pytest.raises(WorkflowValidationError) as exc_info:
        DifferentialDatasetEligibilityValidator().run(
            dataset=_dataset_with_imputed_missing_data_state()
        )

    message = str(exc_info.value)
    assert "Differential analysis" in message
    assert "imputed cells" in message
    assert "observed measurements" in message


def test_differential_imputation_error_is_actionable() -> None:
    with pytest.raises(WorkflowValidationError) as exc_info:
        DifferentialDatasetEligibilityValidator().run(
            dataset=_dataset_with_imputed_missing_data_state()
        )

    message = str(exc_info.value)
    assert "Use a non-imputed dataset" in message
    assert "filter features before imputation" in message
    assert "imputation-aware differential policy" in message


def test_eligibility_validator_rejects_established_linear_phospho_scale() -> None:
    with pytest.raises(
        WorkflowValidationError,
        match="requires established log2-scale phospho intensities",
    ):
        DifferentialDatasetEligibilityValidator().run(
            dataset=_dataset_with_linear_scale()
        )


def test_eligibility_validator_rejects_raw_phospho_scale_if_representable() -> None:
    dataset = _dataset_with_log2_scale()
    object.__setattr__(
        dataset,
        "intensity_scale_state",
        IntensityScaleState.raw(has_total_matrix=False),
    )
    with pytest.raises(
        WorkflowValidationError,
        match="requires established log2-scale phospho intensities",
    ):
        DifferentialDatasetEligibilityValidator().run(dataset=dataset)


def test_eligibility_validator_rejects_unknown_phospho_scale_if_representable() -> None:
    dataset = _dataset_with_log2_scale()
    unknown_state = IntensityScaleState.raw(
        has_total_matrix=False
    ).with_quantitative_meaning(QuantitativeMeaning.UNKNOWN)
    object.__setattr__(dataset, "intensity_scale_state", unknown_state)
    with pytest.raises(
        WorkflowValidationError,
        match="requires established log2-scale phospho intensities",
    ):
        DifferentialDatasetEligibilityValidator().run(dataset=dataset)


def test_eligibility_validator_rejects_declared_but_unestablished_log2_scale() -> None:
    dataset = _dataset_with_log2_scale()
    object.__setattr__(
        dataset,
        "intensity_scale_state",
        IntensityScaleState(
            phospho=MatrixIntensityScaleState.log2(established_by="test.declaration"),
            total=None,
        ),
    )
    with pytest.raises(
        WorkflowValidationError,
        match="requires established log2-scale phospho intensities",
    ):
        DifferentialDatasetEligibilityValidator().run(dataset=dataset)
