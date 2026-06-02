from __future__ import annotations

import pandas as pd
import pytest

from phospy.errors.transformations import TransformationStateEstablishmentError
from phospy.errors.validation import TransformationValidationError
from phospy.science.datasets.builders.preprocessing import (
    build_dataset_processing_state,
)
from phospy.science.datasets.builders.transformation_resolver import (
    DatasetIntensityScaleResolver,
)
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.science.datasets.preprocessing.models import PreprocessingPlan
from phospy.science.transformations.models import (
    IntensityScaleEstablishmentMode,
    IntensityScaleEstablishmentSource,
    IntensityScaleState,
    MatrixIntensityScaleState,
)
from phospy.science.transformations.transformers import IdentityTransformer
from tests.support.site_keys import (
    site_key_context_columns,
    site_key_index_from_display_ids,
)


def _display_ids_for_count(count: int) -> list[str]:
    fallback = [
        "MAPK14;Y182;",
        "GSK3B;S9;",
        "AKT1;T308;",
        "PRKACA;S339;",
        "MAPK1;T185;",
    ]
    display_ids = fallback[:count]
    if len(display_ids) < count:
        for idx in range(len(display_ids), count):
            display_ids.append(f"GENE{idx};S{idx + 1};")
    return display_ids


def _phospho(values: list[float]) -> pd.DataFrame:
    display_ids = _display_ids_for_count(len(values))
    return pd.DataFrame(
        {"sample_a": values},
        index=site_key_index_from_display_ids(display_ids),
    )


def _site_metadata(index: pd.Index) -> pd.DataFrame:
    display_ids = _display_ids_for_count(len(index))
    gene_symbols: list[str] = []
    sites: list[str] = []
    for site_id in display_ids:
        parts = site_id.split(";")
        gene_symbols.append(parts[0])
        sites.append(parts[1] if len(parts) > 1 else "S1")
    rows = len(gene_symbols)
    return pd.DataFrame(
        {
            "site_key": index.astype(str).tolist(),
            "display_id": display_ids,
            **site_key_context_columns(index),
            "gene_symbol": gene_symbols,
            "site": sites,
            "site_sequence": [
                ("A" * 15) + str(site).strip().upper()[0] + ("A" * 15) for site in sites
            ],
            "localisation_confidence": [0.95] * rows,
        },
        index=index.copy(),
    )


def _declared_state(kind: str) -> IntensityScaleState:
    if kind == "log2":
        return IntensityScaleState(
            phospho=MatrixIntensityScaleState.log2(established_by="test.declaration"),
            total=None,
        )
    return IntensityScaleState(
        phospho=MatrixIntensityScaleState.linear(established_by="test.declaration"),
        total=None,
    )


def test_declared_log2_records_declared_establishment_mode() -> None:
    resolver = DatasetIntensityScaleResolver(transformer=IdentityTransformer())
    resolved = resolver.run(
        phospho=_phospho([1.2, 2.1, 3.4]),
        total=None,
        declared_input_scale_state=_declared_state("log2"),
        declared_input_establishment_mode=IntensityScaleEstablishmentMode.DECLARED,
        input_declaration_source="dataset_build_request.input_intensity_scale",
    )

    provenance = resolved.intensity_scale_state.establishment_provenance
    assert (
        resolved.intensity_scale_state.establishment_mode
        is IntensityScaleEstablishmentMode.DECLARED
    )
    assert provenance is not None
    assert (
        provenance.input_declaration_source
        == "dataset_build_request.input_intensity_scale"
    )
    assert provenance.mode is IntensityScaleEstablishmentMode.DECLARED
    assert provenance.source is IntensityScaleEstablishmentSource.DECLARED_BY_USER


def test_transformed_log2_records_transformed_establishment_mode() -> None:
    resolver = DatasetIntensityScaleResolver(transformer=IdentityTransformer())
    resolved = resolver.run(
        phospho=_phospho([4.1, 4.4, 5.2]),
        total=None,
        declared_input_scale_state=_declared_state("log2"),
        declared_input_establishment_mode=IntensityScaleEstablishmentMode.TRANSFORMED,
        establishment_transformer_name=(
            "phospy.science.transformations.transformers.log2.Log2Transformer"
        ),
        scale_establishment_parameters={"operation": "log2", "pseudocount": 1.0},
    )

    provenance = resolved.intensity_scale_state.establishment_provenance
    assert (
        resolved.intensity_scale_state.establishment_mode
        is IntensityScaleEstablishmentMode.TRANSFORMED
    )
    assert provenance is not None
    assert (
        provenance.transformer_name
        == "phospy.science.transformations.transformers.log2.Log2Transformer"
    )
    assert provenance.parameters["operation"] == "log2"
    assert provenance.source is IntensityScaleEstablishmentSource.TRANSFORMED_BY_PHOSPY


def test_identity_pass_through_without_declared_scale_fails_establishment() -> None:
    resolver = DatasetIntensityScaleResolver(transformer=IdentityTransformer())
    with pytest.raises(
        TransformationStateEstablishmentError,
        match="pass-through/identity transformer cannot establish scientific input scale",
    ):
        resolver.run(
            phospho=_phospho([100.0, 101.0]),
            total=None,
        )


def test_suspicious_declared_log2_records_warning() -> None:
    resolver = DatasetIntensityScaleResolver(transformer=IdentityTransformer())
    resolved = resolver.run(
        phospho=_phospho([12000.0, 15000.0, 18000.0]),
        total=None,
        declared_input_scale_state=_declared_state("log2"),
        declared_input_establishment_mode=IntensityScaleEstablishmentMode.DECLARED,
        input_declaration_source="dataset_build_request.input_intensity_scale",
    )
    provenance = resolved.intensity_scale_state.establishment_provenance
    assert provenance is not None
    assert any(
        "declared log2 scale is suspicious" in warning
        for warning in provenance.diagnostic_warnings
    )


def test_suspicious_declared_linear_records_warning() -> None:
    resolver = DatasetIntensityScaleResolver(transformer=IdentityTransformer())
    resolved = resolver.run(
        phospho=_phospho([-2.0, 10.0, 15.0]),
        total=None,
        declared_input_scale_state=_declared_state("linear"),
        declared_input_establishment_mode=IntensityScaleEstablishmentMode.DECLARED,
        input_declaration_source="dataset_build_request.input_intensity_scale",
    )
    provenance = resolved.intensity_scale_state.establishment_provenance
    assert provenance is not None
    assert any(
        "declared linear scale contains negative values" in warning
        for warning in provenance.diagnostic_warnings
    )


def test_declared_log2_impossible_range_records_warning() -> None:
    resolver = DatasetIntensityScaleResolver(transformer=IdentityTransformer())
    resolved = resolver.run(
        phospho=_phospho([-45.0, -5.0, 60.0]),
        total=None,
        declared_input_scale_state=_declared_state("log2"),
        declared_input_establishment_mode=IntensityScaleEstablishmentMode.DECLARED,
        input_declaration_source="dataset_build_request.input_intensity_scale",
    )
    provenance = resolved.intensity_scale_state.establishment_provenance
    assert provenance is not None
    assert any(
        "declared log2 scale has highly suspicious range" in warning
        for warning in provenance.diagnostic_warnings
    )


def test_declared_linear_log_ratio_like_values_record_warning() -> None:
    resolver = DatasetIntensityScaleResolver(transformer=IdentityTransformer())
    resolved = resolver.run(
        phospho=_phospho([-1.2, -0.7, 0.3, 1.1, 2.4]),
        total=None,
        declared_input_scale_state=_declared_state("linear"),
        declared_input_establishment_mode=IntensityScaleEstablishmentMode.DECLARED,
        input_declaration_source="dataset_build_request.input_intensity_scale",
    )
    provenance = resolved.intensity_scale_state.establishment_provenance
    assert provenance is not None
    assert any(
        "consistent with log-ratio style data" in warning
        for warning in provenance.diagnostic_warnings
    )
    assert any(
        "strongly resembles already log-transformed data" in warning
        for warning in provenance.diagnostic_warnings
    )


def test_plausible_declared_log2_has_no_warning_and_no_failure() -> None:
    resolver = DatasetIntensityScaleResolver(transformer=IdentityTransformer())
    resolved = resolver.run(
        phospho=_phospho([0.2, 1.1, 4.9, 7.0]),
        total=None,
        declared_input_scale_state=_declared_state("log2"),
        declared_input_establishment_mode=IntensityScaleEstablishmentMode.DECLARED,
        input_declaration_source="dataset_build_request.input_intensity_scale",
    )
    provenance = resolved.intensity_scale_state.establishment_provenance
    assert provenance is not None
    assert provenance.diagnostic_warnings == ()


def test_analysis_ready_dataset_rejects_unestablished_scale_state() -> None:
    phospho = _phospho([1.0])
    state = IntensityScaleState.raw(has_total_matrix=False)
    processing_state = build_dataset_processing_state(
        plan=PreprocessingPlan.default(),
        intensity_scale_state=state,
    )

    with pytest.raises(
        TransformationValidationError,
        match="must be established through a supported PhosPy path",
    ):
        AnalysisReadyPhosphoDataset(
            phospho=phospho,
            site_metadata=_site_metadata(phospho.index),
            intensity_scale_state=state,
            processing_state=processing_state,
        )
