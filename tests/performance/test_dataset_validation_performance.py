from __future__ import annotations

import math
from typing import cast

import numpy as np
import pandas as pd
import pytest

from phospy.errors.validation import DatasetValidationError
from phospy.science.datasets.models import (
    AnalysisReadyPhosphoDataset,
    ImputationObservationMetadata,
    _analysis_ready_matrix_missing_value_count,
    _require_boolean_observation_mask,
)
from phospy.science.references.models import Organism
from phospy.tables.datasets import SiteMetadataTable
from phospy.validation.datasets.site_metadata import (
    assess_localisation_probability_column,
)
from phospy.validation.identity_contracts import enforce_site_key_matches_metadata
from tests.support.analysis_ready_dataset_factories import (
    trusted_analysis_ready_dataset_from_tables,
)
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
)
from tests.support.performance_contracts import (
    DATASET_VALIDATION_BOUNDED_ALIGNMENT_RUNTIME_SECONDS_MAX,
    DATASET_VALIDATION_BOUNDED_OBSERVATION_MASK_RUNTIME_SECONDS_MAX,
    DATASET_VALIDATION_BOUNDED_SITE_METADATA_RUNTIME_SECONDS_MAX,
    DATASET_VALIDATION_MEDIUM_CONSTRUCTION_RUNTIME_SECONDS_MAX,
    DATASET_VALIDATION_MEDIUM_N_SAMPLES,
    DATASET_VALIDATION_MEDIUM_N_SITES,
    DATASET_VALIDATION_SMALL_N_SAMPLES,
    DATASET_VALIDATION_SMALL_N_SITES,
    deterministic_analysis_ready_dataset_tables,
    measure_runtime_and_peak_mib,
    median_runtime_seconds,
)

pytestmark = [pytest.mark.performance, pytest.mark.release_gate]


@pytest.fixture(scope="module")
def small_dataset_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    return deterministic_analysis_ready_dataset_tables(
        n_sites=DATASET_VALIDATION_SMALL_N_SITES,
        n_samples=DATASET_VALIDATION_SMALL_N_SAMPLES,
        seed=31_001,
        start=10_000,
        gene_prefix="SMALLPERF",
    )


@pytest.fixture(scope="module")
def medium_dataset_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    return deterministic_analysis_ready_dataset_tables(
        n_sites=DATASET_VALIDATION_MEDIUM_N_SITES,
        n_samples=DATASET_VALIDATION_MEDIUM_N_SAMPLES,
        seed=31_002,
        start=20_000,
        gene_prefix="MEDPERF",
    )


def test_realistic_performance_fixtures_have_expected_shapes(
    small_dataset_tables: tuple[pd.DataFrame, pd.DataFrame],
    medium_dataset_tables: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    small_phospho, small_metadata = small_dataset_tables
    medium_phospho, medium_metadata = medium_dataset_tables

    assert small_phospho.shape == (
        DATASET_VALIDATION_SMALL_N_SITES,
        DATASET_VALIDATION_SMALL_N_SAMPLES,
    )
    assert medium_phospho.shape == (
        DATASET_VALIDATION_MEDIUM_N_SITES,
        DATASET_VALIDATION_MEDIUM_N_SAMPLES,
    )
    assert small_metadata.index.equals(small_phospho.index)
    assert medium_metadata.index.equals(medium_phospho.index)
    assert float(small_phospho.min().min()) >= 0.0
    assert float(medium_phospho.min().min()) >= 0.0


def test_medium_dataset_construction_completes_under_generous_threshold(
    medium_dataset_tables: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    phospho, site_metadata = medium_dataset_tables

    dataset, runtime_seconds, _peak_mib = measure_runtime_and_peak_mib(
        lambda: trusted_analysis_ready_dataset_from_tables(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
            intensity_scale_state=supported_linear_intensity_scale_state(
                has_total_matrix=False
            ),
            processing_state=supported_linear_processing_state(has_total_matrix=False),
        ),
        warmup=False,
    )
    dataset = cast(AnalysisReadyPhosphoDataset, dataset)

    assert dataset.phospho.shape == (
        DATASET_VALIDATION_MEDIUM_N_SITES,
        DATASET_VALIDATION_MEDIUM_N_SAMPLES,
    )
    assert dataset.site_metadata.index.equals(dataset.phospho.index)
    assert runtime_seconds < DATASET_VALIDATION_MEDIUM_CONSTRUCTION_RUNTIME_SECONDS_MAX


def test_analysis_ready_missing_value_scan_scales_for_numeric_matrix(
    medium_dataset_tables: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    phospho, _site_metadata = medium_dataset_tables
    object_phospho = phospho.astype(object)

    assert _analysis_ready_matrix_missing_value_count(phospho) == 0
    assert _analysis_ready_matrix_missing_value_count(object_phospho) == 0

    numeric_runtime_seconds = median_runtime_seconds(
        lambda: _analysis_ready_matrix_missing_value_count(phospho),
        repeats=5,
        warmup=True,
    )
    object_runtime_seconds = median_runtime_seconds(
        lambda: _analysis_ready_matrix_missing_value_count(object_phospho),
        repeats=3,
        warmup=True,
    )

    assert numeric_runtime_seconds < object_runtime_seconds


def test_bounded_site_metadata_validation_completes_under_generous_threshold(
    medium_dataset_tables: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    phospho, site_metadata = medium_dataset_tables

    table, runtime_seconds, _peak_mib = measure_runtime_and_peak_mib(
        lambda: SiteMetadataTable(
            frame=site_metadata,
            expected_index=phospho.index,
        ),
        warmup=False,
    )
    table = cast(SiteMetadataTable, table)

    assert table.frame.shape[0] == DATASET_VALIDATION_MEDIUM_N_SITES
    assert table.frame.index.equals(phospho.index)
    assert (
        runtime_seconds < DATASET_VALIDATION_BOUNDED_SITE_METADATA_RUNTIME_SECONDS_MAX
    )


def test_bounded_site_key_metadata_alignment_completes_under_generous_threshold(
    medium_dataset_tables: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    _phospho, site_metadata = medium_dataset_tables

    _result, runtime_seconds, _peak_mib = measure_runtime_and_peak_mib(
        lambda: enforce_site_key_matches_metadata(
            site_metadata=site_metadata,
            field_name="performance.site_metadata",
            error_type=DatasetValidationError,
        ),
        warmup=False,
    )

    assert runtime_seconds < DATASET_VALIDATION_BOUNDED_ALIGNMENT_RUNTIME_SECONDS_MAX


def test_bounded_observation_mask_validation_completes_under_generous_threshold(
    medium_dataset_tables: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    phospho, _site_metadata = medium_dataset_tables
    values = np.ones(phospho.shape, dtype=bool)
    values[::11, 3::7] = False
    observed_mask = pd.DataFrame(
        values,
        index=phospho.index.copy(),
        columns=phospho.columns.copy(),
    )

    metadata, runtime_seconds, _peak_mib = measure_runtime_and_peak_mib(
        lambda: ImputationObservationMetadata(
            observed_mask=observed_mask,
            phospho_index=phospho.index,
            sample_index=phospho.columns,
        ),
        warmup=False,
    )
    metadata = cast(ImputationObservationMetadata, metadata)

    assert metadata.observed_mask.shape == phospho.shape
    assert (
        runtime_seconds
        < DATASET_VALIDATION_BOUNDED_OBSERVATION_MASK_RUNTIME_SECONDS_MAX
    )


def test_optimized_observation_mask_validation_matches_reference_behavior() -> None:
    valid = pd.DataFrame(
        [[True, False], [np.bool_(True), False]],
        index=pd.Index(["site_a", "site_b"]),
        columns=pd.Index(["sample_a", "sample_b"]),
        dtype=object,
    )
    missing = valid.copy(deep=True)
    missing.iloc[0, 0] = pd.NA
    invalid = valid.copy(deep=True)
    invalid.iloc[1, 1] = "true"

    for mask in (valid,):
        _reference_require_boolean_observation_mask(mask)
        _require_boolean_observation_mask(mask)

    for mask, expected in (
        (missing, "missing values"),
        (invalid, "only boolean values"),
    ):
        with pytest.raises(DatasetValidationError, match=expected):
            _reference_require_boolean_observation_mask(mask)
        with pytest.raises(DatasetValidationError, match=expected):
            _require_boolean_observation_mask(mask)


def test_optimized_localisation_probability_assessment_matches_reference_behavior() -> (
    None
):
    site_ids = pd.Index(
        [
            "valid_float",
            "valid_string",
            "blank_string",
            "missing",
            "python_bool",
            "numpy_bool",
            "not_numeric",
            "out_of_range",
            "not_finite",
            "unsupported_type",
        ],
        name="site_key",
    )
    site_metadata = pd.DataFrame(
        {
            "localisation_probability": pd.Series(
                [
                    0.25,
                    "0.95",
                    " ",
                    pd.NA,
                    True,
                    np.bool_(False),
                    "bad",
                    1.1,
                    float("inf"),
                    object(),
                ],
                index=site_ids,
                dtype="object",
            )
        },
        index=site_ids,
    )

    observed = assess_localisation_probability_column(
        site_metadata=site_metadata,
        field_name="performance.site_metadata",
        error_type=DatasetValidationError,
    )
    assert observed is not None
    expected = _reference_assess_localisation_probability_column(site_metadata)

    pd.testing.assert_series_equal(observed.normalized, expected[0])
    pd.testing.assert_series_equal(observed.missing_mask, expected[1])
    pd.testing.assert_series_equal(observed.invalid_mask, expected[2])
    assert observed.invalid_examples == expected[3]


def _reference_require_boolean_observation_mask(mask: pd.DataFrame) -> None:
    values = mask.to_numpy(dtype="object")
    for row_index in range(values.shape[0]):
        for column_index in range(values.shape[1]):
            value = values[row_index, column_index]
            if bool(pd.Series((value,), dtype="object").isna().iat[0]):
                raise DatasetValidationError(
                    "dataset.imputation_observation_mask must not contain "
                    "missing values"
                )
            if isinstance(value, (bool, np.bool_)):
                continue
            raise DatasetValidationError(
                "dataset.imputation_observation_mask must contain only boolean "
                "values; "
                f"invalid_cell=({mask.index[row_index]!r}, "
                f"{mask.columns[column_index]!r})"
            )


def _reference_assess_localisation_probability_column(
    site_metadata: pd.DataFrame,
) -> tuple[pd.Series, pd.Series, pd.Series, tuple[str, ...]]:
    values = site_metadata["localisation_probability"]
    values_index = pd.Index(values.index)
    normalized = pd.Series(pd.NA, index=values_index, dtype="Float64")
    missing_mask = pd.Series(False, index=values_index, dtype="boolean")
    invalid_mask = pd.Series(False, index=values_index, dtype="boolean")
    invalid_examples: list[str] = []

    for site_id, raw_value in values.items():
        parsed = _reference_parse_localisation_probability(raw_value)
        if parsed is None:
            missing_mask.at[site_id] = True
            continue
        if isinstance(parsed, float):
            normalized.at[site_id] = parsed
            continue
        invalid_mask.at[site_id] = True
        if len(invalid_examples) < 5:
            invalid_examples.append(f"{site_id!r}:{raw_value!r}:{parsed}")

    return normalized, missing_mask, invalid_mask, tuple(invalid_examples)


def _reference_parse_localisation_probability(value: object) -> float | str | None:
    if bool(pd.Series((value,), dtype="object").isna().iat[0]):
        return None
    if isinstance(value, bool):
        return "bool_not_allowed"
    if isinstance(value, str):
        token = value.strip()
        if token == "":
            return None
        try:
            numeric_value = float(token)
        except ValueError:
            return "not_numeric"
    elif isinstance(value, (int, float)):
        numeric_value = float(value)
    else:
        return "unsupported_type"
    if not math.isfinite(numeric_value):
        return "not_finite"
    if numeric_value < 0.0 or numeric_value > 1.0:
        return "out_of_range"
    return float(numeric_value)
