from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest

from phospy.errors import PhosPyInputError
from phospy.io.readers.importers import _parse_intensity_column
from phospy.io.readers.maxquant import (
    MaxQuantPhosphositeImporter,
    MaxQuantPhosphositeImportRequest,
)
from tests.support.performance_contracts import (
    DATASET_VALIDATION_MEDIUM_N_SAMPLES,
    DATASET_VALIDATION_MEDIUM_N_SITES,
    IMPORTER_MEDIUM_NORMALISATION_RUNTIME_SECONDS_MAX,
    deterministic_maxquant_source_table,
    measure_runtime_and_peak_mib,
)

pytestmark = [pytest.mark.performance, pytest.mark.release_gate]


@pytest.fixture(scope="module")
def medium_maxquant_source() -> pd.DataFrame:
    return deterministic_maxquant_source_table(
        n_sites=DATASET_VALIDATION_MEDIUM_N_SITES,
        n_samples=DATASET_VALIDATION_MEDIUM_N_SAMPLES,
        seed=41_001,
        start=40_000,
    )


def test_maxquant_importer_medium_fixture_completes_under_generous_threshold(
    medium_maxquant_source: pd.DataFrame,
) -> None:
    result, runtime_seconds, _peak_mib = measure_runtime_and_peak_mib(
        lambda: MaxQuantPhosphositeImporter().run(
            MaxQuantPhosphositeImportRequest(source=medium_maxquant_source)
        ),
        warmup=False,
    )

    assert result.phospho_matrix_candidate.shape == (
        DATASET_VALIDATION_MEDIUM_N_SITES,
        DATASET_VALIDATION_MEDIUM_N_SAMPLES,
    )
    assert result.site_metadata_candidate.shape[0] == DATASET_VALIDATION_MEDIUM_N_SITES
    assert result.peptide_evidence is not None
    assert result.peptide_evidence.shape[0] == DATASET_VALIDATION_MEDIUM_N_SITES
    assert result.localisation_confidence_column == "localisation_confidence"
    assert runtime_seconds < IMPORTER_MEDIUM_NORMALISATION_RUNTIME_SECONDS_MAX


def test_optimized_intensity_parser_matches_reference_behavior() -> None:
    valid = pd.Series(["1.5", "", "NA", pd.NA, 2, True], dtype="object")

    observed = _parse_intensity_column(valid, source_column="Intensity A")
    expected = _reference_parse_intensity_column(valid, source_column="Intensity A")

    np.testing.assert_allclose(observed, expected, equal_nan=True)

    for values, expected_message in (
        (pd.Series(["1.0", "bad"], dtype="object"), "offending_value='bad'"),
        (
            pd.Series(["1.0", float("inf")], dtype="object"),
            "reason='not_finite'",
        ),
    ):
        with pytest.raises(PhosPyInputError, match=expected_message):
            _reference_parse_intensity_column(values, source_column="Intensity A")
        with pytest.raises(PhosPyInputError, match=expected_message):
            _parse_intensity_column(values, source_column="Intensity A")


def _reference_parse_intensity_column(
    series: pd.Series,
    *,
    source_column: str,
) -> list[float]:
    parsed: list[float] = []
    missing_tokens = {"", "na", "n/a", "nan", "null"}
    for position, value in enumerate(series.tolist()):
        if bool(pd.Series((value,), dtype="object").isna().iat[0]) or (
            isinstance(value, str) and value.strip().lower() in missing_tokens
        ):
            parsed.append(float("nan"))
            continue
        try:
            numeric_value = float(value)
        except (TypeError, ValueError) as exc:
            raise PhosPyInputError(
                "failed to parse phosphosite import intensity value: "
                f"source_column='{source_column}', row_position={position}, "
                f"offending_value={value!r}"
            ) from exc
        if not math.isfinite(numeric_value):
            raise PhosPyInputError(
                "failed to parse phosphosite import intensity value: "
                f"source_column='{source_column}', row_position={position}, "
                f"offending_value={value!r}, reason='not_finite'"
            )
        parsed.append(numeric_value)
    return parsed
