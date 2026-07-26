from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pandas.testing as pdt
import pytest

from phospy.api import (
    Contrast,
    DifferentialAnalysisConfig,
    DifferentialAnalysisRequest,
    DifferentialAnalysisWorkflow,
    EmpiricalBayesConfig,
    ExperimentalDesign,
    Organism,
    SampleDesignRecord,
)
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
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

pytestmark = [pytest.mark.parity]

DIFF_TREND_PARITY_DIR = (
    Path(__file__).resolve().parents[1]
    / "fixtures"
    / "rewrite_parity"
    / "differential_limma_trend_large"
)

LOGFC_RTOL = 1.0e-10
LOGFC_ATOL = 1.0e-10
MODERATED_T_CORRELATION_MIN = 0.995
NEG_LOG10_P_CORRELATION_MIN = 0.98
STANDARD_ERROR_CORRELATION_MIN = 0.94
LOG_PRIOR_VARIANCE_CORRELATION_MIN = 0.90
P_VALUE_MEDIAN_ABS_DIFF_MAX = 0.015
P_VALUE_P99_ABS_DIFF_MAX = 0.08
ADJ_P_VALUE_P99_ABS_DIFF_MAX = 0.15
STANDARD_ERROR_P99_ABS_DIFF_MAX = 0.035


def _load_manifest() -> dict[str, Any]:
    return json.loads((DIFF_TREND_PARITY_DIR / "MANIFEST.json").read_text())


def _load_matrix() -> pd.DataFrame:
    return pd.read_csv(DIFF_TREND_PARITY_DIR / "matrix.csv").set_index("site_id")


def _load_design() -> pd.DataFrame:
    return pd.read_csv(DIFF_TREND_PARITY_DIR / "design.csv").set_index("sample")


def _load_expected(result_site_index: pd.Index) -> pd.DataFrame:
    frame = pd.read_csv(DIFF_TREND_PARITY_DIR / "limma_B_vs_A.csv").set_index("site_id")
    frame.index = site_key_index_from_display_ids(
        frame.index.astype(str).tolist(),
        protein_namespace="gene_symbol",
    )
    return frame.reindex(result_site_index)


def _load_simulation_diagnostics(result_site_index: pd.Index) -> pd.DataFrame:
    frame = pd.read_csv(DIFF_TREND_PARITY_DIR / "simulation_diagnostics.csv").set_index(
        "site_id"
    )
    frame.index = site_key_index_from_display_ids(
        frame.index.astype(str).tolist(),
        protein_namespace="gene_symbol",
    )
    return frame.reindex(result_site_index)


def _dataset_from_matrix(matrix: pd.DataFrame) -> AnalysisReadyPhosphoDataset:
    display_ids = matrix.index.astype(str).tolist()
    site_index = site_key_index_from_display_ids(
        display_ids,
        protein_namespace="gene_symbol",
    )
    parsed = [site_id.split(";") for site_id in display_ids]
    phospho = matrix.copy(deep=True)
    phospho.index = site_index
    site_metadata = pd.DataFrame(
        {
            "site_key": site_index.tolist(),
            "display_id": display_ids,
            **site_key_context_columns(site_index),
            "gene_symbol": [parts[0] for parts in parsed],
            "site": [parts[1] for parts in parsed],
            "site_sequence": [
                ("A" * 15) + str(parts[1]).strip().upper()[0] + ("A" * 15)
                for parts in parsed
            ],
            "localisation_confidence": [0.95] * matrix.shape[0],
            "protein_id": [parts[0] for parts in parsed],
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


def _design_from_matrix(design: pd.DataFrame) -> ExperimentalDesign:
    records: list[SampleDesignRecord] = []
    replicate_counts: dict[str, int] = {}
    for sample_id, row in design.iterrows():
        active_terms = [term for term, value in row.items() if float(value) == 1.0]
        if len(active_terms) != 1:
            raise AssertionError(
                "large trend parity fixture design must be one-hot encoded per sample"
            )
        condition = str(active_terms[0])
        replicate_counts.setdefault(condition, 0)
        replicate_counts[condition] += 1
        records.append(
            SampleDesignRecord(
                sample_id=str(sample_id),
                condition=condition,
                biological_replicate_id=f"{condition}_r{replicate_counts[condition]}",
            )
        )
    return ExperimentalDesign(samples=tuple(records))


def _run_fixture_workflow() -> object:
    matrix = _load_matrix()
    request = DifferentialAnalysisRequest(
        dataset=_dataset_from_matrix(matrix),
        design=_design_from_matrix(_load_design()),
        contrasts=(
            Contrast(
                name="B_vs_A",
                numerator_condition="B",
                denominator_condition="A",
            ),
        ),
        config=DifferentialAnalysisConfig(
            empirical_bayes=EmpiricalBayesConfig(method="standard", trend=True)
        ),
    )
    return DifferentialAnalysisWorkflow().run(request)


def _correlation(left: pd.Series | np.ndarray, right: pd.Series | np.ndarray) -> float:
    left_values = np.asarray(left, dtype=float)
    right_values = np.asarray(right, dtype=float)
    finite = np.isfinite(left_values) & np.isfinite(right_values)
    return float(np.corrcoef(left_values[finite], right_values[finite])[0, 1])


def _absolute_diff(
    left: pd.Series | np.ndarray,
    right: pd.Series | np.ndarray,
) -> np.ndarray:
    return np.abs(np.asarray(left, dtype=float) - np.asarray(right, dtype=float))


def test_large_limma_trend_fixture_manifest_hashes_match_checked_in_files() -> None:
    manifest = _load_manifest()

    assert manifest["classification"] == "external_parity"
    assert manifest["external_implementation"]["name"] == "R limma"
    assert manifest["external_implementation"]["limma_version"]
    assert manifest["seed"] == 20260724
    assert manifest["design"]["n_features"] >= 1500
    assert manifest["design"]["sample_counts"] == {"A": 5, "B": 7}

    for file_entry in manifest["files"]:
        path = DIFF_TREND_PARITY_DIR / str(file_entry["relative_path"])
        assert path.is_file(), f"missing fixture file: {path.name}"
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        assert actual == file_entry["sha256"]


def test_large_limma_trend_fixture_records_source_policy_and_versions() -> None:
    provenance = (DIFF_TREND_PARITY_DIR / "PROVENANCE.md").read_text(encoding="utf-8")

    assert "Generated with R version" in provenance
    assert "limma version:" in provenance
    assert "Source policy: deterministic synthetic fixture" in provenance
    assert "Classification: external parity for limma result columns" in provenance
    assert "Design: ~0 + condition with groups A/B and unbalanced 5/7" in provenance
    assert "Contrast: B_vs_A = B - A" in provenance


def test_large_feature_trend_path_matches_limma_coefficients_and_drift_envelope() -> (
    None
):
    result = _run_fixture_workflow()
    observed = result.table_for("B_vs_A")
    expected = _load_expected(observed.index)
    simulation = _load_simulation_diagnostics(observed.index)

    assert observed.shape[0] == 1600
    assert result.empirical_bayes_trend is True
    assert result.mean_variance_trend_diagnostics is not None
    assert result.mean_variance_trend_diagnostics.mean_intensity.shape == (1600,)

    pdt.assert_series_equal(
        observed.loc[:, "logFC"],
        expected.loc[:, "logFC"],
        check_names=False,
        check_exact=False,
        rtol=LOGFC_RTOL,
        atol=LOGFC_ATOL,
    )

    # PhosPy and limma use different trend smoothers. These thresholds lock broad
    # moderated-statistic agreement without pretending that limma's lowess prior
    # and PhosPy's deterministic anchor trend are numerically identical.
    observed_se = np.abs(
        observed.loc[:, "logFC"].to_numpy(dtype=float)
        / observed.loc[:, "t"].to_numpy(dtype=float)
    )
    observed_neg_log10_p = -np.log10(
        np.clip(observed.loc[:, "P.Value"].to_numpy(dtype=float), 1.0e-300, 1.0)
    )
    expected_neg_log10_p = -np.log10(
        np.clip(expected.loc[:, "P.Value"].to_numpy(dtype=float), 1.0e-300, 1.0)
    )
    observed_prior = result.prior_residual_variance_series()
    expected_prior = expected.loc[:, "s2.prior"]

    assert _correlation(observed.loc[:, "t"], expected.loc[:, "t"]) > (
        MODERATED_T_CORRELATION_MIN
    )
    assert _correlation(observed_neg_log10_p, expected_neg_log10_p) > (
        NEG_LOG10_P_CORRELATION_MIN
    )
    assert _correlation(observed_se, expected.loc[:, "SE"]) > (
        STANDARD_ERROR_CORRELATION_MIN
    )
    assert _correlation(np.log(observed_prior), np.log(expected_prior)) > (
        LOG_PRIOR_VARIANCE_CORRELATION_MIN
    )
    np.testing.assert_array_equal(
        np.sign(observed.loc[:, "t"].to_numpy(dtype=float)),
        np.sign(expected.loc[:, "t"].to_numpy(dtype=float)),
    )
    assert (
        float(
            np.nanmedian(
                _absolute_diff(observed.loc[:, "P.Value"], expected.loc[:, "P.Value"])
            )
        )
        < P_VALUE_MEDIAN_ABS_DIFF_MAX
    )
    assert (
        float(
            np.nanquantile(
                _absolute_diff(observed.loc[:, "P.Value"], expected.loc[:, "P.Value"]),
                0.99,
            )
        )
        < P_VALUE_P99_ABS_DIFF_MAX
    )
    assert (
        float(
            np.nanquantile(
                _absolute_diff(
                    observed.loc[:, "adj.P.Val"], expected.loc[:, "adj.P.Val"]
                ),
                0.99,
            )
        )
        < ADJ_P_VALUE_P99_ABS_DIFF_MAX
    )
    assert (
        float(np.nanquantile(_absolute_diff(observed_se, expected.loc[:, "SE"]), 0.99))
        < STANDARD_ERROR_P99_ABS_DIFF_MAX
    )

    shifted = simulation.loc[:, "is_shifted"].astype(bool).to_numpy()
    observed_p = observed.loc[:, "P.Value"].to_numpy(dtype=float)
    expected_p = expected.loc[:, "P.Value"].to_numpy(dtype=float)
    assert float(np.nanmedian(observed_p[shifted])) < float(
        np.nanmedian(observed_p[~shifted])
    )
    assert float(np.nanmedian(expected_p[shifted])) < float(
        np.nanmedian(expected_p[~shifted])
    )
