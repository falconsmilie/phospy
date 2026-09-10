from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, cast

import numpy as np
import pandas as pd
import pandas.testing as pdt
import pytest

from phospy.advanced import (
    DifferentialAnalysisConfig,
    EmpiricalBayesConfig,
)
from phospy.api import (
    Contrast,
    DifferentialAnalysisRequest,
    DifferentialAnalysisWorkflow,
    ExperimentalDesign,
    Organism,
    SampleDesignRecord,
)
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.science.differential.models import (
    EMPIRICAL_BAYES_TREND_COVARIATE_QUANTIFICATION_DEPTH,
    QUANTIFICATION_DEPTH_KIND_PSM_COUNT,
    DifferentialAnalysisResult,
)
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

DEQMS_DEPTH_PARITY_DIR = (
    Path(__file__).resolve().parents[1]
    / "fixtures"
    / "rewrite_parity"
    / "differential_deqms_depth"
)

LOGFC_RTOL = 1.0e-10
LOGFC_ATOL = 1.0e-10
PRIOR_DF_ATOL = 5.0e-2
LOG_PRIOR_VARIANCE_CORRELATION_MIN = 0.998
LOG_PRIOR_VARIANCE_MEDIAN_ABS_DIFF_MAX = 3.0e-2
LOG_PRIOR_VARIANCE_P99_ABS_DIFF_MAX = 7.0e-2
POSTERIOR_VARIANCE_CORRELATION_MIN = 0.999
LOG_POSTERIOR_VARIANCE_P99_ABS_DIFF_MAX = 6.0e-2
MODERATED_T_CORRELATION_MIN = 0.9998
MODERATED_T_P99_ABS_DIFF_MAX = 5.0e-2
P_VALUE_CORRELATION_MIN = 0.9997
P_VALUE_MEDIAN_ABS_DIFF_MAX = 4.0e-3
P_VALUE_P99_ABS_DIFF_MAX = 1.4e-2
ADJ_P_VALUE_CORRELATION_MIN = 0.995
ADJ_P_VALUE_P99_ABS_DIFF_MAX = 6.0e-2
DEPTH_GROUP_LOG_PRIOR_VARIANCE_CORRELATION_MIN = 0.994
DEPTH_GROUP_LOG_PRIOR_VARIANCE_P99_ABS_DIFF_MAX = 7.0e-2
DEPTH_GROUP_POSTERIOR_VARIANCE_CORRELATION_MIN = 0.997
DEPTH_GROUP_LOG_POSTERIOR_VARIANCE_P99_ABS_DIFF_MAX = 6.0e-2
DEPTH_GROUP_MODERATED_T_CORRELATION_MIN = 0.9998
DEPTH_GROUP_P_VALUE_CORRELATION_MIN = 0.9997
DEPTH_GROUP_ADJ_P_VALUE_CORRELATION_MIN = 0.988
DEPTH_GROUP_ADJ_P_VALUE_P99_ABS_DIFF_MAX = 6.0e-2


@pytest.fixture(scope="module")
def workflow_result() -> DifferentialAnalysisResult:
    return _run_fixture_workflow()


def _load_manifest() -> dict[str, Any]:
    return json.loads((DEQMS_DEPTH_PARITY_DIR / "MANIFEST.json").read_text())


def _load_matrix() -> pd.DataFrame:
    return pd.read_csv(DEQMS_DEPTH_PARITY_DIR / "matrix.csv").set_index("site_id")


def _load_design() -> pd.DataFrame:
    return pd.read_csv(DEQMS_DEPTH_PARITY_DIR / "design.csv").set_index("sample")


def _load_depth() -> pd.DataFrame:
    return pd.read_csv(DEQMS_DEPTH_PARITY_DIR / "quantification_depth.csv").set_index(
        "site_id"
    )


def _load_feature_metadata(result_site_index: pd.Index) -> pd.DataFrame:
    frame = pd.read_csv(
        DEQMS_DEPTH_PARITY_DIR / "feature_metadata.csv",
        keep_default_na=False,
    ).set_index("site_id", drop=False)
    frame.index = site_key_index_from_display_ids(
        frame.index.astype(str).tolist(),
        protein_namespace="gene_symbol",
    )
    return frame.reindex(result_site_index)


def _load_expected(result_site_index: pd.Index) -> pd.DataFrame:
    frame = pd.read_csv(DEQMS_DEPTH_PARITY_DIR / "deqms_B_vs_A.csv").set_index(
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
    depth = _load_depth().loc[display_ids, "quantification_depth"].to_numpy(dtype=float)
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
            "quantification_depth": depth,
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
                "DEqMS depth fixture design must be one-hot encoded per sample"
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


def _run_fixture_workflow() -> DifferentialAnalysisResult:
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
            empirical_bayes=EmpiricalBayesConfig(
                method="standard",
                trend=True,
                trend_covariate=(EMPIRICAL_BAYES_TREND_COVARIATE_QUANTIFICATION_DEPTH),
                quantification_depth_kind=QUANTIFICATION_DEPTH_KIND_PSM_COUNT,
            )
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


def _log_absolute_diff(
    left: pd.Series | np.ndarray,
    right: pd.Series | np.ndarray,
) -> np.ndarray:
    return np.abs(np.log(np.asarray(left, dtype=float)) - np.log(np.asarray(right)))


def test_deqms_depth_fixture_manifest_hashes_match_checked_in_files() -> None:
    manifest = _load_manifest()

    assert manifest["manifest_schema_version"] == "fixture-manifest-v1"
    assert manifest["fixture_family"] == "differential_deqms_depth"
    assert manifest["classification"] == "external_parity"
    assert manifest["external_implementation"] == {
        "name": "R DEqMS",
        "r_version": "R version 4.5.2 (2025-10-31 ucrt)",
        "bioconductor_version": "3.22",
        "limma_version": "3.66.0",
        "deqms_version": "1.28.0",
    }
    assert manifest["pinned_environment"] == {
        "r_version": "R version 4.5.2 (2025-10-31 ucrt)",
        "bioconductor_version": "3.22",
        "limma_version": "3.66.0",
        "deqms_version": "1.28.0",
    }
    assert (
        manifest["generator"]
        == "tests/fixtures/rewrite_parity/differential_deqms_depth/generate_fixture.R"
    )
    assert (
        manifest["generator_sha256"]
        == hashlib.sha256(
            (DEQMS_DEPTH_PARITY_DIR / "generate_fixture.R").read_bytes()
        ).hexdigest()
    )
    assert manifest["seed"] == 20260909
    assert manifest["generation_timestamp_utc"] == "2026-09-09T00:00:00Z"
    assert manifest["byte_policy"] == "utf-8 LF with final newline"
    assert (
        "DEqMS and limma outputs are the external scientific authority"
        in (manifest["source_policy"])
    )
    assert "not generated by PhosPy" in manifest["numeric_authority"]
    assert manifest["model_policy"]["deqms_function"] == "DEqMS::spectraCounteBayes"
    assert manifest["model_policy"]["deqms_fit_method"] == "loess"
    assert manifest["model_policy"]["deqms_loess_span"] == 0.75
    assert manifest["model_policy"]["trend_covariate"] == "log2(fit$count)"
    assert manifest["model_policy"]["quantification_depth_kind"] == "psm_count"
    assert manifest["design"]["n_features"] == 144
    assert manifest["design"]["sample_counts"] == {"A": 5, "B": 7}
    assert manifest["design"]["lower_depth_feature_count"] > 40
    assert manifest["design"]["higher_depth_feature_count"] > 40

    declared_paths: set[str] = set()
    for file_entry in manifest["files"]:
        relative_path = str(file_entry["relative_path"])
        assert relative_path not in declared_paths
        declared_paths.add(relative_path)
        path = DEQMS_DEPTH_PARITY_DIR / relative_path
        assert path.is_file(), f"missing fixture file: {relative_path}"
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        assert actual == file_entry["sha256"]


def test_deqms_depth_fixture_records_source_policy_versions_and_terminology() -> None:
    provenance = (DEQMS_DEPTH_PARITY_DIR / "PROVENANCE.md").read_text(encoding="utf-8")
    manifest = _load_manifest()

    assert "Generated with R version" in provenance
    assert "Bioconductor version: 3.22" in provenance
    assert "limma version: 3.66.0" in provenance
    assert "DEqMS version: 1.28.0" in provenance
    assert "DEqMS::spectraCounteBayes" in provenance
    assert "fit.method='loess'" in provenance
    assert "Quantification depth: feature-level synthetic PSM count" in provenance
    assert "PhosPy is not imported or executed by this generator" in provenance
    assert (
        "quantification-depth-aware empirical Bayes moderation inspired by DEqMS"
        in (provenance)
    )
    assert "not exact DEqMS-compatible numerical equivalence" in provenance
    assert manifest["terminology_conclusion"] == (
        "quantification-depth-aware empirical Bayes moderation inspired by DEqMS; "
        "no exact DEqMS-compatible numerical equivalence claim"
    )


def test_depth_aware_moderation_matches_deqms_count_contract_with_tolerances(
    workflow_result: DifferentialAnalysisResult,
) -> None:
    result = workflow_result
    observed = result.table_for("B_vs_A")
    expected = _load_expected(observed.index)

    assert observed.shape[0] == 144
    assert result.empirical_bayes_trend is True
    assert result.mean_variance_trend_diagnostics is None
    assert result.quantification_depth_trend_diagnostics is not None
    assert result.policy_provenance is not None
    assert result.policy_provenance.empirical_bayes.trend_covariate == (
        "quantification_depth"
    )
    assert (
        result.policy_provenance.empirical_bayes.trend_covariate_transformation
        == "log2"
    )
    assert result.policy_provenance.empirical_bayes.quantification_depth_kind == (
        "psm_count"
    )

    pdt.assert_series_equal(
        observed.loc[:, "logFC"],
        expected.loc[:, "logFC"],
        check_names=False,
        check_exact=False,
        rtol=LOGFC_RTOL,
        atol=LOGFC_ATOL,
    )

    # DEqMS scans prior degrees of freedom on a 0.1 grid; PhosPy uses the
    # continuous limma-moment helper already used by its mean-intensity trend.
    np.testing.assert_allclose(
        result.prior_degrees_of_freedom_series().to_numpy(dtype=float),
        expected.loc[:, "sca.dfprior"].to_numpy(dtype=float),
        rtol=0.0,
        atol=PRIOR_DF_ATOL,
    )
    assert (
        _correlation(
            np.log(result.prior_residual_variance_series()),
            np.log(expected.loc[:, "sca.priorvar"]),
        )
        > LOG_PRIOR_VARIANCE_CORRELATION_MIN
    )
    assert (
        float(
            np.nanmedian(
                _log_absolute_diff(
                    result.prior_residual_variance_series(),
                    expected.loc[:, "sca.priorvar"],
                )
            )
        )
        < LOG_PRIOR_VARIANCE_MEDIAN_ABS_DIFF_MAX
    )
    assert (
        float(
            np.nanquantile(
                _log_absolute_diff(
                    result.prior_residual_variance_series(),
                    expected.loc[:, "sca.priorvar"],
                ),
                0.99,
            )
        )
        < LOG_PRIOR_VARIANCE_P99_ABS_DIFF_MAX
    )

    assert (
        _correlation(
            result.posterior_residual_variance_series(),
            expected.loc[:, "sca.postvar"],
        )
        > POSTERIOR_VARIANCE_CORRELATION_MIN
    )
    assert (
        float(
            np.nanquantile(
                _log_absolute_diff(
                    result.posterior_residual_variance_series(),
                    expected.loc[:, "sca.postvar"],
                ),
                0.99,
            )
        )
        < LOG_POSTERIOR_VARIANCE_P99_ABS_DIFF_MAX
    )
    assert _correlation(observed.loc[:, "t"], expected.loc[:, "sca.t"]) > (
        MODERATED_T_CORRELATION_MIN
    )
    assert (
        float(
            np.nanquantile(
                _absolute_diff(observed.loc[:, "t"], expected.loc[:, "sca.t"]), 0.99
            )
        )
        < MODERATED_T_P99_ABS_DIFF_MAX
    )
    assert _correlation(observed.loc[:, "P.Value"], expected.loc[:, "sca.P.Value"]) > (
        P_VALUE_CORRELATION_MIN
    )
    assert (
        float(
            np.nanmedian(
                _absolute_diff(
                    observed.loc[:, "P.Value"], expected.loc[:, "sca.P.Value"]
                )
            )
        )
        < P_VALUE_MEDIAN_ABS_DIFF_MAX
    )
    assert (
        float(
            np.nanquantile(
                _absolute_diff(
                    observed.loc[:, "P.Value"], expected.loc[:, "sca.P.Value"]
                ),
                0.99,
            )
        )
        < P_VALUE_P99_ABS_DIFF_MAX
    )
    assert (
        _correlation(
            observed.loc[:, "adj.P.Val"],
            expected.loc[:, "sca.adj.P.Val"],
        )
        > ADJ_P_VALUE_CORRELATION_MIN
    )
    assert (
        float(
            np.nanquantile(
                _absolute_diff(
                    observed.loc[:, "adj.P.Val"],
                    expected.loc[:, "sca.adj.P.Val"],
                ),
                0.99,
            )
        )
        < ADJ_P_VALUE_P99_ABS_DIFF_MAX
    )


@pytest.mark.parametrize("depth_group", ("lower_depth", "higher_depth"))
def test_lower_and_higher_depth_sites_match_deqms_within_scientific_envelope(
    workflow_result: DifferentialAnalysisResult,
    depth_group: str,
) -> None:
    result = workflow_result
    observed = result.table_for("B_vs_A")
    expected = _load_expected(observed.index)
    metadata = _load_feature_metadata(observed.index)
    group_mask = metadata.loc[:, "depth_group"].eq(depth_group).to_numpy()

    assert int(group_mask.sum()) > 40
    assert (
        _correlation(
            np.log(
                result.prior_residual_variance_series().to_numpy(dtype=float)[
                    group_mask
                ]
            ),
            np.log(expected.loc[:, "sca.priorvar"].to_numpy(dtype=float)[group_mask]),
        )
        > DEPTH_GROUP_LOG_PRIOR_VARIANCE_CORRELATION_MIN
    )
    assert (
        float(
            np.nanquantile(
                _log_absolute_diff(
                    result.prior_residual_variance_series().to_numpy(dtype=float)[
                        group_mask
                    ],
                    expected.loc[:, "sca.priorvar"].to_numpy(dtype=float)[group_mask],
                ),
                0.99,
            )
        )
        < DEPTH_GROUP_LOG_PRIOR_VARIANCE_P99_ABS_DIFF_MAX
    )
    assert (
        _correlation(
            observed.loc[:, "t"].to_numpy(dtype=float)[group_mask],
            expected.loc[:, "sca.t"].to_numpy(dtype=float)[group_mask],
        )
        > DEPTH_GROUP_MODERATED_T_CORRELATION_MIN
    )
    assert (
        _correlation(
            observed.loc[:, "P.Value"].to_numpy(dtype=float)[group_mask],
            expected.loc[:, "sca.P.Value"].to_numpy(dtype=float)[group_mask],
        )
        > DEPTH_GROUP_P_VALUE_CORRELATION_MIN
    )
    assert (
        _correlation(
            result.posterior_residual_variance_series().to_numpy(dtype=float)[
                group_mask
            ],
            expected.loc[:, "sca.postvar"].to_numpy(dtype=float)[group_mask],
        )
        > DEPTH_GROUP_POSTERIOR_VARIANCE_CORRELATION_MIN
    )
    assert (
        float(
            np.nanquantile(
                _log_absolute_diff(
                    result.posterior_residual_variance_series().to_numpy(dtype=float)[
                        group_mask
                    ],
                    expected.loc[:, "sca.postvar"].to_numpy(dtype=float)[group_mask],
                ),
                0.99,
            )
        )
        < DEPTH_GROUP_LOG_POSTERIOR_VARIANCE_P99_ABS_DIFF_MAX
    )
    assert (
        _correlation(
            observed.loc[:, "adj.P.Val"].to_numpy(dtype=float)[group_mask],
            expected.loc[:, "sca.adj.P.Val"].to_numpy(dtype=float)[group_mask],
        )
        > DEPTH_GROUP_ADJ_P_VALUE_CORRELATION_MIN
    )
    assert (
        float(
            np.nanquantile(
                _absolute_diff(
                    observed.loc[:, "adj.P.Val"].to_numpy(dtype=float)[group_mask],
                    expected.loc[:, "sca.adj.P.Val"].to_numpy(dtype=float)[group_mask],
                ),
                0.99,
            )
        )
        < DEPTH_GROUP_ADJ_P_VALUE_P99_ABS_DIFF_MAX
    )


def test_deqms_fixture_shows_depth_changes_prior_for_comparable_sites(
    workflow_result: DifferentialAnalysisResult,
) -> None:
    result = workflow_result
    observed = result.table_for("B_vs_A")
    expected = _load_expected(observed.index)
    metadata = _load_feature_metadata(observed.index)
    matrix = _load_matrix()
    pair_metadata = metadata.loc[
        metadata.loc[:, "comparable_pair_id"].eq("same_intensity_different_depth_pair")
    ]

    assert pair_metadata.index.size == 2
    depth_values = pair_metadata.loc[:, "quantification_depth"].astype(float)
    low_key = cast(str, pair_metadata.index[int(depth_values.argmin())])
    high_key = cast(str, pair_metadata.index[int(depth_values.argmax())])
    low_display = str(pair_metadata.loc[low_key, "site_id"])
    high_display = str(pair_metadata.loc[high_key, "site_id"])

    assert pair_metadata.loc[low_key, "quantification_depth"] == 2
    assert pair_metadata.loc[high_key, "quantification_depth"] == 32
    np.testing.assert_allclose(
        matrix.loc[low_display].to_numpy(dtype=float),
        matrix.loc[high_display].to_numpy(dtype=float),
        rtol=0.0,
        atol=0.0,
    )
    np.testing.assert_allclose(
        observed.loc[[low_key, high_key], "logFC"].to_numpy(dtype=float),
        expected.loc[[low_key, high_key], "logFC"].to_numpy(dtype=float),
        rtol=LOGFC_RTOL,
        atol=LOGFC_ATOL,
    )
    np.testing.assert_allclose(
        result.residual_variance_series()
        .loc[[low_key, high_key]]
        .to_numpy(dtype=float),
        expected.loc[[low_key, high_key], "residual_variance"].to_numpy(dtype=float),
        rtol=1.0e-12,
        atol=1.0e-12,
    )

    deqms_prior = expected.loc[[low_key, high_key], "sca.priorvar"].to_numpy(
        dtype=float
    )
    phospy_prior = (
        result.prior_residual_variance_series()
        .loc[[low_key, high_key]]
        .to_numpy(dtype=float)
    )
    assert deqms_prior[0] > deqms_prior[1] * 4.0
    assert phospy_prior[0] > phospy_prior[1] * 4.0
