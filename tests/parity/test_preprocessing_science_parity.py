from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd
import pandas.testing as pdt
import pytest

from phospy.api.configs import (
    DatasetComparisonBuildingConfig,
    DatasetIntensityTransformConfig,
    DatasetMissingDataConfig,
    DatasetPreprocessingConfig,
    DatasetSiteMatrixConfig,
    DatasetTotalProteinCorrectionConfig,
)
from phospy.science.datasets.builders.preprocessing import DatasetPreprocessor
from phospy.science.datasets.preprocessing.models import PreprocessingPlan
from tests.support.parity_reporting import (
    format_fraction,
    format_shape,
    record_parity_metrics,
)
from tests.support.site_keys import site_key_index_from_display_ids

pytestmark = pytest.mark.parity

ROOT = Path(__file__).resolve().parents[2]
REWRITE_PARITY_ROOT = ROOT / "tests" / "fixtures" / "rewrite_parity"
PROTEIN_CORRECTION_FIXTURES = REWRITE_PARITY_ROOT / "protein_correction"
SITE_MATRIX_FIXTURES = REWRITE_PARITY_ROOT / "site_matrix"
COMPARISON_FIXTURES = REWRITE_PARITY_ROOT / "comparison_building"


def _site_metadata_from_site_ids(site_ids: pd.Index) -> pd.DataFrame:
    rows: list[tuple[str, str, str]] = []
    for site_id in site_ids.astype(str):
        gene_symbol, site, _ = site_id.split(";", 2)
        rows.append((gene_symbol, site, f"SEQ_{gene_symbol}_{site}"))
    return pd.DataFrame(
        {
            "gene_symbol": [row[0] for row in rows],
            "site": [row[1] for row in rows],
            "site_sequence": [row[2] for row in rows],
            "localisation_confidence": [0.95] * len(rows),
        },
        index=site_ids.copy(),
    )


def _plan_without_missing_stage(
    config: DatasetPreprocessingConfig,
) -> PreprocessingPlan:
    plan = PreprocessingPlan.from_config(config)
    return replace(
        plan,
        stage_order=tuple(
            stage for stage in plan.stage_order if stage != "missing_data"
        ),
    )


def test_row_median_missing_data_policy_is_deterministic_and_provenance_backed(
    request: pytest.FixtureRequest,
) -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, float("nan"), float("nan")],
            "sample_b": [2.0, 10.0, float("nan")],
            "sample_c": [3.0, 20.0, 9.0],
            "sample_d": [4.0, float("nan"), float("nan")],
        },
        index=pd.Index(["row_keep", "row_impute", "row_drop"], name="source_row"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "AKT1", "GSK3B"],
            "site": ["Y182", "T308", "S9"],
            "site_sequence": ["SEQ_A", "SEQ_R", "SEQ_C"],
            "localisation_confidence": [0.95] * phospho.shape[0],
        },
        index=phospho.index.copy(),
    )
    config = DatasetPreprocessingConfig(
        missing_data=DatasetMissingDataConfig(
            policy="impute_row_median",
            min_observed_values=2,
        )
    )
    plan = PreprocessingPlan.from_config(config)

    first = DatasetPreprocessor().run(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=None,
        total=None,
        plan=plan,
    )
    second = DatasetPreprocessor().run(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=None,
        total=None,
        plan=plan,
    )

    expected = pd.DataFrame(
        {
            "sample_a": [1.0, 15.0],
            "sample_b": [2.0, 10.0],
            "sample_c": [3.0, 20.0],
            "sample_d": [4.0, 15.0],
        },
        index=pd.Index(["row_keep", "row_impute"], name="source_row"),
    )
    pdt.assert_frame_equal(first.phospho, expected)
    pdt.assert_frame_equal(second.phospho, expected)

    imputed_audit_rows = first.row_audit.loc[
        (first.row_audit.loc[:, "stage"] == "missing_data")
        & (first.row_audit.loc[:, "action"] == "imputed")
    ]
    assert imputed_audit_rows.shape[0] == 1
    imputed_snapshot = imputed_audit_rows.iloc[0]["parameter_snapshot"]
    assert isinstance(imputed_snapshot, dict)
    assert tuple(imputed_snapshot["imputed_columns"]) == ("sample_a", "sample_d")
    assert float(imputed_snapshot["row_median"]) == 15.0

    missing_stage = next(
        stage for stage in first.preprocessing_trace if stage.stage == "missing_data"
    )
    diagnostics = dict(missing_stage.diagnostics)
    assert diagnostics["input_missing_cell_count"] == 5
    assert diagnostics["output_missing_cell_count"] == 0
    assert diagnostics["imputed_cell_count"] == 2
    assert diagnostics["affected_row_count"] == 2
    assert diagnostics["affected_column_count"] == 3
    assert diagnostics["dropped_row_ids"] == ["row_drop"]
    assert diagnostics["stage_order"] == [
        "localisation_confidence",
        "missing_data",
    ]
    assert isinstance(diagnostics["missingness_mask_hash"], str)

    record_parity_metrics(
        request.config,
        family="preprocessing_science",
        metrics=[
            ("row_median output shape", format_shape(*first.phospho.shape)),
            ("row_median imputed cells", int(diagnostics["imputed_cell_count"])),
            ("row_median dropped rows", len(diagnostics["dropped_row_ids"])),
        ],
        notes=(
            "policy lane: missing_data.policy=impute_row_median",
            "deterministic lane: repeated run outputs match exactly",
        ),
    )


def test_minprob_missing_data_policy_is_deterministic_and_seeded(
    request: pytest.FixtureRequest,
) -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [10.0, float("nan"), float("nan"), 4.0],
            "sample_b": [9.0, 8.0, float("nan"), 6.0],
            "sample_c": [11.0, 7.0, 5.0, float("nan")],
        },
        index=pd.Index(["row_keep", "row_impute_a", "row_drop", "row_impute_c"]),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "AKT1", "GSK3B", "PRKACA"],
            "site": ["Y182", "T308", "S9", "S339"],
            "site_sequence": ["SEQ_A", "SEQ_R", "SEQ_C", "SEQ_D"],
            "localisation_confidence": [0.95] * phospho.shape[0],
        },
        index=phospho.index.copy(),
    )
    config = DatasetPreprocessingConfig(
        intensity_transform=DatasetIntensityTransformConfig(
            policy="log2", pseudocount=1.0
        ),
        missing_data=DatasetMissingDataConfig(
            policy="impute_minprob",
            q=0.01,
            width=0.3,
            seed=12345,
            max_missing_fraction_per_row=0.5,
        ),
    )
    plan = PreprocessingPlan.from_config(config)

    first = DatasetPreprocessor().run(
        phospho=phospho.copy(deep=True),
        site_metadata=site_metadata.copy(deep=True),
        sample_metadata=None,
        total=None,
        plan=plan,
    )
    second = DatasetPreprocessor().run(
        phospho=phospho.copy(deep=True),
        site_metadata=site_metadata.copy(deep=True),
        sample_metadata=None,
        total=None,
        plan=plan,
    )

    pdt.assert_frame_equal(first.phospho, second.phospho)
    assert first.phospho.isna().to_numpy().sum() == 0
    assert first.phospho.index.tolist() == ["row_keep", "row_impute_a", "row_impute_c"]

    missing_stage = next(
        stage for stage in first.preprocessing_trace if stage.stage == "missing_data"
    )
    diagnostics = dict(missing_stage.diagnostics)
    assert diagnostics["imputation_method_id"] == "minprob"
    assert diagnostics["left_censored_assumption"] is True
    assert diagnostics["random_seed"] == 12345
    assert diagnostics["dropped_rows_above_max_missing_fraction"] == ["row_drop"]
    assert "sample_a" in diagnostics["per_column_distribution_parameters"]

    record_parity_metrics(
        request.config,
        family="preprocessing_science",
        metrics=[
            ("minprob output shape", format_shape(*first.phospho.shape)),
            ("minprob imputed cells", int(diagnostics["imputed_cell_count"])),
            ("minprob dropped rows", len(diagnostics["dropped_row_ids"])),
        ],
        notes=(
            "policy lane: missing_data.policy=impute_minprob with log2 transform",
            "deterministic lane: fixed-seed repeated run outputs match exactly",
        ),
    )


def test_subtract_log_total_total_protein_correction_matches_rewrite_reference_fixture(
    request: pytest.FixtureRequest,
) -> None:
    phospho = pd.read_csv(
        PROTEIN_CORRECTION_FIXTURES / "reference_input_phospho.csv"
    ).set_index("site_id")
    site_metadata = pd.read_csv(
        PROTEIN_CORRECTION_FIXTURES / "reference_input_site_metadata.csv"
    ).set_index("site_id")
    if "localisation_confidence" not in site_metadata.columns:
        site_metadata.loc[:, "localisation_confidence"] = [0.95] * int(
            site_metadata.shape[0]
        )
    total = pd.read_csv(
        PROTEIN_CORRECTION_FIXTURES / "reference_input_total.csv"
    ).set_index("protein_id")
    expected = pd.read_csv(
        PROTEIN_CORRECTION_FIXTURES / "reference_corrected_matrix.csv"
    ).set_index("site_id")

    preprocessed = DatasetPreprocessor().run(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=None,
        total=total,
        plan=_plan_without_missing_stage(
            DatasetPreprocessingConfig(
                intensity_transform=DatasetIntensityTransformConfig(
                    policy="log2",
                    pseudocount=1.0,
                ),
                total_protein_correction=DatasetTotalProteinCorrectionConfig(
                    policy="subtract_log_total"
                ),
            )
        ),
    )

    observed = preprocessed.phospho.astype(float)
    expected = expected.astype(float)
    pdt.assert_frame_equal(observed, expected)
    assert preprocessed.phospho.columns.equals(phospho.columns)
    assert preprocessed.phospho.index.tolist() == phospho.index.tolist()

    absolute_delta = (observed - expected).abs()
    absolute_delta_values = absolute_delta.to_numpy(dtype=float)
    finite_delta_values = absolute_delta_values[~np.isnan(absolute_delta_values)]
    max_abs_diff = float(finite_delta_values.max()) if finite_delta_values.size else 0.0
    duplicate_site_rows = int(preprocessed.phospho.index.duplicated(keep=False).sum())
    record_parity_metrics(
        request.config,
        family="preprocessing_science",
        metrics=[
            (
                "subtract_log_total output shape",
                format_shape(*preprocessed.phospho.shape),
            ),
            ("subtract_log_total duplicate site rows", duplicate_site_rows),
            ("subtract_log_total max abs diff", max_abs_diff),
        ],
        notes=(
            "fixture lane: tests/fixtures/rewrite_parity/protein_correction/",
            "policy lane: total_protein_correction.policy=subtract_log_total (after log2)",
        ),
    )


def test_site_matrix_build_from_metadata_matches_rewrite_reference_fixture(
    request: pytest.FixtureRequest,
) -> None:
    corrected_fixture = pd.read_csv(
        SITE_MATRIX_FIXTURES / "reference_phospho_corrected.csv"
    )
    expected_matrix_fixture = pd.read_csv(
        SITE_MATRIX_FIXTURES / "reference_expected_matrix.csv"
    )
    expected_input_fixture = pd.read_csv(
        SITE_MATRIX_FIXTURES / "reference_expected_input.csv"
    )

    corrected_cols = tuple(f"phospho_corrected_{position}" for position in range(1, 7))
    phospho = corrected_fixture.loc[:, list(corrected_cols)].astype(float).copy()
    phospho.index = pd.Index(
        corrected_fixture.loc[:, "uid"].astype(str), name="source_uid"
    )

    site_tokens = (
        corrected_fixture.loc[:, "gene_p_site"]
        .astype(str)
        .str.split("_", n=1, expand=True)
    )
    display_ids = [
        f"{gene_symbol};{site};"
        for gene_symbol, site in zip(
            corrected_fixture.loc[:, "gene_names"].astype(str).tolist(),
            site_tokens.loc[:, 1].astype(str).tolist(),
            strict=True,
        )
    ]
    site_keys = site_key_index_from_display_ids(
        display_ids,
        protein_namespace="gene_symbol",
    )
    site_metadata = pd.DataFrame(
        {
            "site_key": site_keys.astype(str).tolist(),
            "display_id": display_ids,
            "gene_symbol": corrected_fixture.loc[:, "gene_names"].astype(str).tolist(),
            "site": site_tokens.loc[:, 1].astype(str).tolist(),
            "site_sequence": corrected_fixture.loc[:, "centralized_sequence"]
            .astype(str)
            .tolist(),
            "localisation_confidence": [0.95] * phospho.shape[0],
        },
        index=phospho.index.copy(),
    )

    preprocessed = DatasetPreprocessor().run(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=None,
        total=None,
        plan=_plan_without_missing_stage(
            DatasetPreprocessingConfig(
                site_matrix=DatasetSiteMatrixConfig(policy="build_from_metadata")
            )
        ),
    )

    expected_phospho = (
        expected_matrix_fixture.set_index(expected_matrix_fixture.columns[0])
        .loc[:, list(corrected_cols)]
        .astype(float)
    )
    expected_display_ids = expected_phospho.index.astype(str).tolist()
    expected_site_keys = site_key_index_from_display_ids(
        expected_display_ids,
        protein_namespace="gene_symbol",
    )
    expected_phospho.index = expected_site_keys.copy()
    expected_site_metadata = pd.DataFrame(
        {
            "gene_symbol": expected_input_fixture.loc[:, "gene_names"]
            .astype(str)
            .tolist(),
            "site": expected_input_fixture.loc[:, "p_site"].astype(str).tolist(),
            "site_sequence": expected_input_fixture.loc[:, "centralized_sequence"]
            .astype(str)
            .tolist(),
            "localisation_confidence": [0.95] * expected_input_fixture.shape[0],
        },
        index=site_key_index_from_display_ids(
            expected_input_fixture.loc[:, "site_id"].astype(str).tolist(),
            protein_namespace="gene_symbol",
        ),
    )

    pdt.assert_frame_equal(preprocessed.phospho, expected_phospho)
    pdt.assert_frame_equal(
        preprocessed.site_metadata.loc[
            :, ["gene_symbol", "site", "site_sequence", "localisation_confidence"]
        ],
        expected_site_metadata.loc[
            :, ["gene_symbol", "site", "site_sequence", "localisation_confidence"]
        ],
    )

    dropped = preprocessed.row_audit.loc[
        (preprocessed.row_audit.loc[:, "stage"] == "site_matrix")
        & (preprocessed.row_audit.loc[:, "action"] == "dropped")
    ]
    dropped_missing_sequence = dropped.loc[
        dropped.loc[:, "reason"]
        == "dropped because site_metadata.site_sequence is missing or blank"
    ]
    dropped_incomplete = dropped.loc[
        dropped.loc[:, "reason"] == "dropped by site_matrix missing-data policy"
    ]
    row_drop_stats = {
        "input_rows": 4,
        "dropped_missing_sequence": int(dropped_missing_sequence.shape[0]),
        "dropped_incomplete_values": int(dropped_incomplete.shape[0]),
        "retained_rows": int(preprocessed.phospho.shape[0]),
    }
    assert row_drop_stats == {
        "input_rows": 4,
        "dropped_missing_sequence": 0,
        "dropped_incomplete_values": 3,
        "retained_rows": 1,
    }

    record_parity_metrics(
        request.config,
        family="preprocessing_science",
        metrics=[
            ("site_matrix output shape", format_shape(*preprocessed.phospho.shape)),
            (
                "site_matrix retained rows",
                format_fraction(
                    int(row_drop_stats["retained_rows"]),
                    int(row_drop_stats["input_rows"]),
                    include_percent=True,
                ),
            ),
            (
                "site_matrix dropped incomplete rows",
                int(row_drop_stats["dropped_incomplete_values"]),
            ),
        ],
        notes=(
            "fixture lane: tests/fixtures/rewrite_parity/site_matrix/",
            "policy lane: site_matrix.policy=build_from_metadata (drop_any_missing + max_mean_signal)",
        ),
    )


def test_comparison_building_explicit_pair_matches_rewrite_reference_fixture(
    request: pytest.FixtureRequest,
) -> None:
    phospho = pd.read_csv(
        COMPARISON_FIXTURES / "reference_pairwise_input_phospho.csv"
    ).set_index("site_id")
    sample_metadata = pd.read_csv(
        COMPARISON_FIXTURES / "reference_pairwise_input_sample_metadata.csv"
    ).set_index("sample_id")
    expected = pd.read_csv(
        COMPARISON_FIXTURES / "reference_pairwise_expected.csv"
    ).set_index("site_id")
    site_metadata = _site_metadata_from_site_ids(phospho.index)

    preprocessed = DatasetPreprocessor().run(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=sample_metadata,
        total=None,
        plan=PreprocessingPlan.from_config(
            DatasetPreprocessingConfig(
                comparisons=DatasetComparisonBuildingConfig(
                    policy="sample_metadata_pairs",
                    pairs=(("sample_a", "sample_b"),),
                )
            )
        ),
    )

    assert preprocessed.comparisons is not None
    pdt.assert_frame_equal(
        preprocessed.comparisons.astype(float), expected.astype(float)
    )

    record_parity_metrics(
        request.config,
        family="preprocessing_science",
        metrics=[
            (
                "comparison explicit output shape",
                format_shape(*preprocessed.comparisons.shape),
            ),
            (
                "comparison explicit columns",
                ",".join(preprocessed.comparisons.columns.astype(str)),
            ),
        ],
        notes=(
            "fixture lane: tests/fixtures/rewrite_parity/comparison_building/",
            "policy lane: comparisons.policy=sample_metadata_pairs with explicit pairs",
        ),
    )


def test_comparison_building_inferred_pairs_match_rewrite_reference_fixture(
    request: pytest.FixtureRequest,
) -> None:
    phospho = pd.read_csv(
        COMPARISON_FIXTURES / "inferred_pairs_input_phospho.csv"
    ).set_index("site_id")
    sample_metadata = pd.read_csv(
        COMPARISON_FIXTURES / "inferred_pairs_input_sample_metadata.csv"
    ).set_index("sample_id")
    expected = pd.read_csv(
        COMPARISON_FIXTURES / "inferred_pairs_expected.csv"
    ).set_index("site_id")
    site_metadata = _site_metadata_from_site_ids(phospho.index)

    preprocessed = DatasetPreprocessor().run(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=sample_metadata,
        total=None,
        plan=PreprocessingPlan.from_config(
            DatasetPreprocessingConfig(
                comparisons=DatasetComparisonBuildingConfig(
                    policy="sample_metadata_pairs"
                )
            )
        ),
    )

    assert preprocessed.comparisons is not None
    assert preprocessed.comparisons.columns.tolist() == expected.columns.tolist()
    pdt.assert_frame_equal(
        preprocessed.comparisons.astype(float), expected.astype(float)
    )

    record_parity_metrics(
        request.config,
        family="preprocessing_science",
        metrics=[
            (
                "comparison inferred output shape",
                format_shape(*preprocessed.comparisons.shape),
            ),
            ("comparison inferred pair count", preprocessed.comparisons.shape[1]),
            (
                "comparison inferred columns",
                ",".join(preprocessed.comparisons.columns.astype(str)),
            ),
        ],
        notes=(
            "fixture lane: tests/fixtures/rewrite_parity/comparison_building/",
            "policy lane: comparisons.policy=sample_metadata_pairs with inferred all-pairs",
        ),
    )
