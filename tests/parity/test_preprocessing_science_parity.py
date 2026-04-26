from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pandas.testing as pdt
import pytest

from phospy.api.configs import (
    DatasetComparisonBuildingConfig,
    DatasetPreprocessingConfig,
    DatasetSiteMatrixConfig,
    DatasetTotalProteinCorrectionConfig,
)
from phospy.datasets.builders.preprocessing import DatasetPreprocessor
from phospy.datasets.preprocessing.models import PreprocessingPlan
from tests.support.parity_reporting import (
    format_fraction,
    format_shape,
    record_parity_metrics,
)

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
        },
        index=site_ids.copy(),
    )


def test_ratio_to_total_total_protein_correction_matches_rewrite_reference_fixture(
    request: pytest.FixtureRequest,
) -> None:
    phospho = pd.read_csv(
        PROTEIN_CORRECTION_FIXTURES / "reference_input_phospho.csv"
    ).set_index("site_id")
    site_metadata = pd.read_csv(
        PROTEIN_CORRECTION_FIXTURES / "reference_input_site_metadata.csv"
    ).set_index("site_id")
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
        plan=PreprocessingPlan.from_config(
            DatasetPreprocessingConfig(
                total_protein_correction=DatasetTotalProteinCorrectionConfig(
                    policy="ratio_to_total"
                )
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
            ("ratio_to_total output shape", format_shape(*preprocessed.phospho.shape)),
            ("ratio_to_total duplicate site rows", duplicate_site_rows),
            ("ratio_to_total max abs diff", max_abs_diff),
        ],
        notes=(
            "fixture lane: tests/fixtures/rewrite_parity/protein_correction/",
            "policy lane: total_protein_correction.policy=ratio_to_total",
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
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": corrected_fixture.loc[:, "gene_names"].astype(str).tolist(),
            "site": site_tokens.loc[:, 1].astype(str).tolist(),
            "site_sequence": corrected_fixture.loc[:, "centralized_sequence"]
            .astype(str)
            .tolist(),
        },
        index=phospho.index.copy(),
    )

    preprocessed = DatasetPreprocessor().run(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan.from_config(
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
    expected_phospho.index = pd.Index(
        expected_phospho.index.astype(str), name="source_uid"
    )
    expected_site_metadata = pd.DataFrame(
        {
            "gene_symbol": expected_input_fixture.loc[:, "gene_names"]
            .astype(str)
            .tolist(),
            "site": expected_input_fixture.loc[:, "p_site"].astype(str).tolist(),
            "site_sequence": expected_input_fixture.loc[:, "centralized_sequence"]
            .astype(str)
            .tolist(),
        },
        index=pd.Index(
            expected_input_fixture.loc[:, "site_id"].astype(str),
            name="source_uid",
        ),
    )

    pdt.assert_frame_equal(preprocessed.phospho, expected_phospho)
    pdt.assert_frame_equal(
        preprocessed.site_metadata.loc[:, ["gene_symbol", "site", "site_sequence"]],
        expected_site_metadata,
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
