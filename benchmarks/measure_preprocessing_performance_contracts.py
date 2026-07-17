#!/usr/bin/env python3
"""Benchmark preprocessing hot paths used in performance contracts.

Targets:
- `phospy.science.datasets.preprocessing.stages.missing_data.MissingDataStage`
- `phospy.science.datasets.preprocessing.stages.normalisation.NormalisationStage`
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from phospy.science.datasets.preprocessing.models import PreprocessingState

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))


def _build_base_matrix(*, n_sites: int, n_samples: int) -> pd.DataFrame:
    from tests.support.performance_contracts import deterministic_matrix

    return deterministic_matrix(n_sites=n_sites, n_samples=n_samples, seed=1439)


def _run_median_center(phospho: pd.DataFrame) -> PreprocessingState:
    from phospy.science.datasets.preprocessing.models import (
        PreprocessingPlan,
        PreprocessingState,
    )
    from phospy.science.datasets.preprocessing.stages.normalisation import (
        NormalisationStage,
    )
    from tests.support.performance_contracts import deterministic_site_metadata

    state = PreprocessingState(
        phospho=phospho,
        site_metadata=deterministic_site_metadata(
            phospho.index, include_protein_id=False
        ),
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan(normalisation_policy="median_center"),
    )
    return NormalisationStage().run(state).state


def _run_quantile_normalisation(phospho: pd.DataFrame) -> PreprocessingState:
    from phospy.science.datasets.preprocessing.models import (
        PreprocessingPlan,
        PreprocessingState,
    )
    from phospy.science.datasets.preprocessing.stages.normalisation import (
        NormalisationStage,
    )
    from tests.support.performance_contracts import deterministic_site_metadata

    state = PreprocessingState(
        phospho=phospho,
        site_metadata=deterministic_site_metadata(
            phospho.index, include_protein_id=False
        ),
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan(normalisation_policy="quantile"),
    )
    return NormalisationStage().run(state).state


def _run_missing_data_imputation(phospho: pd.DataFrame) -> PreprocessingState:
    from phospy.science.datasets.preprocessing.models import (
        PreprocessingPlan,
        PreprocessingState,
    )
    from phospy.science.datasets.preprocessing.stages.missing_data import (
        MissingDataStage,
    )
    from tests.support.performance_contracts import deterministic_site_metadata

    state = PreprocessingState(
        phospho=phospho,
        site_metadata=deterministic_site_metadata(
            phospho.index, include_protein_id=False
        ),
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan(
            missing_data_policy="impute_row_median",
            missing_data_min_observed_values=8,
        ),
    )
    return MissingDataStage().run(state).state


def _run_knn_imputation(phospho: pd.DataFrame) -> object:
    from phospy.science.datasets.preprocessing.models import (
        PreprocessingPlan,
        PreprocessingState,
    )
    from phospy.science.datasets.preprocessing.stages.missing_data.knn import (
        run_knn_policy,
    )
    from tests.support.performance_contracts import deterministic_site_metadata

    state = PreprocessingState(
        phospho=phospho,
        site_metadata=deterministic_site_metadata(
            phospho.index, include_protein_id=False
        ),
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan(
            missing_data_policy="impute_knn",
            missing_data_k=1,
            missing_data_distance="nan_euclidean",
            missing_data_max_missing_fraction_per_row=0.5,
            stage_order=("missing_data",),
        ),
    )
    return run_knn_policy(state)


def _with_knn_target_missingness(
    phospho: pd.DataFrame,
    *,
    missing_target_rows: int,
    missing_cells_per_target_row: int = 1,
) -> pd.DataFrame:
    with_missing = phospho.copy(deep=True)
    target_count = min(int(missing_target_rows), int(with_missing.shape[0]))
    if target_count <= 0:
        return with_missing
    row_positions = np.floor(
        np.linspace(
            0,
            int(with_missing.shape[0]),
            num=target_count,
            endpoint=False,
        )
    ).astype(int, copy=False)
    sample_count = int(with_missing.shape[1])
    cells_per_row = min(int(missing_cells_per_target_row), sample_count)
    for offset, row_position in enumerate(np.unique(row_positions)):
        for missing_offset in range(cells_per_row):
            column_position = (offset + (missing_offset * 7)) % sample_count
            with_missing.iat[int(row_position), int(column_position)] = float("nan")
    return with_missing


def main() -> None:
    from tests.support.performance_contracts import (
        KNN_IMPUTATION_BENCHMARK_TIERS,
        PREPROCESSING_CONTRACT_MISSING_FRACTION,
        PREPROCESSING_CONTRACT_N_SAMPLES,
        PREPROCESSING_CONTRACT_N_SITES,
        median_runtime_and_peak_mib,
        with_missing_fraction,
    )

    repeats = 3
    n_sites = PREPROCESSING_CONTRACT_N_SITES
    n_samples = PREPROCESSING_CONTRACT_N_SAMPLES
    baseline = _build_base_matrix(n_sites=n_sites, n_samples=n_samples)
    with_missing = with_missing_fraction(
        baseline,
        missing_fraction=PREPROCESSING_CONTRACT_MISSING_FRACTION,
        seed=271,
    )

    missing_state, missing_runtime_seconds, missing_peak_mib = (
        median_runtime_and_peak_mib(
            lambda: _run_missing_data_imputation(with_missing),
            repeats=repeats,
            warmup=True,
        )
    )
    median_state, median_runtime_seconds, median_peak_mib = median_runtime_and_peak_mib(
        lambda: _run_median_center(baseline),
        repeats=repeats,
        warmup=True,
    )
    quantile_state, quantile_runtime_seconds, quantile_peak_mib = (
        median_runtime_and_peak_mib(
            lambda: _run_quantile_normalisation(baseline),
            repeats=repeats,
            warmup=True,
        )
    )

    missing_row_audit = missing_state.row_audit
    dropped_rows = 0
    imputed_rows = 0
    if missing_row_audit is not None and not missing_row_audit.empty:
        dropped_rows = int(
            (
                (missing_row_audit.loc[:, "stage"] == "missing_data")
                & (missing_row_audit.loc[:, "action"] == "dropped")
            ).sum()
        )
        imputed_rows = int(
            (
                (missing_row_audit.loc[:, "stage"] == "missing_data")
                & (missing_row_audit.loc[:, "action"] == "imputed")
            ).sum()
        )

    print(f"repeats={repeats}")
    print(f"n_sites={n_sites}")
    print(f"n_samples={n_samples}")
    print(f"missing_fraction={PREPROCESSING_CONTRACT_MISSING_FRACTION:.3f}")
    print(f"missing_runtime_seconds={missing_runtime_seconds:.6f}")
    print(f"missing_peak_mib={missing_peak_mib:.3f}")
    print(f"missing_output_rows={missing_state.phospho.shape[0]}")
    print(f"missing_dropped_rows={dropped_rows}")
    print(f"missing_imputed_rows={imputed_rows}")
    print(f"median_center_runtime_seconds={median_runtime_seconds:.6f}")
    print(f"median_center_peak_mib={median_peak_mib:.3f}")
    print(f"median_center_output_rows={median_state.phospho.shape[0]}")
    print(f"quantile_runtime_seconds={quantile_runtime_seconds:.6f}")
    print(f"quantile_peak_mib={quantile_peak_mib:.3f}")
    print(f"quantile_output_rows={quantile_state.phospho.shape[0]}")

    for knn_tier in KNN_IMPUTATION_BENCHMARK_TIERS:
        knn_baseline = _build_base_matrix(
            n_sites=int(knn_tier.site_count),
            n_samples=int(knn_tier.sample_count),
        )
        knn_input = _with_knn_target_missingness(
            knn_baseline,
            missing_target_rows=int(knn_tier.missing_target_rows),
            missing_cells_per_target_row=int(knn_tier.missing_cells_per_target_row),
        )
        knn_outcome, knn_runtime_seconds, knn_peak_mib = median_runtime_and_peak_mib(
            lambda knn_frame=knn_input: _run_knn_imputation(knn_frame),
            repeats=1,
            warmup=False,
        )
        print(f"knn_{knn_tier.case_id}_n_sites={knn_tier.site_count}")
        print(f"knn_{knn_tier.case_id}_n_samples={knn_tier.sample_count}")
        print(
            f"knn_{knn_tier.case_id}_missing_target_rows={knn_tier.missing_target_rows}"
        )
        print(f"knn_{knn_tier.case_id}_runtime_seconds={knn_runtime_seconds:.6f}")
        print(f"knn_{knn_tier.case_id}_peak_mib={knn_peak_mib:.3f}")
        print(f"knn_{knn_tier.case_id}_imputed_cells={knn_outcome.imputed_cell_count}")


if __name__ == "__main__":
    main()
