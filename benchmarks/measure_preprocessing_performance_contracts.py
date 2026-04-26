#!/usr/bin/env python3
"""Benchmark preprocessing hot paths used in performance contracts.

Targets:
- `phospy.datasets.preprocessing.stages.missing_data.MissingDataStage`
- `phospy.datasets.preprocessing.stages.normalisation.NormalisationStage`
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd

if TYPE_CHECKING:
    from phospy.datasets.preprocessing.models import PreprocessingState

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
    from phospy.datasets.preprocessing.models import (
        PreprocessingPlan,
        PreprocessingState,
    )
    from phospy.datasets.preprocessing.stages.normalisation import NormalisationStage
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
    return NormalisationStage().run(state)


def _run_quantile_normalisation(phospho: pd.DataFrame) -> PreprocessingState:
    from phospy.datasets.preprocessing.models import (
        PreprocessingPlan,
        PreprocessingState,
    )
    from phospy.datasets.preprocessing.stages.normalisation import NormalisationStage
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
    return NormalisationStage().run(state)


def _run_missing_data_imputation(phospho: pd.DataFrame) -> PreprocessingState:
    from phospy.datasets.preprocessing.models import (
        PreprocessingPlan,
        PreprocessingState,
    )
    from phospy.datasets.preprocessing.stages.missing_data import MissingDataStage
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
    return MissingDataStage().run(state)


def main() -> None:
    from tests.support.performance_contracts import (
        median_runtime_and_peak_mib,
        with_missing_fraction,
    )

    repeats = 3
    n_sites = 2_500
    n_samples = 20
    baseline = _build_base_matrix(n_sites=n_sites, n_samples=n_samples)
    with_missing = with_missing_fraction(
        baseline,
        missing_fraction=0.12,
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


if __name__ == "__main__":
    main()
