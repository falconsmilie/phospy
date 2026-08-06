from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from phospy.science.activities.methods.ssgsea_substrate_enrichment import (
    SSGSEA_STATUS_COMPUTED,
    SsgseaSubstrateEnrichmentActivityMethod,
)
from phospy.science.activities.semantics import ActivityInputMatrix
from tests.support.performance_contracts import (
    SSGSEA_ACTIVITY_CONTRACT_N_KINASES,
    SSGSEA_ACTIVITY_CONTRACT_N_PROFILES,
    SSGSEA_ACTIVITY_CONTRACT_N_SITES,
    SSGSEA_ACTIVITY_CONTRACT_PERMUTATIONS,
    SSGSEA_ACTIVITY_CONTRACT_SUBSTRATES_PER_KINASE,
    SSGSEA_ACTIVITY_PEAK_MIB_MAX,
    SSGSEA_ACTIVITY_RUNTIME_SECONDS_MAX,
    deterministic_analysis_ready_site_keys,
    deterministic_kinase_substrate_map,
    measure_runtime_and_peak_mib,
)

pytestmark = [pytest.mark.performance, pytest.mark.release_gate]


def _bounded_ssgsea_inputs(
    *,
    n_sites: int = SSGSEA_ACTIVITY_CONTRACT_N_SITES,
    n_kinases: int = SSGSEA_ACTIVITY_CONTRACT_N_KINASES,
    n_profiles: int = SSGSEA_ACTIVITY_CONTRACT_N_PROFILES,
    substrates_per_kinase: int = SSGSEA_ACTIVITY_CONTRACT_SUBSTRATES_PER_KINASE,
) -> tuple[ActivityInputMatrix, pd.DataFrame]:
    site_ids = deterministic_analysis_ready_site_keys(
        n_sites,
        start=310_000,
        gene_prefix="ACTSITE",
    )
    rng = np.random.default_rng(20260805)
    values = rng.normal(loc=0.0, scale=1.0, size=(int(n_sites), int(n_profiles)))
    values = np.round(values, decimals=2)
    effect_matrix = pd.DataFrame(
        values,
        index=site_ids,
        columns=pd.Index(
            [f"effect_profile_{index + 1:02d}" for index in range(int(n_profiles))],
            name="profile_id",
        ),
        dtype=float,
    )
    membership = deterministic_kinase_substrate_map(
        dataset_site_ids=site_ids,
        eligible_kinase_count=n_kinases,
        substrates_per_kinase=substrates_per_kinase,
    )
    return ActivityInputMatrix.standardised_effect(effect_matrix), membership


def _run_bounded_ssgsea(
    *,
    n_sites: int = SSGSEA_ACTIVITY_CONTRACT_N_SITES,
    n_kinases: int = SSGSEA_ACTIVITY_CONTRACT_N_KINASES,
    n_profiles: int = SSGSEA_ACTIVITY_CONTRACT_N_PROFILES,
    permutation_count: int = SSGSEA_ACTIVITY_CONTRACT_PERMUTATIONS,
):
    activity_input, membership = _bounded_ssgsea_inputs(
        n_sites=n_sites,
        n_kinases=n_kinases,
        n_profiles=n_profiles,
    )
    return SsgseaSubstrateEnrichmentActivityMethod(
        min_substrates=4,
        permutation_count=int(permutation_count),
        random_seed=20260805,
        adjust_p_values=True,
    ).run(
        activity_input=activity_input,
        kinase_substrate_membership=membership,
    )


def test_ssgsea_seeded_permutation_performance_contract_is_bounded() -> None:
    result, runtime_seconds, peak_mib = measure_runtime_and_peak_mib(
        _run_bounded_ssgsea,
        warmup=True,
    )

    assert result.activity_matrix.shape == (
        SSGSEA_ACTIVITY_CONTRACT_N_KINASES,
        SSGSEA_ACTIVITY_CONTRACT_N_PROFILES,
    )
    assert result.substrate_count_matrix.shape == result.activity_matrix.shape
    assert int(result.substrate_count_matrix.min().min()) == (
        SSGSEA_ACTIVITY_CONTRACT_SUBSTRATES_PER_KINASE
    )
    assert result.p_value_matrix is not None
    assert result.q_value_matrix is not None
    assert result.p_value_matrix.stack().between(0.0, 1.0).all()
    assert result.q_value_matrix.stack().between(0.0, 1.0).all()
    assert result.statistics_table is not None
    assert set(result.statistics_table["computability_status"]) == {
        SSGSEA_STATUS_COMPUTED
    }
    assert runtime_seconds < SSGSEA_ACTIVITY_RUNTIME_SECONDS_MAX
    assert peak_mib < SSGSEA_ACTIVITY_PEAK_MIB_MAX


def test_ssgsea_null_calculation_memory_allocation_sanity() -> None:
    result, runtime_seconds, peak_mib = measure_runtime_and_peak_mib(
        lambda: _run_bounded_ssgsea(
            n_sites=600,
            n_kinases=40,
            n_profiles=4,
            permutation_count=32,
        ),
        warmup=True,
    )

    assert result.p_value_matrix is not None
    assert result.p_value_matrix.shape == (40, 4)
    assert runtime_seconds < SSGSEA_ACTIVITY_RUNTIME_SECONDS_MAX
    assert peak_mib < SSGSEA_ACTIVITY_PEAK_MIB_MAX
