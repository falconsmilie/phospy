from __future__ import annotations

import re
from dataclasses import dataclass
from functools import lru_cache

import pandas as pd

from phospy import (
    AnalysisReadyDatasetBuilder,
    KinaseWorkflow,
)
from phospy.api import (
    DatasetBuildRequest,
    KinasePredictionConfig,
    KinaseScoringConfig,
    KinaseWorkflowRequest,
    Organism,
    ReferenceBundle,
)
from tests.support.rewrite_fixture_data import (
    load_public_predmat_input_phospho,
    load_public_predmat_input_site_sequences,
    load_public_predmat_input_substrate_map,
    load_public_predmat_legacy_default_donor,
    load_public_predmat_legacy_r_parity_donor,
    load_public_predmat_rewrite_contract,
    load_public_predmat_rewrite_r_parity,
    load_public_predmat_rewrite_stable,
)

_CANONICAL_SITE_ID_PATTERN = re.compile(r"^\s*[^;]+\s*;\s*[^;]+\s*;\s*$")


@dataclass(frozen=True, slots=True)
class PublicPredmatLaneMetrics:
    adaptive_policy: str
    observed_shape: tuple[int, int]
    expected_shape: tuple[int, int]
    row_identity_match: bool
    column_identity_match: bool
    mean_abs_diff: float
    max_abs_diff: float
    dominant_matches: int
    dominant_total: int
    deterministic_under_seed: bool
    donor_corr: float
    donor_mae: float


@dataclass(frozen=True, slots=True)
class PublicPredmatBenchmarkMetrics:
    stable: PublicPredmatLaneMetrics
    r_parity: PublicPredmatLaneMetrics
    cross_policy_corr: float
    cross_policy_mae: float


@dataclass(frozen=True, slots=True)
class PublicPredmatOrderInvarianceMetrics:
    output_shape: tuple[int, int]
    normalized_equal: bool
    deterministic_under_seed: bool
    dominant_matches: int
    dominant_total: int


def _stack_frame(frame: pd.DataFrame) -> pd.Series:
    try:
        return frame.stack(future_stack=True)
    except TypeError:
        # pandas<2.1 compatibility path (no future_stack argument)
        return frame.stack(dropna=False)


def _canonical_public_site_components(site_id: object) -> tuple[str, str, str]:
    raw_site = str(site_id).strip()
    if _CANONICAL_SITE_ID_PATTERN.fullmatch(raw_site):
        parts = raw_site.split(";")
        gene_symbol = parts[0].strip()
        site = parts[1].strip()
        return f"{gene_symbol};{site};", gene_symbol, site

    gene_symbol = raw_site.split("_", 1)[0].strip()
    site = raw_site
    return f"{gene_symbol};{site};", gene_symbol, site


def _canonicalize_public_site_frame_index(frame: pd.DataFrame) -> pd.DataFrame:
    normalized = frame.copy(deep=True)
    canonical_index = [
        _canonical_public_site_components(site_id)[0] for site_id in normalized.index
    ]
    normalized.index = pd.Index(canonical_index, name=normalized.index.name)
    return normalized


def _build_public_predmat_dataset():
    input_phospho = load_public_predmat_input_phospho()
    raw_site_ids = input_phospho.index.astype(str).tolist()
    phospho = input_phospho.copy(deep=True)
    site_sequences = load_public_predmat_input_site_sequences()
    canonical_components = [
        _canonical_public_site_components(site_id) for site_id in phospho.index
    ]
    phospho.index = pd.Index(
        [canonical_site_id for canonical_site_id, _, _ in canonical_components],
        name=phospho.index.name,
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": [gene_symbol for _, gene_symbol, _ in canonical_components],
            "site": [site for _, _, site in canonical_components],
            "site_sequence": [
                str(site_sequences[str(site_id).strip()]) for site_id in raw_site_ids
            ],
        },
        index=phospho.index.copy(),
    )
    return AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
        )
    )


def _build_public_predmat_references(*, reverse_order: bool) -> ReferenceBundle:
    substrate_map = load_public_predmat_input_substrate_map()
    site_sequences = load_public_predmat_input_site_sequences()

    kinase_items = list(substrate_map.items())
    sequence_items = list(site_sequences.items())
    if reverse_order:
        kinase_items = list(reversed(kinase_items))
        sequence_items = list(reversed(sequence_items))

    substrate_rows = [
        {
            "kinase": str(kinase),
            "substrate_site": _canonical_public_site_components(site_id)[0],
        }
        for kinase, sites in kinase_items
        for site_id in sites
    ]
    sequence_index = [
        _canonical_public_site_components(site_id)[0] for site_id, _ in sequence_items
    ]
    sequence_values = [str(sequence) for _, sequence in sequence_items]
    return ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(substrate_rows),
        site_sequences=pd.DataFrame(
            {"site_sequence": sequence_values},
            index=pd.Index(sequence_index, name="site_id"),
        ),
    )


def _run_public_predmat_lane(
    *, adaptive_policy: str, reverse_reference_order: bool
) -> pd.DataFrame:
    result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=_build_public_predmat_dataset(),
            references=_build_public_predmat_references(
                reverse_order=reverse_reference_order
            ),
            scoring_config=KinaseScoringConfig(min_substrates=2),
            prediction_config=KinasePredictionConfig(
                top_k=4,
                ensemble_size=3,
                mode="adaptive_ensemble",
                adaptive_policy=adaptive_policy,
                n_iterations=2,
                random_state=17,
            ),
            activity_config=None,
        )
    )
    return result.prediction_result.pred_mat


def _collect_lane_metrics(*, adaptive_policy: str) -> PublicPredmatLaneMetrics:
    observed = _run_public_predmat_lane(
        adaptive_policy=adaptive_policy,
        reverse_reference_order=False,
    )
    repeated = _run_public_predmat_lane(
        adaptive_policy=adaptive_policy,
        reverse_reference_order=False,
    )
    expected = (
        load_public_predmat_rewrite_stable()
        if adaptive_policy == "stable"
        else load_public_predmat_rewrite_r_parity()
    )
    donor = (
        load_public_predmat_legacy_default_donor()
        if adaptive_policy == "stable"
        else load_public_predmat_legacy_r_parity_donor()
    )
    expected = _canonicalize_public_site_frame_index(expected)
    donor = _canonicalize_public_site_frame_index(donor)
    observed_aligned = observed.sort_index().sort_index(axis=1)
    expected_aligned = expected.sort_index().sort_index(axis=1)
    delta = (observed_aligned - expected_aligned).abs().pipe(_stack_frame)
    dominant_expected = expected_aligned.idxmax(axis=1)
    dominant_observed = observed_aligned.idxmax(axis=1)
    donor_aligned = donor.sort_index().sort_index(axis=1)
    donor_long = _stack_frame(donor_aligned).rename("donor")
    observed_long = _stack_frame(observed_aligned).rename("observed")
    donor_merged = pd.concat([observed_long, donor_long], axis=1).dropna()
    donor_delta = (donor_merged.loc[:, "observed"] - donor_merged.loc[:, "donor"]).abs()

    return PublicPredmatLaneMetrics(
        adaptive_policy=adaptive_policy,
        observed_shape=(int(observed.shape[0]), int(observed.shape[1])),
        expected_shape=(int(expected.shape[0]), int(expected.shape[1])),
        row_identity_match=sorted(observed.index.astype(str))
        == sorted(expected.index.astype(str)),
        column_identity_match=sorted(observed.columns.astype(str))
        == sorted(expected.columns.astype(str)),
        mean_abs_diff=float(delta.mean()),
        max_abs_diff=float(delta.max()),
        dominant_matches=int((dominant_observed == dominant_expected).sum()),
        dominant_total=int(dominant_expected.shape[0]),
        deterministic_under_seed=observed.equals(repeated),
        donor_corr=float(
            donor_merged.loc[:, "observed"].corr(donor_merged.loc[:, "donor"])
        )
        if not donor_merged.empty
        else 0.0,
        donor_mae=float(donor_delta.mean()) if not donor_delta.empty else 0.0,
    )


@lru_cache(maxsize=1)
def collect_public_predmat_benchmark_metrics() -> PublicPredmatBenchmarkMetrics:
    stable = _collect_lane_metrics(adaptive_policy="stable")
    r_parity = _collect_lane_metrics(adaptive_policy="r_parity")
    stable_frame = (
        _run_public_predmat_lane(
            adaptive_policy="stable",
            reverse_reference_order=False,
        )
        .sort_index()
        .sort_index(axis=1)
    )
    r_parity_frame = (
        _run_public_predmat_lane(
            adaptive_policy="r_parity",
            reverse_reference_order=False,
        )
        .sort_index()
        .sort_index(axis=1)
    )
    merged = pd.concat(
        [
            _stack_frame(stable_frame).rename("stable"),
            _stack_frame(r_parity_frame).rename("r_parity"),
        ],
        axis=1,
    ).dropna()
    cross_delta = (merged.loc[:, "stable"] - merged.loc[:, "r_parity"]).abs()
    return PublicPredmatBenchmarkMetrics(
        stable=stable,
        r_parity=r_parity,
        cross_policy_corr=float(merged.loc[:, "stable"].corr(merged.loc[:, "r_parity"]))
        if not merged.empty
        else 0.0,
        cross_policy_mae=float(cross_delta.mean()) if not cross_delta.empty else 0.0,
    )


@lru_cache(maxsize=1)
def collect_public_predmat_order_invariance_metrics() -> (
    PublicPredmatOrderInvarianceMetrics
):
    reference = _run_public_predmat_lane(
        adaptive_policy="stable",
        reverse_reference_order=False,
    )
    reordered = _run_public_predmat_lane(
        adaptive_policy="stable",
        reverse_reference_order=True,
    )
    repeated = _run_public_predmat_lane(
        adaptive_policy="stable",
        reverse_reference_order=True,
    )
    reference_normalized = reference.sort_index().sort_index(axis=1)
    reordered_normalized = reordered.sort_index().sort_index(axis=1)
    dominant_reference = reference_normalized.idxmax(axis=1)
    dominant_reordered = reordered_normalized.idxmax(axis=1)
    return PublicPredmatOrderInvarianceMetrics(
        output_shape=(int(reference.shape[0]), int(reference.shape[1])),
        normalized_equal=reference_normalized.equals(reordered_normalized),
        deterministic_under_seed=reordered.equals(repeated),
        dominant_matches=int((dominant_reference == dominant_reordered).sum()),
        dominant_total=int(dominant_reference.shape[0]),
    )


@lru_cache(maxsize=1)
def load_public_predmat_contract() -> dict[str, object]:
    return load_public_predmat_rewrite_contract()
