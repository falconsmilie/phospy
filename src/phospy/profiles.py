from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import pandas as pd

from .validation.primitives import validate_positive_int


@dataclass(slots=True)
class KinaseProfileResult:
    """Detached snapshot bundle for kinase substrate-profile outputs.

    The contained tables and lookup structures are produced workflow outputs, not
    live views into the input phosphosite matrix. Mutating them only affects this
    result instance.
    """

    profile_matrix: pd.DataFrame
    substrate_counts: pd.Series
    quantified_substrates: dict[str, list[str]]


def build_kinase_substrate_profiles(
    substrate_map: Mapping[str, Sequence[str]],
    phospho_matrix: pd.DataFrame,
    min_substrates: int = 1,
) -> KinaseProfileResult:
    """Build kinase substrate profiles from quantified phosphosite values.

    This mirrors the core behaviour of PhosR's ``kinaseSubstrateProfile()``:
    for each kinase, intersect its known substrates with the phosphosite matrix,
    use the single quantified row directly when exactly one site is available,
    and otherwise summarise quantified substrates column-wise with the median.
    """

    validate_positive_int(min_substrates, name="min_substrates")

    observed_sites = set(phospho_matrix.index)
    numeric_matrix = phospho_matrix.astype(float)

    profile_rows: dict[str, pd.Series] = {}
    quantified_substrates: dict[str, list[str]] = {}
    substrate_counts: dict[str, int] = {}

    for kinase, substrates in substrate_map.items():
        substrate_sequence = list(substrates)
        substrate_counts[kinase] = sum(
            substrate in observed_sites for substrate in substrate_sequence
        )
        quantified_sites = _quantified_sites(
            substrates=substrate_sequence,
            observed_sites=observed_sites,
        )
        if len(quantified_sites) < min_substrates:
            continue

        profile_rows[kinase] = _aggregate_quantified_sites(
            numeric_matrix.loc[quantified_sites, :],
        )
        quantified_substrates[kinase] = quantified_sites

    if profile_rows:
        profile_matrix = pd.DataFrame.from_dict(profile_rows, orient="index")
        profile_matrix = profile_matrix.loc[:, phospho_matrix.columns.copy()]
    else:
        profile_matrix = pd.DataFrame(
            columns=phospho_matrix.columns.copy(), dtype=float
        )

    profile_matrix.index.name = "kinase"

    count_series = pd.Series(substrate_counts, dtype="int64", name="NumSub")
    count_series.index.name = "kinase"

    return KinaseProfileResult(
        profile_matrix=profile_matrix,
        substrate_counts=count_series,
        quantified_substrates=quantified_substrates,
    )


def _aggregate_quantified_sites(quantified_matrix: pd.DataFrame) -> pd.Series:
    if quantified_matrix.shape[0] == 1:
        return quantified_matrix.iloc[0].astype(float)

    return quantified_matrix.median(axis=0, skipna=False).astype(float)


def _quantified_sites(
    substrates: Sequence[str],
    observed_sites: set[str],
) -> list[str]:
    quantified: list[str] = []
    seen: set[str] = set()

    for substrate in substrates:
        if substrate in observed_sites and substrate not in seen:
            quantified.append(substrate)
            seen.add(substrate)

    return quantified
