from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Literal, TypeAlias

import pandas as pd

from .validation.errors import (
    PhospyValidationError,
)

AggregationMethod: TypeAlias = Literal["median"]


@dataclass(slots=True)
class KinaseProfileResult:
    profile_matrix: pd.DataFrame
    substrate_counts: pd.Series
    quantified_substrates: dict[str, list[str]]


class KinaseProfileBuilder:
    """Construct kinase substrate profiles from known substrate annotations."""

    def __init__(self, aggregation: AggregationMethod = "median") -> None:
        _validate_aggregation(aggregation)
        self.aggregation = aggregation

    def build(
        self,
        substrate_map: Mapping[str, Sequence[str]],
        phospho_matrix: pd.DataFrame,
        min_substrates: int = 1,
    ) -> KinaseProfileResult:
        return build_kinase_substrate_profiles(
            substrate_map=substrate_map,
            phospho_matrix=phospho_matrix,
            min_substrates=min_substrates,
            aggregation=self.aggregation,
        )


def build_kinase_substrate_profiles(
    substrate_map: Mapping[str, Sequence[str]],
    phospho_matrix: pd.DataFrame,
    min_substrates: int = 1,
    aggregation: AggregationMethod = "median",
) -> KinaseProfileResult:
    """Build kinase substrate profiles from quantified phosphosite values.

    This mirrors the core behaviour of PhosR's ``kinaseSubstrateProfile()``:
    for each kinase, intersect its known substrates with the phosphosite matrix,
    use the single quantified row directly when exactly one site is available,
    and otherwise summarise quantified substrates column-wise with the median.
    """

    _validate_positive_int(min_substrates, name="min_substrates")
    _validate_aggregation(aggregation)

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
            aggregation=aggregation,
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


def _aggregate_quantified_sites(
    quantified_matrix: pd.DataFrame,
    aggregation: AggregationMethod,
) -> pd.Series:
    if quantified_matrix.shape[0] == 1:
        return quantified_matrix.iloc[0].astype(float)

    if aggregation == "median":
        return quantified_matrix.apply(
            lambda column: column.median(skipna=False),
            axis=0,
        ).astype(float)

    raise ValueError(f"Unsupported aggregation: {aggregation}")


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


def _validate_positive_int(value: int, name: str) -> None:
    if value < 1:
        raise PhospyValidationError(f"{name} must be at least 1")


def _validate_aggregation(aggregation: AggregationMethod) -> None:
    if aggregation != "median":
        raise PhospyValidationError("aggregation must be 'median'")
