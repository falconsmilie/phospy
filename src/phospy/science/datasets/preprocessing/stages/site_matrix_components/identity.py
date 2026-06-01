from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.science.datasets.preprocessing.models import DuplicateSiteResolutionResult
from phospy.science.datasets.preprocessing.policy_models import (
    SiteMatrixMissingDataPolicy,
)

_SITE_SEQUENCE_COLUMN = "site_sequence"


@dataclass(frozen=True, slots=True)
class SequenceSupportFilterResult:
    phospho: pd.DataFrame
    site_metadata: pd.DataFrame
    scientific_row_key: pd.Series
    dropped_row_count: int
    dropped_rows: tuple[tuple[str, str], ...]

    @property
    def constructed_site_id(self) -> pd.Series:
        return self.scientific_row_key


@dataclass(frozen=True, slots=True)
class MissingDataSiteFilterResult:
    phospho: pd.DataFrame
    dropped_row_count: int
    required_observed_count: int
    dropped_rows: tuple[tuple[str, str, int], ...]


@dataclass(frozen=True, slots=True)
class SiteMatrixAssemblyResult:
    phospho: pd.DataFrame
    site_metadata: pd.DataFrame
    dropped_missing_sequence_row_ids: tuple[str, ...]
    dropped_incomplete_row_ids: tuple[str, ...]
    duplicate_dropped_row_ids: tuple[str, ...]
    dropped_row_ids: tuple[str, ...]


class SequenceSupportFilter:
    """Drop rows lacking usable sequence support and preserve row-level details."""

    def filter(
        self,
        *,
        phospho: pd.DataFrame,
        site_metadata: pd.DataFrame,
        scientific_row_key: pd.Series | None = None,
        constructed_site_id: pd.Series | None = None,
    ) -> SequenceSupportFilterResult:
        if scientific_row_key is None:
            if constructed_site_id is None:
                raise PhosPyInputError(
                    "sequence-support filtering requires scientific_row_key "
                    "(or legacy constructed_site_id)"
                )
            scientific_row_key = constructed_site_id
        if _SITE_SEQUENCE_COLUMN in site_metadata.columns:
            site_sequence = _resolve_optional_string_column(
                site_metadata,
                column_name=_SITE_SEQUENCE_COLUMN,
            )
        else:
            site_sequence = pd.Series(
                pd.NA,
                index=site_metadata.index.copy(),
                dtype="string",
                name=_SITE_SEQUENCE_COLUMN,
            )
        has_sequence = site_sequence.notna()
        dropped_rows = tuple(
            (str(row_id), str(site_id))
            for row_id, site_id in zip(
                phospho.index[~has_sequence].tolist(),
                scientific_row_key.loc[~has_sequence].astype(str).tolist(),
                strict=True,
            )
        )
        return SequenceSupportFilterResult(
            phospho=phospho.loc[has_sequence],
            site_metadata=site_metadata.loc[has_sequence],
            scientific_row_key=scientific_row_key.loc[has_sequence],
            dropped_row_count=int((~has_sequence).sum()),
            dropped_rows=dropped_rows,
        )


class MissingDataSiteFilter:
    """Apply site-matrix missing-data policy and retain row-level drop diagnostics."""

    def filter(
        self,
        *,
        phospho: pd.DataFrame,
        scientific_row_key: pd.Series | None = None,
        constructed_site_id: pd.Series | None = None,
        missing_data_policy: SiteMatrixMissingDataPolicy | str,
        minimum_observed_values: int | None,
    ) -> MissingDataSiteFilterResult:
        if scientific_row_key is None:
            if constructed_site_id is None:
                raise PhosPyInputError(
                    "missing-data site filtering requires scientific_row_key "
                    "(or legacy constructed_site_id)"
                )
            scientific_row_key = constructed_site_id
        resolved_policy = SiteMatrixMissingDataPolicy.parse(
            missing_data_policy,
            field_name="site_matrix.missing_data_policy",
        )
        observed_counts = phospho.notna().sum(axis=1)
        if resolved_policy is SiteMatrixMissingDataPolicy.DROP_ANY_MISSING:
            retained_mask = phospho.notna().all(axis=1)
            required_observed_count = phospho.shape[1]
        elif resolved_policy is SiteMatrixMissingDataPolicy.RETAIN_MISSING:
            retained_mask = pd.Series(True, index=phospho.index)
            required_observed_count = 0
        elif resolved_policy is SiteMatrixMissingDataPolicy.REQUIRE_MIN_OBSERVED_VALUES:
            if minimum_observed_values is None:
                raise PhosPyInputError(
                    "dataset build request preprocessing site-matrix construction "
                    "requires minimum_observed_values when "
                    "site_matrix.missing_data_policy='require_min_observed_values'"
                )
            if minimum_observed_values > phospho.shape[1]:
                raise PhosPyInputError(
                    "dataset build request preprocessing site-matrix construction "
                    "minimum_observed_values cannot exceed phospho sample count "
                    f"({phospho.shape[1]})"
                )
            required_observed_count = minimum_observed_values
            retained_mask = phospho.notna().sum(axis=1) >= required_observed_count
        else:
            raise PhosPyInputError(
                "dataset build request preprocessing_config contains an unsupported "
                "site_matrix.missing_data_policy"
            )

        filtered = phospho.loc[retained_mask]
        dropped_rows = int(len(phospho.index) - len(filtered.index))
        dropped_row_details = tuple(
            (str(row_id), str(site_id), int(observed_value_count))
            for row_id, site_id, observed_value_count in zip(
                phospho.index[~retained_mask].tolist(),
                scientific_row_key.loc[~retained_mask].astype(str).tolist(),
                observed_counts.loc[~retained_mask].tolist(),
                strict=True,
            )
        )
        return MissingDataSiteFilterResult(
            phospho=filtered,
            dropped_row_count=dropped_rows,
            required_observed_count=required_observed_count,
            dropped_rows=dropped_row_details,
        )


class SiteMatrixAssembler:
    """Assemble final site-matrix phospho/metadata tables and dropped row ids."""

    def assemble(
        self,
        *,
        duplicate_site_result: DuplicateSiteResolutionResult,
        output_index_name: str | None,
        dropped_missing_sequence_rows: tuple[tuple[str, str], ...],
        dropped_incomplete_rows: tuple[tuple[str, str, int], ...],
    ) -> SiteMatrixAssemblyResult:
        final_phospho = duplicate_site_result.phospho.sort_index(kind="stable")
        final_site_metadata = duplicate_site_result.site_metadata.reindex(
            final_phospho.index
        )
        if "display_id" in final_site_metadata.columns:
            display_ids = (
                final_site_metadata.loc[:, "display_id"].astype("string").str.strip()
            )
            if bool(display_ids.notna().all()) and bool((display_ids != "").all()):
                if not bool(display_ids.duplicated().any()):
                    canonical_display_index = pd.Index(
                        display_ids.astype(str).tolist(),
                        name=output_index_name,
                    )
                    final_phospho.index = canonical_display_index
                    final_site_metadata.index = canonical_display_index.copy()
                    if "site_key" in final_site_metadata.columns:
                        final_site_metadata.loc[:, "site_key"] = (
                            canonical_display_index.astype(str).tolist()
                        )
        final_site_index = pd.Index(
            final_phospho.index.tolist(), name=output_index_name
        )
        final_phospho.index = final_site_index
        final_site_metadata.index = final_site_index.copy()

        duplicate_dropped_row_ids = tuple(
            str(row_id)
            for row_id in duplicate_site_result.duplicate_site_resolution.loc[
                ~duplicate_site_result.duplicate_site_resolution.loc[:, "retained"],
                "source_row_id",
            ]
            .astype(str)
            .tolist()
        )
        dropped_missing_sequence_row_ids = tuple(
            row_id for row_id, _ in dropped_missing_sequence_rows
        )
        dropped_incomplete_row_ids = tuple(
            row_id for row_id, _, _ in dropped_incomplete_rows
        )
        dropped_row_ids = _unique_strings_preserve_order(
            (
                *dropped_missing_sequence_row_ids,
                *dropped_incomplete_row_ids,
                *duplicate_dropped_row_ids,
            )
        )
        return SiteMatrixAssemblyResult(
            phospho=final_phospho,
            site_metadata=final_site_metadata,
            dropped_missing_sequence_row_ids=dropped_missing_sequence_row_ids,
            dropped_incomplete_row_ids=dropped_incomplete_row_ids,
            duplicate_dropped_row_ids=duplicate_dropped_row_ids,
            dropped_row_ids=dropped_row_ids,
        )


def _resolve_optional_string_column(
    site_metadata: pd.DataFrame,
    *,
    column_name: str,
) -> pd.Series:
    column = site_metadata.loc[:, column_name]
    normalized = column.astype("string").str.strip()
    missing_mask = column.isna() | normalized.isna() | (normalized == "")
    return normalized.where(~missing_mask, other=pd.NA)


def _unique_strings_preserve_order(values: tuple[str, ...]) -> tuple[str, ...]:
    seen: set[str] = set()
    ordered: list[str] = []
    for value in values:
        if value in seen:
            continue
        seen.add(value)
        ordered.append(value)
    return tuple(ordered)
