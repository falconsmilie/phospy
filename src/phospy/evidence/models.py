"""Peptide-level evidence models for phosphosite workflows."""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass

import pandas as pd

from phospy._frame_ownership import (
    export_dataframe,
    export_series,
    own_dataframe,
    own_series,
)
from phospy.errors.input import PhosPyInputError
from phospy.evidence.multi_site import (
    MultiSiteHandlingConfig,
    MultiSiteObservation,
    build_multi_site_observation,
    resolve_site_mapping_frame,
)
from phospy.sites.identifiers import canonicalize_site_identifier
from phospy.validation.common.dataframes import (
    require_columns,
    require_dataframe,
    require_finite_numeric_dataframe,
    require_no_duplicate_labels,
    require_non_empty_dataframe,
    require_unique_columns,
    require_unique_row_pairs,
)

_REQUIRED_EVIDENCE_COLUMNS: tuple[str, ...] = (
    "peptide_row_id",
    "site_id",
    "unique_feature_id",
    "gene_symbol",
    "protein_accession",
    "site_string",
    "peptide_sequence",
    "modified_peptide_sequence",
    "multi_site",
    "provenance_source",
)
_OPTIONAL_EVIDENCE_COLUMNS: tuple[str, ...] = (
    "site_sequence",
    "localisation_confidence",
    "missingness_flags",
    "imputation_flags",
)


@dataclass(frozen=True, slots=True, init=False)
class PeptideEvidenceRecord:
    """One peptide-level evidence record aligned to one peptide-table row."""

    peptide_row_id: str
    site_id: str | None
    unique_feature_id: str
    gene_symbol: str
    protein_accession: str
    site_string: str
    sample_intensities: pd.Series
    peptide_sequence: str
    modified_peptide_sequence: str
    site_sequence: str | None
    localisation_confidence: float | None
    missingness_flags: tuple[str, ...]
    imputation_flags: tuple[str, ...]
    multi_site: bool
    provenance_source: str

    def __init__(
        self,
        *,
        peptide_row_id: str,
        site_id: str | None,
        unique_feature_id: str,
        gene_symbol: str,
        protein_accession: str,
        site_string: str,
        sample_intensities: pd.Series,
        peptide_sequence: str,
        modified_peptide_sequence: str,
        site_sequence: str | None = None,
        localisation_confidence: float | None = None,
        missingness_flags: Sequence[str] = (),
        imputation_flags: Sequence[str] = (),
        multi_site: bool,
        provenance_source: str,
        _assume_owned: bool = False,
    ) -> None:
        peptide_row_id_value = _canonical_non_empty_string(
            peptide_row_id,
            field_name="peptide_evidence_record.peptide_row_id",
        )
        site_id_value = _canonical_optional_site_id(
            site_id,
            field_name="peptide_evidence_record.site_id",
        )
        unique_feature_id_value = _canonical_non_empty_string(
            unique_feature_id,
            field_name="peptide_evidence_record.unique_feature_id",
        )
        gene_symbol_value = _canonical_non_empty_string(
            gene_symbol,
            field_name="peptide_evidence_record.gene_symbol",
        )
        protein_accession_value = _canonical_non_empty_string(
            protein_accession,
            field_name="peptide_evidence_record.protein_accession",
        )
        site_string_value = _canonical_non_empty_string(
            site_string,
            field_name="peptide_evidence_record.site_string",
        )
        peptide_sequence_value = _canonical_non_empty_string(
            peptide_sequence,
            field_name="peptide_evidence_record.peptide_sequence",
        )
        modified_peptide_sequence_value = _canonical_non_empty_string(
            modified_peptide_sequence,
            field_name="peptide_evidence_record.modified_peptide_sequence",
        )
        provenance_source_value = _canonical_non_empty_string(
            provenance_source,
            field_name="peptide_evidence_record.provenance_source",
        )
        site_sequence_value = _canonical_optional_string(
            site_sequence,
            field_name="peptide_evidence_record.site_sequence",
        )
        localisation_confidence_value = _canonical_optional_float(
            localisation_confidence,
            field_name="peptide_evidence_record.localisation_confidence",
            minimum=0.0,
            maximum=1.0,
        )
        missingness_flags_value = _canonical_flag_values(
            missingness_flags,
            field_name="peptide_evidence_record.missingness_flags",
        )
        imputation_flags_value = _canonical_flag_values(
            imputation_flags,
            field_name="peptide_evidence_record.imputation_flags",
        )
        if not isinstance(multi_site, bool):
            raise PhosPyInputError("peptide_evidence_record.multi_site must be a bool")

        sample_intensity_series = own_series(
            sample_intensities,
            field_name="peptide_evidence_record.sample_intensities",
            error_type=PhosPyInputError,
            assume_owned=_assume_owned,
        )
        _validate_sample_intensity_series(
            sample_intensity_series,
            field_name="peptide_evidence_record.sample_intensities",
        )

        object.__setattr__(self, "peptide_row_id", peptide_row_id_value)
        object.__setattr__(self, "site_id", site_id_value)
        object.__setattr__(self, "unique_feature_id", unique_feature_id_value)
        object.__setattr__(self, "gene_symbol", gene_symbol_value)
        object.__setattr__(self, "protein_accession", protein_accession_value)
        object.__setattr__(self, "site_string", site_string_value)
        object.__setattr__(self, "sample_intensities", sample_intensity_series)
        object.__setattr__(self, "peptide_sequence", peptide_sequence_value)
        object.__setattr__(
            self,
            "modified_peptide_sequence",
            modified_peptide_sequence_value,
        )
        object.__setattr__(self, "site_sequence", site_sequence_value)
        object.__setattr__(
            self, "localisation_confidence", localisation_confidence_value
        )
        object.__setattr__(self, "missingness_flags", missingness_flags_value)
        object.__setattr__(self, "imputation_flags", imputation_flags_value)
        object.__setattr__(self, "multi_site", multi_site)
        object.__setattr__(self, "provenance_source", provenance_source_value)

    def sample_intensity_series(self) -> pd.Series:
        return export_series(self.sample_intensities)


@dataclass(frozen=True, slots=True, init=False)
class SiteEvidenceMapping:
    """Validated peptide-row to canonical-site mapping table."""

    _frame: pd.DataFrame

    def __init__(self, frame: pd.DataFrame, *, _assume_owned: bool = False) -> None:
        frame = own_dataframe(
            frame,
            field_name="site_evidence_mapping",
            error_type=PhosPyInputError,
            assume_owned=_assume_owned,
        )
        require_dataframe(
            frame,
            field_name="site_evidence_mapping",
            allow_empty=True,
            error_type=PhosPyInputError,
        )
        require_unique_columns(
            frame,
            field_name="site_evidence_mapping",
            error_type=PhosPyInputError,
        )
        require_columns(
            frame,
            field_name="site_evidence_mapping",
            required_columns=("peptide_row_id", "site_id"),
            error_type=PhosPyInputError,
        )
        canonical = frame.loc[:, ["peptide_row_id", "site_id"]].copy(deep=True)
        if not canonical.empty:
            canonical.loc[:, "peptide_row_id"] = canonical.loc[:, "peptide_row_id"].map(
                lambda value: _canonical_non_empty_string(
                    value,
                    field_name="site_evidence_mapping.peptide_row_id",
                )
            )
            canonical.loc[:, "site_id"] = canonical.loc[:, "site_id"].map(
                lambda value: _canonical_required_site_id(
                    value,
                    field_name="site_evidence_mapping.site_id",
                )
            )
            require_unique_row_pairs(
                canonical,
                field_name="site_evidence_mapping",
                column_names=("peptide_row_id", "site_id"),
                error_type=PhosPyInputError,
            )
        object.__setattr__(self, "_frame", canonical)

    def to_dataframe(self) -> pd.DataFrame:
        return export_dataframe(self._frame)

    def validate_peptide_row_ids(
        self,
        peptide_row_ids: Iterable[str],
        *,
        field_name: str,
    ) -> None:
        allowed = {
            _canonical_non_empty_string(value, field_name=field_name)
            for value in peptide_row_ids
        }
        mapped = set(self._frame.loc[:, "peptide_row_id"].tolist())
        unknown = sorted(mapped.difference(allowed))
        if unknown:
            preview = ", ".join(repr(value) for value in unknown[:5])
            suffix = "" if len(unknown) <= 5 else " ..."
            raise PhosPyInputError(
                f"site_evidence_mapping contains peptide_row_id values not present in "
                f"{field_name}: {preview}{suffix}"
            )

    def site_ids_for_peptide(self, peptide_row_id: str) -> tuple[str, ...]:
        peptide_row_id_value = _canonical_non_empty_string(
            peptide_row_id,
            field_name="site_evidence_mapping.peptide_row_id",
        )
        site_ids = self._frame.loc[
            self._frame.loc[:, "peptide_row_id"] == peptide_row_id_value,
            "site_id",
        ].tolist()
        return tuple(site_ids)


@dataclass(frozen=True, slots=True, init=False)
class PeptideEvidenceTable:
    """Owned peptide-level evidence table with validated site mapping."""

    _frame: pd.DataFrame
    _sample_intensity_columns: tuple[str, ...]
    _site_mapping: SiteEvidenceMapping
    _multi_site_observations: tuple[MultiSiteObservation, ...]
    _multi_site_handling_config: MultiSiteHandlingConfig

    def __init__(
        self,
        *,
        frame: pd.DataFrame,
        sample_intensity_columns: Sequence[str],
        site_mapping: SiteEvidenceMapping | pd.DataFrame | None = None,
        multi_site_handling_config: MultiSiteHandlingConfig | None = None,
        _assume_owned: bool = False,
    ) -> None:
        frame = own_dataframe(
            frame,
            field_name="peptide_evidence_table",
            error_type=PhosPyInputError,
            assume_owned=_assume_owned,
        )
        require_dataframe(
            frame,
            field_name="peptide_evidence_table",
            allow_empty=False,
            error_type=PhosPyInputError,
        )
        require_non_empty_dataframe(
            frame,
            field_name="peptide_evidence_table",
            error_type=PhosPyInputError,
        )
        require_unique_columns(
            frame,
            field_name="peptide_evidence_table",
            error_type=PhosPyInputError,
        )
        require_columns(
            frame,
            field_name="peptide_evidence_table",
            required_columns=_REQUIRED_EVIDENCE_COLUMNS,
            error_type=PhosPyInputError,
        )

        sample_columns = _resolve_sample_intensity_columns(
            frame=frame,
            sample_intensity_columns=sample_intensity_columns,
        )
        canonical = frame.copy(deep=True)
        for column_name in _REQUIRED_EVIDENCE_COLUMNS:
            if column_name == "site_id":
                site_field_name = f"peptide_evidence_table.{column_name}"
                canonical[column_name] = canonical.loc[:, column_name].map(
                    lambda value, _site_field_name=site_field_name: (
                        _canonical_optional_site_id(
                            value,
                            field_name=_site_field_name,
                        )
                    )
                )
                continue
            if column_name == "multi_site":
                canonical[column_name] = canonical.loc[:, column_name].map(
                    _validate_bool_value
                )
                continue
            field_name = f"peptide_evidence_table.{column_name}"
            canonical[column_name] = canonical.loc[:, column_name].map(
                lambda value, _field_name=field_name: _canonical_non_empty_string(
                    value,
                    field_name=_field_name,
                )
            )
        if canonical.loc[:, "peptide_row_id"].duplicated().any():
            duplicate_ids = (
                canonical.loc[canonical.loc[:, "peptide_row_id"].duplicated(keep=False)]
                .loc[:, "peptide_row_id"]
                .drop_duplicates()
                .tolist()
            )
            preview = ", ".join(repr(value) for value in duplicate_ids[:5])
            suffix = "" if len(duplicate_ids) <= 5 else " ..."
            raise PhosPyInputError(
                "peptide_evidence_table.peptide_row_id must be unique; "
                f"duplicates={preview}{suffix}"
            )
        for optional_column in _OPTIONAL_EVIDENCE_COLUMNS:
            if optional_column not in canonical.columns:
                continue
            if optional_column == "localisation_confidence":
                canonical[optional_column] = canonical.loc[:, optional_column].map(
                    lambda value: _canonical_optional_float(
                        value,
                        field_name="peptide_evidence_table.localisation_confidence",
                        minimum=0.0,
                        maximum=1.0,
                    )
                )
                continue
            if optional_column in ("missingness_flags", "imputation_flags"):
                field_name = f"peptide_evidence_table.{optional_column}"
                canonical[optional_column] = canonical.loc[:, optional_column].map(
                    lambda value, _field_name=field_name: _parse_optional_flag_cell(
                        value,
                        field_name=_field_name,
                    )
                )
                continue
            field_name = f"peptide_evidence_table.{optional_column}"
            canonical[optional_column] = canonical.loc[:, optional_column].map(
                lambda value, _field_name=field_name: _canonical_optional_string(
                    value,
                    field_name=_field_name,
                )
            )

        sample_view = canonical.loc[:, list(sample_columns)]
        _validate_sample_intensity_frame(sample_view)

        multi_site_config = (
            MultiSiteHandlingConfig()
            if multi_site_handling_config is None
            else multi_site_handling_config
        )
        if not isinstance(multi_site_config, MultiSiteHandlingConfig):
            raise PhosPyInputError(
                "peptide_evidence_table.multi_site_handling_config must be "
                "a MultiSiteHandlingConfig or None"
            )
        observations = _build_multi_site_observations(canonical)

        mapping = _resolve_mapping(
            canonical,
            observations=observations,
            site_mapping=site_mapping,
            multi_site_handling_config=multi_site_config,
        )
        mapping.validate_peptide_row_ids(
            canonical.loc[:, "peptide_row_id"].tolist(),
            field_name="peptide_evidence_table.peptide_row_id",
        )
        if site_mapping is not None:
            _validate_default_site_links(
                canonical=canonical,
                mapping=mapping,
            )
        object.__setattr__(self, "_frame", canonical)
        object.__setattr__(self, "_sample_intensity_columns", sample_columns)
        object.__setattr__(self, "_site_mapping", mapping)
        object.__setattr__(self, "_multi_site_observations", observations)
        object.__setattr__(self, "_multi_site_handling_config", multi_site_config)

    @property
    def sample_intensity_columns(self) -> tuple[str, ...]:
        return self._sample_intensity_columns

    @property
    def site_mapping(self) -> SiteEvidenceMapping:
        return self._site_mapping

    @property
    def multi_site_handling_config(self) -> MultiSiteHandlingConfig:
        return self._multi_site_handling_config

    def multi_site_observations(self) -> tuple[MultiSiteObservation, ...]:
        return self._multi_site_observations

    def kinase_sequence_site_mapping(self) -> pd.DataFrame:
        resolved = resolve_site_mapping_frame(
            observations=self._multi_site_observations,
            policy=self._multi_site_handling_config.kinase_sequence_scoring_policy,
        )
        return export_dataframe(resolved)

    def multi_site_policy_provenance(self) -> dict[str, object]:
        return {
            "statistical_modeling_policy": (
                self._multi_site_handling_config.statistical_modeling_policy
            ),
            "kinase_sequence_scoring_policy": (
                self._multi_site_handling_config.kinase_sequence_scoring_policy
            ),
            "multi_site_rows": int(
                sum(1 for row in self._multi_site_observations if row.is_multi_site)
            ),
            "single_site_rows": int(
                sum(1 for row in self._multi_site_observations if not row.is_multi_site)
            ),
        }

    def to_dataframe(self) -> pd.DataFrame:
        return export_dataframe(self._frame)

    def records(self) -> tuple[PeptideEvidenceRecord, ...]:
        records: list[PeptideEvidenceRecord] = []
        for _, row in self._frame.iterrows():
            sample_intensities = row.loc[list(self._sample_intensity_columns)].copy()
            missingness_flags = _parse_optional_flag_cell(
                row.get("missingness_flags", None),
                field_name="peptide_evidence_table.missingness_flags",
            )
            imputation_flags = _parse_optional_flag_cell(
                row.get("imputation_flags", None),
                field_name="peptide_evidence_table.imputation_flags",
            )
            records.append(
                PeptideEvidenceRecord(
                    peptide_row_id=str(row.loc["peptide_row_id"]),
                    site_id=(
                        None if pd.isna(row.loc["site_id"]) else str(row.loc["site_id"])
                    ),
                    unique_feature_id=str(row.loc["unique_feature_id"]),
                    gene_symbol=str(row.loc["gene_symbol"]),
                    protein_accession=str(row.loc["protein_accession"]),
                    site_string=str(row.loc["site_string"]),
                    sample_intensities=sample_intensities,
                    peptide_sequence=str(row.loc["peptide_sequence"]),
                    modified_peptide_sequence=str(row.loc["modified_peptide_sequence"]),
                    site_sequence=(
                        None
                        if "site_sequence" not in self._frame.columns
                        or pd.isna(row.get("site_sequence", None))
                        else str(row.loc["site_sequence"])
                    ),
                    localisation_confidence=(
                        None
                        if "localisation_confidence" not in self._frame.columns
                        or pd.isna(row.get("localisation_confidence", None))
                        else float(row.loc["localisation_confidence"])
                    ),
                    missingness_flags=missingness_flags,
                    imputation_flags=imputation_flags,
                    multi_site=bool(row.loc["multi_site"]),
                    provenance_source=str(row.loc["provenance_source"]),
                    _assume_owned=False,
                )
            )
        return tuple(records)


def _resolve_sample_intensity_columns(
    *,
    frame: pd.DataFrame,
    sample_intensity_columns: Sequence[str],
) -> tuple[str, ...]:
    if isinstance(sample_intensity_columns, str):
        raise PhosPyInputError(
            "peptide_evidence_table.sample_intensity_columns must be a sequence of "
            "column names"
        )
    columns = tuple(sample_intensity_columns)
    if not columns:
        raise PhosPyInputError(
            "peptide_evidence_table.sample_intensity_columns must contain at least "
            "one column name"
        )
    require_no_duplicate_labels(
        pd.Index(columns),
        field_name="peptide_evidence_table.sample_intensity_columns",
        error_type=PhosPyInputError,
    )
    missing = [column for column in columns if column not in frame.columns]
    if missing:
        joined = ", ".join(missing)
        raise PhosPyInputError(
            "peptide_evidence_table.sample_intensity_columns includes missing "
            f"columns: {joined}"
        )
    reserved = set(_REQUIRED_EVIDENCE_COLUMNS + _OPTIONAL_EVIDENCE_COLUMNS)
    colliding = [column for column in columns if column in reserved]
    if colliding:
        joined = ", ".join(colliding)
        raise PhosPyInputError(
            "peptide_evidence_table.sample_intensity_columns cannot include reserved "
            f"metadata columns: {joined}"
        )
    return columns


def _validate_sample_intensity_frame(frame: pd.DataFrame) -> None:
    if any(pd.api.types.is_bool_dtype(frame[column]) for column in frame.columns):
        boolean_columns = [
            str(column)
            for column in frame.columns
            if pd.api.types.is_bool_dtype(frame[column])
        ]
        joined = ", ".join(boolean_columns)
        raise PhosPyInputError(
            "peptide_evidence_table sample intensity columns must be numeric; "
            f"boolean columns are invalid: {joined}"
        )
    non_numeric = [
        str(column)
        for column in frame.columns
        if not pd.api.types.is_numeric_dtype(frame[column])
    ]
    if non_numeric:
        joined = ", ".join(non_numeric)
        raise PhosPyInputError(
            "peptide_evidence_table sample intensity columns must be numeric; "
            f"non-numeric columns: {joined}"
        )
    require_finite_numeric_dataframe(
        frame,
        field_name="peptide_evidence_table.sample_intensities",
        error_type=PhosPyInputError,
        allow_missing=True,
    )


def _validate_sample_intensity_series(series: pd.Series, *, field_name: str) -> None:
    require_no_duplicate_labels(
        pd.Index(series.index),
        field_name=f"{field_name}.index",
        error_type=PhosPyInputError,
    )
    if any(
        not isinstance(label, str) or not label.strip()
        for label in series.index.tolist()
    ):
        raise PhosPyInputError(
            f"{field_name}.index must contain non-empty string sample labels"
        )
    sample_frame = pd.DataFrame([series.to_dict()])
    _validate_sample_intensity_frame(sample_frame)


def _resolve_mapping(
    canonical: pd.DataFrame,
    *,
    observations: tuple[MultiSiteObservation, ...],
    site_mapping: SiteEvidenceMapping | pd.DataFrame | None,
    multi_site_handling_config: MultiSiteHandlingConfig,
) -> SiteEvidenceMapping:
    if site_mapping is None:
        resolved = resolve_site_mapping_frame(
            observations=observations,
            policy=multi_site_handling_config.statistical_modeling_policy,
        )
        default_linkable = set(
            canonical.loc[canonical.loc[:, "site_id"].notna(), "peptide_row_id"].astype(
                str
            )
        )
        filtered = resolved.loc[
            resolved.loc[:, "peptide_row_id"].astype(str).isin(default_linkable),
            :,
        ].copy(deep=True)
        resolved_mapping = filtered.loc[:, ["peptide_row_id", "site_id"]].copy(
            deep=True
        )
        return SiteEvidenceMapping(resolved_mapping, _assume_owned=True)
    if isinstance(site_mapping, pd.DataFrame):
        return SiteEvidenceMapping(site_mapping)
    if isinstance(site_mapping, SiteEvidenceMapping):
        return site_mapping
    raise PhosPyInputError(
        "peptide_evidence_table.site_mapping must be a SiteEvidenceMapping, "
        "pandas DataFrame, or None"
    )


def _validate_default_site_links(
    *,
    canonical: pd.DataFrame,
    mapping: SiteEvidenceMapping,
) -> None:
    expected = canonical.loc[
        canonical.loc[:, "site_id"].notna(),
        ["peptide_row_id", "site_id"],
    ].drop_duplicates()
    if expected.empty:
        return
    merged = expected.merge(
        mapping._frame,
        how="left",
        on=["peptide_row_id", "site_id"],
        indicator=True,
    )
    missing_pairs = merged.loc[merged.loc[:, "_merge"] != "both", :]
    if missing_pairs.empty:
        return
    preview_rows = (
        missing_pairs.loc[:, ["peptide_row_id", "site_id"]]
        .drop_duplicates()
        .itertuples(index=False, name=None)
    )
    preview_list = list(preview_rows)
    preview = ", ".join(repr(pair) for pair in preview_list[:5])
    suffix = "" if len(preview_list) <= 5 else " ..."
    raise PhosPyInputError(
        "peptide_evidence_table.site_mapping must include each non-missing default "
        f"(peptide_row_id, site_id) pair from peptide_evidence_table: {preview}{suffix}"
    )


def _build_multi_site_observations(
    canonical: pd.DataFrame,
) -> tuple[MultiSiteObservation, ...]:
    observations: list[MultiSiteObservation] = []
    for _, row in canonical.iterrows():
        observations.append(
            build_multi_site_observation(
                peptide_row_id=str(row.loc["peptide_row_id"]),
                gene_symbol=str(row.loc["gene_symbol"]),
                site_string=str(row.loc["site_string"]),
                declared_multi_site=bool(row.loc["multi_site"]),
                field_name=("peptide_evidence_table.site_string"),
            )
        )
    return tuple(observations)


def _validate_bool_value(value: object) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)) and value in (0, 1):
        raise PhosPyInputError(
            "peptide_evidence_table.multi_site must use boolean values, not numeric "
            "0/1 encodings"
        )
    raise PhosPyInputError("peptide_evidence_table.multi_site must be a bool column")


def _canonical_non_empty_string(value: object, *, field_name: str) -> str:
    if _is_missing(value):
        raise PhosPyInputError(f"{field_name} must not contain missing values")
    if not isinstance(value, str):
        raise PhosPyInputError(f"{field_name} must contain non-empty string values")
    stripped = value.strip()
    if stripped == "":
        raise PhosPyInputError(f"{field_name} must contain non-empty string values")
    return stripped


def _canonical_optional_string(value: object, *, field_name: str) -> str | None:
    if _is_missing(value):
        return None
    if not isinstance(value, str):
        raise PhosPyInputError(f"{field_name} must be a string when provided")
    stripped = value.strip()
    if stripped == "":
        return None
    return stripped


def _canonical_optional_float(
    value: object,
    *,
    field_name: str,
    minimum: float | None = None,
    maximum: float | None = None,
) -> float | None:
    if _is_missing(value):
        return None
    if isinstance(value, bool) or not isinstance(value, int | float):
        raise PhosPyInputError(f"{field_name} must be numeric when provided")
    resolved = float(value)
    if resolved in (float("inf"), float("-inf")):
        raise PhosPyInputError(f"{field_name} must be finite when provided")
    if minimum is not None and resolved < minimum:
        raise PhosPyInputError(
            f"{field_name} must be greater than or equal to {minimum} when provided"
        )
    if maximum is not None and resolved > maximum:
        raise PhosPyInputError(
            f"{field_name} must be less than or equal to {maximum} when provided"
        )
    return resolved


def _canonical_optional_site_id(value: object, *, field_name: str) -> str | None:
    if _is_missing(value):
        return None
    if isinstance(value, str) and value.strip() == "":
        return None
    return canonicalize_site_identifier(
        value,
        field_name=field_name,
        error_type=PhosPyInputError,
    )


def _canonical_required_site_id(value: object, *, field_name: str) -> str:
    canonical = _canonical_optional_site_id(value, field_name=field_name)
    if canonical is None:
        raise PhosPyInputError(
            f"{field_name} must not contain missing site identifiers"
        )
    return canonical


def _canonical_flag_values(
    values: Sequence[str],
    *,
    field_name: str,
) -> tuple[str, ...]:
    return tuple(
        _canonical_non_empty_string(value, field_name=field_name) for value in values
    )


def _parse_optional_flag_cell(value: object, *, field_name: str) -> tuple[str, ...]:
    if _is_missing(value):
        return ()
    if isinstance(value, str):
        if value.strip() == "":
            return ()
        return tuple(
            _canonical_non_empty_string(token, field_name=field_name)
            for token in value.split("|")
            if token.strip()
        )
    if isinstance(value, Sequence) and not isinstance(value, str):
        return _canonical_flag_values(tuple(value), field_name=field_name)
    raise PhosPyInputError(
        f"{field_name} must be either a '|' delimited string, a sequence of strings, "
        "or missing"
    )


def _is_missing(value: object) -> bool:
    try:
        return bool(pd.isna(value))
    except (TypeError, ValueError):
        return False
