"""Foundation classes for upstream phosphosite table importers."""

from __future__ import annotations

import math
from os import PathLike
from pathlib import Path

import numpy as np
import pandas as pd

from phospy.contracts.requests import PhosphositeImportRequest
from phospy.contracts.results import (
    IMPORTER_QUALITY_STATUS_NOT_APPLICABLE,
    IMPORTER_QUALITY_STATUS_REPORTED,
    ImporterDetectedIntensityColumn,
    ImporterDuplicateKeySummary,
    ImporterFlaggedRowSummary,
    ImporterLocalisationConfidenceSummary,
    ImporterMissingIntensitySummary,
    ImporterQualityCount,
    ImporterQualityReport,
    PhosphositeImportResult,
)
from phospy.errors.input import PhosPyInputError, UnsupportedInputFormatError
from phospy.io.readers.tables import (
    supported_table_input_formats,
    table_format_from_path,
)
from phospy.science.evidence.localisation import (
    LOCALISATION_CONFIDENCE_OUTPUT_COLUMN,
    normalise_localisation_confidence_series,
)
from phospy.science.evidence.multi_site import parse_phospho_site_tokens
from phospy.science.sites.identifiers import canonicalize_site_components
from phospy.validation.common.dataframes import require_dataframe
from phospy.validation.datasets.importers import (
    PhosphositeImportRequestValidator,
    normalise_sample_column_mapping,
    peptide_evidence_requested,
    require_import_source_columns,
    required_import_source_columns,
)

_CSV = "csv"
_TSV = "tsv"
_PARQUET = "parquet"
_MISSING_NUMERIC_TOKENS = frozenset({"", "na", "n/a", "nan", "null"})


class MappedPhosphositeTableImporter:
    """Translate an explicitly column-mapped upstream table into PhosPy candidates.

    This importer is a framework component, not a MaxQuant/FragPipe/Spectronaut/
    DIA-NN parser. Tool-specific importers should stay small by adapting their
    known column names into ``PhosphositeImportRequest`` and delegating shared
    normalisation here.
    """

    def __init__(
        self,
        *,
        request_validator: PhosphositeImportRequestValidator | None = None,
    ) -> None:
        self._request_validator = request_validator or (
            PhosphositeImportRequestValidator()
        )

    def run(self, request: PhosphositeImportRequest) -> PhosphositeImportResult:
        validated = self._request_validator.run(request)
        sample_mapping = normalise_sample_column_mapping(
            validated.sample_intensity_columns
        )
        source = _read_upstream_table(validated.source)
        require_dataframe(
            source,
            field_name="phosphosite import source",
            allow_empty=False,
            error_type=PhosPyInputError,
        )
        require_import_source_columns(
            source,
            required_columns=required_import_source_columns(
                validated,
                sample_mapping,
            ),
        )

        row_index = _resolve_row_index(source, validated)
        phospho = _build_phospho_candidate(
            source,
            row_index=row_index,
            sample_mapping=sample_mapping,
        )
        (
            site_metadata,
            localisation_column,
            localisation_diagnostics,
            localisation_warnings,
        ) = _build_site_metadata_candidate(
            source,
            request=validated,
            row_index=row_index,
        )
        peptide_evidence = _build_peptide_evidence_candidate(
            source,
            request=validated,
            row_index=row_index,
            phospho=phospho,
            localisation_column=localisation_column,
        )
        diagnostics, warnings = _build_import_diagnostics(
            source=source,
            site_metadata=site_metadata,
            sample_mapping=sample_mapping,
            peptide_evidence=peptide_evidence,
            localisation_diagnostics=localisation_diagnostics,
            localisation_warnings=localisation_warnings,
            source_name=validated.source_name,
        )
        quality_report = _build_import_quality_report(
            source=source,
            phospho=phospho,
            site_metadata=site_metadata,
            sample_mapping=sample_mapping,
            localisation_diagnostics=localisation_diagnostics,
            warnings=warnings,
            source_name=validated.source_name,
        )
        return PhosphositeImportResult(
            phospho_matrix_candidate=phospho,
            site_metadata_candidate=site_metadata,
            peptide_evidence=peptide_evidence,
            sample_column_mapping=sample_mapping,
            localisation_confidence_column=localisation_column,
            warnings=warnings,
            diagnostics=diagnostics,
            source_name=validated.source_name,
            quality_report=quality_report,
            _assume_owned=True,
        )


ColumnMappedPhosphositeImporter = MappedPhosphositeTableImporter


def _read_upstream_table(source: object) -> pd.DataFrame:
    if isinstance(source, pd.DataFrame):
        return source.copy(deep=True)
    if isinstance(source, str):
        if source.strip() == "":
            raise PhosPyInputError("phosphosite import source path cannot be empty")
        path = Path(source.strip())
    elif isinstance(source, (Path, PathLike)):
        path = Path(source)
    else:
        raise UnsupportedInputFormatError(
            "phosphosite import source must be a pandas DataFrame or a file path "
            "(str/pathlib.Path)"
        )
    table_format = table_format_from_path(path)
    try:
        if table_format == _CSV:
            return pd.read_csv(
                path,
                dtype=object,
                keep_default_na=False,
                na_values=[],
            )
        if table_format == _TSV:
            return pd.read_csv(
                path,
                sep="\t",
                dtype=object,
                keep_default_na=False,
                na_values=[],
            )
        if table_format == _PARQUET:
            return pd.read_parquet(path)
    except FileNotFoundError as exc:
        raise PhosPyInputError(f"input file does not exist: {path}") from exc
    except PermissionError as exc:
        raise PhosPyInputError(
            f"permission denied while reading phosphosite import source: {path}"
        ) from exc
    except ImportError as exc:
        raise UnsupportedInputFormatError(
            "parquet input requires optional parquet dependencies (for example pyarrow)"
        ) from exc
    except (
        OSError,
        UnicodeDecodeError,
        ValueError,
        pd.errors.EmptyDataError,
        pd.errors.ParserError,
    ) as exc:
        raise PhosPyInputError(
            f"failed to parse phosphosite import source '{path}': {exc}"
        ) from exc
    raise UnsupportedInputFormatError(
        f"unsupported phosphosite import source format for '{path}'. supported "
        f"formats: {supported_table_input_formats()}"
    )


def _resolve_row_index(
    source: pd.DataFrame,
    request: PhosphositeImportRequest,
) -> pd.Index:
    if request.row_id_column is None:
        return pd.Index(
            [
                f"phosphosite_import_row_{position + 1}"
                for position in range(source.shape[0])
            ],
            name="source_row_id",
        )
    row_ids = [
        _required_text(
            value,
            field_name=f"phosphosite import source.{request.row_id_column}",
            row_position=position,
        )
        for position, value in enumerate(source.loc[:, request.row_id_column].tolist())
    ]
    index = pd.Index(row_ids, name=str(request.row_id_column))
    if index.has_duplicates:
        duplicates = list(dict.fromkeys(index[index.duplicated(keep=False)].tolist()))
        preview = ", ".join(repr(value) for value in duplicates[:5])
        suffix = "" if len(duplicates) <= 5 else " ..."
        raise PhosPyInputError(
            "phosphosite import request row_id_column must identify unique source "
            f"rows; duplicate_row_ids={preview}{suffix}"
        )
    return index


def _build_phospho_candidate(
    source: pd.DataFrame,
    *,
    row_index: pd.Index,
    sample_mapping: dict[str, str],
) -> pd.DataFrame:
    columns: dict[str, list[float]] = {}
    for source_column, sample_id in sample_mapping.items():
        columns[sample_id] = _parse_intensity_column(
            source.loc[:, source_column],
            source_column=source_column,
        )
    return pd.DataFrame(columns, index=row_index.copy(), dtype=float)


def _build_site_metadata_candidate(
    source: pd.DataFrame,
    *,
    request: PhosphositeImportRequest,
    row_index: pd.Index,
) -> tuple[pd.DataFrame, str | None, dict[str, object] | None, tuple[str, ...]]:
    metadata = pd.DataFrame(index=row_index.copy())
    metadata.loc[:, "gene_symbol"] = _required_text_column(
        source.loc[:, request.gene_symbol_column],
        field_name=f"phosphosite import source.{request.gene_symbol_column}",
    )
    metadata.loc[:, "site"] = _required_text_column(
        source.loc[:, request.site_column],
        field_name=f"phosphosite import source.{request.site_column}",
    )
    optional_column_map = {
        "protein_id": request.protein_id_column,
        "protein_accession": request.protein_accession_column,
        "protein_identifier": request.protein_identifier_column,
        "protein_namespace": request.protein_namespace_column,
        "organism": request.organism_column,
        "isoform_id": request.isoform_id_column,
        "site_sequence": request.site_sequence_column,
        "display_id": request.display_id_column,
        "site_key": request.site_key_column,
    }
    for output_column, source_column in optional_column_map.items():
        if source_column is None:
            continue
        metadata.loc[:, output_column] = _optional_text_column(
            source.loc[:, source_column],
            field_name=f"phosphosite import source.{source_column}",
        )

    if request.localisation_confidence_column is None:
        return metadata, None, None, ()

    normalised, report, warnings = normalise_localisation_confidence_series(
        source.loc[:, request.localisation_confidence_column],
        source_column=request.localisation_confidence_column,
        scale=request.localisation_confidence_scale,
        output_column=LOCALISATION_CONFIDENCE_OUTPUT_COLUMN,
    )
    metadata.loc[:, LOCALISATION_CONFIDENCE_OUTPUT_COLUMN] = normalised.tolist()
    return (
        metadata,
        LOCALISATION_CONFIDENCE_OUTPUT_COLUMN,
        report.to_payload(),
        warnings,
    )


def _build_peptide_evidence_candidate(
    source: pd.DataFrame,
    *,
    request: PhosphositeImportRequest,
    row_index: pd.Index,
    phospho: pd.DataFrame,
    localisation_column: str | None,
) -> pd.DataFrame | None:
    if not peptide_evidence_requested(request):
        return None
    peptide_row_ids = _resolve_peptide_row_ids(source, request, row_index=row_index)
    evidence = pd.DataFrame(
        {
            "peptide_row_id": peptide_row_ids,
            "site_id": _resolve_peptide_site_ids(source, request),
            "unique_feature_id": _required_text_column(
                source.loc[:, str(request.unique_feature_id_column)],
                field_name=(
                    f"phosphosite import source.{request.unique_feature_id_column}"
                ),
            ),
            "gene_symbol": _required_text_column(
                source.loc[:, request.gene_symbol_column],
                field_name=f"phosphosite import source.{request.gene_symbol_column}",
            ),
            "protein_accession": _resolve_evidence_protein_accession(source, request),
            "site_string": _required_text_column(
                source.loc[:, str(request.peptide_site_string_column)],
                field_name=(
                    f"phosphosite import source.{request.peptide_site_string_column}"
                ),
            ),
            "peptide_sequence": _required_text_column(
                source.loc[:, str(request.peptide_sequence_column)],
                field_name=(
                    f"phosphosite import source.{request.peptide_sequence_column}"
                ),
            ),
            "modified_peptide_sequence": _required_text_column(
                source.loc[:, str(request.modified_peptide_sequence_column)],
                field_name=(
                    "phosphosite import source."
                    f"{request.modified_peptide_sequence_column}"
                ),
            ),
            "multi_site": _resolve_multi_site_flags(source, request),
            "provenance_source": [request.source_name] * int(source.shape[0]),
        },
        index=row_index.copy(),
    )
    for column_name in phospho.columns.astype(str).tolist():
        evidence.loc[:, column_name] = phospho.loc[:, column_name].tolist()
    if request.site_sequence_column is not None:
        evidence.loc[:, "site_sequence"] = _optional_text_column(
            source.loc[:, request.site_sequence_column],
            field_name=f"phosphosite import source.{request.site_sequence_column}",
        )
    if localisation_column is not None:
        evidence.loc[:, localisation_column] = _optional_float_column(
            phospho.index,
            source.loc[:, str(request.localisation_confidence_column)],
            source_column=str(request.localisation_confidence_column),
            scale=request.localisation_confidence_scale,
        )
    return evidence


def _resolve_peptide_row_ids(
    source: pd.DataFrame,
    request: PhosphositeImportRequest,
    *,
    row_index: pd.Index,
) -> list[str]:
    if request.peptide_row_id_column is None:
        return row_index.astype(str).tolist()
    return _required_text_column(
        source.loc[:, request.peptide_row_id_column],
        field_name=f"phosphosite import source.{request.peptide_row_id_column}",
    )


def _resolve_peptide_site_ids(
    source: pd.DataFrame,
    request: PhosphositeImportRequest,
) -> list[str | None]:
    if request.peptide_site_id_column is not None:
        return _optional_site_id_column(
            source.loc[:, request.peptide_site_id_column],
            field_name=f"phosphosite import source.{request.peptide_site_id_column}",
        )
    gene_values = source.loc[:, request.gene_symbol_column].tolist()
    site_values = source.loc[:, str(request.peptide_site_string_column)].tolist()
    site_ids: list[str | None] = []
    for position, (gene_symbol, site_string) in enumerate(
        zip(gene_values, site_values, strict=True)
    ):
        site_ids.append(
            _site_id_from_gene_and_site_string(
                gene_symbol,
                site_string,
                row_position=position,
            )
        )
    return site_ids


def _resolve_evidence_protein_accession(
    source: pd.DataFrame,
    request: PhosphositeImportRequest,
) -> list[str]:
    source_column = (
        request.protein_accession_column
        or request.protein_identifier_column
        or request.protein_id_column
    )
    if source_column is None:  # pragma: no cover - validator owns this branch.
        raise PhosPyInputError(
            "phosphosite import peptide evidence requires protein context"
        )
    return _required_text_column(
        source.loc[:, source_column],
        field_name=f"phosphosite import source.{source_column}",
    )


def _resolve_multi_site_flags(
    source: pd.DataFrame,
    request: PhosphositeImportRequest,
) -> list[bool]:
    values = source.loc[:, str(request.peptide_site_string_column)].tolist()
    flags: list[bool] = []
    for position, value in enumerate(values):
        tokens = _parse_site_string_tokens(value, row_position=position)
        flags.append(len(tokens) > 1)
    return flags


def _build_import_diagnostics(
    *,
    source: pd.DataFrame,
    site_metadata: pd.DataFrame,
    sample_mapping: dict[str, str],
    peptide_evidence: pd.DataFrame | None,
    localisation_diagnostics: dict[str, object] | None,
    localisation_warnings: tuple[str, ...],
    source_name: str,
) -> tuple[dict[str, object], tuple[str, ...]]:
    site_labels = _diagnostic_site_labels(site_metadata)
    duplicate_site_rows = int(pd.Series(site_labels, dtype="object").duplicated().sum())
    multi_site_rows = _diagnostic_multi_site_rows(site_metadata)
    diagnostics: dict[str, object] = {
        "source_name": source_name,
        "input_row_count": int(source.shape[0]),
        "phospho_candidate_shape": (
            int(site_metadata.shape[0]),
            int(len(sample_mapping)),
        ),
        "site_metadata_candidate_shape": (
            int(site_metadata.shape[0]),
            int(site_metadata.shape[1]),
        ),
        "sample_column_mapping": dict(sample_mapping),
        "duplicate_site_candidate_rows": duplicate_site_rows,
        "multi_site_candidate_rows": multi_site_rows,
        "peptide_evidence_rows": (
            0 if peptide_evidence is None else int(peptide_evidence.shape[0])
        ),
    }
    if localisation_diagnostics is not None:
        diagnostics["localisation_confidence"] = localisation_diagnostics
    warnings = list(localisation_warnings)
    if duplicate_site_rows:
        warnings.append(
            "duplicate site candidates were retained for builder-owned duplicate "
            "site handling"
        )
    if multi_site_rows:
        warnings.append(
            "multi-site candidates were retained and reported; use peptide-evidence "
            "handoff with an explicit multi_site_policy when site-level resolution "
            "is required"
        )
    return diagnostics, tuple(dict.fromkeys(warnings))


def _build_import_quality_report(
    *,
    source: pd.DataFrame,
    phospho: pd.DataFrame,
    site_metadata: pd.DataFrame,
    sample_mapping: dict[str, str],
    localisation_diagnostics: dict[str, object] | None,
    warnings: tuple[str, ...],
    source_name: str,
) -> ImporterQualityReport:
    rows_read = int(source.shape[0])
    rows_retained = int(site_metadata.shape[0])
    return ImporterQualityReport(
        source_name=source_name,
        row_count_status=IMPORTER_QUALITY_STATUS_REPORTED,
        rows_read=rows_read,
        rows_retained=rows_retained,
        rows_dropped=max(rows_read - rows_retained, 0),
        intensity_column_status=IMPORTER_QUALITY_STATUS_REPORTED,
        detected_intensity_columns=tuple(
            ImporterDetectedIntensityColumn(
                source_column=source_column,
                sample_id=sample_id,
            )
            for source_column, sample_id in sample_mapping.items()
        ),
        missing_intensity=_build_missing_intensity_summary(
            phospho=phospho,
            sample_mapping=sample_mapping,
        ),
        localisation_confidence=_build_localisation_quality_summary(
            localisation_diagnostics,
        ),
        flagged_rows=ImporterFlaggedRowSummary(
            contaminant=ImporterQualityCount(
                status=IMPORTER_QUALITY_STATUS_NOT_APPLICABLE,
                reason="generic mapped importer does not parse contaminant flags",
            ),
            reverse=ImporterQualityCount(
                status=IMPORTER_QUALITY_STATUS_NOT_APPLICABLE,
                reason="generic mapped importer does not parse reverse flags",
            ),
            decoy=ImporterQualityCount(
                status=IMPORTER_QUALITY_STATUS_NOT_APPLICABLE,
                reason="generic mapped importer does not parse decoy flags",
            ),
        ),
        duplicate_keys=_build_duplicate_key_quality_summary(site_metadata),
        warnings=warnings,
    )


def _build_missing_intensity_summary(
    *,
    phospho: pd.DataFrame,
    sample_mapping: dict[str, str],
) -> ImporterMissingIntensitySummary:
    missing_by_sample_id: dict[str, int] = {}
    missing_by_source_column: dict[str, int] = {}
    for source_column, sample_id in sample_mapping.items():
        missing_count = int(phospho.loc[:, sample_id].isna().sum())
        missing_by_sample_id[sample_id] = missing_count
        missing_by_source_column[source_column] = missing_count
    return ImporterMissingIntensitySummary(
        status=IMPORTER_QUALITY_STATUS_REPORTED,
        total_missing_values=sum(missing_by_sample_id.values()),
        rows_with_any_missing_intensity=int(phospho.isna().any(axis=1).sum()),
        missing_values_by_sample_id=missing_by_sample_id,
        missing_values_by_source_column=missing_by_source_column,
    )


def _build_localisation_quality_summary(
    localisation_diagnostics: dict[str, object] | None,
) -> ImporterLocalisationConfidenceSummary:
    if localisation_diagnostics is None:
        return ImporterLocalisationConfidenceSummary(
            status=IMPORTER_QUALITY_STATUS_NOT_APPLICABLE,
            reason="localisation confidence column was not mapped",
        )
    invalid_examples = localisation_diagnostics.get("invalid_examples", ())
    if not isinstance(invalid_examples, (list, tuple)):
        invalid_examples = ()
    return ImporterLocalisationConfidenceSummary(
        status=IMPORTER_QUALITY_STATUS_REPORTED,
        source_column=str(localisation_diagnostics.get("source_column", "")),
        output_column=str(localisation_diagnostics.get("output_column", "")),
        scale=str(localisation_diagnostics.get("scale", "")),
        row_count=_quality_diagnostic_int(localisation_diagnostics, "row_count"),
        missing_count=_quality_diagnostic_int(
            localisation_diagnostics,
            "missing_count",
        ),
        invalid_count=_quality_diagnostic_int(
            localisation_diagnostics,
            "invalid_count",
        ),
        invalid_examples=tuple(str(value) for value in invalid_examples),
    )


def _quality_diagnostic_int(payload: dict[str, object], field_name: str) -> int:
    value = payload.get(field_name, 0)
    if isinstance(value, bool) or not isinstance(value, (int, float, str)):
        raise PhosPyInputError(
            f"importer localisation diagnostics {field_name} must be int-compatible"
        )
    return int(value)


def _build_duplicate_key_quality_summary(
    site_metadata: pd.DataFrame,
) -> ImporterDuplicateKeySummary:
    return ImporterDuplicateKeySummary(
        site_key=_quality_duplicate_count(
            site_metadata,
            column_name="site_key",
            missing_reason="site_key column was not mapped",
        ),
        display_key=_quality_duplicate_count(
            site_metadata,
            column_name="display_id",
            missing_reason="display_id column was not mapped",
        ),
        duplicate_site_candidate_rows=int(
            pd.Series(_diagnostic_site_labels(site_metadata), dtype="object")
            .duplicated()
            .sum()
        ),
    )


def _quality_duplicate_count(
    site_metadata: pd.DataFrame,
    *,
    column_name: str,
    missing_reason: str,
) -> ImporterQualityCount:
    if column_name not in site_metadata.columns:
        return ImporterQualityCount(
            status=IMPORTER_QUALITY_STATUS_NOT_APPLICABLE,
            reason=missing_reason,
        )
    return ImporterQualityCount(
        status=IMPORTER_QUALITY_STATUS_REPORTED,
        count=_diagnostic_duplicate_non_missing_rows(site_metadata.loc[:, column_name]),
        source_column=column_name,
    )


def _diagnostic_duplicate_non_missing_rows(values: pd.Series) -> int:
    tokens: list[str] = []
    for value in values.tolist():
        if _is_missing(value):
            continue
        token = str(value).strip()
        if token == "":
            continue
        tokens.append(token)
    return int(pd.Series(tokens, dtype="object").duplicated().sum())


def _parse_intensity_column(series: pd.Series, *, source_column: str) -> list[float]:
    values = pd.Series(series, dtype="object")
    missing_mask = values.isna() | values.map(
        lambda value: (
            isinstance(value, str) and value.strip().lower() in _MISSING_NUMERIC_TOKENS
        )
    )
    numeric = pd.to_numeric(values.mask(missing_mask), errors="coerce")
    parse_failures = (~missing_mask) & numeric.isna()
    if bool(parse_failures.any()):
        position = int(np.flatnonzero(parse_failures.to_numpy(dtype=bool))[0])
        value = values.iloc[position]
        raise PhosPyInputError(
            "failed to parse phosphosite import intensity value: "
            f"source_column='{source_column}', row_position={position}, "
            f"offending_value={value!r}"
        )
    numeric_values = numeric.to_numpy(dtype=float, copy=False, na_value=np.nan)
    finite_values = np.isfinite(numeric_values)
    not_finite = (~missing_mask.to_numpy(dtype=bool, copy=False)) & (~finite_values)
    if bool(not_finite.any()):
        position = int(np.flatnonzero(not_finite)[0])
        value = values.iloc[position]
        raise PhosPyInputError(
            "failed to parse phosphosite import intensity value: "
            f"source_column='{source_column}', row_position={position}, "
            f"offending_value={value!r}, reason='not_finite'"
        )
    return numeric_values.tolist()


def _required_text_column(series: pd.Series, *, field_name: str) -> list[str]:
    return [
        _required_text(value, field_name=field_name, row_position=position)
        for position, value in enumerate(series.tolist())
    ]


def _optional_text_column(series: pd.Series, *, field_name: str) -> list[str | None]:
    return [
        _optional_text(value, field_name=field_name, row_position=position)
        for position, value in enumerate(series.tolist())
    ]


def _optional_site_id_column(series: pd.Series, *, field_name: str) -> list[str | None]:
    values: list[str | None] = []
    for position, value in enumerate(series.tolist()):
        token = _optional_text(value, field_name=field_name, row_position=position)
        if token is None:
            values.append(None)
            continue
        try:
            gene_symbol, site = token.split(";", maxsplit=2)[:2]
            values.append(
                canonicalize_site_components(
                    gene_symbol,
                    site,
                    field_name=field_name,
                    error_type=PhosPyInputError,
                )
            )
        except (ValueError, PhosPyInputError) as exc:
            raise PhosPyInputError(
                f"{field_name} row_position={position} must contain site IDs in "
                "'GENE;SITE;' format when provided"
            ) from exc
    return values


def _optional_float_column(
    index: pd.Index,
    values: pd.Series,
    *,
    source_column: str,
    scale: str,
) -> list[float | None]:
    normalised, _, _ = normalise_localisation_confidence_series(
        pd.Series(values.tolist(), index=index.copy()),
        source_column=source_column,
        scale=scale,
    )
    return [
        None if _is_missing(value) else float(value)
        for value in normalised.astype("object").tolist()
    ]


def _site_id_from_gene_and_site_string(
    gene_symbol: object,
    site_string: object,
    *,
    row_position: int,
) -> str:
    gene = _required_text(
        gene_symbol,
        field_name="phosphosite import source.gene_symbol",
        row_position=row_position,
    )
    tokens = _parse_site_string_tokens(site_string, row_position=row_position)
    joint_site = ",".join(token.token for token in tokens)
    return canonicalize_site_components(
        gene,
        joint_site,
        field_name="phosphosite import peptide evidence site_id",
        error_type=PhosPyInputError,
    )


def _parse_site_string_tokens(
    value: object,
    *,
    row_position: int,
):
    try:
        return parse_phospho_site_tokens(
            value,
            field_name=(
                "phosphosite import source.peptide_site_string_column "
                f"row_position={row_position}"
            ),
        )
    except PhosPyInputError:
        raise


def _diagnostic_site_labels(site_metadata: pd.DataFrame) -> list[str]:
    labels: list[str] = []
    for gene_symbol, site in zip(
        site_metadata.loc[:, "gene_symbol"].tolist(),
        site_metadata.loc[:, "site"].tolist(),
        strict=True,
    ):
        labels.append(
            f"{str(gene_symbol).strip().upper()};{str(site).strip().upper()};"
        )
    return labels


def _diagnostic_multi_site_rows(site_metadata: pd.DataFrame) -> int:
    count = 0
    for site in site_metadata.loc[:, "site"].tolist():
        site_text = str(site)
        if "," in site_text or ";" in site_text:
            count += 1
    return count


def _required_text(
    value: object,
    *,
    field_name: str,
    row_position: int,
) -> str:
    if _is_missing(value):
        raise PhosPyInputError(
            f"{field_name} must not contain missing values; row_position={row_position}"
        )
    token = str(value).strip()
    if token == "":
        raise PhosPyInputError(
            f"{field_name} must contain non-empty values; row_position={row_position}"
        )
    return token


def _optional_text(
    value: object,
    *,
    field_name: str,
    row_position: int,
) -> str | None:
    if _is_missing(value):
        return None
    token = str(value).strip()
    if token == "":
        return None
    return token


def _is_missing(value: object) -> bool:
    if value is None or value is pd.NA or value is pd.NaT:
        return True
    if isinstance(value, float):
        return math.isnan(value)
    if isinstance(value, np.floating):
        scalar_value: object = value
        return str(scalar_value).lower() == "nan"
    if isinstance(value, (np.datetime64, np.timedelta64)):
        temporal_value: object = value
        return str(temporal_value) == "NaT"
    return False


__all__ = [
    "ColumnMappedPhosphositeImporter",
    "MappedPhosphositeTableImporter",
]
