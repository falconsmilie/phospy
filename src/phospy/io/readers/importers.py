"""Foundation classes for upstream phosphosite table importers."""

from __future__ import annotations

import math
from os import PathLike
from pathlib import Path

import pandas as pd

from phospy.contracts.requests import PhosphositeImportRequest
from phospy.contracts.results import PhosphositeImportResult
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
        return PhosphositeImportResult(
            phospho_matrix_candidate=phospho,
            site_metadata_candidate=site_metadata,
            peptide_evidence=peptide_evidence,
            sample_column_mapping=sample_mapping,
            localisation_confidence_column=localisation_column,
            warnings=warnings,
            diagnostics=diagnostics,
            source_name=validated.source_name,
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


def _parse_intensity_column(series: pd.Series, *, source_column: str) -> list[float]:
    parsed: list[float] = []
    for position, value in enumerate(series.tolist()):
        if _is_missing_numeric(value):
            parsed.append(float("nan"))
            continue
        try:
            numeric_value = float(value)
        except (TypeError, ValueError) as exc:
            raise PhosPyInputError(
                "failed to parse phosphosite import intensity value: "
                f"source_column='{source_column}', row_position={position}, "
                f"offending_value={value!r}"
            ) from exc
        if not math.isfinite(numeric_value):
            raise PhosPyInputError(
                "failed to parse phosphosite import intensity value: "
                f"source_column='{source_column}', row_position={position}, "
                f"offending_value={value!r}, reason='not_finite'"
            )
        parsed.append(numeric_value)
    return parsed


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
        None if bool(pd.isna(value)) else float(value)
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


def _is_missing_numeric(value: object) -> bool:
    if _is_missing(value):
        return True
    if isinstance(value, str) and value.strip().lower() in _MISSING_NUMERIC_TOKENS:
        return True
    return False


def _is_missing(value: object) -> bool:
    try:
        return bool(pd.Series((value,), dtype="object").isna().iat[0])
    except (TypeError, ValueError):
        return False


__all__ = [
    "ColumnMappedPhosphositeImporter",
    "MappedPhosphositeTableImporter",
]
