"""FragPipe importer result augmentation and quality reporting."""
# pyright: reportMissingTypeStubs=false

from __future__ import annotations

from dataclasses import replace

import pandas as pd

from phospy.contracts.results import (
    IMPORTER_QUALITY_STATUS_NOT_APPLICABLE,
    IMPORTER_QUALITY_STATUS_REPORTED,
    ImporterFlaggedRowSummary,
    ImporterQualityCount,
    ImporterQualityReport,
    PhosphositeImportResult,
)
from phospy.errors.input import PhosPyInputError
from phospy.io.readers.fragpipe.columns import _resolved_columns_payload
from phospy.io.readers.fragpipe.constants import (
    _ADAPTED_AMBIGUOUS_COLUMN,
    _ADAPTED_CANDIDATE_SITES_COLUMN,
    _ADAPTED_MODIFIED_PHOSPHO_COUNT_COLUMN,
    _ADAPTED_ROW_ID_COLUMN,
    _ADAPTED_SITE_PROBABILITIES_COLUMN,
    _FRAGPIPE_CONTAMINANT_OUTPUT_COLUMN,
    _FRAGPIPE_DECOY_OUTPUT_COLUMN,
)
from phospy.io.readers.fragpipe.models import _ResolvedFragPipeColumns
from phospy.validation.datasets.fragpipe import FRAGPIPE_FLAG_POLICY_FLAG


def _augment_mapped_result(
    mapped_result: PhosphositeImportResult,
    *,
    adapted: pd.DataFrame,
    flags: pd.DataFrame,
    contaminant_policy: str,
    decoy_policy: str,
    resolved: _ResolvedFragPipeColumns,
    filter_diagnostics: dict[str, object],
    adapter_diagnostics: dict[str, object],
    warnings: tuple[str, ...],
) -> PhosphositeImportResult:
    site_metadata = mapped_result.site_metadata_candidate
    peptide_evidence = mapped_result.peptide_evidence
    row_ids = adapted.loc[:, _ADAPTED_ROW_ID_COLUMN].astype(str).tolist()
    sidecars = adapted.loc[
        :,
        [
            _ADAPTED_CANDIDATE_SITES_COLUMN,
            _ADAPTED_SITE_PROBABILITIES_COLUMN,
            _ADAPTED_AMBIGUOUS_COLUMN,
            _ADAPTED_MODIFIED_PHOSPHO_COUNT_COLUMN,
        ],
    ].copy(deep=True)
    sidecars.index = pd.Index(row_ids, name=_ADAPTED_ROW_ID_COLUMN)
    for column_name in sidecars.columns.astype(str).tolist():
        site_metadata[column_name] = sidecars.loc[
            site_metadata.index,
            column_name,
        ].tolist()
        if peptide_evidence is not None:
            peptide_evidence[column_name] = sidecars.loc[
                peptide_evidence.index,
                column_name,
            ].tolist()

    flag_values = flags.copy(deep=True)
    flag_values.index = pd.Index(row_ids, name=_ADAPTED_ROW_ID_COLUMN)
    if contaminant_policy == FRAGPIPE_FLAG_POLICY_FLAG:
        site_metadata[_FRAGPIPE_CONTAMINANT_OUTPUT_COLUMN] = (
            flag_values.loc[site_metadata.index, _FRAGPIPE_CONTAMINANT_OUTPUT_COLUMN]
            .astype(bool)
            .tolist()
        )
        if peptide_evidence is not None:
            peptide_evidence[_FRAGPIPE_CONTAMINANT_OUTPUT_COLUMN] = (
                flag_values.loc[
                    peptide_evidence.index,
                    _FRAGPIPE_CONTAMINANT_OUTPUT_COLUMN,
                ]
                .astype(bool)
                .tolist()
            )
    if decoy_policy == FRAGPIPE_FLAG_POLICY_FLAG:
        site_metadata[_FRAGPIPE_DECOY_OUTPUT_COLUMN] = (
            flag_values.loc[site_metadata.index, _FRAGPIPE_DECOY_OUTPUT_COLUMN]
            .astype(bool)
            .tolist()
        )
        if peptide_evidence is not None:
            peptide_evidence[_FRAGPIPE_DECOY_OUTPUT_COLUMN] = (
                flag_values.loc[peptide_evidence.index, _FRAGPIPE_DECOY_OUTPUT_COLUMN]
                .astype(bool)
                .tolist()
            )

    diagnostics = dict(mapped_result.diagnostics)
    diagnostics["fragpipe"] = {
        "source_type": "fragpipe_ptmprophet_phosphosite",
        "resolved_columns": _resolved_columns_payload(resolved),
        "filtering": filter_diagnostics,
        "adaptation": adapter_diagnostics,
    }
    combined_warnings = tuple(dict.fromkeys((*mapped_result.warnings, *warnings)))
    quality_report = _augment_quality_report(
        mapped_result.quality_report,
        resolved=resolved,
        filter_diagnostics=filter_diagnostics,
        adapter_diagnostics=adapter_diagnostics,
        contaminant_policy=contaminant_policy,
        decoy_policy=decoy_policy,
        warnings=combined_warnings,
    )
    return PhosphositeImportResult(
        phospho_matrix_candidate=mapped_result.phospho_matrix_candidate,
        site_metadata_candidate=site_metadata,
        peptide_evidence=peptide_evidence,
        sample_column_mapping=mapped_result.sample_column_mapping,
        localisation_confidence_column=mapped_result.localisation_confidence_column,
        warnings=combined_warnings,
        diagnostics=diagnostics,
        source_name=mapped_result.source_name,
        quality_report=quality_report,
        _assume_owned=True,
    )


def _augment_quality_report(
    mapped_report: ImporterQualityReport,
    *,
    resolved: _ResolvedFragPipeColumns,
    filter_diagnostics: dict[str, object],
    adapter_diagnostics: dict[str, object],
    contaminant_policy: str,
    decoy_policy: str,
    warnings: tuple[str, ...],
) -> ImporterQualityReport:
    format_specific = dict(mapped_report.format_specific)
    format_specific["fragpipe_ptmprophet"] = {
        "resolved_columns": _resolved_columns_payload(resolved),
        "filtering": dict(filter_diagnostics),
        "adaptation": dict(adapter_diagnostics),
    }
    return replace(
        mapped_report,
        rows_read=_diagnostic_int(filter_diagnostics, "input_row_count"),
        rows_retained=_diagnostic_int(filter_diagnostics, "retained_row_count"),
        rows_dropped=_diagnostic_int(filter_diagnostics, "removed_rows"),
        localisation_confidence=replace(
            mapped_report.localisation_confidence,
            source_column=resolved.ptmprophet_probabilities,
        ),
        flagged_rows=ImporterFlaggedRowSummary(
            contaminant=_fragpipe_flag_quality_count(
                count=_diagnostic_int(filter_diagnostics, "contaminant_rows"),
                explicit_column=resolved.contaminant,
                fallback_column=resolved.protein_accession,
                policy=contaminant_policy,
                prefix_count=_diagnostic_int(
                    filter_diagnostics,
                    "contaminant_prefix_rows",
                ),
                label="contaminant",
            ),
            reverse=ImporterQualityCount(
                status=IMPORTER_QUALITY_STATUS_NOT_APPLICABLE,
                reason="FragPipe importer reports decoy flags instead of reverse flags",
            ),
            decoy=_fragpipe_flag_quality_count(
                count=_diagnostic_int(filter_diagnostics, "decoy_rows"),
                explicit_column=resolved.decoy,
                fallback_column=resolved.protein_accession,
                policy=decoy_policy,
                prefix_count=_diagnostic_int(filter_diagnostics, "decoy_prefix_rows"),
                label="decoy",
            ),
        ),
        format_specific=format_specific,
        warnings=warnings,
    )


def _fragpipe_flag_quality_count(
    *,
    count: int,
    explicit_column: str | None,
    fallback_column: str,
    policy: str,
    prefix_count: int,
    label: str,
) -> ImporterQualityCount:
    if explicit_column is None:
        return ImporterQualityCount(
            status=IMPORTER_QUALITY_STATUS_REPORTED,
            count=count,
            source_column=fallback_column,
            policy=policy,
            reason=f"{label} count derived from protein accession prefixes",
        )
    reason = None
    if prefix_count:
        reason = f"{label} count includes protein accession prefix matches"
    return ImporterQualityCount(
        status=IMPORTER_QUALITY_STATUS_REPORTED,
        count=count,
        source_column=explicit_column,
        policy=policy,
        reason=reason,
    )


def _diagnostic_int(payload: dict[str, object], field_name: str) -> int:
    value = payload[field_name]
    if isinstance(value, bool) or not isinstance(value, int):
        raise PhosPyInputError(f"FragPipe diagnostic {field_name} must be an int")
    return int(value)


__all__ = ["_augment_mapped_result"]
