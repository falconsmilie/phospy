"""Table-fingerprint invariants for analysis-ready dataset construction."""

from __future__ import annotations

import pandas as pd

from phospy.errors.validation import DatasetValidationError
from phospy.provenance.hashing import fingerprint_optional_table_strict
from phospy.provenance.models import RunProvenance, TableFingerprint
from phospy.science.datasets.direct_construction import (
    DIRECT_CONSTRUCTION_WORKFLOW_NAME,
)


def _fingerprints_for_analysis_ready_tables(
    *,
    phospho: pd.DataFrame,
    site_metadata: pd.DataFrame,
    sample_metadata: pd.DataFrame | None,
    total: pd.DataFrame | None,
    comparisons: pd.DataFrame | None,
    imputation_observation_mask: pd.DataFrame | None,
) -> tuple[TableFingerprint, ...]:
    entries: tuple[tuple[str, pd.DataFrame | None], ...] = (
        ("dataset.phospho", phospho),
        ("dataset.site_metadata", site_metadata),
        ("dataset.sample_metadata", sample_metadata),
        ("dataset.total", total),
        ("dataset.comparisons", comparisons),
        ("dataset.imputation_observation_mask", imputation_observation_mask),
    )
    fingerprints: list[TableFingerprint] = []
    for name, table in entries:
        fingerprint = fingerprint_optional_table_strict(table, name=name)
        if fingerprint is None:
            continue
        fingerprints.append(fingerprint)
    return tuple(fingerprints)


def _require_trusted_provenance_table_fingerprints(
    *,
    provenance: RunProvenance,
    actual_fingerprints: tuple[TableFingerprint, ...],
) -> None:
    _require_fingerprint_sets_match(
        expected=provenance.output_tables,
        actual=actual_fingerprints,
        field_name="run_provenance.output_tables",
        expected_source="actual analysis-ready dataset tables",
    )
    if provenance.workflow_name == DIRECT_CONSTRUCTION_WORKFLOW_NAME:
        _require_fingerprint_sets_match(
            expected=provenance.input_tables,
            actual=actual_fingerprints,
            field_name="run_provenance.input_tables",
            expected_source="actual analysis-ready dataset tables",
        )


def _require_fingerprint_sets_match(
    *,
    expected: tuple[TableFingerprint, ...],
    actual: tuple[TableFingerprint, ...],
    field_name: str,
    expected_source: str,
) -> None:
    expected_by_name = _fingerprint_map(expected, field_name=field_name)
    actual_by_name = _fingerprint_map(actual, field_name=field_name)
    missing = sorted(set(expected_by_name) - set(actual_by_name))
    unexpected = sorted(set(actual_by_name) - set(expected_by_name))
    if missing or unexpected:
        detail_parts: list[str] = []
        if missing:
            detail_parts.append("missing fingerprints: " + ", ".join(missing))
        if unexpected:
            detail_parts.append("unexpected fingerprints: " + ", ".join(unexpected))
        raise DatasetValidationError(
            f"{field_name} must match {expected_source}; " + "; ".join(detail_parts)
        )
    for name in expected_by_name:
        _require_fingerprint_matches(
            expected=expected_by_name[name],
            actual=actual_by_name[name],
            field_name=f"{field_name}.{name}",
            expected_source=expected_source,
        )


def _fingerprint_map(
    fingerprints: tuple[TableFingerprint, ...],
    *,
    field_name: str,
) -> dict[str, TableFingerprint]:
    result: dict[str, TableFingerprint] = {}
    for fingerprint in fingerprints:
        if fingerprint.name in result:
            raise DatasetValidationError(
                f"{field_name} contains duplicate table fingerprint "
                f"{fingerprint.name!r}"
            )
        result[fingerprint.name] = fingerprint
    return result


def _require_fingerprint_matches(
    *,
    expected: TableFingerprint,
    actual: TableFingerprint,
    field_name: str,
    expected_source: str,
) -> None:
    checks: tuple[tuple[str, object, object], ...] = (
        ("rows", expected.rows, actual.rows),
        ("columns", expected.columns, actual.columns),
        ("index_name", expected.index_name, actual.index_name),
        ("column_names", expected.column_names, actual.column_names),
        ("dtypes", expected.dtypes, actual.dtypes),
        ("index_structure", expected.index_structure, actual.index_structure),
        (
            "column_index_structure",
            expected.column_index_structure,
            actual.column_index_structure,
        ),
        (
            "exact_hash_algorithm",
            expected.exact_hash_algorithm,
            actual.exact_hash_algorithm,
        ),
        ("exact_hash_value", expected.exact_hash_value, actual.exact_hash_value),
        (
            "tolerance_hash_algorithm",
            expected.tolerance_hash_algorithm,
            actual.tolerance_hash_algorithm,
        ),
        (
            "tolerance_hash_value",
            expected.tolerance_hash_value,
            actual.tolerance_hash_value,
        ),
    )
    mismatched = [
        name
        for name, expected_value, actual_value in checks
        if expected_value != actual_value
    ]
    if mismatched:
        raise DatasetValidationError(
            f"{field_name} table fingerprint mismatch for {expected.name!r}; "
            f"does not match {expected_source}; mismatched fields: "
            + ", ".join(mismatched)
            + f"; expected exact digest {expected.exact_hash_value}; "
            f"actual exact digest {actual.exact_hash_value}; "
            f"expected tolerance digest {expected.tolerance_hash_value}; "
            f"actual tolerance digest {actual.tolerance_hash_value}"
        )
