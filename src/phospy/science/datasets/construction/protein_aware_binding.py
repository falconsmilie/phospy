"""Protein-aware preparation sidecar binding checks for dataset construction."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import TYPE_CHECKING

import pandas as pd

from phospy.errors.base import PhosPyError
from phospy.errors.validation import DatasetValidationError
from phospy.frames.comparison import dataframe_equals
from phospy.frames.validation import format_label_examples, require_exact_index_match
from phospy.provenance.hashing import fingerprint_optional_table_strict
from phospy.provenance.serialization.tables import table_fingerprints_from_payload
from phospy.science.configs.preprocessing.total_protein import (
    DATASET_PROTEIN_AWARE_PREPARATION_POLICY_PREPARE_MODEL_INPUTS,
)
from phospy.science.datasets.construction.fingerprints import (
    _require_fingerprint_sets_match,
)
from phospy.science.datasets.preprocessing.protein_aware_models import (
    PROTEIN_AWARE_PREPARATION_SCHEMA_VERSION,
    ProteinAwarePreparationReport,
    ProteinAwarePreparationResult,
)

if TYPE_CHECKING:
    from phospy.science.datasets.processing_state import DatasetPreprocessingReport

_BINDING_FINGERPRINTS_KEY = "dataset_binding_table_fingerprints"


def validate_protein_aware_preparation_binding(
    *,
    phospho: pd.DataFrame,
    site_metadata: pd.DataFrame,
    total: pd.DataFrame | None,
    preprocessing_report: DatasetPreprocessingReport | None,
    protein_aware_preparation: ProteinAwarePreparationResult | None,
) -> None:
    """Validate that an attached protein-aware sidecar matches dataset tables."""

    if protein_aware_preparation is None:
        return
    if total is None:
        raise DatasetValidationError(
            "dataset.protein_aware_preparation requires dataset.total because "
            "prepared protein covariates must be bound to owned total-protein rows"
        )

    report = protein_aware_preparation.report
    _validate_report_contract(report)
    _validate_report_site_order(phospho=phospho, report=report)
    _validate_preprocessing_report_binding(
        preprocessing_report=preprocessing_report,
        protein_aware_preparation=protein_aware_preparation,
    )

    matched_pairs = protein_aware_preparation.matched_pairs_dataframe()
    protein_covariates = protein_aware_preparation.protein_covariate_matrix_dataframe()
    _validate_matched_pair_order(matched_pairs=matched_pairs, report=report)
    _validate_matched_pair_report_binding(matched_pairs=matched_pairs, report=report)
    _validate_matched_pair_dataset_binding(
        matched_pairs=matched_pairs,
        phospho=phospho,
        site_metadata=site_metadata,
        total=total,
    )
    _validate_sample_binding(
        phospho=phospho,
        total=total,
        report=report,
        protein_covariates=protein_covariates,
    )
    _validate_covariate_binding(
        matched_pairs=matched_pairs,
        total=total,
        protein_covariates=protein_covariates,
    )
    _validate_table_fingerprint_binding(
        phospho=phospho,
        site_metadata=site_metadata,
        total=total,
        policy_parameters=report.policy_parameters,
    )


def _validate_report_contract(report: ProteinAwarePreparationReport) -> None:
    schema_version = report.schema_version
    if schema_version != PROTEIN_AWARE_PREPARATION_SCHEMA_VERSION:
        raise DatasetValidationError(
            "dataset.protein_aware_preparation.report.schema_version must be "
            f"{PROTEIN_AWARE_PREPARATION_SCHEMA_VERSION}; got {schema_version!r}"
        )
    preparation_policy = report.preparation_policy
    if (
        preparation_policy
        != DATASET_PROTEIN_AWARE_PREPARATION_POLICY_PREPARE_MODEL_INPUTS
    ):
        raise DatasetValidationError(
            "dataset.protein_aware_preparation.report.preparation_policy must be "
            "'prepare_model_inputs' when a sidecar is attached"
        )


def _validate_report_site_order(
    *,
    phospho: pd.DataFrame,
    report: ProteinAwarePreparationReport,
) -> None:
    site_keys = _string_index(phospho.index)
    report_site_keys = tuple(row.site_key for row in report.site_eligibility)
    if report_site_keys == site_keys:
        return
    raise DatasetValidationError(
        "dataset.protein_aware_preparation.report.site_eligibility must cover "
        "dataset.phospho.index exactly once and in dataset order; "
        f"dataset_site_count={len(site_keys)}, "
        f"report_site_count={len(report_site_keys)}; "
        + _site_key_mismatch_summary(
            left=report_site_keys,
            right=site_keys,
            left_name="report.site_eligibility.site_key",
            right_name="dataset.phospho.index",
        )
    )


def _validate_preprocessing_report_binding(
    *,
    preprocessing_report: DatasetPreprocessingReport | None,
    protein_aware_preparation: ProteinAwarePreparationResult,
) -> None:
    if preprocessing_report is None:
        return
    report = preprocessing_report.protein_aware_preparation
    if report is None:
        return
    if report.scientifically_equals(protein_aware_preparation.report):
        return
    raise DatasetValidationError(
        "dataset.preprocessing_report.protein_aware_preparation must be "
        "scientifically equal to dataset.protein_aware_preparation.report"
    )


def _validate_matched_pair_order(
    *,
    matched_pairs: pd.DataFrame,
    report: ProteinAwarePreparationReport,
) -> None:
    matched_site_keys = tuple(matched_pairs.loc[:, "site_key"].astype(str).tolist())
    if matched_site_keys == report.eligible_site_keys:
        return
    raise DatasetValidationError(
        "dataset.protein_aware_preparation.matched_pairs.site_key values must "
        "match report eligible_site_keys in order"
    )


def _validate_matched_pair_report_binding(
    *,
    matched_pairs: pd.DataFrame,
    report: ProteinAwarePreparationReport,
) -> None:
    report_rows_by_site = {row.site_key: row for row in report.site_eligibility}
    for row in matched_pairs.itertuples(index=False):
        site_key = str(row.site_key)
        report_row = report_rows_by_site[site_key]
        _validate_matched_pair_report_value(
            field_name="protein_identifier",
            site_key=site_key,
            matched_value=str(row.protein_identifier),
            report_value=report_row.protein_identifier,
        )
        _validate_matched_pair_report_value(
            field_name="total_protein_row_key",
            site_key=site_key,
            matched_value=str(row.total_protein_row_key),
            report_value=report_row.total_protein_row_key,
        )


def _validate_matched_pair_report_value(
    *,
    field_name: str,
    site_key: str,
    matched_value: str,
    report_value: str | None,
) -> None:
    if report_value is None:
        raise DatasetValidationError(
            "dataset.protein_aware_preparation.report.site_eligibility "
            f"{field_name} must be present for eligible site_key={site_key!r}"
        )
    if matched_value == report_value:
        return
    raise DatasetValidationError(
        "dataset.protein_aware_preparation.matched_pairs and "
        "report.site_eligibility must agree on "
        f"{field_name} for site_key={site_key!r}; "
        f"matched_pairs={matched_value!r}, report={report_value!r}"
    )


def _validate_matched_pair_dataset_binding(
    *,
    matched_pairs: pd.DataFrame,
    phospho: pd.DataFrame,
    site_metadata: pd.DataFrame,
    total: pd.DataFrame,
) -> None:
    site_keys = _string_index(phospho.index)
    site_key_set = set(site_keys)
    matched_site_keys = tuple(matched_pairs.loc[:, "site_key"].astype(str).tolist())
    missing_sites = [
        site_key for site_key in matched_site_keys if site_key not in site_key_set
    ]
    if missing_sites:
        raise DatasetValidationError(
            "dataset.protein_aware_preparation.matched_pairs.site_key values must "
            "exist in dataset.phospho.index; missing="
            + format_label_examples(missing_sites)
        )

    for row in matched_pairs.itertuples(index=False):
        site_key = str(row.site_key)
        expected_protein_identifier = str(
            site_metadata.loc[site_key, "protein_identifier"]
        )
        actual_protein_identifier = str(row.protein_identifier)
        if actual_protein_identifier != expected_protein_identifier:
            raise DatasetValidationError(
                "dataset.protein_aware_preparation.matched_pairs."
                "protein_identifier must match dataset.site_metadata."
                f"protein_identifier for site_key={site_key!r}; "
                f"expected={expected_protein_identifier!r}, "
                f"actual={actual_protein_identifier!r}"
            )

    total_row_key_set = set(_string_index(total.index))
    matched_total_row_keys = tuple(
        matched_pairs.loc[:, "total_protein_row_key"].astype(str).tolist()
    )
    missing_total_rows = [
        row_key
        for row_key in matched_total_row_keys
        if row_key not in total_row_key_set
    ]
    if missing_total_rows:
        raise DatasetValidationError(
            "dataset.protein_aware_preparation.matched_pairs."
            "total_protein_row_key values must exist in dataset.total.index; "
            "missing=" + format_label_examples(missing_total_rows)
        )


def _validate_sample_binding(
    *,
    phospho: pd.DataFrame,
    total: pd.DataFrame,
    report: ProteinAwarePreparationReport,
    protein_covariates: pd.DataFrame,
) -> None:
    phospho_columns = pd.Index(_string_index(phospho.columns))
    total_columns = pd.Index(_string_index(total.columns))
    report_phospho_columns = pd.Index(report.sample_alignment.phospho_sample_columns)
    report_total_columns = pd.Index(
        report.sample_alignment.total_protein_sample_columns
    )
    require_exact_index_match(
        left=report_phospho_columns,
        right=phospho_columns,
        left_name=(
            "dataset.protein_aware_preparation.report.sample_alignment."
            "phospho_sample_columns"
        ),
        right_name="dataset.phospho.columns",
        error_type=DatasetValidationError,
    )
    require_exact_index_match(
        left=report_total_columns,
        right=total_columns,
        left_name=(
            "dataset.protein_aware_preparation.report.sample_alignment."
            "total_protein_sample_columns"
        ),
        right_name="dataset.total.columns",
        error_type=DatasetValidationError,
    )
    require_exact_index_match(
        left=pd.Index(protein_covariates.columns.astype(str).tolist()),
        right=phospho_columns,
        left_name="dataset.protein_aware_preparation.protein_covariate_matrix.columns",
        right_name="dataset.phospho.columns",
        error_type=DatasetValidationError,
    )


def _validate_covariate_binding(
    *,
    matched_pairs: pd.DataFrame,
    total: pd.DataFrame,
    protein_covariates: pd.DataFrame,
) -> None:
    expected_row_keys = _dedupe_preserving_order(
        tuple(matched_pairs.loc[:, "total_protein_row_key"].astype(str).tolist())
    )
    actual_row_keys = tuple(protein_covariates.index.astype(str).tolist())
    if actual_row_keys != expected_row_keys:
        raise DatasetValidationError(
            "dataset.protein_aware_preparation.protein_covariate_matrix.index "
            "must match distinct matched total_protein_row_key values in order; "
            f"expected={list(expected_row_keys)!r}, actual={list(actual_row_keys)!r}"
        )
    if not expected_row_keys:
        return
    expected_covariates = total.loc[
        list(expected_row_keys),
        protein_covariates.columns.tolist(),
    ]
    if dataframe_equals(protein_covariates, expected_covariates):
        return
    raise DatasetValidationError(
        "dataset.protein_aware_preparation.protein_covariate_matrix must equal "
        "the corresponding owned dataset.total rows and sample columns exactly"
    )


def _validate_table_fingerprint_binding(
    *,
    phospho: pd.DataFrame,
    site_metadata: pd.DataFrame,
    total: pd.DataFrame,
    policy_parameters: Mapping[str, object],
) -> None:
    value = policy_parameters.get(_BINDING_FINGERPRINTS_KEY)
    if value is None:
        raise DatasetValidationError(
            "dataset.protein_aware_preparation.report.policy_parameters."
            f"{_BINDING_FINGERPRINTS_KEY} is required to prove the preparation "
            "sidecar is bound to the owning analysis-ready dataset tables"
        )
    try:
        expected_fingerprints = table_fingerprints_from_payload(
            value,
            field_name=(
                "dataset.protein_aware_preparation.report.policy_parameters."
                f"{_BINDING_FINGERPRINTS_KEY}"
            ),
        )
    except PhosPyError as exc:
        raise DatasetValidationError(
            "dataset.protein_aware_preparation.report.policy_parameters."
            f"{_BINDING_FINGERPRINTS_KEY} must contain table fingerprint payloads"
        ) from exc
    actual_fingerprints = tuple(
        fingerprint
        for fingerprint in (
            fingerprint_optional_table_strict(phospho, name="dataset.phospho"),
            fingerprint_optional_table_strict(
                site_metadata,
                name="dataset.site_metadata",
            ),
            fingerprint_optional_table_strict(total, name="dataset.total"),
        )
        if fingerprint is not None
    )
    _require_fingerprint_sets_match(
        expected=expected_fingerprints,
        actual=actual_fingerprints,
        field_name=(
            "dataset.protein_aware_preparation.report.policy_parameters."
            f"{_BINDING_FINGERPRINTS_KEY}"
        ),
        expected_source="actual analysis-ready dataset tables",
    )


def _string_index(index: pd.Index) -> tuple[str, ...]:
    return tuple(str(value) for value in index.tolist())


def _site_key_mismatch_summary(
    *,
    left: Sequence[str],
    right: Sequence[str],
    left_name: str,
    right_name: str,
) -> str:
    left_index = pd.Index(left)
    right_index = pd.Index(right)
    left_only = [value for value in left if value not in set(right)]
    right_only = [value for value in right if value not in set(left)]
    first_mismatch = _first_mismatch(left=left, right=right)
    labels_match_set_order_differs = (
        not left_only and not right_only and not left_index.equals(right_index)
    )
    if first_mismatch is None:
        mismatch_text = "none"
    else:
        position, left_value, right_value = first_mismatch
        mismatch_text = (
            f"position {position}, {left_name}={left_value!r}, "
            f"{right_name}={right_value!r}"
        )
    return "; ".join(
        (
            f"Only in {left_name}: {format_label_examples(left_only)}",
            f"Only in {right_name}: {format_label_examples(right_only)}",
            f"First positional mismatch: {mismatch_text}",
            (
                "Labels match as a set but order differs: "
                f"{str(labels_match_set_order_differs).lower()}"
            ),
        )
    )


def _first_mismatch(
    *,
    left: Sequence[str],
    right: Sequence[str],
) -> tuple[int, str, str] | None:
    shared_length = min(len(left), len(right))
    for position in range(shared_length):
        if left[position] != right[position]:
            return position, left[position], right[position]
    if len(left) == len(right):
        return None
    if len(left) > len(right):
        return shared_length, left[shared_length], "<missing>"
    return shared_length, "<missing>", right[shared_length]


def _dedupe_preserving_order(values: Sequence[str]) -> tuple[str, ...]:
    seen: set[str] = set()
    output: list[str] = []
    for value in values:
        if value in seen:
            continue
        seen.add(value)
        output.append(value)
    return tuple(output)


__all__: list[str] = []
