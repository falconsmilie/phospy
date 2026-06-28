"""Batch-correction metadata and design adequacy validation."""
# pyright: reportUnnecessaryIsInstance=false, reportUnknownMemberType=false
# Runtime boundary guards are intentionally retained for untyped external callers.

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import cast

import numpy as np
import pandas as pd

from phospy.contracts.configs.preprocessing.internal_batch_correction import (
    SUPPORTED_INTERNAL_BATCH_CORRECTION_EXECUTED_STAGE_ORDER,
)
from phospy.errors.input import PhosPyInputError
from phospy.provenance.environment import BATCH_CORRECTION_ENVIRONMENT_DEPENDENCIES
from phospy.provenance.models import (
    BATCH_CORRECTION_SELECTED_SITE_KEY_ROW_SENTINELS,
    BatchCorrectionProvenance,
)

_APPLIED_NATIVE_SPS_RUV_METHODS = frozenset({"sps_ruv_style"})
_UNSUPPORTED_SPS_RUV_METHODS = frozenset({"ruv_iii_style"})
_MISSING_PROVENANCE_MESSAGE = (
    "corrected_preprocessing_output with applied native SPS/RUV-style correction "
    "requires typed BatchCorrectionProvenance"
)
_NOT_PROVIDED_VALUES = frozenset({"not_provided", "not provided"})
_MISSING_ENVIRONMENT_VALUES = _NOT_PROVIDED_VALUES | frozenset({"unknown"})
_SELECTED_SITE_KEY_ROW_SENTINELS = (
    BATCH_CORRECTION_SELECTED_SITE_KEY_ROW_SENTINELS | _NOT_PROVIDED_VALUES
)
_STRICT_CONTROL_SOURCE_TYPE_MARKERS = frozenset({"packaged", "reference", "external"})
_CALLER_CONTROL_SOURCE_AUDIT_FIELDS = (
    "organism",
    "identifier_namespace",
    "source_version",
    "license",
    "redistribution",
)
_STRICT_CONTROL_SOURCE_REQUIRED_FIELDS = (
    "organism",
    "identifier_namespace",
    "source_name",
    "source_version",
    "license",
    "redistribution",
    "selection_method",
)
_MISSING = object()


@dataclass(frozen=True, slots=True)
class ResolvedBatchDesignMetadata:
    """Batch, condition, and optional replicate labels in matrix sample order."""

    batch_by_sample: Mapping[str, str]
    condition_by_sample: Mapping[str, str]
    sample_order: tuple[str, ...]
    replicate_by_sample: Mapping[str, str] | None = None
    replicate_structure_diagnostics: ReplicateStructureDiagnostics | None = None

    @property
    def batch_labels(self) -> tuple[str, ...]:
        return tuple(self.batch_by_sample[sample] for sample in self.sample_order)

    @property
    def condition_labels(self) -> tuple[str, ...]:
        return tuple(self.condition_by_sample[sample] for sample in self.sample_order)

    @property
    def replicate_labels(self) -> tuple[str, ...] | None:
        if self.replicate_by_sample is None:
            return None
        return tuple(self.replicate_by_sample[sample] for sample in self.sample_order)


@dataclass(frozen=True, slots=True)
class ReplicateStructureDiagnostics:
    """Structural diagnostics for provenance-only replicate metadata."""

    replicate_column: str
    sample_count: int
    replicate_count: int
    singleton_count: int
    singleton_replicates: tuple[str, ...]
    all_same: bool
    all_unique: bool
    perfectly_confounded_with_batch: bool
    perfectly_confounded_with_condition: bool
    diagnostic_flags: tuple[str, ...]

    def to_payload(self) -> dict[str, object]:
        return {
            "replicate_column": self.replicate_column,
            "sample_count": self.sample_count,
            "replicate_count": self.replicate_count,
            "singleton_count": self.singleton_count,
            "singleton_replicates": list(self.singleton_replicates),
            "all_same": self.all_same,
            "all_unique": self.all_unique,
            "perfectly_confounded_with_batch": (self.perfectly_confounded_with_batch),
            "perfectly_confounded_with_condition": (
                self.perfectly_confounded_with_condition
            ),
            "diagnostic_flags": list(self.diagnostic_flags),
            "policy": "provenance_only_structural_issues_are_rejected",
            "used_for_numerical_factor_estimation": False,
            "ruv_iii_semantics_enabled": False,
        }


class SampleMetadataAlignmentValidator:
    """Validate sample metadata without repairing or dropping rows."""

    def run(
        self,
        *,
        phospho: pd.DataFrame,
        sample_metadata: pd.DataFrame | None,
        required_columns: Sequence[str],
        context: str = "batch-correction",
    ) -> tuple[str, ...]:
        if sample_metadata is None:
            raise PhosPyInputError(
                f"{context} validation requires sample_metadata input data"
            )
        if not isinstance(phospho, pd.DataFrame):
            raise PhosPyInputError(f"{context} validation requires phospho DataFrame")
        if not isinstance(sample_metadata, pd.DataFrame):
            raise PhosPyInputError(
                f"{context} validation requires sample_metadata DataFrame"
            )

        sample_order = _normalize_sample_index(
            phospho.columns,
            field_name=f"{context} phospho.columns",
        )
        metadata_index = _normalize_sample_index(
            sample_metadata.index,
            field_name=f"{context} sample_metadata.index",
        )
        _require_unique_sample_index(
            sample_order,
            field_name=f"{context} phospho.columns",
        )
        _require_unique_sample_index(
            metadata_index,
            field_name=f"{context} sample_metadata.index",
        )
        _require_unique_metadata_columns(sample_metadata, context=context)

        metadata_samples = set(metadata_index.tolist())
        matrix_samples = set(sample_order.tolist())
        missing_samples = [
            sample for sample in sample_order.tolist() if sample not in metadata_samples
        ]
        extra_samples = [
            sample for sample in metadata_index.tolist() if sample not in matrix_samples
        ]
        if missing_samples or extra_samples:
            details: list[str] = []
            if missing_samples:
                details.append(
                    "missing sample_metadata rows for matrix columns: "
                    f"{_format_labels(missing_samples)}"
                )
            if extra_samples:
                details.append(
                    "sample_metadata rows not present in matrix columns: "
                    f"{_format_labels(extra_samples)}"
                )
            raise PhosPyInputError(
                f"{context} sample_metadata is misaligned with matrix columns; "
                + "; ".join(details)
            )

        for column in required_columns:
            _require_metadata_column(
                sample_metadata,
                column=str(column).strip(),
                context=context,
            )
        return tuple(str(sample) for sample in sample_order.tolist())


class BatchStructureValidator:
    """Validate batch labels in aligned sample metadata."""

    def run(
        self,
        *,
        sample_metadata: pd.DataFrame,
        sample_order: Sequence[str],
        batch_column: str,
        context: str = "batch-correction",
    ) -> Mapping[str, str]:
        column = str(batch_column).strip()
        _require_metadata_column(sample_metadata, column=column, context=context)
        return _resolve_column_labels(
            sample_metadata,
            sample_order=tuple(sample_order),
            column=column,
            label_kind="batch",
            context=context,
        )


class ConditionStructureValidator:
    """Validate condition columns and resolve protected joint condition strata.

    Multiple condition columns are intentionally combined as observed joint
    strata, not as additive protected-condition terms.
    """

    def run(
        self,
        *,
        sample_metadata: pd.DataFrame,
        sample_order: Sequence[str],
        condition_columns: Sequence[str],
        context: str = "batch-correction",
    ) -> Mapping[str, str]:
        if isinstance(condition_columns, str) or not isinstance(
            condition_columns, Sequence
        ):
            raise PhosPyInputError(
                f"{context} validation condition_columns must be a non-empty "
                "sequence of column names"
            )
        columns = tuple(str(column).strip() for column in condition_columns)
        if not columns or any(column == "" for column in columns):
            raise PhosPyInputError(
                f"{context} validation condition_columns must contain non-empty "
                "column names"
            )
        if len(set(columns)) != len(columns):
            raise PhosPyInputError(
                f"{context} validation condition_columns must not contain duplicates"
            )
        for column in columns:
            _require_metadata_column(sample_metadata, column=column, context=context)

        sample_tuple = tuple(sample_order)
        resolved_by_column = tuple(
            _resolve_column_labels(
                sample_metadata,
                sample_order=sample_tuple,
                column=column,
                label_kind="condition",
                context=context,
            )
            for column in columns
        )
        if len(columns) == 1:
            return resolved_by_column[0]
        return {
            sample: "|".join(
                f"{column}={labels[sample]}"
                for column, labels in zip(columns, resolved_by_column, strict=True)
            )
            for sample in sample_tuple
        }


class ReplicateStructureValidator:
    """Validate optional replicate labels when a request requires them."""

    def run(
        self,
        *,
        sample_metadata: pd.DataFrame,
        sample_order: Sequence[str],
        replicate_column: str | None,
        required: bool,
        context: str = "batch-correction",
    ) -> Mapping[str, str] | None:
        if replicate_column is None:
            if required:
                raise PhosPyInputError(
                    f"{context} validation requires replicate_column metadata"
                )
            return None
        column = str(replicate_column).strip()
        if column == "":
            raise PhosPyInputError(
                f"{context} validation replicate_column must be non-empty when provided"
            )
        _require_metadata_column(sample_metadata, column=column, context=context)
        return _resolve_column_labels(
            sample_metadata,
            sample_order=tuple(sample_order),
            column=column,
            label_kind="replicate",
            context=context,
        )


class ReplicateStructureDiagnosticHelper:
    """Compute diagnostics for supplied replicate labels."""

    def run(
        self,
        *,
        replicate_column: str | None,
        replicate_by_sample: Mapping[str, str] | None,
        batch_by_sample: Mapping[str, str],
        condition_by_sample: Mapping[str, str],
        sample_order: Sequence[str],
    ) -> ReplicateStructureDiagnostics | None:
        if replicate_column is None or replicate_by_sample is None:
            return None
        column = str(replicate_column).strip()
        samples = tuple(str(sample) for sample in sample_order)
        replicate_labels = tuple(replicate_by_sample[sample] for sample in samples)
        batch_labels = tuple(batch_by_sample[sample] for sample in samples)
        condition_labels = tuple(condition_by_sample[sample] for sample in samples)
        replicate_levels = _levels_in_order(replicate_labels)
        replicate_counts = {
            label: replicate_labels.count(label) for label in replicate_levels
        }
        singleton_replicates = tuple(
            label for label in replicate_levels if replicate_counts[label] == 1
        )
        replicate_count = len(replicate_levels)
        sample_count = len(samples)
        all_same = replicate_count == 1 and sample_count > 0
        all_unique = replicate_count == sample_count and sample_count > 0
        batch_confounded = _same_label_partition(replicate_labels, batch_labels)
        condition_confounded = _same_label_partition(
            replicate_labels,
            condition_labels,
        )
        diagnostic_flags: list[str] = []
        if all_same:
            diagnostic_flags.append("all_same_replicate_labels")
        if all_unique:
            diagnostic_flags.append("all_unique_replicate_labels")
        if batch_confounded:
            diagnostic_flags.append("batch_confounded_replicate_labels")
        if condition_confounded:
            diagnostic_flags.append("condition_confounded_replicate_labels")
        return ReplicateStructureDiagnostics(
            replicate_column=column,
            sample_count=sample_count,
            replicate_count=replicate_count,
            singleton_count=len(singleton_replicates),
            singleton_replicates=singleton_replicates,
            all_same=all_same,
            all_unique=all_unique,
            perfectly_confounded_with_batch=batch_confounded,
            perfectly_confounded_with_condition=condition_confounded,
            diagnostic_flags=tuple(diagnostic_flags),
        )


class BatchDesignMetadataValidator:
    """Coordinate sample, batch, condition, and replicate metadata validation."""

    def __init__(
        self,
        *,
        sample_alignment_validator: SampleMetadataAlignmentValidator | None = None,
        batch_structure_validator: BatchStructureValidator | None = None,
        condition_structure_validator: ConditionStructureValidator | None = None,
        replicate_structure_validator: ReplicateStructureValidator | None = None,
        replicate_structure_diagnostic_helper: (
            ReplicateStructureDiagnosticHelper | None
        ) = None,
    ) -> None:
        self._sample_alignment_validator = (
            sample_alignment_validator or SampleMetadataAlignmentValidator()
        )
        self._batch_structure_validator = (
            batch_structure_validator or BatchStructureValidator()
        )
        self._condition_structure_validator = (
            condition_structure_validator or ConditionStructureValidator()
        )
        self._replicate_structure_validator = (
            replicate_structure_validator or ReplicateStructureValidator()
        )
        self._replicate_structure_diagnostic_helper = (
            replicate_structure_diagnostic_helper
            or ReplicateStructureDiagnosticHelper()
        )

    def run(
        self,
        *,
        phospho: pd.DataFrame,
        sample_metadata: pd.DataFrame | None,
        batch_column: str,
        condition_columns: Sequence[str],
        replicate_column: str | None = None,
        require_replicate_column: bool = False,
        context: str = "batch-correction",
    ) -> ResolvedBatchDesignMetadata:
        required_columns: list[str] = [str(batch_column).strip()]
        required_columns.extend(str(column).strip() for column in condition_columns)
        if replicate_column is not None:
            required_columns.append(str(replicate_column).strip())
        sample_order = self._sample_alignment_validator.run(
            phospho=phospho,
            sample_metadata=sample_metadata,
            required_columns=tuple(required_columns),
            context=context,
        )
        if sample_metadata is None:
            raise PhosPyInputError(
                f"{context} validation requires sample_metadata input data"
            )
        metadata_frame = sample_metadata
        aligned = cast(pd.DataFrame, metadata_frame.copy(deep=False))
        aligned.index = _normalize_sample_index(
            metadata_frame.index,
            field_name=f"{context} sample_metadata.index",
        )
        aligned = cast(pd.DataFrame, aligned.reindex(sample_order))
        batch_by_sample = self._batch_structure_validator.run(
            sample_metadata=aligned,
            sample_order=sample_order,
            batch_column=batch_column,
            context=context,
        )
        condition_by_sample = self._condition_structure_validator.run(
            sample_metadata=aligned,
            sample_order=sample_order,
            condition_columns=condition_columns,
            context=context,
        )
        replicate_by_sample = self._replicate_structure_validator.run(
            sample_metadata=aligned,
            sample_order=sample_order,
            replicate_column=replicate_column,
            required=require_replicate_column,
            context=context,
        )
        replicate_structure_diagnostics = (
            self._replicate_structure_diagnostic_helper.run(
                replicate_column=replicate_column,
                replicate_by_sample=replicate_by_sample,
                batch_by_sample=batch_by_sample,
                condition_by_sample=condition_by_sample,
                sample_order=sample_order,
            )
        )
        return ResolvedBatchDesignMetadata(
            batch_by_sample=batch_by_sample,
            condition_by_sample=condition_by_sample,
            replicate_by_sample=replicate_by_sample,
            replicate_structure_diagnostics=replicate_structure_diagnostics,
            sample_order=sample_order,
        )


class DesignRankValidator:
    """Validate full column rank for numeric design matrices."""

    def run(self, design_matrix: pd.DataFrame, *, context: str) -> int:
        if not isinstance(design_matrix, pd.DataFrame):
            raise PhosPyInputError(f"{context} design matrix must be a DataFrame")
        if not design_matrix.index.is_unique:
            raise PhosPyInputError(
                f"{context} design matrix sample labels must be unique"
            )
        if not design_matrix.columns.is_unique:
            raise PhosPyInputError(
                f"{context} design matrix coefficient labels must be unique"
            )
        try:
            values = np.asarray(design_matrix.to_numpy(dtype=float), dtype=float)
        except (TypeError, ValueError) as exc:
            raise PhosPyInputError(
                f"{context} design matrix must contain numeric values"
            ) from exc
        if not np.isfinite(values).all():
            raise PhosPyInputError(
                f"{context} design matrix must contain finite numeric values"
            )
        rank = _matrix_rank(values)
        column_count = int(values.shape[1])
        if rank < column_count:
            coefficients = ", ".join(str(label) for label in design_matrix.columns)
            raise PhosPyInputError(
                f"{context} design matrix is rank-deficient; condition, batch, "
                "or replicate terms are collinear or confounded "
                f"(rank={rank}, columns={column_count}, coefficients={coefficients})"
            )
        return rank


class BatchCorrectionAdequacyValidator:
    """Validate fixed-effect batch-correction design adequacy.

    The validator only inspects sample labels and categorical design rank. It does
    not correct matrices, estimate coefficients, or mutate metadata.
    """

    def __init__(
        self,
        *,
        design_rank_validator: DesignRankValidator | None = None,
    ) -> None:
        self._design_rank_validator = design_rank_validator or DesignRankValidator()

    def run(
        self,
        *,
        batch_by_sample: Mapping[str, object],
        condition_by_sample: Mapping[str, object],
        sample_order: Sequence[str],
        preserve_condition_effects: bool,
        context: str = "linear_residualize_batch",
    ) -> None:
        samples = _normalize_sample_order(sample_order, context=context)
        if preserve_condition_effects is not True:
            raise PhosPyInputError(
                f"{context} requires "
                "preprocessing_config.batch_correction.preserve_condition_effects=True; "
                "refusing batch correction because condition effects would not be "
                "explicitly preserved"
            )
        if len(samples) < 2:
            raise PhosPyInputError(
                f"{context} requires at least two samples to "
                "estimate batch effects while preserving condition effects"
            )

        batch_labels = _resolve_labels(
            batch_by_sample,
            sample_order=samples,
            label_kind="batch",
            context=context,
        )
        condition_labels = _resolve_labels(
            condition_by_sample,
            sample_order=samples,
            label_kind="condition",
            context=context,
        )
        batch_levels = _levels_in_order(batch_labels)
        condition_levels = _levels_in_order(condition_labels)
        if len(batch_levels) < 2:
            raise PhosPyInputError(
                f"{context} requires at least two batch levels; "
                f"observed {len(batch_levels)}"
            )

        singleton_batches = _singleton_levels(batch_labels)
        if singleton_batches:
            raise PhosPyInputError(
                f"{context} requires at least two samples in each "
                "batch level to estimate batch effects; singleton batch levels: "
                f"{_format_labels(singleton_batches)}"
            )

        preservation_design = _treatment_coded_design(
            condition_labels,
            include_intercept=True,
        )
        preservation_columns = int(preservation_design.shape[1])
        preservation_rank = _matrix_rank(preservation_design)
        if preservation_rank < preservation_columns:
            raise PhosPyInputError(
                f"{context} condition preservation design is "
                "rank-deficient; condition effects cannot be explicitly preserved "
                f"(rank={preservation_rank}, columns={preservation_columns})"
            )
        if len(samples) <= preservation_rank:
            raise PhosPyInputError(
                f"{context} condition preservation design is "
                "saturated; batch effects cannot be estimated while preserving "
                "condition effects "
                f"(samples={len(samples)}, condition_design_rank={preservation_rank})"
            )

        batch_terms = _treatment_coded_design(batch_labels, include_intercept=False)
        full_design = np.concatenate((preservation_design, batch_terms), axis=1)
        full_columns = int(full_design.shape[1])
        if len(samples) <= full_columns:
            raise PhosPyInputError(
                f"{context} requires more samples than estimable "
                "condition-plus-batch design parameters; "
                f"samples={len(samples)}, design_columns={full_columns}. Add "
                "replicate samples or reduce batch/condition levels."
            )

        if len(condition_levels) > 1 and _condition_is_determined_by_batch(
            batch_labels=batch_labels,
            condition_labels=condition_labels,
        ):
            raise PhosPyInputError(
                f"{context} cannot run because batch and condition "
                "are perfectly confounded: each batch level contains only one "
                "condition level, so removing batch would remove biological "
                "condition signal"
            )
        if len(condition_levels) > 1 and _batch_is_determined_by_condition(
            batch_labels=batch_labels,
            condition_labels=condition_labels,
        ):
            raise PhosPyInputError(
                f"{context} cannot run because batch and condition "
                "are perfectly confounded: each condition level contains only one "
                "batch level, so batch cannot be estimated while preserving "
                "condition effects"
            )

        full_frame = pd.DataFrame(
            full_design,
            index=pd.Index(samples, name="sample"),
            columns=pd.Index(
                (
                    *(f"condition[{index}]" for index in range(preservation_columns)),
                    *(f"batch[{index}]" for index in range(batch_terms.shape[1])),
                ),
                name="coefficient",
            ),
        )
        try:
            full_rank = self._design_rank_validator.run(
                full_frame,
                context=f"{context} batch/condition",
            )
        except PhosPyInputError as exc:
            raise PhosPyInputError(
                f"{context} batch/condition design is "
                "rank-deficient; batch effects are not estimable while preserving "
                "condition effects "
                f"(rank={_matrix_rank(full_design)}, columns={full_columns})"
            ) from exc
        batch_degrees = int(batch_terms.shape[1])
        batch_rank_after_condition = full_rank - preservation_rank
        if full_rank < full_columns or batch_rank_after_condition < batch_degrees:
            raise PhosPyInputError(
                f"{context} batch/condition design is "
                "rank-deficient; batch effects are not estimable while preserving "
                "condition effects "
                f"(rank={full_rank}, columns={full_columns}, "
                f"estimable_batch_degrees={batch_rank_after_condition}, "
                f"batch_degrees={batch_degrees})"
            )


def validate_applied_native_sps_ruv_correction_provenance(
    *,
    method: object,
    status: object,
    provenance: object,
) -> None:
    """Require complete typed provenance for applied SPS/RUV-style outputs."""

    if str(status).strip() != "applied":
        return
    normalized_method = _normalise_method(method)
    if normalized_method in _UNSUPPORTED_SPS_RUV_METHODS:
        raise PhosPyInputError(
            "corrected_preprocessing_output declares unsupported SPS/RUV-style "
            f"batch correction method {normalized_method!r}; regenerate with a "
            "supported native method and complete BatchCorrectionProvenance"
        )
    if not _is_sps_ruv_style_label(normalized_method):
        return
    if normalized_method not in _APPLIED_NATIVE_SPS_RUV_METHODS:
        raise PhosPyInputError(
            "corrected_preprocessing_output declares ambiguous or unsupported "
            f"SPS/RUV-style batch correction method {normalized_method!r}; "
            "applied corrected outputs require a supported native method and "
            "typed BatchCorrectionProvenance"
        )
    if provenance is None:
        raise PhosPyInputError(_MISSING_PROVENANCE_MESSAGE)
    if not isinstance(provenance, BatchCorrectionProvenance):
        raise PhosPyInputError(
            _MISSING_PROVENANCE_MESSAGE
            + "; untyped provenance payloads are not accepted for applied "
            "SPS/RUV-style corrected outputs"
        )
    _validate_complete_sps_ruv_provenance(
        provenance,
        expected_method=normalized_method,
    )


def _validate_complete_sps_ruv_provenance(
    provenance: BatchCorrectionProvenance,
    *,
    expected_method: str,
) -> None:
    requested_method = _normalise_method(provenance.requested_method)
    if requested_method != expected_method:
        raise PhosPyInputError(
            "corrected_preprocessing_output BatchCorrectionProvenance method must "
            f"match the applied correction method; expected {expected_method!r}, "
            f"observed {requested_method!r}"
        )
    if requested_method in _UNSUPPORTED_SPS_RUV_METHODS:
        raise PhosPyInputError(
            "corrected_preprocessing_output BatchCorrectionProvenance declares "
            f"unsupported method {requested_method!r}"
        )
    _require_environment_provenance(provenance)
    selected_site_key_rows = normalize_applied_selected_site_key_rows(
        provenance.selected_site_key_rows
    )
    n_unwanted_factors = _extract_sps_ruv_n_unwanted_factors(provenance)
    _require_selected_control_count_for_unwanted_factors(
        selected_site_key_rows=selected_site_key_rows,
        n_unwanted_factors=n_unwanted_factors,
    )
    _require_unique_selected_control_site_rows(selected_site_key_rows)
    _require_control_site_source_metadata(
        provenance.control_site_source,
        selected_site_key_rows=selected_site_key_rows,
    )
    _require_non_empty_mapping(
        provenance.batch_metadata,
        field_name="BatchCorrectionProvenance.batch_metadata",
    )
    _require_non_empty_mapping(
        provenance.design_metadata,
        field_name="BatchCorrectionProvenance.design_metadata",
    )
    _require_non_empty_mapping(
        provenance.missing_value_policy,
        field_name="BatchCorrectionProvenance.missing_value_policy",
    )
    if not provenance.observation_masks:
        raise PhosPyInputError(
            "corrected_preprocessing_output BatchCorrectionProvenance must include "
            "observation mask fingerprints for SPS/RUV-style missingness provenance"
        )
    input_matrix_fingerprint = getattr(provenance, "input_matrix_fingerprint", None)
    if input_matrix_fingerprint is None:
        raise PhosPyInputError(
            "corrected_preprocessing_output BatchCorrectionProvenance is missing "
            "input/output matrix fingerprints: input_matrix_fingerprint is required"
        )
    if provenance.output_matrix_fingerprint is None:
        raise PhosPyInputError(
            "corrected_preprocessing_output BatchCorrectionProvenance is missing "
            "input/output matrix fingerprints: output_matrix_fingerprint is required"
        )
    _require_supported_stage_order(provenance.preprocessing_stage_order)
    _require_non_empty_mapping(
        provenance.diagnostics,
        field_name="BatchCorrectionProvenance.diagnostics",
    )
    _reject_not_provided_required_mapping(
        provenance.resolved_parameters,
        field_name="BatchCorrectionProvenance.resolved_parameters",
    )
    _reject_not_provided_required_mapping(
        provenance.batch_metadata,
        field_name="BatchCorrectionProvenance.batch_metadata",
    )
    _reject_not_provided_required_mapping(
        provenance.design_metadata,
        field_name="BatchCorrectionProvenance.design_metadata",
    )
    _reject_not_provided_required_mapping(
        provenance.missing_value_policy,
        field_name="BatchCorrectionProvenance.missing_value_policy",
    )


def _require_environment_provenance(provenance: BatchCorrectionProvenance) -> None:
    if _is_missing_environment_text(provenance.phospy_version):
        raise PhosPyInputError(
            "corrected_preprocessing_output BatchCorrectionProvenance must include "
            "non-empty phospy_version for applied SPS/RUV-style corrected output"
        )
    if _is_missing_environment_text(provenance.python_version):
        raise PhosPyInputError(
            "corrected_preprocessing_output BatchCorrectionProvenance must include "
            "non-empty python_version for applied SPS/RUV-style corrected output"
        )
    dependency_versions = provenance.dependency_versions
    if not isinstance(dependency_versions, Mapping) or not dependency_versions:
        raise PhosPyInputError(
            "corrected_preprocessing_output BatchCorrectionProvenance must include "
            "non-empty dependency_versions for applied SPS/RUV-style corrected output"
        )

    missing_dependencies = tuple(
        dependency
        for dependency in BATCH_CORRECTION_ENVIRONMENT_DEPENDENCIES
        if _is_missing_environment_text(dependency_versions.get(dependency))
    )
    if missing_dependencies:
        raise PhosPyInputError(
            "corrected_preprocessing_output BatchCorrectionProvenance "
            "dependency_versions must include versions for native "
            "batch-correction dependencies: "
            f"{_format_labels(missing_dependencies)}"
        )


def _extract_sps_ruv_n_unwanted_factors(
    provenance: BatchCorrectionProvenance,
) -> int:
    resolved_parameters = provenance.resolved_parameters
    value = resolved_parameters.get("n_unwanted_factors", _MISSING)
    if value is _MISSING:
        config = resolved_parameters.get("config")
        if isinstance(config, Mapping):
            value = cast(Mapping[str, object], config).get(
                "n_unwanted_factors",
                _MISSING,
            )
    if value is _MISSING:
        raise PhosPyInputError(
            "corrected_preprocessing_output BatchCorrectionProvenance selected "
            "controls cannot be validated against unwanted-factor count because "
            "resolved_parameters.n_unwanted_factors is missing"
        )
    if isinstance(value, bool) or not isinstance(value, int):
        raise PhosPyInputError(
            "corrected_preprocessing_output BatchCorrectionProvenance selected "
            "controls require a positive integer unwanted-factor count; "
            f"observed n_unwanted_factors={value!r}"
        )
    if value < 1:
        raise PhosPyInputError(
            "corrected_preprocessing_output BatchCorrectionProvenance selected "
            "controls require unwanted-factor count n_unwanted_factors >= 1; "
            f"observed n_unwanted_factors={value}"
        )
    return value


def _require_selected_control_count_for_unwanted_factors(
    *,
    selected_site_key_rows: Sequence[str],
    n_unwanted_factors: int,
) -> None:
    duplicates = _duplicates_in_order(selected_site_key_rows)
    selected_count = len(_unique_in_order(selected_site_key_rows))
    required_count = n_unwanted_factors + 1
    if selected_count < required_count:
        duplicate_detail = (
            ""
            if not duplicates
            else "; selected_site_key_rows duplicate identifiers: "
            f"{_format_labels(duplicates)}"
        )
        raise PhosPyInputError(
            "corrected_preprocessing_output BatchCorrectionProvenance selected "
            "controls are too few for unwanted-factor count; "
            f"unique_selected_controls={selected_count}, "
            f"n_unwanted_factors={n_unwanted_factors}, "
            f"required_selected_controls={required_count}"
            f"{duplicate_detail}"
        )


def _normalise_method(method: object) -> str:
    normalized = str(method).strip().lower()
    if not normalized:
        raise PhosPyInputError(
            "corrected_preprocessing_output applied batch correction method is "
            "missing or empty"
        )
    return normalized


def _is_sps_ruv_style_label(method: str) -> bool:
    return "ruv" in method or method.startswith("sps_") or "_sps_" in method


def normalize_applied_selected_site_key_rows(rows: Sequence[object]) -> tuple[str, ...]:
    """Normalize and validate applied selected control row identifiers."""

    if not rows:
        raise PhosPyInputError(
            "corrected_preprocessing_output BatchCorrectionProvenance selected "
            "controls/control provenance must include selected_site_key_rows"
        )

    normalized_rows: list[str] = []
    missing_rows: list[int] = []
    blank_rows: list[int] = []
    sentinel_rows: list[int] = []
    for position, row in enumerate(tuple(rows)):
        if _is_missing_value(row):
            missing_rows.append(position)
            continue
        normalized = str(row).strip()
        if normalized == "":
            blank_rows.append(position)
            continue
        if _is_selected_site_key_row_sentinel(normalized):
            sentinel_rows.append(position)
            continue
        normalized_rows.append(normalized)

    if missing_rows:
        raise PhosPyInputError(
            "corrected_preprocessing_output BatchCorrectionProvenance selected "
            "controls/control provenance selected_site_key_rows contains missing "
            f"site_key rows at positions {_format_positions(missing_rows)}"
        )
    if blank_rows:
        raise PhosPyInputError(
            "corrected_preprocessing_output BatchCorrectionProvenance selected "
            "controls/control provenance selected_site_key_rows contains blank "
            f"site_key rows at positions {_format_positions(blank_rows)}"
        )
    if sentinel_rows:
        raise PhosPyInputError(
            "corrected_preprocessing_output BatchCorrectionProvenance selected "
            "controls/control provenance selected_site_key_rows contains sentinel "
            f"site_key rows at positions {_format_positions(sentinel_rows)}"
        )

    return tuple(normalized_rows)


def _is_selected_site_key_row_sentinel(value: object) -> bool:
    return str(value).strip().lower() in _SELECTED_SITE_KEY_ROW_SENTINELS


def _require_unique_selected_control_site_rows(rows: Sequence[str]) -> None:
    duplicates = _duplicates_in_order(rows)
    if duplicates:
        raise PhosPyInputError(
            "corrected_preprocessing_output BatchCorrectionProvenance "
            "selected_site_key_rows contains duplicate selected control row "
            f"identifiers: {_format_labels(duplicates)}"
        )


def _require_control_site_source_metadata(
    source: Mapping[str, object],
    *,
    selected_site_key_rows: Sequence[str],
) -> None:
    _require_non_empty_mapping(
        source,
        field_name="BatchCorrectionProvenance.control_site_source",
    )
    _reject_not_provided_required_mapping(
        source,
        field_name="BatchCorrectionProvenance.control_site_source",
    )
    source_type = _source_type(source)
    if source_type is None or _is_not_provided(source_type):
        raise PhosPyInputError(
            "corrected_preprocessing_output BatchCorrectionProvenance selected "
            "controls/control provenance must include control source metadata"
        )

    if _has_strict_control_source_type(source):
        missing = tuple(
            field_name
            for field_name in _STRICT_CONTROL_SOURCE_REQUIRED_FIELDS
            if not _has_non_missing_text(source.get(field_name))
        )
        if missing:
            raise PhosPyInputError(
                "corrected_preprocessing_output BatchCorrectionProvenance "
                "packaged/reference/external control-source metadata is "
                f"incomplete; missing {_format_labels(missing)}"
            )
        return

    missing_without_reason = tuple(
        field_name
        for field_name in _CALLER_CONTROL_SOURCE_AUDIT_FIELDS
        if not _has_non_missing_text(source.get(field_name))
        and not _has_metadata_missing_reason(
            source,
            field_name,
            selected_site_key_rows=selected_site_key_rows,
        )
    )
    if missing_without_reason:
        raise PhosPyInputError(
            "corrected_preprocessing_output BatchCorrectionProvenance control "
            "source is missing "
            f"{_format_labels(missing_without_reason)} without explicit rationale"
        )

    has_source_name = _has_non_missing_text(source.get("source_name"))
    has_source_version = _has_non_missing_text(source.get("source_version"))
    has_unavailable_reason = _has_non_missing_text(
        source.get("source_version_unavailable_reason")
    )
    has_missing_reason = _has_metadata_missing_reason(
        source,
        "source_version",
        selected_site_key_rows=selected_site_key_rows,
    )
    if has_source_name and not (
        has_source_version or has_unavailable_reason or has_missing_reason
    ):
        raise PhosPyInputError(
            "corrected_preprocessing_output BatchCorrectionProvenance control "
            "source declares source_name without source_version or explicit "
            "source_version_unavailable_reason"
        )
    if source_type == "caller_supplied" and not (
        has_source_version or has_unavailable_reason or has_missing_reason
    ):
        raise PhosPyInputError(
            "corrected_preprocessing_output BatchCorrectionProvenance "
            "caller_supplied control source must record source_version or "
            "source_version_unavailable_reason"
        )


def _source_type(source: Mapping[str, object]) -> str | None:
    for key in ("source_type", "source"):
        value = source.get(key)
        if _has_non_missing_text(value):
            return str(value).strip().lower()
    return None


def _has_strict_control_source_type(source: Mapping[str, object]) -> bool:
    for key in (
        "source_type",
        "source",
        "control_site_set_source_type",
        "source_name",
    ):
        value = source.get(key)
        if _has_non_missing_text(value) and _is_strict_control_source_type(
            str(value).strip().lower()
        ):
            return True
    return False


def _is_strict_control_source_type(source_type: str | None) -> bool:
    if source_type is None:
        return False
    tokens = frozenset(source_type.replace("-", "_").split("_"))
    return bool(tokens & _STRICT_CONTROL_SOURCE_TYPE_MARKERS)


def _has_metadata_missing_reason(
    source: Mapping[str, object],
    field_name: str,
    *,
    selected_site_key_rows: Sequence[str],
) -> bool:
    reasons = source.get("metadata_missing_reason")
    if isinstance(reasons, Mapping) and _has_non_missing_text(
        cast(Mapping[str, object], reasons).get(field_name)
    ):
        return True
    if _has_non_missing_text(source.get(f"{field_name}_missing_reason")):
        return True
    if field_name == "source_version" and _has_non_missing_text(
        source.get("source_version_unavailable_reason")
    ):
        return True
    return _has_metadata_missing_reason_by_site_key(
        source,
        field_name,
        selected_site_key_rows=selected_site_key_rows,
    )


def _has_metadata_missing_reason_by_site_key(
    source: Mapping[str, object],
    field_name: str,
    *,
    selected_site_key_rows: Sequence[str],
) -> bool:
    reasons_by_site_key = source.get("metadata_missing_reason_by_site_key")
    if not isinstance(reasons_by_site_key, Mapping):
        return False
    selected = tuple(str(site_key) for site_key in selected_site_key_rows)
    if not selected:
        return False
    by_site_key = cast(Mapping[str, object], reasons_by_site_key)
    for site_key in selected:
        site_reasons = by_site_key.get(site_key)
        if not isinstance(site_reasons, Mapping):
            return False
        if not _has_non_missing_text(
            cast(Mapping[str, object], site_reasons).get(field_name)
        ):
            return False
    return True


def _require_supported_stage_order(stage_order: Sequence[str]) -> None:
    normalized = tuple(str(stage).strip() for stage in tuple(stage_order))
    if not normalized:
        raise PhosPyInputError(
            "corrected_preprocessing_output BatchCorrectionProvenance stage order "
            "is missing"
        )
    if normalized != SUPPORTED_INTERNAL_BATCH_CORRECTION_EXECUTED_STAGE_ORDER:
        supported = " -> ".join(
            SUPPORTED_INTERNAL_BATCH_CORRECTION_EXECUTED_STAGE_ORDER
        )
        observed = " -> ".join(normalized)
        raise PhosPyInputError(
            "corrected_preprocessing_output BatchCorrectionProvenance stage order "
            f"is unsupported; observed {observed!r}; supported stage order is "
            f"{supported}"
        )


def _require_non_empty_mapping(
    value: Mapping[str, object],
    *,
    field_name: str,
) -> None:
    if not isinstance(value, Mapping) or not value:
        raise PhosPyInputError(f"{field_name} must be a non-empty object")


def _reject_not_provided_required_mapping(
    value: Mapping[str, object],
    *,
    field_name: str,
) -> None:
    for key, item in value.items():
        if _is_not_provided(item):
            raise PhosPyInputError(
                f"{field_name}.{key} must not be recorded as not_provided for "
                "applied SPS/RUV-style corrected output"
            )


def _has_non_missing_text(value: object) -> bool:
    return not _is_missing_required_text(value) and not _is_not_provided(value)


def _is_missing_environment_text(value: object) -> bool:
    return (
        _is_missing_required_text(value)
        or str(value).strip().lower() in _MISSING_ENVIRONMENT_VALUES
    )


def _is_missing_required_text(value: object) -> bool:
    return value is None or str(value).strip() == ""


def _is_not_provided(value: object) -> bool:
    return str(value).strip().lower() in _NOT_PROVIDED_VALUES


def _normalize_sample_order(
    sample_order: Sequence[str],
    *,
    context: str = "linear_residualize_batch",
) -> tuple[str, ...]:
    samples = tuple(str(sample).strip() for sample in sample_order)
    blank_positions = [
        position for position, sample in enumerate(samples) if sample == ""
    ]
    if blank_positions:
        raise PhosPyInputError(
            f"{context} sample_order contains blank sample labels at "
            f"positions {_format_positions(blank_positions)}"
        )
    if len(set(samples)) != len(samples):
        duplicates = list(
            dict.fromkeys(sample for sample in samples if samples.count(sample) > 1)
        )
        raise PhosPyInputError(
            f"{context} sample_order contains duplicate sample "
            f"labels: {_format_labels(duplicates)}"
        )
    return samples


def _normalize_sample_index(index: pd.Index, *, field_name: str) -> pd.Index:
    normalized: list[str] = []
    missing_positions: list[int] = []
    blank_positions: list[int] = []
    for position, value in enumerate(index.tolist()):
        if _is_missing_value(value):
            missing_positions.append(position)
            continue
        label = str(value).strip()
        if label == "":
            blank_positions.append(position)
            continue
        normalized.append(label)
    if missing_positions:
        raise PhosPyInputError(
            f"{field_name} must not contain missing sample labels; positions: "
            f"{_format_positions(missing_positions)}"
        )
    if blank_positions:
        raise PhosPyInputError(
            f"{field_name} must not contain blank sample labels; positions: "
            f"{_format_positions(blank_positions)}"
        )
    return pd.Index(normalized, name=index.name)


def _require_unique_sample_index(index: pd.Index, *, field_name: str) -> None:
    if index.is_unique:
        return
    seen: set[str] = set()
    duplicates: list[str] = []
    for value in list(index):
        label = str(value)
        if label in seen and label not in duplicates:
            duplicates.append(label)
        seen.add(label)
    raise PhosPyInputError(
        f"{field_name} contains duplicate sample labels: {_format_labels(duplicates)}"
    )


def _require_unique_metadata_columns(
    sample_metadata: pd.DataFrame,
    *,
    context: str,
) -> None:
    if sample_metadata.columns.is_unique:
        return
    seen: set[str] = set()
    duplicate_columns: list[str] = []
    for value in list(sample_metadata.columns):
        column = str(value)
        if column in seen and column not in duplicate_columns:
            duplicate_columns.append(column)
        seen.add(column)
    raise PhosPyInputError(
        f"{context} sample_metadata contains duplicated columns: "
        f"{_format_labels(duplicate_columns)}"
    )


def _require_metadata_column(
    sample_metadata: pd.DataFrame,
    *,
    column: str,
    context: str,
) -> None:
    if column == "":
        raise PhosPyInputError(f"{context} metadata column names must be non-empty")
    matches = [
        candidate for candidate in sample_metadata.columns if candidate == column
    ]
    if len(matches) == 1:
        return
    if len(matches) > 1:
        raise PhosPyInputError(
            f"{context} sample_metadata contains duplicated column {column!r}"
        )
    raise PhosPyInputError(
        f"{context} sample_metadata is missing required column {column!r}"
    )


def _resolve_column_labels(
    sample_metadata: pd.DataFrame,
    *,
    sample_order: tuple[str, ...],
    column: str,
    label_kind: str,
    context: str,
) -> dict[str, str]:
    labels: dict[str, str] = {}
    missing_samples: list[str] = []
    blank_samples: list[str] = []
    for sample in sample_order:
        value = sample_metadata.at[sample, column]
        if _is_missing_value(value):
            missing_samples.append(sample)
            continue
        label = str(value).strip()
        if label == "":
            blank_samples.append(sample)
            continue
        labels[sample] = label
    if missing_samples:
        raise PhosPyInputError(
            f"{context} sample_metadata column {column!r} contains missing "
            f"{label_kind} labels for samples: {_format_labels(missing_samples)}"
        )
    if blank_samples:
        raise PhosPyInputError(
            f"{context} sample_metadata column {column!r} contains blank "
            f"{label_kind} labels for samples: {_format_labels(blank_samples)}"
        )
    return labels


def _resolve_labels(
    labels_by_sample: Mapping[str, object],
    *,
    sample_order: tuple[str, ...],
    label_kind: str,
    context: str = "linear_residualize_batch",
) -> tuple[str, ...]:
    labels: list[str] = []
    missing_samples: list[str] = []
    blank_samples: list[str] = []
    for sample in sample_order:
        if sample not in labels_by_sample:
            missing_samples.append(sample)
            continue
        value = labels_by_sample[sample]
        if _is_missing_value(value):
            missing_samples.append(sample)
            continue
        label = str(value).strip()
        if label == "":
            blank_samples.append(sample)
            continue
        labels.append(label)

    if missing_samples:
        raise PhosPyInputError(
            f"{context} requires {label_kind} labels for every "
            f"sample; missing {label_kind} labels for samples: "
            f"{_format_labels(missing_samples)}"
        )
    if blank_samples:
        raise PhosPyInputError(
            f"{context} requires {label_kind} labels for every "
            f"sample; blank {label_kind} labels for samples: "
            f"{_format_labels(blank_samples)}"
        )
    return tuple(labels)


def _levels_in_order(labels: Sequence[str]) -> tuple[str, ...]:
    levels: list[str] = []
    seen: set[str] = set()
    for label in labels:
        if label in seen:
            continue
        seen.add(label)
        levels.append(label)
    return tuple(levels)


def _unique_in_order(labels: Sequence[str]) -> tuple[str, ...]:
    return _levels_in_order(labels)


def _duplicates_in_order(labels: Sequence[str]) -> tuple[str, ...]:
    counts: dict[str, int] = {}
    for label in labels:
        counts[label] = counts.get(label, 0) + 1
    return tuple(label for label in _levels_in_order(labels) if counts[label] > 1)


def _singleton_levels(labels: Sequence[str]) -> tuple[str, ...]:
    counts: dict[str, int] = {}
    for label in labels:
        counts[label] = counts.get(label, 0) + 1
    return tuple(label for label in _levels_in_order(labels) if counts[label] == 1)


def _treatment_coded_design(
    labels: Sequence[str],
    *,
    include_intercept: bool,
) -> np.ndarray:
    levels = _levels_in_order(labels)
    row_width = (1 if include_intercept else 0) + max(len(levels) - 1, 0)
    if row_width == 0:
        return np.empty((len(labels), 0), dtype=float)

    rows: list[list[float]] = []
    for label in labels:
        row: list[float] = []
        if include_intercept:
            row.append(1.0)
        row.extend(1.0 if label == level else 0.0 for level in levels[1:])
        rows.append(row)
    return np.asarray(rows, dtype=float)


def _matrix_rank(matrix: np.ndarray) -> int:
    if matrix.size == 0:
        return 0
    return int(np.linalg.matrix_rank(matrix))


def _condition_is_determined_by_batch(
    *,
    batch_labels: Sequence[str],
    condition_labels: Sequence[str],
) -> bool:
    conditions_by_batch: dict[str, set[str]] = {}
    for batch, condition in zip(batch_labels, condition_labels, strict=True):
        conditions_by_batch.setdefault(batch, set()).add(condition)
    return all(len(conditions) == 1 for conditions in conditions_by_batch.values())


def _batch_is_determined_by_condition(
    *,
    batch_labels: Sequence[str],
    condition_labels: Sequence[str],
) -> bool:
    batches_by_condition: dict[str, set[str]] = {}
    for batch, condition in zip(batch_labels, condition_labels, strict=True):
        batches_by_condition.setdefault(condition, set()).add(batch)
    return all(len(batches) == 1 for batches in batches_by_condition.values())


def _same_label_partition(
    left_labels: Sequence[str],
    right_labels: Sequence[str],
) -> bool:
    left = tuple(left_labels)
    right = tuple(right_labels)
    if len(left) != len(right) or not left:
        return False
    for first_index in range(len(left)):
        for second_index in range(len(left)):
            same_left = left[first_index] == left[second_index]
            same_right = right[first_index] == right[second_index]
            if same_left != same_right:
                return False
    return True


def _is_missing_value(value: object) -> bool:
    return bool(pd.Series((value,), dtype="object").isna().iat[0])


def _format_positions(positions: Sequence[int]) -> str:
    preview = ", ".join(str(position) for position in tuple(positions)[:8])
    suffix = "" if len(positions) <= 8 else ", ..."
    return f"[{preview}{suffix}]"


def _format_labels(labels: Sequence[str]) -> str:
    preview = ", ".join(repr(value) for value in tuple(labels)[:5])
    suffix = "" if len(labels) <= 5 else " ..."
    return f"{preview}{suffix}"


__all__ = [
    "BatchCorrectionAdequacyValidator",
    "BatchDesignMetadataValidator",
    "BatchStructureValidator",
    "ConditionStructureValidator",
    "DesignRankValidator",
    "ReplicateStructureDiagnosticHelper",
    "ReplicateStructureDiagnostics",
    "ReplicateStructureValidator",
    "ResolvedBatchDesignMetadata",
    "SampleMetadataAlignmentValidator",
    "normalize_applied_selected_site_key_rows",
    "validate_applied_native_sps_ruv_correction_provenance",
]
