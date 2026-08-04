"""Design metadata and adequacy validation for batch correction."""
# pyright: reportUnnecessaryIsInstance=false, reportUnknownMemberType=false
# Runtime boundary guards are intentionally retained for untyped external callers.

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import cast

import numpy as np
import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.validation.datasets._batch_correction_helpers import (
    batch_is_determined_by_condition,
    condition_is_determined_by_batch,
    format_labels,
    levels_in_order,
    matrix_rank,
    normalize_sample_index,
    normalize_sample_order,
    require_metadata_column,
    require_unique_metadata_columns,
    require_unique_sample_index,
    resolve_column_labels,
    resolve_labels,
    same_label_partition,
    singleton_levels,
    treatment_coded_design,
)


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

        sample_order = normalize_sample_index(
            phospho.columns,
            field_name=f"{context} phospho.columns",
        )
        metadata_index = normalize_sample_index(
            sample_metadata.index,
            field_name=f"{context} sample_metadata.index",
        )
        require_unique_sample_index(
            sample_order,
            field_name=f"{context} phospho.columns",
        )
        require_unique_sample_index(
            metadata_index,
            field_name=f"{context} sample_metadata.index",
        )
        require_unique_metadata_columns(sample_metadata, context=context)

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
                    f"{format_labels(missing_samples)}"
                )
            if extra_samples:
                details.append(
                    "sample_metadata rows not present in matrix columns: "
                    f"{format_labels(extra_samples)}"
                )
            raise PhosPyInputError(
                f"{context} sample_metadata is misaligned with matrix columns; "
                + "; ".join(details)
            )

        for column in required_columns:
            require_metadata_column(
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
        require_metadata_column(sample_metadata, column=column, context=context)
        return resolve_column_labels(
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
            require_metadata_column(sample_metadata, column=column, context=context)

        sample_tuple = tuple(sample_order)
        resolved_by_column = tuple(
            resolve_column_labels(
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
        require_metadata_column(sample_metadata, column=column, context=context)
        return resolve_column_labels(
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
        replicate_levels = levels_in_order(replicate_labels)
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
        batch_confounded = same_label_partition(replicate_labels, batch_labels)
        condition_confounded = same_label_partition(
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
        aligned = cast(pd.DataFrame, metadata_frame.copy(deep=False))  # pyright: ignore[reportUnnecessaryCast] - retained for pandas-stubs compatibility across supported targets.
        aligned.index = normalize_sample_index(
            metadata_frame.index,
            field_name=f"{context} sample_metadata.index",
        )
        aligned = cast(pd.DataFrame, aligned.reindex(sample_order))  # pyright: ignore[reportUnnecessaryCast] - retained for pandas-stubs compatibility across supported targets.
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
        rank = matrix_rank(values)
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
        samples = normalize_sample_order(sample_order, context=context)
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

        batch_labels = resolve_labels(
            batch_by_sample,
            sample_order=samples,
            label_kind="batch",
            context=context,
        )
        condition_labels = resolve_labels(
            condition_by_sample,
            sample_order=samples,
            label_kind="condition",
            context=context,
        )
        batch_levels = levels_in_order(batch_labels)
        condition_levels = levels_in_order(condition_labels)
        if len(batch_levels) < 2:
            raise PhosPyInputError(
                f"{context} requires at least two batch levels; "
                f"observed {len(batch_levels)}"
            )

        singleton_batches = singleton_levels(batch_labels)
        if singleton_batches:
            raise PhosPyInputError(
                f"{context} requires at least two samples in each "
                "batch level to estimate batch effects; singleton batch levels: "
                f"{format_labels(singleton_batches)}"
            )

        preservation_design = treatment_coded_design(
            condition_labels,
            include_intercept=True,
        )
        preservation_columns = int(preservation_design.shape[1])
        preservation_rank = matrix_rank(preservation_design)
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

        batch_terms = treatment_coded_design(batch_labels, include_intercept=False)
        full_design = np.concatenate((preservation_design, batch_terms), axis=1)
        full_columns = int(full_design.shape[1])
        if len(samples) <= full_columns:
            raise PhosPyInputError(
                f"{context} requires more samples than estimable "
                "condition-plus-batch design parameters; "
                f"samples={len(samples)}, design_columns={full_columns}. Add "
                "replicate samples or reduce batch/condition levels."
            )

        if len(condition_levels) > 1 and condition_is_determined_by_batch(
            batch_labels=batch_labels,
            condition_labels=condition_labels,
        ):
            raise PhosPyInputError(
                f"{context} cannot run because batch and condition "
                "are perfectly confounded: each batch level contains only one "
                "condition level, so removing batch would remove biological "
                "condition signal"
            )
        if len(condition_levels) > 1 and batch_is_determined_by_condition(
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
                f"(rank={matrix_rank(full_design)}, columns={full_columns})"
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
]
