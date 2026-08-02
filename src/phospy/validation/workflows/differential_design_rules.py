"""Rule-family validators for differential experimental-design contracts."""

from __future__ import annotations

import math
from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import cast

import numpy as np
import numpy.typing as npt
import pandas as pd

from phospy.contracts.configs.differential import (
    PAIRED_DESIGN_POLICY_REJECT,
    SUPPORTED_PAIRED_DESIGN_POLICIES,
    PairedDesignPolicy,
)
from phospy.errors.validation import WorkflowValidationError
from phospy.science.datasets.internal_view import DatasetInternalView
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.science.design.models import (
    FIXED_EFFECT_COVARIATE_KIND_BATCH,
    FIXED_EFFECT_COVARIATE_KIND_CATEGORICAL,
    FIXED_EFFECT_COVARIATE_KIND_CONTINUOUS,
    Contrast,
    ExperimentalDesign,
    SampleDesignRecord,
)
from phospy.science.differential.linear_model import (
    DifferentialDesignDecomposition,
    DifferentialDesignDecompositionError,
    decompose_differential_design,
)


@dataclass(frozen=True, slots=True)
class ExperimentalDesignInputValidation:
    """Validated scalar request policy for experimental-design validation."""

    fixed_block_requested: bool


@dataclass(frozen=True, slots=True)
class ExperimentalDesignSampleAlignment:
    """Dataset/design sample alignment facts consumed by the coordinator."""

    dataset_sample_ids: tuple[str, ...]
    design_sample_ids: tuple[str, ...]
    missing_samples: tuple[str, ...]
    extra_samples: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class ExperimentalDesignConditionReplicateValidation:
    """Condition grouping and replicate-count validation facts."""

    condition_to_records: dict[str, tuple[SampleDesignRecord, ...]]
    condition_labels: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class FixedBlockDesignValidation:
    """Fixed-block coverage facts for a validated paired design."""

    block_ids: tuple[str, ...]


class ExperimentalDesignInputValidator:
    """Validate request-level experimental-design input types and policy."""

    def run(
        self,
        *,
        dataset: AnalysisReadyPhosphoDataset,
        design: ExperimentalDesign,
        allow_design_subset: bool,
        minimum_condition_replicates: int,
        paired_design_policy: PairedDesignPolicy,
    ) -> ExperimentalDesignInputValidation:
        if not isinstance(cast(object, dataset), AnalysisReadyPhosphoDataset):
            raise WorkflowValidationError(
                "differential workflow request dataset must be AnalysisReadyPhosphoDataset"
            )
        if not isinstance(cast(object, design), ExperimentalDesign):
            raise WorkflowValidationError(
                "differential workflow request design must be ExperimentalDesign"
            )
        if not isinstance(cast(object, allow_design_subset), bool):
            raise WorkflowValidationError(
                "differential workflow request allow_design_subset must be a bool"
            )
        if not isinstance(cast(object, minimum_condition_replicates), int):
            raise WorkflowValidationError(
                "differential workflow request minimum_condition_replicates must be an int"
            )
        if minimum_condition_replicates < 1:
            raise WorkflowValidationError(
                "differential workflow request minimum_condition_replicates must be >= 1"
            )
        return ExperimentalDesignInputValidation(
            fixed_block_requested=self._validate_paired_design_policy(
                records=design.samples,
                paired_design_policy=paired_design_policy,
            )
        )

    @staticmethod
    def _validate_paired_design_policy(
        *,
        records: tuple[SampleDesignRecord, ...],
        paired_design_policy: PairedDesignPolicy,
    ) -> bool:
        if paired_design_policy not in SUPPORTED_PAIRED_DESIGN_POLICIES:
            supported = ", ".join(
                repr(value) for value in SUPPORTED_PAIRED_DESIGN_POLICIES
            )
            raise WorkflowValidationError(
                "differential.paired_design_policy must be one of: " + supported
            )

        samples_with_block_id = [
            record.sample_id for record in records if record.block_id is not None
        ]
        samples_missing_block_id = [
            record.sample_id for record in records if record.block_id is None
        ]

        if paired_design_policy == PAIRED_DESIGN_POLICY_REJECT:
            if samples_with_block_id:
                raise WorkflowValidationError(
                    "experimental design includes block_id values while "
                    "differential.paired_design_policy='reject'. Set "
                    "differential.paired_design_policy='fixed_block' to request "
                    "fixed-block validation and execution. Samples with block_id: "
                    + ", ".join(samples_with_block_id)
                )
            return False

        if samples_missing_block_id:
            raise WorkflowValidationError(
                "differential.paired_design_policy='fixed_block' requires block_id "
                "for every design sample; missing block_id for samples: "
                + ", ".join(samples_missing_block_id)
            )
        return True


class ExperimentalDesignContrastSetValidator:
    """Validate contrast container shape and stable contrast names."""

    def run(self, contrasts: tuple[Contrast, ...]) -> tuple[Contrast, ...]:
        normalized = tuple(contrasts)
        for contrast in normalized:
            if not isinstance(cast(object, contrast), Contrast):
                raise WorkflowValidationError(
                    "differential workflow request contrasts must contain Contrast values"
                )
        duplicate_names = [
            name
            for name, count in Counter(contrast.name for contrast in normalized).items()
            if count > 1
        ]
        if duplicate_names:
            raise WorkflowValidationError(
                "experimental design contains duplicate contrast names: "
                + ", ".join(sorted(duplicate_names))
            )
        if not normalized:
            raise WorkflowValidationError(
                "experimental design contrasts must include at least one contrast"
            )
        return normalized


class ExperimentalDesignSampleAlignmentValidator:
    """Validate dataset/design sample-set compatibility."""

    def run(
        self,
        *,
        dataset: AnalysisReadyPhosphoDataset,
        design: ExperimentalDesign,
        allow_design_subset: bool,
        fixed_block_requested: bool,
        dataset_view: DatasetInternalView | None = None,
    ) -> ExperimentalDesignSampleAlignment:
        resolved_dataset_view = dataset_view or DatasetInternalView(dataset)
        dataset_sample_ids = tuple(
            str(label) for label in resolved_dataset_view.phospho.columns
        )
        design_sample_ids = design.sample_ids()
        dataset_sample_set = set(dataset_sample_ids)
        design_sample_set = set(design_sample_ids)

        missing_samples = tuple(sorted(dataset_sample_set - design_sample_set))
        extra_samples = tuple(sorted(design_sample_set - dataset_sample_set))
        if extra_samples:
            raise WorkflowValidationError(
                "experimental design contains samples not present in dataset: "
                + ", ".join(extra_samples)
            )
        if fixed_block_requested and missing_samples:
            raise WorkflowValidationError(
                "differential.paired_design_policy='fixed_block' requires every "
                "dataset sample to have an explicit design row; no samples are "
                "silently dropped from fixed-block designs. Missing design rows "
                "for dataset samples: " + ", ".join(missing_samples)
            )
        if missing_samples and not allow_design_subset:
            raise WorkflowValidationError(
                "experimental design is missing required dataset samples: "
                + ", ".join(missing_samples)
            )
        return ExperimentalDesignSampleAlignment(
            dataset_sample_ids=dataset_sample_ids,
            design_sample_ids=design_sample_ids,
            missing_samples=missing_samples,
            extra_samples=extra_samples,
        )


class ExperimentalDesignFixedEffectValidator:
    """Validate optional field alignment and declared fixed-effect covariates."""

    def run(self, *, design: ExperimentalDesign) -> None:
        declares_batch = any(
            covariate.kind == FIXED_EFFECT_COVARIATE_KIND_BATCH
            for covariate in design.fixed_effects
        )
        if not declares_batch:
            self._validate_optional_field_alignment(
                design.samples,
                field_name="batch",
            )
        self._validate_optional_field_alignment(
            design.samples,
            field_name="block_id",
        )
        self._validate_fixed_effect_inputs(design)

    @staticmethod
    def _validate_optional_field_alignment(
        records: tuple[SampleDesignRecord, ...],
        *,
        field_name: str,
    ) -> None:
        values = [getattr(record, field_name) for record in records]
        has_value = [value is not None for value in values]
        if any(has_value) and not all(has_value):
            raise WorkflowValidationError(
                f"experimental design optional field {field_name!r} must be set for "
                "all samples or for no samples"
            )

    @staticmethod
    def _validate_fixed_effect_inputs(design: ExperimentalDesign) -> None:
        declared_covariate_names = {
            covariate.name
            for covariate in design.fixed_effects
            if covariate.kind != FIXED_EFFECT_COVARIATE_KIND_BATCH
        }
        declares_batch = any(
            covariate.kind == FIXED_EFFECT_COVARIATE_KIND_BATCH
            for covariate in design.fixed_effects
        )
        for record in design.samples:
            extra_names = sorted(set(record.covariates) - declared_covariate_names)
            if extra_names:
                raise WorkflowValidationError(
                    "experimental design sample covariates must be explicitly "
                    "declared as fixed effects; "
                    f"sample={record.sample_id!r}, undeclared=" + ", ".join(extra_names)
                )
        if (
            any(record.batch is not None for record in design.samples)
            and not declares_batch
        ):
            raise WorkflowValidationError(
                "experimental design batch values must be explicitly declared as a "
                "batch fixed effect"
            )

        for covariate in design.fixed_effects:
            missing_samples: list[str] = []
            non_numeric_samples: list[str] = []
            non_finite_samples: list[str] = []
            invalid_level_samples: list[str] = []
            observed_levels: set[str] = set()
            for record in design.samples:
                if covariate.kind == FIXED_EFFECT_COVARIATE_KIND_BATCH:
                    value = record.batch
                    missing = value is None
                else:
                    missing = covariate.name not in record.covariates
                    value = None if missing else record.covariates[covariate.name]
                if missing:
                    missing_samples.append(record.sample_id)
                    continue
                if covariate.kind == FIXED_EFFECT_COVARIATE_KIND_CONTINUOUS:
                    if isinstance(value, bool) or not isinstance(value, int | float):
                        non_numeric_samples.append(record.sample_id)
                        continue
                    if not math.isfinite(float(value)):
                        non_finite_samples.append(record.sample_id)
                    continue
                if covariate.kind in {
                    FIXED_EFFECT_COVARIATE_KIND_BATCH,
                    FIXED_EFFECT_COVARIATE_KIND_CATEGORICAL,
                }:
                    level = _normalize_categorical_level_for_validation(value)
                    if level is None:
                        invalid_level_samples.append(record.sample_id)
                    else:
                        observed_levels.add(level)
            if missing_samples and (covariate.required or covariate.include_in_model):
                requirement = (
                    "is required for modelling"
                    if covariate.include_in_model
                    else "is required"
                )
                raise WorkflowValidationError(
                    "experimental design fixed-effect covariate "
                    f"{covariate.name!r} {requirement} but missing for samples: "
                    + ", ".join(missing_samples)
                )
            if non_numeric_samples:
                raise WorkflowValidationError(
                    "experimental design continuous covariate "
                    f"{covariate.name!r} must be numeric for samples: "
                    + ", ".join(non_numeric_samples)
                )
            if non_finite_samples:
                raise WorkflowValidationError(
                    "experimental design continuous covariate "
                    f"{covariate.name!r} must contain finite values for samples: "
                    + ", ".join(non_finite_samples)
                )
            if invalid_level_samples:
                raise WorkflowValidationError(
                    "experimental design categorical covariate "
                    f"{covariate.name!r} must have non-empty finite string or "
                    "numeric levels for samples: " + ", ".join(invalid_level_samples)
                )
            if (
                covariate.include_in_model
                and covariate.kind
                in {
                    FIXED_EFFECT_COVARIATE_KIND_BATCH,
                    FIXED_EFFECT_COVARIATE_KIND_CATEGORICAL,
                }
                and len(observed_levels) < 2
            ):
                observed = ", ".join(sorted(observed_levels)) or "<none>"
                raise WorkflowValidationError(
                    "experimental design categorical fixed-effect covariate "
                    f"{covariate.name!r} must have at least two observed levels "
                    "when included in the differential model; observed levels: "
                    f"{observed}"
                )


class ExperimentalDesignConditionReplicateValidator:
    """Validate condition presence, contrast references, and replicate counts."""

    def run(
        self,
        *,
        design: ExperimentalDesign,
        contrasts: tuple[Contrast, ...],
        minimum_condition_replicates: int,
    ) -> ExperimentalDesignConditionReplicateValidation:
        condition_to_records = self._records_by_condition(design.samples)
        known_conditions = tuple(condition_to_records.keys())
        if len(known_conditions) < 2:
            raise WorkflowValidationError(
                "experimental design must contain at least two conditions"
            )

        for contrast in contrasts:
            if contrast.numerator_condition not in condition_to_records:
                raise WorkflowValidationError(
                    "contrast references unknown numerator condition: "
                    f"{contrast.numerator_condition!r}"
                )
            if contrast.denominator_condition not in condition_to_records:
                raise WorkflowValidationError(
                    "contrast references unknown denominator condition: "
                    f"{contrast.denominator_condition!r}"
                )
            for condition in (
                contrast.numerator_condition,
                contrast.denominator_condition,
            ):
                replicate_count = self._count_replicates_for_condition(
                    condition_to_records[condition]
                )
                if replicate_count < minimum_condition_replicates:
                    raise WorkflowValidationError(
                        "insufficient replicate counts for contrast "
                        f"{contrast.name!r}: condition={condition!r}, "
                        f"replicates={replicate_count}, required="
                        f"{minimum_condition_replicates}"
                    )

        return ExperimentalDesignConditionReplicateValidation(
            condition_to_records=condition_to_records,
            condition_labels=known_conditions,
        )

    @staticmethod
    def _records_by_condition(
        records: tuple[SampleDesignRecord, ...],
    ) -> dict[str, tuple[SampleDesignRecord, ...]]:
        grouped: defaultdict[str, list[SampleDesignRecord]] = defaultdict(list)
        for record in records:
            grouped[record.condition].append(record)
        return {condition: tuple(values) for condition, values in grouped.items()}

    @staticmethod
    def _count_replicates_for_condition(
        records: tuple[SampleDesignRecord, ...],
    ) -> int:
        biological_ids = [record.biological_replicate_id for record in records]
        if any(value is not None for value in biological_ids):
            if not all(value is not None for value in biological_ids):
                raise WorkflowValidationError(
                    "biological_replicate_id must be set for all samples within each "
                    "condition when replicate-aware counting is requested"
                )
            return len({str(value) for value in biological_ids if value is not None})
        return len(records)


class FixedBlockDesignValidator:
    """Validate fixed-block differential design scientific eligibility."""

    def run(
        self,
        *,
        records: tuple[SampleDesignRecord, ...],
        contrasts: tuple[Contrast, ...],
    ) -> FixedBlockDesignValidation:
        records_by_block: defaultdict[str, list[SampleDesignRecord]] = defaultdict(list)
        for record in records:
            if record.block_id is None:
                raise WorkflowValidationError(
                    "differential.paired_design_policy='fixed_block' requires "
                    "block_id for every design sample; missing block_id for "
                    f"samples: {record.sample_id}"
                )
            records_by_block[record.block_id].append(record)

        incomplete_blocks: list[str] = []
        for block_id in sorted(records_by_block):
            block_records = records_by_block[block_id]
            if len(block_records) >= 2:
                continue
            sample_ids = ", ".join(record.sample_id for record in block_records)
            incomplete_blocks.append(f"{block_id} (samples: {sample_ids})")
        if incomplete_blocks:
            raise WorkflowValidationError(
                "differential.paired_design_policy='fixed_block' requires every "
                "block to contain at least 2 samples; incomplete blocks: "
                + "; ".join(incomplete_blocks)
            )

        condition_sets_by_block = {
            block_id: {record.condition for record in block_records}
            for block_id, block_records in records_by_block.items()
        }
        blocks_by_condition: defaultdict[str, set[str]] = defaultdict(set)
        for block_id, condition_set in condition_sets_by_block.items():
            for condition in condition_set:
                blocks_by_condition[condition].add(block_id)
        confounded_conditions = [
            condition
            for condition in sorted(blocks_by_condition)
            if all(
                condition_sets_by_block[block_id] == {condition}
                for block_id in blocks_by_condition[condition]
            )
        ]
        if confounded_conditions:
            raise WorkflowValidationError(
                "differential.paired_design_policy='fixed_block' found condition "
                "perfectly confounded with block; each condition must occur in at "
                "least one block that also contains another condition. Confounded "
                "conditions: " + ", ".join(confounded_conditions)
            )

        for contrast in contrasts:
            numerator = contrast.numerator_condition
            denominator = contrast.denominator_condition
            numerator_blocks = {
                block_id
                for block_id, condition_set in condition_sets_by_block.items()
                if numerator in condition_set
            }
            denominator_blocks = {
                block_id
                for block_id, condition_set in condition_sets_by_block.items()
                if denominator in condition_set
            }
            shared_blocks = sorted(numerator_blocks & denominator_blocks)
            if not shared_blocks:
                raise WorkflowValidationError(
                    "differential.paired_design_policy='fixed_block' contrast "
                    f"{contrast.name!r} is non-estimable because condition is "
                    "perfectly confounded with block; no block contains both "
                    f"numerator condition {numerator!r} and denominator condition "
                    f"{denominator!r}"
                )

            invalid_blocks: list[str] = []
            for block_id in sorted(condition_sets_by_block):
                condition_set = condition_sets_by_block[block_id]
                missing_terms: list[str] = []
                if numerator not in condition_set:
                    missing_terms.append(f"numerator condition {numerator!r}")
                if denominator not in condition_set:
                    missing_terms.append(f"denominator condition {denominator!r}")
                if not missing_terms:
                    continue
                observed = ", ".join(sorted(condition_set)) or "<none>"
                missing = " and ".join(missing_terms)
                invalid_blocks.append(
                    f"{block_id} missing {missing} (observed conditions: {observed})"
                )
            if invalid_blocks:
                raise WorkflowValidationError(
                    "differential.paired_design_policy='fixed_block' contrast "
                    f"{contrast.name!r} has invalid block condition coverage; "
                    "every block must contain both "
                    f"numerator condition {numerator!r} and denominator condition "
                    f"{denominator!r}. Blocks: " + "; ".join(invalid_blocks)
                )

        return FixedBlockDesignValidation(block_ids=tuple(sorted(records_by_block)))


class ContrastFrameBuilder:
    """Build contrast vectors aligned to resolved design coefficients."""

    def run(
        self,
        *,
        coefficient_labels: tuple[str, ...],
        condition_labels: tuple[str, ...],
        contrasts: tuple[Contrast, ...],
    ) -> pd.DataFrame:
        missing_condition_coefficients = [
            condition
            for condition in condition_labels
            if condition not in set(coefficient_labels)
        ]
        if missing_condition_coefficients:
            raise WorkflowValidationError(
                "experimental design contrast vectors are invalid because the "
                "design matrix is missing condition coefficients: "
                + ", ".join(missing_condition_coefficients)
            )
        frame = pd.DataFrame(
            0.0,
            index=pd.Index(coefficient_labels, name="coefficient"),
            columns=pd.Index(
                [contrast.name for contrast in contrasts], name="contrast"
            ),
        )
        for contrast in contrasts:
            frame.at[contrast.denominator_condition, contrast.name] = -1.0
            frame.at[contrast.numerator_condition, contrast.name] = 1.0
        return frame


class ResolvedDifferentialDesignMatrixValidator:
    """Validate matrix numeric domain, rank, conditioning, and estimability."""

    def run(
        self,
        *,
        design_frame: pd.DataFrame,
        contrast_frame: pd.DataFrame,
    ) -> DifferentialDesignDecomposition:
        if not design_frame.index.is_unique:
            raise WorkflowValidationError(
                "experimental design matrix sample labels must be unique"
            )
        if not design_frame.columns.is_unique:
            raise WorkflowValidationError(
                "experimental design matrix coefficient labels must be unique"
            )
        if not contrast_frame.index.equals(design_frame.columns):
            raise WorkflowValidationError(
                "experimental design contrast vectors must be aligned to every "
                "design coefficient, including fixed-effect coefficients"
            )
        try:
            design_values: npt.NDArray[np.float64] = np.asarray(
                design_frame.to_numpy(dtype=float),
                dtype=np.float64,
            )
            contrast_values: npt.NDArray[np.float64] = np.asarray(
                contrast_frame.to_numpy(dtype=float),
                dtype=np.float64,
            )
        except (TypeError, ValueError) as exc:
            raise WorkflowValidationError(
                "experimental design and contrast matrices must contain numeric values"
            ) from exc
        if not np.isfinite(design_values).all():
            raise WorkflowValidationError(
                "experimental design matrix must contain only finite numeric values"
            )
        if not np.isfinite(contrast_values).all():
            raise WorkflowValidationError(
                "experimental design contrast vectors must contain only finite "
                "numeric values"
            )

        try:
            design_decomposition = decompose_differential_design(design_values)
        except DifferentialDesignDecompositionError as exc:
            coefficients = ", ".join(str(label) for label in design_frame.columns)
            raise WorkflowValidationError(
                "experimental design matrix is rank deficient or too "
                "ill-conditioned for stable differential linear-model fitting; "
                "condition and fixed-effect terms may be collinear or confounded. "
                "Remove redundant/confounded covariates or respecify the design "
                "before running differential analysis; coefficients="
                f"{coefficients}; {exc}"
            ) from exc

        invalid_positions = design_decomposition.invalid_contrast_positions(
            contrast_values
        )
        if invalid_positions:
            invalid_contrasts = [
                str(contrast_frame.columns[int(position)])
                for position in invalid_positions
            ]
            raise WorkflowValidationError(
                "experimental design contrast vectors are not estimable under the "
                "resolved fixed-effect design matrix; update contrasts or remove "
                "collinear design terms. Invalid contrasts: "
                + ", ".join(invalid_contrasts)
            )
        return design_decomposition


def _normalize_categorical_level_for_validation(value: object) -> str | None:
    if isinstance(value, bool) or not isinstance(value, str | int | float):
        return None
    if isinstance(value, str):
        level = value.strip()
        return level or None
    if not math.isfinite(float(value)):
        return None
    return str(value).strip() or None


__all__ = [
    "ContrastFrameBuilder",
    "ExperimentalDesignConditionReplicateValidation",
    "ExperimentalDesignConditionReplicateValidator",
    "ExperimentalDesignContrastSetValidator",
    "ExperimentalDesignFixedEffectValidator",
    "ExperimentalDesignInputValidation",
    "ExperimentalDesignInputValidator",
    "ExperimentalDesignSampleAlignment",
    "ExperimentalDesignSampleAlignmentValidator",
    "FixedBlockDesignValidation",
    "FixedBlockDesignValidator",
    "ResolvedDifferentialDesignMatrixValidator",
]
