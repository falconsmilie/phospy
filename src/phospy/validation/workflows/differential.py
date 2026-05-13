"""Experimental-design validation for differential workflows."""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import cast

import pandas as pd

from phospy.errors.validation import WorkflowValidationError
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.science.design.models import (
    Contrast,
    ExperimentalDesign,
    SampleDesignRecord,
)
from phospy.science.transformations.models import IntensityScaleKind

_DIFFERENTIAL_LOGFC_SCALE_ERROR_MESSAGE = (
    "Differential analysis reports logFC and therefore requires established "
    "log2-scale phospho intensities. Build the dataset with log2 preprocessing "
    "enabled, or provide an analysis-ready dataset with validated log2 "
    "intensity-scale state."
)


@dataclass(frozen=True, slots=True)
class ValidatedExperimentalDesignContract:
    """Validated design contract resolved into matrix-ready DataFrames."""

    design: ExperimentalDesign
    contrasts: tuple[Contrast, ...]
    analysis_sample_ids: tuple[str, ...]
    condition_labels: tuple[str, ...]
    design_frame: pd.DataFrame
    contrast_frame: pd.DataFrame


class ExperimentalDesignContractValidator:
    """Validate typed experimental design and contrast definitions."""

    def run(
        self,
        *,
        dataset: AnalysisReadyPhosphoDataset,
        design: ExperimentalDesign,
        contrasts: tuple[Contrast, ...],
        allow_design_subset: bool,
        minimum_condition_replicates: int,
    ) -> ValidatedExperimentalDesignContract:
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

        normalized_contrasts = self._validate_contrasts(contrasts)
        if not normalized_contrasts:
            raise WorkflowValidationError(
                "experimental design contrasts must include at least one contrast"
            )

        dataset_sample_ids = tuple(
            str(label)
            for label in dataset._borrow_phospho_frame().columns  # pyright: ignore[reportPrivateUsage] - workflow boundary reads trusted internal dataset snapshots
        )
        design_sample_ids = design.sample_ids()
        dataset_sample_set = set(dataset_sample_ids)
        design_sample_set = set(design_sample_ids)

        missing_samples = sorted(dataset_sample_set - design_sample_set)
        extra_samples = sorted(design_sample_set - dataset_sample_set)
        if extra_samples:
            raise WorkflowValidationError(
                "experimental design contains samples not present in dataset: "
                + ", ".join(extra_samples)
            )
        if missing_samples and not allow_design_subset:
            raise WorkflowValidationError(
                "experimental design is missing required dataset samples: "
                + ", ".join(missing_samples)
            )

        self._validate_optional_field_alignment(
            design.samples,
            field_name="batch",
        )
        self._validate_optional_field_alignment(
            design.samples,
            field_name="block",
        )

        if any(record.batch is not None for record in design.samples):
            raise WorkflowValidationError(
                "unsupported design features: batch-aware differential modelling is "
                "not available in this release"
            )
        if any(record.block is not None for record in design.samples):
            raise WorkflowValidationError(
                "unsupported design features: blocking/paired differential modelling "
                "is not available in this release"
            )

        condition_to_records = self._records_by_condition(design.samples)
        known_conditions = tuple(condition_to_records.keys())
        if len(known_conditions) < 2:
            raise WorkflowValidationError(
                "experimental design must contain at least two conditions"
            )

        for contrast in normalized_contrasts:
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

        design_frame = self._build_design_frame(
            design=design,
            condition_labels=known_conditions,
        )
        contrast_frame = self._build_contrast_frame(
            condition_labels=known_conditions,
            contrasts=normalized_contrasts,
        )
        return ValidatedExperimentalDesignContract(
            design=design,
            contrasts=normalized_contrasts,
            analysis_sample_ids=design_sample_ids,
            condition_labels=known_conditions,
            design_frame=design_frame,
            contrast_frame=contrast_frame,
        )

    @staticmethod
    def _validate_contrasts(contrasts: tuple[Contrast, ...]) -> tuple[Contrast, ...]:
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
        return normalized

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

    @staticmethod
    def _build_design_frame(
        *,
        design: ExperimentalDesign,
        condition_labels: tuple[str, ...],
    ) -> pd.DataFrame:
        sample_ids = design.sample_ids()
        frame = pd.DataFrame(
            0.0,
            index=pd.Index(sample_ids, name="sample"),
            columns=pd.Index(condition_labels, name="coefficient"),
        )
        for record in design.samples:
            frame.at[record.sample_id, record.condition] = 1.0
        return frame

    @staticmethod
    def _build_contrast_frame(
        *,
        condition_labels: tuple[str, ...],
        contrasts: tuple[Contrast, ...],
    ) -> pd.DataFrame:
        frame = pd.DataFrame(
            0.0,
            index=pd.Index(condition_labels, name="coefficient"),
            columns=pd.Index(
                [contrast.name for contrast in contrasts], name="contrast"
            ),
        )
        for contrast in contrasts:
            frame.at[contrast.denominator_condition, contrast.name] = -1.0
            frame.at[contrast.numerator_condition, contrast.name] = 1.0
        return frame


class DifferentialDatasetEligibilityValidator:
    """Validate dataset quantitative-scale eligibility for differential logFC."""

    def run(self, *, dataset: AnalysisReadyPhosphoDataset) -> None:
        if not isinstance(cast(object, dataset), AnalysisReadyPhosphoDataset):
            raise WorkflowValidationError(
                "differential workflow request dataset must be AnalysisReadyPhosphoDataset"
            )
        phospho_scale = dataset.intensity_scale_state.phospho
        if (
            not dataset.intensity_scale_state.is_established
            or phospho_scale.kind is not IntensityScaleKind.LOG2
        ):
            raise WorkflowValidationError(_DIFFERENTIAL_LOGFC_SCALE_ERROR_MESSAGE)


__all__ = [
    "DifferentialDatasetEligibilityValidator",
    "ExperimentalDesignContractValidator",
    "ValidatedExperimentalDesignContract",
]
