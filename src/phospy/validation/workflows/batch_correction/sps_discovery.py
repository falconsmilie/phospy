"""Private validation for SPS discovery reference inputs."""

from __future__ import annotations

from collections import Counter
from collections.abc import Mapping, Sequence

import pandas as pd

from phospy.errors.validation import ReferenceValidationError
from phospy.frames.validation import (
    require_finite_numeric_dataframe,
    require_numeric_dataframe,
)
from phospy.science.batch_correction.sps_discovery import (
    SpsDiscoveryConfig,
    SpsDiscoveryValidationIssue,
    SpsDiscoveryValidationResult,
    SpsReferenceDataset,
)
from phospy.science.sites.validation import require_site_key_index


class SpsDiscoveryRequestValidator:
    """Aggregate reference failures without exposing a public validator API."""

    def run(
        self,
        *,
        reference_datasets: object,
        config: object,
    ) -> SpsDiscoveryValidationResult:
        issues: list[SpsDiscoveryValidationIssue] = []
        if not isinstance(config, SpsDiscoveryConfig):
            return SpsDiscoveryValidationResult(
                issues=(
                    SpsDiscoveryValidationIssue(
                        code="invalid_discovery_config",
                        message="config must be an SpsDiscoveryConfig",
                        field_name="config",
                    ),
                ),
                reference_dataset_count=0,
                potential_overlap_site_count=0,
            )
        if not isinstance(reference_datasets, Sequence) or isinstance(
            reference_datasets, str | bytes | bytearray
        ):
            return SpsDiscoveryValidationResult(
                issues=(
                    SpsDiscoveryValidationIssue(
                        code="invalid_reference_collection",
                        message=(
                            "reference_datasets must be a sequence of "
                            "SpsReferenceDataset values"
                        ),
                        field_name="reference_datasets",
                    ),
                ),
                reference_dataset_count=0,
                potential_overlap_site_count=0,
            )
        try:
            datasets = tuple(reference_datasets)
        except Exception as exc:
            return SpsDiscoveryValidationResult(
                issues=(
                    SpsDiscoveryValidationIssue(
                        code="invalid_reference_collection",
                        message=(
                            "reference_datasets could not be read as a stable "
                            f"sequence: {exc}"
                        ),
                        field_name="reference_datasets",
                    ),
                ),
                reference_dataset_count=0,
                potential_overlap_site_count=0,
            )

        if len(datasets) < config.minimum_reference_datasets:
            issues.append(
                SpsDiscoveryValidationIssue(
                    code="insufficient_reference_datasets",
                    message=(
                        "SPS discovery requires at least "
                        f"{config.minimum_reference_datasets} reference datasets; "
                        f"received {len(datasets)}"
                    ),
                    field_name="reference_datasets",
                )
            )

        dataset_ids: list[str] = []
        site_occurrences: Counter[str] = Counter()
        structurally_valid_dataset_count = 0
        for position, dataset in enumerate(datasets):
            if not isinstance(dataset, SpsReferenceDataset):
                issues.append(
                    SpsDiscoveryValidationIssue(
                        code="invalid_reference_dataset",
                        message=(
                            "reference dataset entries must be "
                            "SpsReferenceDataset values"
                        ),
                        field_name=f"reference_datasets[{position}]",
                    )
                )
                continue
            dataset_ids.append(dataset.dataset_id)
            dataset_issues, site_keys = self._validate_dataset(dataset)
            issues.extend(dataset_issues)
            if not dataset_issues:
                structurally_valid_dataset_count += 1
                site_occurrences.update(site_keys)

        duplicate_ids = sorted(
            dataset_id
            for dataset_id, count in Counter(dataset_ids).items()
            if count > 1
        )
        if duplicate_ids:
            issues.append(
                SpsDiscoveryValidationIssue(
                    code="duplicate_reference_dataset_id",
                    message=(
                        "reference dataset_id values must be unique; duplicates: "
                        + ", ".join(duplicate_ids)
                    ),
                    field_name="reference_datasets",
                )
            )

        potential_overlap = sum(
            count >= config.minimum_datasets_per_site
            for count in site_occurrences.values()
        )
        if (
            structurally_valid_dataset_count >= config.minimum_reference_datasets
            and potential_overlap < config.minimum_shared_sites
        ):
            issues.append(
                SpsDiscoveryValidationIssue(
                    code="insufficient_potential_overlap",
                    message=(
                        "SPS discovery found "
                        f"{potential_overlap} sites present in at least "
                        f"{config.minimum_datasets_per_site} reference datasets; "
                        f"at least {config.minimum_shared_sites} are required"
                    ),
                    field_name="reference_datasets",
                )
            )

        return SpsDiscoveryValidationResult(
            issues=tuple(issues),
            reference_dataset_count=len(datasets),
            potential_overlap_site_count=potential_overlap,
        )

    def _validate_dataset(
        self,
        dataset: SpsReferenceDataset,
    ) -> tuple[list[SpsDiscoveryValidationIssue], tuple[str, ...]]:
        issues: list[SpsDiscoveryValidationIssue] = []
        dataset_id = dataset.dataset_id
        intensities = dataset._intensities_snapshot()
        if not isinstance(intensities, pd.DataFrame):
            return (
                [
                    SpsDiscoveryValidationIssue(
                        code="missing_intensity_data",
                        message="reference intensities must be a pandas DataFrame",
                        dataset_id=dataset_id,
                        field_name="intensities",
                    )
                ],
                (),
            )
        if intensities.empty or intensities.shape[1] == 0:
            issues.append(
                SpsDiscoveryValidationIssue(
                    code="missing_intensity_data",
                    message="reference intensity data must contain sites and samples",
                    dataset_id=dataset_id,
                    field_name="intensities",
                )
            )

        if intensities.index.name != "site_key":
            issues.append(
                SpsDiscoveryValidationIssue(
                    code="missing_site_key_identity",
                    message=(
                        "reference intensity rows must use an index named 'site_key'; "
                        "display identifiers are not cross-dataset identity"
                    ),
                    dataset_id=dataset_id,
                    field_name="intensities.index",
                )
            )
        duplicate_site_keys = tuple(
            str(value)
            for value in intensities.index[intensities.index.duplicated()].unique()
        )
        if duplicate_site_keys:
            issues.append(
                SpsDiscoveryValidationIssue(
                    code="duplicate_site_key",
                    message=(
                        "reference intensity site_key values must be unique; "
                        "duplicates: " + ", ".join(duplicate_site_keys[:5])
                    ),
                    dataset_id=dataset_id,
                    field_name="intensities.index",
                )
            )

        site_keys: tuple[str, ...] = ()
        if intensities.index.name == "site_key" and not duplicate_site_keys:
            try:
                require_site_key_index(
                    intensities.index,
                    field_name=(
                        f"sps_reference_dataset[{dataset_id!r}].intensities.index"
                    ),
                    error_type=ReferenceValidationError,
                )
            except ReferenceValidationError as exc:
                issues.append(
                    SpsDiscoveryValidationIssue(
                        code="invalid_site_key",
                        message=str(exc),
                        dataset_id=dataset_id,
                        field_name="intensities.index",
                    )
                )
            else:
                site_keys = tuple(str(value) for value in intensities.index)

        try:
            require_numeric_dataframe(
                intensities,
                field_name=(f"sps_reference_dataset[{dataset_id!r}].intensities"),
                error_type=ReferenceValidationError,
            )
        except ReferenceValidationError as exc:
            issues.append(
                SpsDiscoveryValidationIssue(
                    code="non_numeric_intensity_data",
                    message=str(exc),
                    dataset_id=dataset_id,
                    field_name="intensities",
                )
            )
        else:
            complex_columns = tuple(
                str(column)
                for column in intensities.columns
                if pd.api.types.is_complex_dtype(intensities[column].dtype)
            )
            if complex_columns:
                issues.append(
                    SpsDiscoveryValidationIssue(
                        code="non_numeric_intensity_data",
                        message=(
                            "reference intensities must be real-valued; complex "
                            "columns are invalid: " + ", ".join(complex_columns)
                        ),
                        dataset_id=dataset_id,
                        field_name="intensities",
                    )
                )
            else:
                try:
                    require_finite_numeric_dataframe(
                        intensities,
                        field_name=(
                            f"sps_reference_dataset[{dataset_id!r}].intensities"
                        ),
                        error_type=ReferenceValidationError,
                        allow_missing=True,
                    )
                except ReferenceValidationError as exc:
                    issues.append(
                        SpsDiscoveryValidationIssue(
                            code="non_finite_intensity_data",
                            message=str(exc),
                            dataset_id=dataset_id,
                            field_name="intensities",
                        )
                    )

        conditions = dataset.condition_by_sample
        if not isinstance(conditions, Mapping):
            issues.append(
                SpsDiscoveryValidationIssue(
                    code="missing_condition_metadata",
                    message="condition_by_sample must be a sample-to-group mapping",
                    dataset_id=dataset_id,
                    field_name="condition_by_sample",
                )
            )
            return issues, site_keys

        invalid_sample_ids = tuple(
            column
            for column in intensities.columns
            if not isinstance(column, str) or not column.strip()
        )
        if invalid_sample_ids:
            issues.append(
                SpsDiscoveryValidationIssue(
                    code="invalid_sample_id",
                    message=(
                        "reference intensity sample columns must be non-empty strings"
                    ),
                    dataset_id=dataset_id,
                    field_name="intensities.columns",
                )
            )
        sample_ids = tuple(str(column) for column in intensities.columns)
        if len(set(sample_ids)) != len(sample_ids):
            issues.append(
                SpsDiscoveryValidationIssue(
                    code="duplicate_sample_id",
                    message="reference intensity sample columns must be unique",
                    dataset_id=dataset_id,
                    field_name="intensities.columns",
                )
            )
        condition_keys = tuple(conditions.keys())
        invalid_keys = tuple(key for key in condition_keys if not isinstance(key, str))
        invalid_values = tuple(
            key
            for key, value in conditions.items()
            if not isinstance(value, str) or not value.strip()
        )
        if invalid_keys or invalid_values:
            issues.append(
                SpsDiscoveryValidationIssue(
                    code="invalid_condition_metadata",
                    message=(
                        "condition_by_sample keys and group labels must be "
                        "non-empty strings"
                    ),
                    dataset_id=dataset_id,
                    field_name="condition_by_sample",
                )
            )
        else:
            missing_samples = sorted(set(sample_ids).difference(condition_keys))
            extra_samples = sorted(set(condition_keys).difference(sample_ids))
            if missing_samples or extra_samples:
                issues.append(
                    SpsDiscoveryValidationIssue(
                        code="unresolvable_condition_metadata",
                        message=(
                            "condition_by_sample must match intensity sample columns "
                            f"exactly; missing={missing_samples}, extra={extra_samples}"
                        ),
                        dataset_id=dataset_id,
                        field_name="condition_by_sample",
                    )
                )
            if len({str(value).strip() for value in conditions.values()}) < 2:
                issues.append(
                    SpsDiscoveryValidationIssue(
                        code="insufficient_condition_groups",
                        message=(
                            "condition_by_sample must resolve at least two distinct "
                            "condition/group labels"
                        ),
                        dataset_id=dataset_id,
                        field_name="condition_by_sample",
                    )
                )

        return issues, site_keys


__all__: list[str] = []
