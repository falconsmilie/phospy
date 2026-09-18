"""Private validation for SPS discovery reference inputs."""

from __future__ import annotations

from collections import Counter
from collections.abc import Mapping, Sequence
from typing import cast

import pandas as pd

from phospy.errors.validation import ReferenceValidationError
from phospy.frames.validation import (
    require_finite_numeric_dataframe,
    require_numeric_dataframe,
)
from phospy.science.batch_correction.sps_discovery import (
    SPS_REFERENCE_MATRIX_FINGERPRINT_PARAMETER,
    SpsDiscoveryConfig,
    SpsDiscoveryValidationIssue,
    SpsDiscoveryValidationResult,
    SpsReferenceDataset,
    _is_accepted_sps_meaning_evidence,
    _sps_meaning_evidence_matches_matrix,
    _sps_reference_intensity_fingerprint,
    _sps_scale_evidence_matches_matrix,
)
from phospy.science.sites.site_keys import decode_site_key
from phospy.science.sites.validation import require_site_key_index
from phospy.science.transformations.models import (
    IntensityScaleEvidenceLevel,
    IntensityScaleKind,
    IntensityScaleState,
    QuantitativeMeaning,
    QuantitativeMeaningEvidenceMode,
)

_SPS_QUANTITATIVE_REQUIREMENT = (
    "The PhosR-style SPS statistic requires pre-established condition-relative "
    "log2 measurements centred on the reference/control baseline."
)


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
        organisms_by_dataset: list[tuple[str, str]] = []
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
            if dataset.organism is not None:
                organisms_by_dataset.append(
                    (dataset.dataset_id, dataset.organism.value)
                )
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

        distinct_organisms = sorted({organism for _, organism in organisms_by_dataset})
        if len(distinct_organisms) > 1:
            declarations = ", ".join(
                f"{dataset_id}={organism}"
                for dataset_id, organism in organisms_by_dataset
            )
            issues.append(
                SpsDiscoveryValidationIssue(
                    code="incompatible_reference_organisms",
                    message=(
                        "SPS discovery requires biologically coherent references "
                        "from one organism; observed " + declarations
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
        required_context = (
            (
                "organism",
                dataset.organism,
                "missing_reference_organism",
                "reference organism must be declared with a supported Organism",
            ),
            (
                "baseline_context",
                dataset.baseline_context,
                "missing_baseline_context",
                (
                    "reference baseline_context must identify the biological "
                    "reference/control baseline represented by zero"
                ),
            ),
            (
                "reference_context",
                dataset.reference_context,
                "missing_reference_context",
                (
                    "reference_context must describe the biological context in "
                    "which the reference evidence was produced"
                ),
            ),
        )
        for field_name, value, code, message in required_context:
            if value is None:
                issues.append(
                    SpsDiscoveryValidationIssue(
                        code=code,
                        message=message,
                        dataset_id=dataset_id,
                        field_name=field_name,
                    )
                )
        missing_source_fields = tuple(
            field_name
            for field_name in ("source_name", "source_version", "source_uri")
            if getattr(dataset, field_name) is None
        )
        if missing_source_fields:
            issues.append(
                SpsDiscoveryValidationIssue(
                    code="missing_source_identity",
                    message=(
                        "reference source identity must be reconstructable; missing "
                        + ", ".join(missing_source_fields)
                    ),
                    dataset_id=dataset_id,
                    field_name="source_identity",
                )
            )
        intensities = dataset._intensities_snapshot()
        issues.extend(
            self._validate_quantitative_state(dataset, intensities=intensities)
        )
        if not isinstance(intensities, pd.DataFrame):
            issues.append(
                SpsDiscoveryValidationIssue(
                    code="missing_intensity_data",
                    message="reference intensities must be a pandas DataFrame",
                    dataset_id=dataset_id,
                    field_name="intensities",
                )
            )
            return issues, ()
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
                encoded_organisms = {
                    decode_site_key(
                        site_key,
                        field_name=(
                            f"sps_reference_dataset[{dataset_id!r}].intensities.index"
                        ),
                        error_type=ReferenceValidationError,
                    ).organism
                    for site_key in site_keys
                }
                if dataset.organism is not None and encoded_organisms != {
                    dataset.organism
                }:
                    observed = ", ".join(
                        sorted(organism.value for organism in encoded_organisms)
                    )
                    issues.append(
                        SpsDiscoveryValidationIssue(
                            code="reference_organism_site_key_mismatch",
                            message=(
                                "reference organism must match every encoded "
                                "site_key organism; declared "
                                f"{dataset.organism.value!r}, observed {observed!r}"
                            ),
                            dataset_id=dataset_id,
                            field_name="intensities.index",
                        )
                    )

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

    def _validate_quantitative_state(
        self,
        dataset: SpsReferenceDataset,
        *,
        intensities: object,
    ) -> list[SpsDiscoveryValidationIssue]:
        state = dataset.intensity_scale_state
        dataset_id = dataset.dataset_id
        field_name = "intensity_scale_state"
        if not isinstance(state, IntensityScaleState):
            return [
                SpsDiscoveryValidationIssue(
                    code="unknown_quantitative_state",
                    message=(
                        f"{_SPS_QUANTITATIVE_REQUIREMENT} Supply a governed "
                        "IntensityScaleState; numeric values are never used to "
                        "infer SPS quantitative semantics."
                    ),
                    dataset_id=dataset_id,
                    field_name=field_name,
                )
            ]

        issues: list[SpsDiscoveryValidationIssue] = []
        establishment = state.establishment_provenance
        if not state.is_established or establishment is None:
            issues.append(
                SpsDiscoveryValidationIssue(
                    code="missing_scale_establishment_evidence",
                    message=(
                        f"{_SPS_QUANTITATIVE_REQUIREMENT} The supplied scale state "
                        "does not carry supported establishment provenance."
                    ),
                    dataset_id=dataset_id,
                    field_name=field_name,
                )
            )
        elif establishment.evidence_level is IntensityScaleEvidenceLevel.UNKNOWN:
            issues.append(
                SpsDiscoveryValidationIssue(
                    code="missing_scale_establishment_evidence",
                    message=(
                        f"{_SPS_QUANTITATIVE_REQUIREMENT} Unknown scale evidence "
                        "does not establish the required log2 state."
                    ),
                    dataset_id=dataset_id,
                    field_name=field_name,
                )
            )

        if state.kind is not IntensityScaleKind.LOG2:
            issues.append(
                SpsDiscoveryValidationIssue(
                    code="incompatible_sps_intensity_scale",
                    message=(
                        f"{_SPS_QUANTITATIVE_REQUIREMENT} Received governed scale "
                        f"{state.kind.value!r}."
                    ),
                    dataset_id=dataset_id,
                    field_name=field_name,
                )
            )

        if state.quantity is not QuantitativeMeaning.CONTRAST_LOG2_FOLD_CHANGE:
            meaning = None if state.quantity is None else state.quantity.value
            code = (
                "absolute_abundance_quantitative_state"
                if state.quantity
                in {
                    QuantitativeMeaning.PHOSPHOSITE_ABUNDANCE,
                    QuantitativeMeaning.PHOSPHOSITE_LOG_ABUNDANCE,
                }
                else "incompatible_sps_quantitative_meaning"
            )
            issues.append(
                SpsDiscoveryValidationIssue(
                    code=code,
                    message=(
                        f"{_SPS_QUANTITATIVE_REQUIREMENT} Received quantitative "
                        f"meaning {meaning!r}; SPS discovery does not centre or "
                        "transform absolute-abundance input."
                    ),
                    dataset_id=dataset_id,
                    field_name=field_name,
                )
            )

        meaning_provenance = state.quantitative_meaning_provenance
        if meaning_provenance is None:
            issues.append(
                SpsDiscoveryValidationIssue(
                    code="missing_baseline_establishment_evidence",
                    message=(
                        f"{_SPS_QUANTITATIVE_REQUIREMENT} Quantitative-meaning "
                        "provenance must establish the reference-centred contrast."
                    ),
                    dataset_id=dataset_id,
                    field_name=field_name,
                )
            )
        elif not _is_accepted_sps_meaning_evidence(meaning_provenance):
            evidence_mode = cast(
                QuantitativeMeaningEvidenceMode,
                meaning_provenance.evidence_mode,
            )
            issues.append(
                SpsDiscoveryValidationIssue(
                    code="missing_baseline_establishment_evidence",
                    message=(
                        f"{_SPS_QUANTITATIVE_REQUIREMENT} Evidence mode "
                        f"{evidence_mode.value!r} does not "
                        "establish that reference/control centring occurred."
                    ),
                    dataset_id=dataset_id,
                    field_name=field_name,
                )
            )
        elif meaning_provenance.target_quantity is not state.quantity:
            issues.append(
                SpsDiscoveryValidationIssue(
                    code="inconsistent_quantitative_state_evidence",
                    message=(
                        f"{_SPS_QUANTITATIVE_REQUIREMENT} Quantitative-meaning "
                        "provenance does not match the resolved matrix meaning."
                    ),
                    dataset_id=dataset_id,
                    field_name=field_name,
                )
            )
        if (
            isinstance(intensities, pd.DataFrame)
            and not intensities.index.has_duplicates
            and not intensities.columns.has_duplicates
        ):
            intensity_fingerprint = _sps_reference_intensity_fingerprint(
                dataset_id=dataset_id,
                intensities=intensities,
            )
            if establishment is not None and not _sps_scale_evidence_matches_matrix(
                establishment,
                intensity_fingerprint=intensity_fingerprint,
            ):
                binding = establishment.parameters.get(
                    SPS_REFERENCE_MATRIX_FINGERPRINT_PARAMETER
                )
                code = (
                    "missing_scale_matrix_binding_evidence"
                    if binding is None
                    else "mismatched_scale_matrix_binding_evidence"
                )
                issues.append(
                    SpsDiscoveryValidationIssue(
                        code=code,
                        message=(
                            f"{_SPS_QUANTITATIVE_REQUIREMENT} Log2 scale evidence "
                            "must be fingerprint-bound to this submitted SPS "
                            "reference matrix; construct the reference with "
                            "SpsReferenceDataset.from_condition_relative_log2(...)."
                        ),
                        dataset_id=dataset_id,
                        field_name=field_name,
                    )
                )
            if meaning_provenance is not None and not (
                _sps_meaning_evidence_matches_matrix(
                    meaning_provenance,
                    intensity_fingerprint=intensity_fingerprint,
                )
            ):
                code = (
                    "missing_baseline_matrix_binding_evidence"
                    if meaning_provenance.output_table_fingerprint is None
                    else "mismatched_baseline_matrix_binding_evidence"
                )
                issues.append(
                    SpsDiscoveryValidationIssue(
                        code=code,
                        message=(
                            f"{_SPS_QUANTITATIVE_REQUIREMENT} Reference/control "
                            "centering evidence must be fingerprint-bound to this "
                            "submitted SPS reference matrix; evidence from another "
                            "matrix cannot be reused."
                        ),
                        dataset_id=dataset_id,
                        field_name=field_name,
                    )
                )
        return issues


__all__: list[str] = []
