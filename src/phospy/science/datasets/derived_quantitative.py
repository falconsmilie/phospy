"""Derived analysis-ready quantitative dataset models."""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.errors.validation import DatasetValidationError
from phospy.provenance.derived_quantitative import (
    DerivedQuantitativeDataProvenance,
)
from phospy.provenance.hashing import fingerprint_optional_table
from phospy.provenance.immutability import thaw_json_value
from phospy.provenance.models import RunProvenance, TableFingerprint
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.science.datasets.processing_state import DatasetProcessingState
from phospy.science.references.models import Organism
from phospy.science.transformations.models import IntensityScaleState


class DerivedAnalysisReadyPhosphoDataset(AnalysisReadyPhosphoDataset):
    """Analysis-ready dataset whose quantitative matrices were derived from a parent."""

    __slots__ = ("derived_lineage",)

    def __init__(
        self,
        *,
        phospho: pd.DataFrame,
        site_metadata: pd.DataFrame,
        intensity_scale_state: IntensityScaleState,
        processing_state: DatasetProcessingState,
        derived_lineage: DerivedQuantitativeDataProvenance,
        provenance: RunProvenance,
        sample_metadata: pd.DataFrame | None = None,
        total: pd.DataFrame | None = None,
        comparisons: pd.DataFrame | None = None,
        imputation_observation_mask: pd.DataFrame | None = None,
        organism: Organism | None = None,
        allow_opaque_site_values: bool = False,
    ) -> None:
        self._init_from_derived_tables(
            phospho=phospho,
            site_metadata=site_metadata,
            intensity_scale_state=intensity_scale_state,
            processing_state=processing_state,
            derived_lineage=derived_lineage,
            provenance=provenance,
            sample_metadata=sample_metadata,
            total=total,
            comparisons=comparisons,
            imputation_observation_mask=imputation_observation_mask,
            organism=organism,
            allow_opaque_site_values=allow_opaque_site_values,
            _assume_owned=False,
        )

    def _init_from_derived_tables(
        self,
        *,
        phospho: pd.DataFrame,
        site_metadata: pd.DataFrame,
        intensity_scale_state: IntensityScaleState,
        processing_state: DatasetProcessingState,
        derived_lineage: DerivedQuantitativeDataProvenance,
        provenance: RunProvenance,
        sample_metadata: pd.DataFrame | None = None,
        total: pd.DataFrame | None = None,
        comparisons: pd.DataFrame | None = None,
        imputation_observation_mask: pd.DataFrame | None = None,
        organism: Organism | None = None,
        allow_opaque_site_values: bool = False,
        _assume_owned: bool,
    ) -> None:
        if not isinstance(derived_lineage, DerivedQuantitativeDataProvenance):
            raise PhosPyInputError(
                "derived dataset requires DerivedQuantitativeDataProvenance"
            )
        if not isinstance(provenance, RunProvenance):
            raise PhosPyInputError("derived dataset requires RunProvenance")
        super().__init__(
            phospho=phospho,
            site_metadata=site_metadata,
            intensity_scale_state=intensity_scale_state,
            processing_state=processing_state,
            sample_metadata=sample_metadata,
            total=total,
            comparisons=comparisons,
            imputation_observation_mask=imputation_observation_mask,
            organism=organism,
            preprocessing_report=None,
            protein_aware_preparation=None,
            provenance=provenance,
            allow_opaque_site_values=allow_opaque_site_values,
            _assume_owned=_assume_owned,
        )
        _require_lineage_matches_owned_state(
            dataset=self,
            derived_lineage=derived_lineage,
            provenance=provenance,
        )
        object.__setattr__(self, "derived_lineage", derived_lineage)

    @classmethod
    def _from_owned_derived_tables(
        cls,
        *,
        phospho: pd.DataFrame,
        site_metadata: pd.DataFrame,
        intensity_scale_state: IntensityScaleState,
        processing_state: DatasetProcessingState,
        derived_lineage: DerivedQuantitativeDataProvenance,
        provenance: RunProvenance,
        sample_metadata: pd.DataFrame | None = None,
        total: pd.DataFrame | None = None,
        comparisons: pd.DataFrame | None = None,
        imputation_observation_mask: pd.DataFrame | None = None,
        organism: Organism | None = None,
        allow_opaque_site_values: bool = False,
    ) -> DerivedAnalysisReadyPhosphoDataset:
        """Package-private owned transfer for internally produced derived tables."""

        dataset = cls.__new__(cls)
        dataset._init_from_derived_tables(
            phospho=phospho,
            site_metadata=site_metadata,
            intensity_scale_state=intensity_scale_state,
            processing_state=processing_state,
            sample_metadata=sample_metadata,
            total=total,
            comparisons=comparisons,
            imputation_observation_mask=imputation_observation_mask,
            organism=organism,
            provenance=provenance,
            derived_lineage=derived_lineage,
            allow_opaque_site_values=allow_opaque_site_values,
            _assume_owned=True,
        )
        return dataset


def _require_lineage_matches_owned_state(
    *,
    dataset: DerivedAnalysisReadyPhosphoDataset,
    derived_lineage: DerivedQuantitativeDataProvenance,
    provenance: RunProvenance,
) -> None:
    actual_fingerprints = _fingerprints_for_owned_derived_dataset(dataset)
    _require_sample_mapping_matches_actual_phospho(
        dataset=dataset,
        derived_lineage=derived_lineage,
    )
    _require_fingerprint_sets_match(
        expected=derived_lineage.derived_dataset_fingerprints,
        actual=actual_fingerprints,
        field_name="derived_quantitative_data.derived_dataset_fingerprints",
        expected_source="actual derived dataset tables",
    )
    _require_fingerprint_sets_match(
        expected=derived_lineage.parent_dataset_fingerprints,
        actual=provenance.input_tables,
        field_name="run_provenance.input_tables",
        expected_source="derived_lineage.parent_dataset_fingerprints",
    )
    _require_fingerprint_sets_match(
        expected=derived_lineage.derived_dataset_fingerprints,
        actual=provenance.output_tables,
        field_name="run_provenance.output_tables",
        expected_source="derived_lineage.derived_dataset_fingerprints",
    )
    if provenance.workflow_name != derived_lineage.derivation_type:
        raise DatasetValidationError(
            "derived dataset run_provenance.workflow_name must match "
            "derived_lineage.derivation_type"
        )
    _require_run_provenance_embeds_lineage(
        provenance=provenance,
        derived_lineage=derived_lineage,
    )


def _fingerprints_for_owned_derived_dataset(
    dataset: DerivedAnalysisReadyPhosphoDataset,
) -> tuple[TableFingerprint, ...]:
    entries: tuple[tuple[str, pd.DataFrame | None], ...] = (
        ("dataset.phospho", dataset._phospho),
        ("dataset.site_metadata", dataset._site_metadata),
        ("dataset.sample_metadata", dataset._sample_metadata),
        ("dataset.total", dataset._total),
        ("dataset.comparisons", dataset._comparisons),
        (
            "dataset.imputation_observation_mask",
            _owned_imputation_observation_mask(dataset),
        ),
    )
    return _collect_table_fingerprints(entries)


def _owned_imputation_observation_mask(
    dataset: DerivedAnalysisReadyPhosphoDataset,
) -> pd.DataFrame | None:
    metadata = dataset._imputation_observation_metadata
    if metadata is None:
        return None
    return metadata.observed_mask_dataframe()


def _collect_table_fingerprints(
    entries: Sequence[tuple[str, pd.DataFrame | None]],
) -> tuple[TableFingerprint, ...]:
    fingerprints: list[TableFingerprint] = []
    for name, table in entries:
        fingerprint = fingerprint_optional_table(table, name=name)
        if fingerprint is None:
            continue
        fingerprints.append(fingerprint)
    return tuple(fingerprints)


def _require_sample_mapping_matches_actual_phospho(
    *,
    dataset: DerivedAnalysisReadyPhosphoDataset,
    derived_lineage: DerivedQuantitativeDataProvenance,
) -> None:
    expected_columns = tuple(
        item.output_sample_id for item in derived_lineage.sample_mapping
    )
    actual_columns = tuple(str(label) for label in dataset._phospho.columns.tolist())
    if expected_columns != actual_columns:
        raise DatasetValidationError(
            "derived_quantitative_data.sample_mapping output_sample_id values "
            "must match actual derived dataset.phospho columns"
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
            detail_parts.append("missing actual tables: " + ", ".join(missing))
        if unexpected:
            detail_parts.append("unexpected actual tables: " + ", ".join(unexpected))
        raise DatasetValidationError(
            f"{field_name} table set must match {expected_source}; "
            + "; ".join(detail_parts)
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
            f"{field_name} does not match {expected_source}; "
            "mismatched fields: " + ", ".join(mismatched)
        )


def _require_run_provenance_embeds_lineage(
    *,
    provenance: RunProvenance,
    derived_lineage: DerivedQuantitativeDataProvenance,
) -> None:
    payload = provenance.workflow_parameters.get("derived_quantitative_data")
    if not isinstance(payload, Mapping):
        raise DatasetValidationError(
            "derived dataset run_provenance.workflow_parameters must include "
            "typed derived_quantitative_data lineage"
        )
    try:
        restored = DerivedQuantitativeDataProvenance.from_payload(
            _thaw_mapping(payload)
        )
    except PhosPyInputError as exc:
        raise DatasetValidationError(
            "derived dataset run_provenance.workflow_parameters "
            "derived_quantitative_data must decode as typed lineage"
        ) from exc
    if restored.to_payload() != derived_lineage.to_payload():
        raise DatasetValidationError(
            "derived dataset run_provenance.workflow_parameters "
            "derived_quantitative_data must match derived_lineage"
        )


def _thaw_mapping(payload: Mapping[str, object]) -> Mapping[str, object]:
    thawed = thaw_json_value(payload, field_name="derived_dataset.lineage")
    if not isinstance(thawed, Mapping):
        raise DatasetValidationError(
            "derived dataset lineage payload must thaw to a JSON object"
        )
    return thawed


__all__ = ["DerivedAnalysisReadyPhosphoDataset"]
