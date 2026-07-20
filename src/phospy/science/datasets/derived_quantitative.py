"""Derived analysis-ready quantitative dataset models."""

from __future__ import annotations

from collections import Counter
from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.errors.validation import DatasetValidationError
from phospy.provenance.derived_quantitative import (
    DerivedQuantitativeDataProvenance,
)
from phospy.provenance.hashing import fingerprint_optional_table_strict
from phospy.provenance.immutability import thaw_json_value
from phospy.provenance.models import RunProvenance, TableFingerprint
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.science.datasets.processing_state import DatasetProcessingState
from phospy.science.references.models import Organism
from phospy.science.transformations.models import IntensityScaleState

_CERTIFIED_PARENT_STATE_AUTHORITY = object()
_MATRIX_TRANSFORMATION_TABLES: Mapping[str, str] = {
    "phospho": "dataset.phospho",
    "total_protein": "dataset.total",
    "sample_metadata": "dataset.sample_metadata",
    "imputation_observation_mask": "dataset.imputation_observation_mask",
    "comparisons": "dataset.comparisons",
}
_TECHNICAL_REPLICATE_AGGREGATION = "technical_replicate_aggregation"
_TECHNICAL_REPLICATE_AGGREGATION_METHODS = frozenset(("mean", "median"))


@dataclass(frozen=True, slots=True, init=False)
class CertifiedDerivedQuantitativeParentState:
    """Certified parent state used to validate derived quantitative lineage."""

    parent_dataset_type: str
    parent_dataset_fingerprints: tuple[TableFingerprint, ...]
    phospho_sample_ids: tuple[str, ...]
    intensity_scale: str
    quantitative_meaning: str

    def __init__(
        self,
        *,
        parent_dataset_type: str,
        parent_dataset_fingerprints: tuple[TableFingerprint, ...],
        phospho_sample_ids: tuple[str, ...],
        intensity_scale: str,
        quantitative_meaning: str,
        _authority: object,
    ) -> None:
        if _authority is not _CERTIFIED_PARENT_STATE_AUTHORITY:
            raise PhosPyInputError(
                "CertifiedDerivedQuantitativeParentState must be created from an "
                "actual parent dataset or a certified internal snapshot"
            )
        object.__setattr__(
            self,
            "parent_dataset_type",
            _required_text(
                parent_dataset_type,
                field_name="derived_quantitative_parent_state.parent_dataset_type",
            ),
        )
        object.__setattr__(
            self,
            "parent_dataset_fingerprints",
            _require_table_fingerprint_tuple(
                parent_dataset_fingerprints,
                field_name=(
                    "derived_quantitative_parent_state.parent_dataset_fingerprints"
                ),
            ),
        )
        object.__setattr__(
            self,
            "phospho_sample_ids",
            _required_text_tuple(
                phospho_sample_ids,
                field_name="derived_quantitative_parent_state.phospho_sample_ids",
            ),
        )
        object.__setattr__(
            self,
            "intensity_scale",
            _required_text(
                intensity_scale,
                field_name="derived_quantitative_parent_state.intensity_scale",
            ),
        )
        object.__setattr__(
            self,
            "quantitative_meaning",
            _required_text(
                quantitative_meaning,
                field_name="derived_quantitative_parent_state.quantitative_meaning",
            ),
        )

    @classmethod
    def from_dataset(
        cls,
        dataset: AnalysisReadyPhosphoDataset,
    ) -> CertifiedDerivedQuantitativeParentState:
        """Capture certified state from an already-validated parent dataset."""

        if not isinstance(dataset, AnalysisReadyPhosphoDataset):
            raise PhosPyInputError(
                "derived quantitative parent_state must be captured from an "
                "AnalysisReadyPhosphoDataset"
            )
        phospho = dataset.phospho
        return cls._from_owned_parent_tables(
            parent_dataset_type=type(dataset).__name__,
            phospho=phospho,
            site_metadata=dataset.site_metadata,
            sample_metadata=dataset.sample_metadata,
            total=dataset.total,
            comparisons=dataset.comparisons,
            imputation_observation_mask=dataset.imputation_observed_mask_dataframe(),
            intensity_scale_state=dataset.intensity_scale_state,
        )

    @classmethod
    def _from_owned_parent_tables(
        cls,
        *,
        parent_dataset_type: str,
        phospho: pd.DataFrame,
        site_metadata: pd.DataFrame,
        intensity_scale_state: IntensityScaleState,
        sample_metadata: pd.DataFrame | None = None,
        total: pd.DataFrame | None = None,
        comparisons: pd.DataFrame | None = None,
        imputation_observation_mask: pd.DataFrame | None = None,
    ) -> CertifiedDerivedQuantitativeParentState:
        """Package-private certified snapshot path for internally owned tables."""

        return cls(
            parent_dataset_type=parent_dataset_type,
            parent_dataset_fingerprints=_fingerprints_for_table_state(
                phospho=phospho,
                site_metadata=site_metadata,
                sample_metadata=sample_metadata,
                total=total,
                comparisons=comparisons,
                imputation_observation_mask=imputation_observation_mask,
            ),
            phospho_sample_ids=tuple(str(label) for label in phospho.columns.tolist()),
            intensity_scale=str(intensity_scale_state.label),
            quantitative_meaning=_quantitative_meaning_for_intensity_state(
                intensity_scale_state
            ),
            _authority=_CERTIFIED_PARENT_STATE_AUTHORITY,
        )


class DerivedAnalysisReadyPhosphoDataset(AnalysisReadyPhosphoDataset):
    """Analysis-ready dataset whose quantitative matrices were derived from a parent."""

    __slots__ = ("derived_lineage", "derived_parent_state")

    def __init__(
        self,
        *,
        phospho: pd.DataFrame,
        site_metadata: pd.DataFrame,
        intensity_scale_state: IntensityScaleState,
        processing_state: DatasetProcessingState,
        parent_state: CertifiedDerivedQuantitativeParentState,
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
            parent_state=parent_state,
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
        parent_state: CertifiedDerivedQuantitativeParentState,
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
        if not isinstance(parent_state, CertifiedDerivedQuantitativeParentState):
            raise PhosPyInputError(
                "derived dataset requires CertifiedDerivedQuantitativeParentState"
            )
        if not isinstance(provenance, RunProvenance):
            raise PhosPyInputError("derived dataset requires RunProvenance")
        self._init_analysis_ready_tables(
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
            assume_owned=_assume_owned,
        )
        _require_lineage_matches_owned_state(
            dataset=self,
            parent_state=parent_state,
            derived_lineage=derived_lineage,
            provenance=provenance,
        )
        object.__setattr__(self, "derived_lineage", derived_lineage)
        object.__setattr__(self, "derived_parent_state", parent_state)

    @classmethod
    def _from_owned_derived_tables(
        cls,
        *,
        phospho: pd.DataFrame,
        site_metadata: pd.DataFrame,
        intensity_scale_state: IntensityScaleState,
        processing_state: DatasetProcessingState,
        parent_state: CertifiedDerivedQuantitativeParentState,
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
            parent_state=parent_state,
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
    parent_state: CertifiedDerivedQuantitativeParentState,
    derived_lineage: DerivedQuantitativeDataProvenance,
    provenance: RunProvenance,
) -> None:
    actual_fingerprints = _fingerprints_for_owned_derived_dataset(dataset)
    _require_sample_mapping_matches_actual_phospho(
        dataset=dataset,
        derived_lineage=derived_lineage,
    )
    _require_parent_identity_matches_certified_state(
        parent_state=parent_state,
        derived_lineage=derived_lineage,
    )
    _require_sample_mapping_matches_certified_parent(
        parent_state=parent_state,
        derived_lineage=derived_lineage,
    )
    _require_quantitative_state_matches_certified_boundaries(
        parent_state=parent_state,
        dataset=dataset,
        derived_lineage=derived_lineage,
    )
    _require_fingerprint_sets_match(
        expected=derived_lineage.derived_dataset_fingerprints,
        actual=actual_fingerprints,
        field_name="derived_quantitative_data.derived_dataset_fingerprints",
        expected_source="actual derived dataset tables",
    )
    _require_matrices_transformed_match_table_state(
        parent_state=parent_state,
        actual_derived_fingerprints=actual_fingerprints,
        derived_lineage=derived_lineage,
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
    return _fingerprints_for_table_state(
        phospho=dataset._phospho,
        site_metadata=dataset._site_metadata,
        sample_metadata=dataset._sample_metadata,
        total=dataset._total,
        comparisons=dataset._comparisons,
        imputation_observation_mask=_owned_imputation_observation_mask(dataset),
    )


def _fingerprints_for_table_state(
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
        fingerprint = fingerprint_optional_table_strict(table, name=name)
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


def _require_parent_identity_matches_certified_state(
    *,
    parent_state: CertifiedDerivedQuantitativeParentState,
    derived_lineage: DerivedQuantitativeDataProvenance,
) -> None:
    if derived_lineage.parent_dataset_type != parent_state.parent_dataset_type:
        raise DatasetValidationError(
            "derived_quantitative_data.parent_dataset_type must match certified "
            "parent dataset state"
        )
    _require_fingerprint_sets_match(
        expected=derived_lineage.parent_dataset_fingerprints,
        actual=parent_state.parent_dataset_fingerprints,
        field_name="derived_quantitative_data.parent_dataset_fingerprints",
        expected_source="certified parent dataset state",
    )


def _require_sample_mapping_matches_certified_parent(
    *,
    parent_state: CertifiedDerivedQuantitativeParentState,
    derived_lineage: DerivedQuantitativeDataProvenance,
) -> None:
    parent_sample_ids = parent_state.phospho_sample_ids
    declared_input_sample_ids = tuple(
        sample_id
        for item in derived_lineage.sample_mapping
        for sample_id in item.input_sample_ids
    )
    missing_from_parent = sorted(
        set(declared_input_sample_ids) - set(parent_sample_ids)
    )
    if missing_from_parent:
        raise DatasetValidationError(
            "derived_quantitative_data.sample_mapping input_sample_ids must exist "
            "in certified parent dataset.phospho columns; missing: "
            + ", ".join(missing_from_parent)
        )
    duplicate_input_ids = sorted(
        sample_id
        for sample_id, count in Counter(declared_input_sample_ids).items()
        if count > 1
    )
    if duplicate_input_ids:
        raise DatasetValidationError(
            "derived_quantitative_data.sample_mapping input_sample_ids must not "
            "duplicate parent dataset.phospho columns for this operation; "
            "duplicates: " + ", ".join(duplicate_input_ids)
        )
    omitted_parent_ids = [
        sample_id
        for sample_id in parent_sample_ids
        if sample_id not in set(declared_input_sample_ids)
    ]
    if omitted_parent_ids:
        raise DatasetValidationError(
            "derived_quantitative_data.sample_mapping input_sample_ids must cover "
            "every certified parent dataset.phospho column for this operation; "
            "omitted: " + ", ".join(omitted_parent_ids)
        )
    _require_operation_specific_sample_mapping(derived_lineage=derived_lineage)


def _require_operation_specific_sample_mapping(
    *,
    derived_lineage: DerivedQuantitativeDataProvenance,
) -> None:
    if derived_lineage.derivation_type != _TECHNICAL_REPLICATE_AGGREGATION:
        return
    if (
        derived_lineage.aggregation_method
        not in _TECHNICAL_REPLICATE_AGGREGATION_METHODS
    ):
        raise DatasetValidationError(
            "technical_replicate_aggregation derived_quantitative_data."
            "aggregation_method must be 'mean' or 'median'"
        )
    group_keys = tuple(
        (item.condition, item.biological_replicate_id)
        for item in derived_lineage.sample_mapping
    )
    duplicate_group_keys = sorted(
        f"{condition}:{biological_replicate_id}"
        for (condition, biological_replicate_id), count in Counter(group_keys).items()
        if count > 1
    )
    if duplicate_group_keys:
        raise DatasetValidationError(
            "technical_replicate_aggregation sample_mapping must have one output "
            "per condition+biological_replicate_id group; duplicate groups: "
            + ", ".join(duplicate_group_keys)
        )
    biological_replicate_counts = Counter(
        item.biological_replicate_id for item in derived_lineage.sample_mapping
    )
    for item in derived_lineage.sample_mapping:
        expected_output_sample_id = item.biological_replicate_id
        if biological_replicate_counts[item.biological_replicate_id] > 1:
            expected_output_sample_id = (
                f"{item.condition}__{item.biological_replicate_id}"
            )
        if item.output_sample_id != expected_output_sample_id:
            raise DatasetValidationError(
                "technical_replicate_aggregation sample_mapping output_sample_id "
                "must match the condition+biological_replicate_id grouping rule; "
                f"expected {expected_output_sample_id!r} for "
                f"condition={item.condition!r}, "
                f"biological_replicate_id={item.biological_replicate_id!r}"
            )
        duplicate_technical_ids = sorted(
            technical_id
            for technical_id, count in Counter(item.technical_replicate_ids).items()
            if count > 1
        )
        if duplicate_technical_ids:
            raise DatasetValidationError(
                "technical_replicate_aggregation sample_mapping "
                "technical_replicate_ids must be unique within each output "
                "sample; duplicates: " + ", ".join(duplicate_technical_ids)
            )
        if len(item.technical_replicate_ids) > len(item.input_sample_ids):
            raise DatasetValidationError(
                "technical_replicate_aggregation sample_mapping "
                "technical_replicate_ids cannot outnumber input_sample_ids"
            )


def _require_quantitative_state_matches_certified_boundaries(
    *,
    parent_state: CertifiedDerivedQuantitativeParentState,
    dataset: DerivedAnalysisReadyPhosphoDataset,
    derived_lineage: DerivedQuantitativeDataProvenance,
) -> None:
    output_intensity_scale = str(dataset.intensity_scale_state.label)
    output_quantitative_meaning = _quantitative_meaning_for_intensity_state(
        dataset.intensity_scale_state
    )
    if derived_lineage.input_intensity_scale != parent_state.intensity_scale:
        raise DatasetValidationError(
            "derived_quantitative_data.input_intensity_scale must match certified "
            "parent dataset intensity scale"
        )
    if derived_lineage.output_intensity_scale != output_intensity_scale:
        raise DatasetValidationError(
            "derived_quantitative_data.output_intensity_scale must match actual "
            "derived dataset intensity scale"
        )
    if derived_lineage.quantitative_meaning != parent_state.quantitative_meaning:
        raise DatasetValidationError(
            "derived_quantitative_data.quantitative_meaning must match certified "
            "parent dataset quantitative meaning"
        )
    if derived_lineage.quantitative_meaning != output_quantitative_meaning:
        raise DatasetValidationError(
            "derived_quantitative_data.quantitative_meaning must match actual "
            "derived dataset quantitative meaning"
        )


def _require_matrices_transformed_match_table_state(
    *,
    parent_state: CertifiedDerivedQuantitativeParentState,
    actual_derived_fingerprints: tuple[TableFingerprint, ...],
    derived_lineage: DerivedQuantitativeDataProvenance,
) -> None:
    expected = _actual_matrices_transformed(
        parent_fingerprints=parent_state.parent_dataset_fingerprints,
        derived_fingerprints=actual_derived_fingerprints,
    )
    declared = dict(derived_lineage.matrices_transformed)
    missing_keys = sorted(set(expected) - set(declared))
    unexpected_keys = sorted(set(declared) - set(expected))
    if missing_keys or unexpected_keys:
        detail_parts: list[str] = []
        if missing_keys:
            detail_parts.append("missing keys: " + ", ".join(missing_keys))
        if unexpected_keys:
            detail_parts.append("unexpected keys: " + ", ".join(unexpected_keys))
        raise DatasetValidationError(
            "derived_quantitative_data.matrices_transformed keys must match "
            "certified parent and actual derived table state; "
            + "; ".join(detail_parts)
        )
    mismatched = sorted(
        key
        for key, expected_value in expected.items()
        if declared[key] != expected_value
    )
    if mismatched:
        details = ", ".join(f"{key} expected {expected[key]!r}" for key in mismatched)
        raise DatasetValidationError(
            "derived_quantitative_data.matrices_transformed must match certified "
            "parent and actual derived table state; mismatched flags: " + details
        )


def _actual_matrices_transformed(
    *,
    parent_fingerprints: tuple[TableFingerprint, ...],
    derived_fingerprints: tuple[TableFingerprint, ...],
) -> dict[str, bool]:
    return {
        key: (
            _fingerprint_by_name(parent_fingerprints, table_name)
            != _fingerprint_by_name(derived_fingerprints, table_name)
        )
        for key, table_name in _MATRIX_TRANSFORMATION_TABLES.items()
    }


def _fingerprint_by_name(
    fingerprints: tuple[TableFingerprint, ...],
    name: str,
) -> TableFingerprint | None:
    for fingerprint in fingerprints:
        if fingerprint.name == name:
            return fingerprint
    return None


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


def _quantitative_meaning_for_intensity_state(
    intensity_scale_state: IntensityScaleState,
) -> str:
    quantity = intensity_scale_state.quantity
    if quantity is None:
        return "unknown"
    return quantity.value


def _require_table_fingerprint_tuple(
    values: object,
    *,
    field_name: str,
) -> tuple[TableFingerprint, ...]:
    if isinstance(values, (str, bytes, bytearray)) or not isinstance(values, Sequence):
        raise PhosPyInputError(f"{field_name} must be a sequence")
    fingerprints = tuple(values)
    if not fingerprints:
        raise PhosPyInputError(f"{field_name} must not be empty")
    for fingerprint in fingerprints:
        if not isinstance(fingerprint, TableFingerprint):
            raise PhosPyInputError(
                f"{field_name} must contain only TableFingerprint values"
            )
    return fingerprints


def _required_text(value: object, *, field_name: str) -> str:
    if not isinstance(value, str):
        raise PhosPyInputError(f"{field_name} must be a non-empty string")
    normalized = value.strip()
    if not normalized:
        raise PhosPyInputError(f"{field_name} must be a non-empty string")
    return normalized


def _required_text_tuple(values: object, *, field_name: str) -> tuple[str, ...]:
    if isinstance(values, (str, bytes, bytearray)) or not isinstance(values, Sequence):
        raise PhosPyInputError(f"{field_name} must be a sequence of strings")
    result = tuple(
        _required_text(value, field_name=f"{field_name}[]") for value in values
    )
    if not result:
        raise PhosPyInputError(f"{field_name} must not be empty")
    return result


__all__ = [
    "CertifiedDerivedQuantitativeParentState",
    "DerivedAnalysisReadyPhosphoDataset",
]
