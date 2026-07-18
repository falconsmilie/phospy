"""Technical-replicate planning and aggregation for differential workflows."""

from __future__ import annotations

import warnings
from collections import Counter
from collections.abc import Mapping
from dataclasses import dataclass

import pandas as pd

from phospy.errors.validation import WorkflowValidationError
from phospy.provenance.derived_quantitative import (
    TECHNICAL_REPLICATE_AGGREGATION_DERIVATION_TYPE,
    TECHNICAL_REPLICATE_AGGREGATOR_IMPLEMENTATION,
    DerivedQuantitativeDataProvenance,
    DerivedSampleMapping,
    build_derived_quantitative_run_provenance,
)
from phospy.provenance.environment import collect_environment_provenance
from phospy.provenance.hashing import fingerprint_optional_table
from phospy.provenance.models import JsonValue, TableFingerprint
from phospy.science.datasets.derived_quantitative import (
    CertifiedDerivedQuantitativeParentState,
    DerivedAnalysisReadyPhosphoDataset,
)
from phospy.science.datasets.internal_view import DatasetInternalView
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.science.design.models import ExperimentalDesign, SampleDesignRecord
from phospy.science.differential.policy_models import TechnicalReplicatePolicy

_TECHNICAL_REPLICATE_RESOLVER_DEPRECATION_MESSAGE = (
    "TechnicalReplicateResolver is deprecated and will be removed in a future "
    "release. Use TechnicalReplicateAggregationPlanner to create an explicit "
    "TechnicalReplicateAggregationPlan, then apply it with "
    "TechnicalReplicateAggregator."
)


@dataclass(frozen=True, slots=True)
class TechnicalReplicateAggregationGroup:
    condition: str
    biological_replicate_id: str
    output_sample_id: str
    input_sample_ids: tuple[str, ...]
    technical_replicate_ids: tuple[str, ...]
    batch: str | None
    block_id: str | None
    covariates: Mapping[str, str | float]


@dataclass(frozen=True, slots=True)
class TechnicalReplicateAggregationPlan:
    """Explicit technical-replicate aggregation plan for workflow preparation."""

    technical_replicate_policy: TechnicalReplicatePolicy
    groups: tuple[TechnicalReplicateAggregationGroup, ...]
    aggregate_phospho: bool
    aggregate_total_protein: bool

    @property
    def requires_aggregation(self) -> bool:
        return bool(self.groups)


@dataclass(frozen=True, slots=True)
class TechnicalReplicateResolution:
    """Resolved dataset/design pair after technical-replicate aggregation."""

    dataset: AnalysisReadyPhosphoDataset
    design: ExperimentalDesign
    workflow_provenance: Mapping[str, object] | None = None


class TechnicalReplicateAggregationPlanner:
    """Build an explicit aggregation plan from technical-replicate policy."""

    def run(
        self,
        *,
        dataset: AnalysisReadyPhosphoDataset,
        design: ExperimentalDesign,
        technical_replicate_policy: TechnicalReplicatePolicy,
    ) -> TechnicalReplicateAggregationPlan:
        if not isinstance(technical_replicate_policy, TechnicalReplicatePolicy):
            raise WorkflowValidationError(
                "differential workflow request technical_replicate_policy must be "
                "TechnicalReplicatePolicy"
            )

        repeated_groups = self._find_repeated_biological_groups(design)
        if not repeated_groups:
            return TechnicalReplicateAggregationPlan(
                technical_replicate_policy=technical_replicate_policy,
                groups=(),
                aggregate_phospho=False,
                aggregate_total_protein=False,
            )

        if technical_replicate_policy == TechnicalReplicatePolicy.REJECT:
            rendered_groups = ", ".join(
                f"{condition}:{biological_replicate_id}"
                for condition, biological_replicate_id in repeated_groups
            )
            raise WorkflowValidationError(
                "experimental design contains repeated biological_replicate_id values "
                "within condition groups (condition:biological_replicate_id="
                f"{rendered_groups}). Technical replicates require explicit aggregation; "
                "set technical_replicate_policy='mean' or "
                "technical_replicate_policy='median'."
            )

        groups = self._build_resolved_groups(design=design)
        self._validate_required_sample_ids_exist(
            dataset=dataset,
            groups=groups,
        )
        dataset_view = DatasetInternalView(dataset)
        return TechnicalReplicateAggregationPlan(
            technical_replicate_policy=technical_replicate_policy,
            groups=groups,
            aggregate_phospho=True,
            aggregate_total_protein=dataset_view.total is not None,
        )

    @staticmethod
    def _validate_required_sample_ids_exist(
        *,
        dataset: AnalysisReadyPhosphoDataset,
        groups: tuple[TechnicalReplicateAggregationGroup, ...],
    ) -> None:
        phospho = DatasetInternalView(dataset).phospho
        required_sample_ids = [
            sample_id for group in groups for sample_id in group.input_sample_ids
        ]
        missing = sorted(set(required_sample_ids) - set(phospho.columns.astype(str)))
        if missing:
            raise WorkflowValidationError(
                "technical replicate aggregation design samples are missing from "
                "dataset.phospho columns: " + ", ".join(missing)
            )

    @staticmethod
    def _find_repeated_biological_groups(
        design: ExperimentalDesign,
    ) -> tuple[tuple[str, str], ...]:
        group_counts: Counter[tuple[str, str]] = Counter()
        for record in design.samples:
            if record.biological_replicate_id is None:
                continue
            group_counts[(record.condition, record.biological_replicate_id)] += 1
        repeated = sorted(key for key, count in group_counts.items() if count > 1)
        return tuple(repeated)

    @staticmethod
    def _build_resolved_groups(
        *, design: ExperimentalDesign
    ) -> tuple[TechnicalReplicateAggregationGroup, ...]:
        missing_biological_id_samples = [
            record.sample_id
            for record in design.samples
            if record.biological_replicate_id is None
        ]
        if missing_biological_id_samples:
            raise WorkflowValidationError(
                "technical_replicate_policy aggregation requires "
                "biological_replicate_id on every design sample; missing for: "
                + ", ".join(missing_biological_id_samples)
            )

        grouped_records: dict[tuple[str, str], list[SampleDesignRecord]] = {}
        for record in design.samples:
            biological_replicate_id = record.biological_replicate_id
            if biological_replicate_id is None:
                raise WorkflowValidationError(
                    "technical_replicate_policy aggregation requires "
                    "biological_replicate_id on every design sample"
                )
            group_key = (record.condition, biological_replicate_id)
            if group_key not in grouped_records:
                grouped_records[group_key] = []
            grouped_records[group_key].append(record)

        biological_id_counts: Counter[str] = Counter(
            biological_replicate_id
            for _, biological_replicate_id in grouped_records.keys()
        )
        resolved_groups: list[TechnicalReplicateAggregationGroup] = []
        for (condition, biological_replicate_id), records in grouped_records.items():
            output_sample_id = biological_replicate_id
            if biological_id_counts[biological_replicate_id] > 1:
                output_sample_id = f"{condition}__{biological_replicate_id}"
            technical_ids = tuple(
                technical_id
                for technical_id in dict.fromkeys(
                    record.technical_replicate_id for record in records
                )
                if technical_id is not None
            )
            batch = (
                TechnicalReplicateAggregationPlanner._require_consistent_optional_field(
                    records=tuple(records),
                    field_name="batch",
                )
            )
            block_id = (
                TechnicalReplicateAggregationPlanner._require_consistent_optional_field(
                    records=tuple(records),
                    field_name="block_id",
                )
            )
            covariates = (
                TechnicalReplicateAggregationPlanner._require_consistent_covariates(
                    records=tuple(records)
                )
            )
            resolved_groups.append(
                TechnicalReplicateAggregationGroup(
                    condition=condition,
                    biological_replicate_id=biological_replicate_id,
                    output_sample_id=output_sample_id,
                    input_sample_ids=tuple(record.sample_id for record in records),
                    technical_replicate_ids=technical_ids,
                    batch=batch,
                    block_id=block_id,
                    covariates=covariates,
                )
            )
        return tuple(resolved_groups)

    @staticmethod
    def _require_consistent_optional_field(
        *,
        records: tuple[SampleDesignRecord, ...],
        field_name: str,
    ) -> str | None:
        values = [getattr(record, field_name) for record in records]
        first = values[0]
        if any(value != first for value in values):
            sample_ids = ", ".join(record.sample_id for record in records)
            biological_replicate_id = records[0].biological_replicate_id
            raise WorkflowValidationError(
                "technical replicate aggregation requires consistent "
                f"{field_name!r} within each condition+biological_replicate_id group; "
                f"condition={records[0].condition!r}, "
                f"biological_replicate_id={biological_replicate_id!r}, "
                f"samples={sample_ids}"
            )
        return first

    @staticmethod
    def _require_consistent_covariates(
        *,
        records: tuple[SampleDesignRecord, ...],
    ) -> Mapping[str, str | float]:
        first = dict(records[0].covariates)
        if any(dict(record.covariates) != first for record in records[1:]):
            sample_ids = ", ".join(record.sample_id for record in records)
            biological_replicate_id = records[0].biological_replicate_id
            raise WorkflowValidationError(
                "technical replicate aggregation requires consistent covariates "
                "within each condition+biological_replicate_id group; "
                f"condition={records[0].condition!r}, "
                f"biological_replicate_id={biological_replicate_id!r}, "
                f"samples={sample_ids}"
            )
        return first


class TechnicalReplicateAggregator:
    """Apply an explicit technical-replicate aggregation plan."""

    def run(
        self,
        *,
        dataset: AnalysisReadyPhosphoDataset,
        design: ExperimentalDesign,
        aggregation_plan: TechnicalReplicateAggregationPlan,
    ) -> TechnicalReplicateResolution:
        if not isinstance(aggregation_plan, TechnicalReplicateAggregationPlan):
            raise WorkflowValidationError(
                "technical_replicate_aggregation_plan must be "
                "TechnicalReplicateAggregationPlan"
            )
        if not aggregation_plan.requires_aggregation:
            return TechnicalReplicateResolution(
                dataset=dataset,
                design=design,
                workflow_provenance=None,
            )

        groups = aggregation_plan.groups
        technical_replicate_policy = aggregation_plan.technical_replicate_policy
        aggregated_dataset = self._aggregate_dataset(
            dataset=dataset,
            groups=groups,
            technical_replicate_policy=technical_replicate_policy,
        )
        aggregated_design = ExperimentalDesign(
            samples=tuple(
                SampleDesignRecord(
                    sample_id=group.output_sample_id,
                    condition=group.condition,
                    biological_replicate_id=group.biological_replicate_id,
                    technical_replicate_id=None,
                    batch=group.batch,
                    block_id=group.block_id,
                    covariates=group.covariates,
                )
                for group in groups
            ),
            fixed_effects=design.fixed_effects,
        )
        provenance_groups: list[dict[str, object]] = []
        for group in groups:
            provenance_groups.append(
                {
                    "condition": group.condition,
                    "biological_replicate_id": group.biological_replicate_id,
                    "output_sample_id": group.output_sample_id,
                    "input_sample_ids": list(group.input_sample_ids),
                    "source_sample_ids": list(group.input_sample_ids),
                    "technical_replicate_ids": list(group.technical_replicate_ids),
                    "n_technical_replicates": len(group.input_sample_ids),
                    "aggregation_method": technical_replicate_policy.value,
                }
            )
        grouped_samples = [group.output_sample_id for group in groups]
        source_samples = sorted(
            {sample_id for group in groups for sample_id in group.input_sample_ids}
        )
        aggregate_total = aggregation_plan.aggregate_total_protein
        derived_lineage = (
            aggregated_dataset.derived_lineage
            if isinstance(aggregated_dataset, DerivedAnalysisReadyPhosphoDataset)
            else None
        )
        workflow_provenance: dict[str, object] = {
            "technical_replicate_policy": technical_replicate_policy.value,
            "aggregation_policy": technical_replicate_policy.value,
            "aggregation_method": technical_replicate_policy.value,
            "grouped_samples": grouped_samples,
            "source_samples": source_samples,
            "matrices_aggregated": {
                "phospho": True,
                "total_protein": aggregate_total,
            },
            "both_phospho_and_total_aggregated": bool(aggregate_total),
            "groups": provenance_groups,
        }
        if derived_lineage is not None:
            workflow_provenance["derived_quantitative_data"] = (
                derived_lineage.to_payload()
            )
            workflow_provenance["derived_dataset_type"] = type(
                aggregated_dataset
            ).__name__
        return TechnicalReplicateResolution(
            dataset=aggregated_dataset,
            design=aggregated_design,
            workflow_provenance=workflow_provenance,
        )

    @staticmethod
    def _aggregate_dataset(
        *,
        dataset: AnalysisReadyPhosphoDataset,
        groups: tuple[TechnicalReplicateAggregationGroup, ...],
        technical_replicate_policy: TechnicalReplicatePolicy,
    ) -> AnalysisReadyPhosphoDataset:
        dataset_view = DatasetInternalView(dataset)
        phospho = dataset_view.phospho
        parent_state = CertifiedDerivedQuantitativeParentState.from_dataset(dataset)
        aggregated_phospho = TechnicalReplicateAggregator._aggregate_numeric_matrix(
            matrix=phospho,
            groups=groups,
            technical_replicate_policy=technical_replicate_policy,
        )
        total = dataset_view.total
        aggregated_total = None
        if total is not None:
            aggregated_total = TechnicalReplicateAggregator._aggregate_numeric_matrix(
                matrix=total,
                groups=groups,
                technical_replicate_policy=technical_replicate_policy,
            )
        sample_metadata = dataset_view.sample_metadata
        aggregated_sample_metadata = None
        if sample_metadata is not None:
            aggregated_sample_metadata = (
                TechnicalReplicateAggregator._aggregate_metadata(
                    metadata=sample_metadata,
                    groups=groups,
                )
            )
        comparisons = dataset_view.comparisons
        aggregated_imputation_observation_mask = (
            dataset_view.aggregate_imputation_observation_mask(
                sample_groups=tuple(
                    (group.output_sample_id, group.input_sample_ids) for group in groups
                )
            )
        )
        parent_fingerprints = _collect_table_fingerprints(
            (
                ("dataset.phospho", phospho),
                ("dataset.site_metadata", dataset_view.site_metadata),
                ("dataset.sample_metadata", sample_metadata),
                ("dataset.total", total),
                ("dataset.comparisons", comparisons),
                (
                    "dataset.imputation_observation_mask",
                    dataset.imputation_observed_mask_dataframe(),
                ),
            )
        )
        derived_fingerprints = _collect_table_fingerprints(
            (
                ("dataset.phospho", aggregated_phospho),
                ("dataset.site_metadata", dataset_view.site_metadata),
                ("dataset.sample_metadata", aggregated_sample_metadata),
                ("dataset.total", aggregated_total),
                ("dataset.comparisons", comparisons),
                (
                    "dataset.imputation_observation_mask",
                    aggregated_imputation_observation_mask,
                ),
            )
        )
        environment = collect_environment_provenance()
        quantitative_meaning = dataset.intensity_scale_state.quantity
        lineage = DerivedQuantitativeDataProvenance(
            derivation_type=TECHNICAL_REPLICATE_AGGREGATION_DERIVATION_TYPE,
            parent_dataset_type=type(dataset).__name__,
            derived_dataset_type="DerivedAnalysisReadyPhosphoDataset",
            parent_dataset_fingerprints=parent_fingerprints,
            derived_dataset_fingerprints=derived_fingerprints,
            sample_mapping=tuple(
                DerivedSampleMapping(
                    output_sample_id=group.output_sample_id,
                    input_sample_ids=group.input_sample_ids,
                    condition=group.condition,
                    biological_replicate_id=group.biological_replicate_id,
                    technical_replicate_ids=group.technical_replicate_ids,
                )
                for group in groups
            ),
            aggregation_method=technical_replicate_policy.value,
            input_intensity_scale=str(dataset.intensity_scale_state.label),
            output_intensity_scale=str(dataset.intensity_scale_state.label),
            quantitative_meaning=(
                "unknown"
                if quantitative_meaning is None
                else quantitative_meaning.value
            ),
            missingness_policy=_missingness_policy_payload(dataset),
            matrices_transformed={
                "phospho": True,
                "total_protein": aggregated_total is not None,
                "sample_metadata": aggregated_sample_metadata is not None,
                "imputation_observation_mask": (
                    aggregated_imputation_observation_mask is not None
                ),
                "comparisons": False,
            },
            implementation=TECHNICAL_REPLICATE_AGGREGATOR_IMPLEMENTATION,
            implementation_version=environment.package_version,
            parameters={
                "aggregation_axis": "samples",
                "source_grouping": "condition+biological_replicate_id",
            },
        )
        provenance = build_derived_quantitative_run_provenance(
            lineage=lineage,
            environment=environment,
            reference_context=dataset.reference_context,
        )
        return DerivedAnalysisReadyPhosphoDataset._from_owned_derived_tables(  # pyright: ignore[reportPrivateUsage] - workflow owns freshly derived tables and lineage before domain certification
            phospho=aggregated_phospho,
            site_metadata=dataset_view.site_metadata.copy(deep=True),
            intensity_scale_state=dataset.intensity_scale_state,
            processing_state=dataset.processing_state,
            sample_metadata=aggregated_sample_metadata,
            total=None if aggregated_total is None else aggregated_total,
            comparisons=None if comparisons is None else comparisons.copy(deep=True),
            imputation_observation_mask=aggregated_imputation_observation_mask,
            organism=dataset.organism,
            parent_state=parent_state,
            provenance=provenance,
            derived_lineage=lineage,
            allow_opaque_site_values=dataset.opaque_site_values_allowed,
        )

    @staticmethod
    def _aggregate_numeric_matrix(
        *,
        matrix: pd.DataFrame,
        groups: tuple[TechnicalReplicateAggregationGroup, ...],
        technical_replicate_policy: TechnicalReplicatePolicy,
    ) -> pd.DataFrame:
        aggregated_columns: list[pd.Series] = []
        output_labels: list[str] = []
        for group in groups:
            source = matrix.loc[:, list(group.input_sample_ids)]
            if technical_replicate_policy == TechnicalReplicatePolicy.MEAN:
                collapsed = source.mean(axis=1)
            elif technical_replicate_policy == TechnicalReplicatePolicy.MEDIAN:
                collapsed = source.median(axis=1)
            else:
                raise WorkflowValidationError(
                    "technical_replicate_policy must be 'mean' or 'median' for "
                    "aggregation"
                )
            collapsed.name = group.output_sample_id
            aggregated_columns.append(collapsed.astype(float))
            output_labels.append(group.output_sample_id)
        aggregated = pd.concat(aggregated_columns, axis=1)
        aggregated.columns = pd.Index(output_labels)
        aggregated.index = matrix.index.copy()
        return aggregated

    @staticmethod
    def _aggregate_metadata(
        *,
        metadata: pd.DataFrame,
        groups: tuple[TechnicalReplicateAggregationGroup, ...],
    ) -> pd.DataFrame:
        rows: list[pd.Series] = []
        for group in groups:
            rows.append(metadata.loc[group.input_sample_ids[0], :].copy(deep=True))
        aggregated = pd.DataFrame(
            rows,
            index=pd.Index(
                [group.output_sample_id for group in groups],
                name=metadata.index.name,
            ),
            columns=metadata.columns.copy(),
        )
        return aggregated


def _collect_table_fingerprints(
    entries: tuple[tuple[str, pd.DataFrame | None], ...],
) -> tuple[TableFingerprint, ...]:
    fingerprints: list[TableFingerprint] = []
    for name, table in entries:
        fingerprint = fingerprint_optional_table(table, name=name)
        if fingerprint is None:
            continue
        fingerprints.append(fingerprint)
    return tuple(fingerprints)


def _missingness_policy_payload(
    dataset: AnalysisReadyPhosphoDataset,
) -> dict[str, JsonValue]:
    missing_data = dataset.processing_state.missing_data
    diagnostics = missing_data.diagnostics
    imputation_method_id: str | None = None
    if diagnostics is not None:
        raw_imputation_method_id = diagnostics.get("imputation_method_id")
        imputation_method_id = (
            None if raw_imputation_method_id is None else str(raw_imputation_method_id)
        )
    return {
        "policy": missing_data.policy.value,
        "complete_matrix": bool(missing_data.complete_matrix),
        "imputed": bool(missing_data.imputed),
        "has_missing_values": missing_data.has_missing_values,
        "missing_value_count": missing_data.missing_value_count,
        "imputation_method_id": imputation_method_id,
        "numeric_aggregation_skip_missing": True,
        "imputation_observation_mask_aggregation": "all_source_cells_observed",
    }


class TechnicalReplicateResolver:
    """Deprecated wrapper that plans and then applies aggregation."""

    def __init__(
        self,
        *,
        planner: TechnicalReplicateAggregationPlanner | None = None,
        aggregator: TechnicalReplicateAggregator | None = None,
    ) -> None:
        warnings.warn(
            _TECHNICAL_REPLICATE_RESOLVER_DEPRECATION_MESSAGE,
            DeprecationWarning,
            stacklevel=2,
        )
        self._planner = planner or TechnicalReplicateAggregationPlanner()
        self._aggregator = aggregator or TechnicalReplicateAggregator()

    def run(
        self,
        *,
        dataset: AnalysisReadyPhosphoDataset,
        design: ExperimentalDesign,
        technical_replicate_policy: TechnicalReplicatePolicy,
    ) -> TechnicalReplicateResolution:
        warnings.warn(
            _TECHNICAL_REPLICATE_RESOLVER_DEPRECATION_MESSAGE,
            DeprecationWarning,
            stacklevel=2,
        )
        plan = self._planner.run(
            dataset=dataset,
            design=design,
            technical_replicate_policy=technical_replicate_policy,
        )
        return self._aggregator.run(
            dataset=dataset,
            design=design,
            aggregation_plan=plan,
        )


__all__ = [
    "TechnicalReplicateAggregationGroup",
    "TechnicalReplicateAggregationPlan",
    "TechnicalReplicateAggregationPlanner",
    "TechnicalReplicateAggregator",
    "TechnicalReplicateResolution",
    "TechnicalReplicateResolver",
]
